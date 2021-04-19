library(presto)
library(dplyr)
library(Seurat)
source("Script_functions.R")

# extract epithelial stem cells from each organ 
## cells are within the age range comparable among organs, i.e. age.range <- c(11,14,15,17,18,19)
multi.epi <- readRDS("Res_combined_age_range_cell_type_selected_epi_all_organ.rds")
selected.ct <- c("Stem_cell","Antrum stem cell", "Corpus stem cell","Bud tip progenitor","Basal cell","KRT15 high/SERPINF+ basal cell")
sc.combined <- subset(multi.epi, cells=colnames(multi.epi)[which(multi.epi$Cell_type%in%selected.ct)])
sc.combined <- FindVariableFeatures(object = sc.combined, selection.method = "vst", nfeatures = 2000)
sc.combined <- ScaleData(object = sc.combined, verbose = T)
sc.combined <- RunPCA(sc.combined)
usefulPCs <- 1:10
sc.combined <- RunUMAP(sc.combined, dims = usefulPCs)
DimPlot(sc.combined, reduction = "umap", group.by = "Tissue", label=T)
saveRDS(sc.combined, file="Res_stem_cells_from_various_organs.rds")

# identify DE genes across organs and tissues
## DE genes between organs
de.res <- wilcoxauc(sc.combined, group_by = "Organ")
sig.res <- de.res[which(de.res$logFC>0.25 & de.res$auc>0.6 & de.res$pct_in>25 & de.res$padj<0.01),]
top.res <- sig.res %>% group_by(group) %>% top_n(50, wt=auc)
top.genes <- unique(top.res$feature)
de.res.list <- list("sig.res"=sig.res, "top.res"=top.res, "top.genes"=top.genes)
saveRDS(de.res.list, file="Res_DE_between_stem_cell_of_various_organs.rds")

## DE genes between different small intestine regions
cells <- colnames(sc.combined)[which(sc.combined@meta.data$Organ=="SI")]
X <- sc.combined@assays$RNA@data[,cells]
y <- sc.combined@meta.data[cells, "Tissue"]
de.res <- wilcoxauc(X = X, y = y)
sig.res <- de.res[which(de.res$logFC>0.25 & de.res$auc>0.6 & de.res$pct_in>25 & de.res$padj<0.01),]
table(sig.res$group)
top.res <- sig.res %>% group_by(group) %>% top_n(50, wt=auc)
top.genes <- unique(top.res$feature)
de.res.list <- list("sig.res"=sig.res, "top.res"=top.res, "top.genes"=top.genes)
saveRDS(de.res.list, file="Res_DE_between_stem_cell_of_various_SI_regions.rds")

## DE genes between different stomach regions
cells <- colnames(sc.combined)[which(sc.combined@meta.data$Organ=="Stomach")]
X <- sc.combined@assays$RNA@data[,cells]
y <- sc.combined@meta.data[cells, "Tissue"]
de.res <- wilcoxauc(X = X, y = y)
sig.res <- de.res[which(de.res$logFC>0.25 & de.res$auc>0.6 & de.res$pct_in>25 & de.res$padj<0.01),]
table(sig.res$group)
top.res <- sig.res %>% group_by(group) %>% top_n(50, wt=auc)
top.genes <- unique(top.res$feature)
de.res.list <- list("sig.res"=sig.res, "top.res"=top.res, "top.genes"=top.genes)
saveRDS(de.res.list, file="Res_DE_between_stem_cell_of_various_stomach_regions.rds")

stomach.res <- de.res.list
g.stomach <- stomach.res$top.genes

si.res <- readRDS("Res_DE_between_stem_cell_of_various_SI_regions.rds")
g.si <- si.res$top.genes

organ.res <- readRDS("Res_DE_between_stem_cell_of_various_organs.rds")
g.organ <- organ.res$top.genes 

selected.hvg <- unique(c(g.organ, g.stomach, g.si))
length(selected.hvg)
de.list <- list("all"=selected.hvg, "between-organ"=g.organ, "between-si-region"=g.si, "between-stomach-region"=g.stomach)
saveRDS(de.list, file="Res_categorical_comparison_DE_genes_between_tissues.rds")

de.res <- readRDS("Res_categorical_comparison_DE_genes_between_tissues.rds")
selected.hvg <- de.res$all
write.table(selected.hvg, file="List_DEG_for_stem_cell_spatial_deconvolution.txt", quote=F, row.names = F,
            col.names = F)
de.res <- readRDS("Res_DE_between_stem_cell_of_various_SI_regions.rds")
top.res <- de.res$top.res
write.table(top.res, file="Table_DE_between_stem_cell_of_various_SI_regions.txt", sep="\t",
            quote=F, row.names = F)
de.res <- readRDS("Res_DE_between_stem_cell_of_various_stomach_regions.rds")
top.res <- de.res$top.res
write.table(top.res, file="Table_DE_between_stem_cell_of_various_stomach_regions.txt", sep="\t",
            quote=F, row.names = F)
de.res <- readRDS("Res_DE_between_stem_cell_of_various_organs.rds")
top.res <- de.res$top.res
write.table(top.res, file="Table_DE_between_stem_cell_of_various_organs.txt", sep="\t",
            quote=F, row.names = F)

# use quadratic programming to resolve the stem cell state transition across tissues
## get reference expression pattern
## use the lung bud tip progenitor and colon stem cells to represent the extreme of stem cell heterogeneity across tissues  
cells.1 <- colnames(sc.combined)[which(sc.combined@meta.data$Cell_type=="Bud tip progenitor")]
cells.2 <- colnames(sc.combined)[which(sc.combined@meta.data$Organ=="Colon")]
expr.1 <- apply(as.matrix(sc.combined@assays$RNA@data[selected.hvg, cells.1]), 1, mean)
expr.2 <- apply(as.matrix(sc.combined@assays$RNA@data[selected.hvg, cells.2]), 1, mean)
X <- cbind(expr.1, expr.2)
saveRDS(X,file="Dat_lung_and_colon_stem_cell_expr_as_ref_for_qp.rds")
que.expr <- as.matrix(sc.combined@assays$RNA@data[selected.hvg, ])
## run quadratic programming
library(quadprog)
identity.que <- matrix(NA, nrow=ncol(que.expr), ncol=4)
rownames(identity.que) <- colnames(que.expr)
colnames(identity.que) <- c("frxn_lung", "frxn_colon", "Lagrangian", "Error")
for (j in seq(ncol(que.expr))){
  Y <- as.matrix(que.expr[,j])
  Rinv <- solve(chol(t(X) %*% X));
  C <- cbind(rep(1,2), diag(2))
  b <- c(1,rep(0,2))
  d <- t(Y) %*% X  
  QP<-solve.QP(Dmat = Rinv, factorized = TRUE, dvec = d, Amat = C, bvec = b, meq = 1)
  Error<-sum(abs(Y-X%*%QP$solution))
  identity.que[j,] <- c(QP$solution[1],QP$solution[2],QP$Lagrangian[1],Error)
}

# use colon identity fraction as a pseudospace representation, and identify genes showing variable expression along the pseudodtime
colon.prop <- identity.que[,2]
pt.vec=colon.prop
expr.mat=as.matrix(sc.combined@assays$RNA@data)
## order cells according to the pseudospace score
## group every 20 cells into a cell bin
cell.num.per.bin=20
tissue.orders <- c("Lung-distal", "Lung-tracheal-epi", "Lung-airway-trachea", "Lung-airway", "Esophagus", "Stomach-corpus", "Stomach-antrum",
                   "Duodenum", "Jejunum", "Ileum", "Colon")
tissue.vec <- sc.combined@meta.data$Tissue
start.idx <- seq(from=1, to=length(pt.vec), by=cell.num.per.bin)
end.idx <- c(start.idx[-1]-1,length(pt.vec))
idx <- cbind(start.idx, end.idx)
pt.rank <- rank(pt.vec, ties.method = "first")
tissue.num.by.pt.bin <- sapply(seq(nrow(idx)), function(i){
  idx.start <- idx[i,1]
  idx.end <- idx[i,2]
  idx.cell <- which(pt.rank>=idx.start & pt.rank<=idx.end)
  sapply(tissue.orders, function(x){ sum(tissue.vec[idx.cell]==x) })
})
prop.by.bin <- t(t(tissue.num.by.pt.bin)/apply(tissue.num.by.pt.bin, 2, sum))
saveRDS(prop.by.bin, file="Res_tissue_proportion_of_each_pseudocell_bin.rds")

## use stacked bar plot to represent tissue composition of each cell bin 
pdf("Plot_stacked_barplot_tissue_proportion_by_pt_bin-33.pdf")
barplot(prop.by.bin, col=tissue.cols[rownames(prop.by.bin)], border = NA, space = 0)
dev.off()

## get cell bin average expression levels
expr.by.pt.bin <- getExprByPt(pt.vec=pt.vec, expr.mat=expr.mat, cell.num.per.bin=cell.num.per.bin)
saveRDS(expr.by.pt.bin, file="Res_expr_bin_ordered_by_colon_identity_fractino_20cells_per_bin.rds")

## identify genes that show variable expression levels across the pseudospace trajectory
age.test.p <- splineBasedAgeTest(pseudotime.vec=seq(ncol(expr.by.pt.bin)), expr.mat=expr.by.pt.bin, df=5, mode="stringent")
padj <- p.adjust(age.test.p, method="bonferroni")
logfc <- apply(expr.by.pt.bin, 1, function(vec){
  max(vec)-min(vec)
})
age.test.res <- data.frame("Residual_comparison_P"=age.test.p, "Bonferroni_corrected_P"=padj, "Log(max-min)"=logfc, stringsAsFactors = F)
saveRDS(age.test.res, file="Res_age_test_along_colon_identity_fraction_pseudospace.rds")  
sc.genes <- rownames(age.test.res)[which(age.test.res$Bonferroni_corrected_P<0.01 & age.test.res$Log.max.min.>1)]
saveRDS(sc.genes, file="Res_pseudospace_variable_genes.rds")
write.table(sc.genes, file="List_pseudospace_variable_genes.txt", quote=F, row.names = F, col.names = F)

# focus on annotated transcription factors
human.tf <- read.table("~/Work/Annotation/HumanTFDB/List_HumanTFDB_TF.txt", sep="\t", stringsAsFactors = F, head=T)
sc.tf <- intersect(sc.genes, human.tf$Symbol)
length(sc.genes)
length(sc.tf)
sc.pt.genes <- list("all"=sc.genes, "tf"=sc.tf)
sapply(sc.pt.genes, length)
saveRDS(sc.pt.genes, file="Res_lung_to_colon_axis_patterning_related_genes.rds")

# visualize expression of selected genes in heatmap
selected.genes <- union(sc.tf, intersect(selected.hvg, human.tf$Symbol))
saveRDS(selected.genes, file="Res_DE_TF.rds")
sc.expr <- expr.by.pt.bin[selected.genes,]
hc.row <- hclust(as.dist(1-cor(t(sc.expr))), method="ward.D2")
mat <- sc.expr[hc.row$labels[hc.row$order],]
max.idx <- apply(mat, 1, which.max)
gene.order.list <- lapply(sort(unique(max.idx)), function(i){
  rownames(mat)[which(max.idx==i)]
})
gene.orders <- unlist(gene.order.list)
input <- mat[gene.orders,]
input.norm <- t(apply(input, 1, function(vec){
  (vec-min(vec))/(max(vec)-min(vec))
}))

colors <- colorRampPalette(c("#0571b0","#92c5de", "#f7f7f7", "#f4a582", "#ca0020"))(50)
library(gplots)
pdf("Plot_heatmap_basal_cell_markers.pdf")
heatmap.2(input.norm, trace="none",main="", density.info="none",dendrogram="none",scale="none", 
          col=colors,Rowv=FALSE, Colv=FALSE, cexRow=0.5, cexCol = 0.5, co)
dev.off()
