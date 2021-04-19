library(Seurat)
library(destiny)
library(presto)
source("Script_functions.R")

# intestinal epithelial stem cell maturation analysis
## load data
hio <- readRDS("Res_HIO_before_transplantation_with_DE.rds")
thio.epi <- readRDS("Res_tHIO_epi_on_intestine_lineage_with_RSS_based_UMAP.rds")
adult.epi <- readRDS("Res_adult_epi.rds")

dir.create("stem_cell_maturation")
setwd("stem_cell_maturation/")
## combine the intestinal epithelial cells of HIO, stem cells of tHIO, fetal and adult duodenum
## exclude mitochondrial genes, ribosomal genes and sex chromosome genes before normalization
### get counts
hio.cells <- colnames(hio)[which(hio@meta.data$Selected_for_epi_comparison_with_tHIO_and_vivo)]
hio.epi <- subset(hio, cells=hio.cells)
hio.counts <- hio.epi@assays$RNA@counts
colnames(hio.counts) <- paste(colnames(hio.counts), "hio", sep="_")

thio.cells <- colnames(thio.epi)[which(thio.epi@meta.data$Cell_type=="Stem cell")]
thio.counts <- thio.epi@assays$RNA@counts[,thio.cells]
colnames(thio.counts) <- paste(colnames(thio.counts), "thio", sep="_")

fetal.cells <- colnames(duo.epi)[which(duo.epi@meta.data$Cell_type=="Stem_cell")]
fetal.counts <- duo.epi@assays$RNA@counts[,fetal.cells]
colnames(fetal.counts) <- paste(colnames(fetal.counts), "fetal", sep="_")

adult.cells <- colnames(adult.epi)[which(adult.epi@meta.data$Cell_type=="Stem-cell")]
adult.counts <- adult.epi@assays$RNA@counts[,adult.cells]
colnames(adult.counts) <- paste(colnames(adult.counts), "adult", sep="_")

### get meta.data
hio.meta <- hio.epi@meta.data[,c("Age", "Tissue", "RNA_snn_res.0.6")]
hio.meta[,3] <- paste0("Cluster", hio.meta[,3])
rownames(hio.meta) <- paste(rownames(hio.meta),"hio",sep="_")
colnames(hio.meta)[3] <- "Cell_type" 

thio.meta <- thio.epi@meta.data[thio.cells,c("Age", "Tissue", "Cell_type")]
thio.meta$Age <- paste0("d", (as.numeric(sub("w", "", thio.meta$Age))+4)*7)
rownames(thio.meta) <- paste(rownames(thio.meta),"thio",sep="_")

fetal.meta <- duo.epi@meta.data[fetal.cells,c( "Age", "Tissue", "Cell_type")]
fetal.meta$Tissue <- paste("Fetal", fetal.meta$Tissue, sep="_")
rownames(fetal.meta) <- paste(rownames(fetal.meta),"fetal",sep="_")

adult.meta <- adult.epi@meta.data[adult.cells,c("Age", "Tissue", "Cell_type")]
adult.meta$Tissue <- paste("Adult", adult.meta$Tissue, sep="_")
rownames(adult.meta) <- paste(rownames(adult.meta),"adult",sep="_")

### get list confounding genes
dir <- "~/Work/Annotation/confound_genes/" 
all.files <- list.files(dir, pattern=".txt")
confound.genes <- lapply(seq(length(all.files)), function(i){
  file <- paste0(dir, all.files[i])
  g <- readLines(file)
  return(g)
})
names(confound.genes) <- c("Cell_cycle", "Experiment_induced", "HB", "MT", "Ribosome", "Sex_chr") 
genes.to.remove <- unique(c(confound.genes[["MT"]], confound.genes[["Ribosome"]], confound.genes[["Sex_chr"]]))

seu.list <- list()
for(x in c("hio", "thio", "fetal","adult")){
  count.mat <- get(paste0(x,".counts"))
  count.mat <- count.mat[setdiff(rownames(count.mat), genes.to.remove),]
  meta.mat <- get(paste0(x,".meta"))
  print(paste(dim(count.mat)[1], dim(meta.mat)[1]))
  seu.obj <- CreateSeuratObject(counts=count.mat, meta=meta.mat)
  seu.obj <- NormalizeData(object = seu.obj, normalization.method = "LogNormalize", scale.factor = 1e4)
  seu.obj@meta.data$orig.ident <- paste(seu.obj@meta.data$orig.ident,x,sep="_")
  seu.list[[x]] <- seu.obj
}
combined <- merge(x=seu.list[[1]], y=seu.list[-1])
stage.vec <- sapply(colnames(combined), function(x){strsplit(x, "_")[[1]][3]})
unique(stage.vec)
combined@meta.data$Stage <- stage.vec
saveRDS(combined, file="Res_HIO_tHIO_fetal_adult_duo_stem_cell_LOG_RNA_data.rds")

## pseudotime analysis: resolve the heterogeneity of stem cells of fetal and adult duodenum using DiffusionMap
pt.vec <- rep(NA, ncol(combined))
names(pt.vec) <- colnames(combined)
### for fetal cells
fetal.cells <- colnames(combined)[which(combined@meta.data$Stage=="fetal")]
input.mat <- Embeddings(combined, reduction = "css")[fetal.cells,]
dm <- DiffusionMap(data=input.mat, distance="euclidean", k=20)
dcs <- dm@eigenvectors
dc.values <- rank(dcs[,1])
names(dc.values) <- fetal.cells
pt.vec[names(dc.values)] <- dc.values
age.values <- combined@meta.data[fetal.cells, "Age"]
## visually examine the relationship between sample ages and rank of diffusion component 1 (DC1)
par(mfrow=c(1,2))
plotFeature2(coor=Embeddings(combined, reduction = "umap_css")[fetal.cells, ], values=dc.values, main="DC1 rank")
plotFeature2(coor=Embeddings(combined, reduction = "umap_css")[fetal.cells,], values=age.values, main="Age", add.label=T)

### for adult
adult.cells <- colnames(combined)[which(combined@meta.data$Stage=="adult")]
cells.to.remove <- adult.cells[which(Embeddings(combined, reduction = "umap_css")[adult.cells, 2]>10)]
adult.cells <- setdiff(adult.cells, cells.to.remove)
input.mat <- Embeddings(combined, reduction = "css")[adult.cells,]
dm <- DiffusionMap(data=input.mat, distance="euclidean", k=20)
dcs <- dm@eigenvectors
dc.values <- rank(-dcs[,1])
names(dc.values) <- adult.cells
pt.vec[names(dc.values)] <- dc.values
combined@meta.data$Pt_by_group <- pt.vec[colnames(combined)]
saveRDS(combined, file="Res_HIO_tHIO_fetal_adult_duo_stem_cell_LOG_RNA_data.rds")

## identify genes showing variable expression along fetal pseudotime based on spline regression
## every 50 cells were grouped into a bin, and get the average expression of each cell bin
## run the test based on cell bin average expression levels
idx <- which(combined@meta.data$Stage=="fetal")
pt.vec <- combined@meta.data$Pt_by_group[idx]
expr.mat <- as.matrix(combined@assays$RNA@data[,idx])
expr.pt <- getExprByPt(pt.vec=pt.vec, expr.mat=expr.mat, cell.num.per.bin=50)
age.test.p <- splineBasedAgeTest(pseudotime.vec=seq(ncol(expr.pt)), expr.mat=expr.pt, df=5, mode="stringent")
padj <- p.adjust(age.test.p, method="BH")
logfc <- apply(expr.pt, 1, function(vec){
  max(vec)-min(vec)
})
res <- data.frame("Residual_comparison_P"=age.test.p, "BH_corrected_P"=padj, "Log(max-min)"=logfc, stringsAsFactors = F)
## only the genes with BH-corrected P < 0.01 and difference between maximum and minimum of spline interpolated values > 0.5 are considered as genes with significant expression changes along pseudotime, pseudotime-dependent genes 
fetal.pt.genes <- rownames(res)[which(res$BH_corrected_P<0.01 & res$Log.max.min.>0.5)]
## classify the fetal stem cell pseudotime-dependent genes into modules based on their expression pattern 
fetal.pt.genes <- pt.gene.list[["fetal"]]
fetal.bin.expr <- expr.by.pt.bin[["fetal"]][fetal.pt.genes,]
hc.fetal.pt.genes <- hclust(as.dist(1-cor(t(fetal.bin.expr), method="spearman")), method="ward.D2")
plot(hc.fetal.pt.genes, hang=-1, cex=0.1)
g=cutree(hc.fetal.pt.genes, k=8)
expr.pt <- fetal.bin.expr
pt.gene.modules <- plotClusterExprProfile(expr=expr.pt, time.vec=seq(ncol(expr.pt)), group.vec=rep("fetal", ncol(expr.pt)), cluster.vec=g, group.cols="#31a354", return.value=T, to.plot=T, plot.name="Plot_fetal_pt_gene_module_expr_profile.pdf", add.legend=F, legend.pos="topleft", cex.legend=2)
saveRDS(pt.gene.modules, file="Res_pseudotime_dependent_gene_module_expression_along_pseudotime.rds")

## only take those genes belonging to cluster 2 and cluster 4 which shows monotomous increase/decrease expression patterns alon pseudotime for following quadratic programming for HIO, tHIO and fetal stem cells, so as to characterize the stem cell state transition among these datasets
fetal.pt.genes <- names(g)[which(g%in%c(2,4))]
## prepare reference expression pattern for quadratic programming
## use HIO cluster 12 as the most immature state, and the day 132 as the most mature state in fetal
cells.1 <- colnames(combined)[which(combined@meta.data$Cell_type=="Cluster12")]
cells.2 <- colnames(combined)[which(combined@meta.data$Cell_type=="Stem_cell" & combined@meta.data$Age=="d132")]
expr.1 <- apply(as.matrix(combined@assays$RNA@data[fetal.pt.genes, cells.1]), 1, mean)
expr.2 <- apply(as.matrix(combined@assays$RNA@data[fetal.pt.genes, cells.2]), 1, mean)
X <- cbind(expr.1, expr.2)
## prepare expression matrix of query cells
que.cells <- colnames(combined)[which(combined@meta.data$Stage%in%c("hio", "thio", "fetal"))]
que.expr <- as.matrix(combined@assays$RNA@data[fetal.pt.genes, que.cells])
## run quadratic programming, and use the obtained day-132 fetal stem cell component fraction as a quantification of stem cell maturity
##library(quadprog)
identity.que <- matrix(NA, nrow=ncol(que.expr), ncol=4)
rownames(identity.que) <- colnames(que.expr)
colnames(identity.que) <- c("frxn_cl12", "frxn_fetLatSC", "Lagrangian", "Error")
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

## some estimated fraction proportion is reported to larger than 1 or less than 0, with negligible difference. this happened because the values are float, and it needs to be adjusted  
hio.fetal.qp <- identity.que
sc.prop <- hio.fetal.qp[,2]
sc.prop[which(sc.prop>1)] <- 1
sc.prop[which(sc.prop<0)] <- 0
hio.fetal.qp[["Float_adjusted_SC_prop"]] <- sc.prop
summary(sc.prop)
zero.cells <- rownames(hio.fetal.qp)[which(sc.prop==0)]
length(zero.cells)
zero.cells.2 <- sub("_hio", "", zero.cells)
table(hio@meta.data[zero.cells.2, "Age"])
age.vec <- hio@meta.data[zero.cells.2, "Age"]
names(age.vec) <- zero.cells
## for the cells with maturity=0, order them in the ascending order of organoid age
sorted.zero.cells <- unlist(lapply(c("d0", "d3" ,"d7"), function(x){
  names(age.vec)[which(age.vec==x)]
}))
non.zero.cells <- rownames(hio.fetal.qp)[which(sc.prop>0)]
length(sorted.zero.cells)+length(non.zero.cells)
dim(hio.fetal.qp)
hio.fetal.qp <- hio.fetal.qp[c(sorted.zero.cells, non.zero.cells),]
saveRDS(hio.fetal.qp, file="Res_organoid_fetal_QP_result.rds")

## use stem cell identity fraction as a pseudotime representing stem cell maturation, and identify genes showing variable expression along the pseudodtime
#### only take out fetal stem cells and HIO selected intestinal epithelial cells
cells <- colnames(combined)[which(combined$Stage%in%c("fetal", "hio"))]
age.vec <- combined@meta.data[cells, "Age"]
names(age.vec) <- cells

## order cells in the ascending order of stem cell maturity, and group every 100 cells into a cell bin 
ages <- unique(combined@meta.data[cells, "Age"])
age.order <- paste0("d", sort(as.numeric(sub("d","",ages))))
cells <- unlist(lapply(age.order, function(x){
  names(age.vec)[which(age.vec==x)]
}))
expr.mat <- as.matrix(combined@assays$RNA@data[,cells])
sc.prop <- hio.fetal.qp[cells,"Float_adjusted_SC_prop"]
names(sc.prop) <- cells
saveRDS(sc.prop, file="Res_fetal_and_HIO_stem_cell_ISC_maturity.rds")
cell.num.per.bin=100
expr.by.pt.bin <- getExprByPt(pt.vec=sc.prop, expr.mat=expr.mat, cell.num.per.bin=cell.num.per.bin)
saveRDS(expr.by.pt.bin, file="Res_HIO_and_fetal_stem_cell_expr_by_stem_cell_identity_fraction_bin.rds")

## get dominating age of each bin
### get cell index of each bin
pt.vec=sc.prop
start.idx <- seq(from=1, to=length(pt.vec), by=cell.num.per.bin)
end.idx <- c(start.idx[-1]-1,length(pt.vec))
idx <- cbind(start.idx, end.idx)
pt.rank <- rank(pt.vec, ties.method = "first")
bin.vec <- rep(NA, length(cells))
names(bin.vec) <- cells
for (i in seq(nrow(idx))){
  idx.start <- idx[i,1]
  idx.end <- idx[i,2]
  idx.cell <- which(pt.rank>=idx.start & pt.rank<=idx.end)
  bin.vec[idx.cell] <- i
}
saveRDS(bin.vec, file="Res_bin_assignment_for_HIO_and_fetal_SC.rds")

age.vec <- combined@meta.data[cells, "Age"]
freq.age.by.bin <- getMostFrequentGroupByPt(pt.vec=pt.vec, group.vec=age.vec, cell.num.per.bin=cell.num.per.bin)
names(freq.age.by.bin) <- paste0("bin_", seq(length(freq.age.by.bin)))
saveRDS(freq.age.by.bin, file="Res_most_frequent_age_of_each_cell_bin.rds")

## identify genes with variable expression levels along stem cell state transition from HIO to fetal, i.e. ISC development associated genes
age.test.p <- splineBasedAgeTest(pseudotime.vec=seq(ncol(expr.by.pt.bin)), expr.mat=expr.by.pt.bin, df=5, mode="stringent")
padj <- p.adjust(age.test.p, method="BH")
logfc <- apply(expr.by.pt.bin, 1, function(vec){
  max(vec)-min(vec)
})
age.test.res <- data.frame("Residual_comparison_P"=age.test.p, "BH_corrected_P"=padj, "Log(max-min)"=logfc, stringsAsFactors = F)
saveRDS(age.test.res, file="Res_age_test_along_HIO_and_fetal_stem_cell_identity_fraction_pseudotime.rds")  
## obtain ISC development associated genes based on filtering on BH-corrected P values and magnitude of expression changes
sc.genes <- rownames(age.test.res)[which(age.test.res$BH_corrected_P<0.01 & age.test.res$Log.max.min.>0.5)]

## add adult stem cells
selected.cells <- c(cells, colnames(combined)[which(combined$Stage=="adult")])
hfa <- subset(combined, cells=selected.cells)

sc.bin.vec <- rep(73, ncol(hfa))
names(sc.bin.vec) <- colnames(hfa)
sc.bin.vec[names(bin.vec)] <- bin.vec
hfa$Stem_cell_maturation_bin <- sc.bin.vec

## resolve the heterogeneity within adult stem cells
idx <- which(hfa$Stage=="adult")
expr.mat <- as.matrix(hfa@assays$RNA@data[,idx])
pt.vec <- hfa@meta.data$Pt_by_group[idx]
expr.by.pt.bin <- getExprByPt(pt.vec=pt.vec, expr.mat=expr.mat, cell.num.per.bin=50)
saveRDS(expr.by.pt.bin, file="Res_adult_stem_cell_expr_by_adult_stem_cell_pseudotime_bin.rds")

## get cell bin within adult stem cells
cells <- colnames(hfa)[idx]
cell.num.per.bin=50
start.idx <- seq(from=1, to=length(pt.vec), by=cell.num.per.bin)
end.idx <- c(start.idx[-1]-1,length(pt.vec))
idx <- cbind(start.idx, end.idx)
pt.rank <- rank(pt.vec, ties.method = "first")
bin.vec <- rep(NA, length(cells))
names(bin.vec) <- cells
for (i in seq(nrow(idx))){
  idx.start <- idx[i,1]
  idx.end <- idx[i,2]
  idx.cell <- which(pt.rank>=idx.start & pt.rank<=idx.end)
  bin.vec[idx.cell] <- i+72
}
hfa@meta.data[cells, "Stem_cell_maturation_bin"] <- bin.vec
saveRDS(hfa, file="Res_HIO_fetal_and_adult_stem_cell_combined.rds")

## get cell bin average expression levels from HIO to adult
expr.by.pt.bin <- getAveExpr(seu.obj = hfa, feature.to.calc = "Stem_cell_maturation_bin", colname.prefix = NULL)
saveRDS(expr.by.pt.bin, file="Res_HIO_fetal_and_adult_stem_cell_expr_by_bin.rds")

## create the dendrogram of cell bin from HIO to adult based on ISC development associated genes
## classify bins into 6 groups. i.e. classify the ISC state transition into 6 phases
a <- expr.by.pt.bin[sc.genes,]
hc.hfa <- hclust(as.dist(1-cor(a, method="spearman")), method="ward.D2") 
plot(hc.hfa, hang=-1)
hc.hfa$order <- seq(ncol(a))
saveRDS(hc.hfa, file="Res_sorted_HIO_fetal_ISC_maturity_dependent_genes_based_HIO_to_adult_cell_bin_hc_sorted.rds")
pdf("Plot_hc_HIO_to_adult_cell_bin.pdf")
plot(hc.hfa, hang=-1)
dev.off()

## assign phase index to each bin 
hfa <- subset(sc.combined, cells=selected.cells)
sc.phase.vec <- rep("Phase-6", ncol(hfa))
names(sc.phase.vec) <- colnames(hfa)
bin.vec <- hfa$Stem_cell_maturation_bin
sc.phase.vec[bin.vec%in%c(1:32)] <- "Phase-1"
sc.phase.vec[bin.vec%in%c(33:42)] <- "Phase-2"
sc.phase.vec[bin.vec%in%c(43:57)] <- "Phase-3"
sc.phase.vec[bin.vec%in%c(58:64)] <- "Phase-4"
sc.phase.vec[bin.vec%in%c(65:72)] <- "Phase-5"
hfa$Stem_cell_maturation_phase <- sc.stage.vec

## check expression pattern of curated signaling genes along stem cell state transition
signal.gene.list <- readRDS("/nas/groups/treutlein/USERS/Qianhui_Yu/Endoderm/include_new_SI/human_endoderm_map/remove_sex_MT_ribo_genes/LR_and_TF/with_Han_list/Res_Nat_comm_Han_curated_mouse_signaling_gene_list.rds")
signal.gene.ave.expr <- c()
rownames.vec <- c()
for(pathway in names(signal.gene.list)){
  for(gene.type in setdiff(names(signal.gene.list[[pathway]]), "Antagonists")){
    mouse.genes <- signal.gene.list[[pathway]][[gene.type]]
    human.genes <- unique(mouse2human[mouse2human[,"Mouse_symbol"]%in%mouse.genes, "Human_symbol"])
    human.genes <- intersect(human.genes, rownames(expr.by.pt.bin))
    if(length(human.genes)<1){
      next
    } 
    if(length(human.genes)>1){
      expr.mat <- expr.by.pt.bin[human.genes,]
      vec <- paste(pathway, gene.type, rownames(expr.mat), sep="_")
      
    }else if(length(human.genes)==1){
      expr.mat <- expr.by.pt.bin[human.genes,]
      vec <- paste(pathway, gene.type, human.genes, sep="_")
    }
    signal.gene.ave.expr <- rbind(signal.gene.ave.expr, expr.mat)
    rownames.vec <- c(rownames.vec, vec)
  }
}
rownames(signal.gene.ave.expr) <- rownames.vec
saveRDS(signal.gene.ave.expr, file="Res_development_signal_gene_expr_by_stem_cell_maturation_Pt_bin.rds")

## filter out non-expressed genes (sd==0)
sd.vec <- apply(signal.gene.ave.expr, 1, sd)
sum(sd.vec>0)
signal.gene.ave.expr <- signal.gene.ave.expr[names(sd.vec)[which(sd.vec>0)],]
dim(signal.gene.ave.expr)

## check the maximum expresion level of each gene, only keep those with max.expr>0.1
max.expr <- apply(signal.gene.ave.expr, 1, max)
summary(max.expr)
selected.signal.gene.ave.expr <- signal.gene.ave.expr[names(max.expr)[which(max.expr > 0.1)],]
dim(selected.signal.gene.ave.expr)
saveRDS(selected.signal.gene.ave.expr, file="Res_after_filtering_development_signal_gene_expr_by_stem_cell_maturation_pt_bin.rds")

########################################
# choose the genes to generate heatmap #
########################################
## plot heatmap for all pathway together
expr.mat <- selected.signal.gene.ave.expr
## plot stage-dependent genes (combine all adult stem cells as a group)
res <- readRDS("Res_HIO_fetal_and_adult_stem_cell_maturation_stage_markers.rds")
top.res <- res$top.res
selected.genes <- unique(top.res$feature)
length(selected.genes)
expr.mat <- expr.by.pt.bin[selected.genes,]
dim(expr.mat)
#############################################################################
## for both ###
hc.row <- hclust(as.dist(1-cor(t(expr.mat), method="pearson")), method="ward.D2")
row.orders <- hc.row$labels[hc.row$order]
expr.mat <- expr.mat[row.orders, ]
###
## for curated pathway genes ############
mat1 <- do.call('rbind',strsplit(rownames(expr.mat), "_"))
max.ct <- apply(expr.mat, 1, which.max)
idx.combined <- c()
for(pathway in unique(mat1[,1])){
  for(ct in seq(ncol(expr.mat))){
    idx <- which(max.ct==ct & mat1[,1]==pathway)
    idx.combined <- c(idx.combined, idx)
  }
  idx.combined <- c(idx.combined, NA)
}
# for others genes ########################
max.ct <- apply(expr.mat, 1, which.max)
idx.combined <- c()
for(ct in seq(ncol(expr.mat))){
  idx <- which(max.ct==ct)
  idx.combined <- c(idx.combined, idx)
}
###########################

## for both ###
expr.mat <- expr.mat[idx.combined,]
expr.mat <- rbind(expr.by.pt.bin[c("SOX11","SHH","LGR5","OLFM4"),], rep(NA, ncol(expr.by.pt.bin)), expr.mat)
rownames(expr.mat)[5] <- NA
input.norm <- t(apply(expr.mat, 1, function(vec){(vec-min(vec))/(max(vec)-min(vec))}))

## assign side bar colors for curated signaling genes
pathway.cols <- setNames(c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#ffffff", "#ffffff"), 
                         c("Bmp", "Fgf", "Hedgehog", "Notch", "RA", "Wnt", "LGR5", "OLFM4"))
gene.type.cols <- setNames(c("#52BE80", "#F4D03F", "#8E44AD", "#3498DB", "#ffffff", "#ffffff"),
                           c("Target genes", "Ligands", "Receptors",  "Antagonists", "LGR5", "OLFM4"))
saveRDS(pathway.cols, file="Res_signaling_pathway_cols.rds")
saveRDS(gene.type.cols, file="Res_gene_type_cols.rds")
mat1 <- do.call('rbind',strsplit(rownames(input.norm), "_"))
row.pathway.col <- pathway.cols[mat1[,1]]
row.gene.type.col <- gene.type.cols[mat1[,2]]
pdf("Plot_heatmap_stem_cell_maturation_stage_marker_gene_expr_across_fetal_and_adult_sc_pt_bin-pathway-3.pdf", height=14, width=10); par(cex=1)
heatmap.2(input.norm, trace="none",main="", density.info="none",
          dendrogram="none", Rowv=FALSE, Colv=FALSE, scale="none", keysize = 0.8, key=FALSE, 
          #          col=blue2Red,cexRow=0.8, cexCol = 0.8, margins=c(20,15))
          col=blue2Red,cexRow=0.7, cexCol = 1, RowSideColors = row.pathway.col, margins=c(20,15))
dev.off()
pdf("Plot_heatmap_stem_cell_maturation_stage_marker_gene_expr_across_fetal_and_adult_sc_pt_bin-gene_type-2.pdf", height=14, width=10); par(cex=1)
heatmap.2(input.norm, trace="none",main="", density.info="none",
          dendrogram="none", Rowv=FALSE, Colv=FALSE, scale="none", keysize = 0.8, key=FALSE, 
          #col=blue2Red,cexRow=1, cexCol = 1, margins=c(20,15))
          col=blue2Red,cexRow=0.7, cexCol = 1, RowSideColors = row.gene.type.col, margins=c(20,15))
dev.off()

## identify stem cell phase-dependent genes
de.res <- wilcoxauc(X=hfa, group_by = "Stem_cell_maturation_phase")
de.res$pct_diff <- de.res$pct_in - de.res$pct_out
sig.res <- de.res[which(de.res$logFC>0.25 & de.res$auc>0.6 & de.res$pct_in>20 & de.res$padj<0.01 & de.res$pct_diff>20),]
top.res <- sig.res %>% group_by(group) %>% top_n(50, wt=auc)
res <- list("sig.res"=sig.res, "top.res"=top.res)
saveRDS(res, file="Res_HIO_fetal_and_adult_stem_cell_maturation_phase_markers.rds")

## exclude the genes that are significantly enriched in enterocyte
### for fetal
cell.1 <- colnames(duo.epi)[which(duo.epi$Cell_type=="Stem_cell")]
cell.2 <- colnames(duo.epi)[which(duo.epi$Cell_type=="Enterocyte")]
cells <- c(cell.1, cell.2)
X <- duo.epi@assays$RNA@data[,cells]
y <- duo.epi@meta.data[cells, "Cell_type"]
de.res <- wilcoxauc(X=X, y = y)
de.res$pct_diff <- de.res$pct_in - de.res$pct_out
sig.res <- de.res[which(de.res$logFC>0.2 & de.res$auc>0.7 & de.res$pct_diff>20 & de.res$pct_in>25 & de.res$padj<0.01),]
g.e <- sig.res$feature[which(sig.res$group=="Enterocyte")]
g.sc <- sig.res$feature[which(sig.res$group=="Stem_cell")]
g.list <- list("Stem_cell"=g.sc, "Enterocyte"=g.e)
saveRDS(g.list, file="Res_stem_cell_vs_enterocyte_feature_genes_from_fetal_duo.rds")
fetal.e <- g.list$Enterocyte

### for adult
cell.1 <- colnames(adult.epi)[which(adult.epi$Cell_type=="Stem-cell")]
cell.2 <- colnames(adult.epi)[which(adult.epi$Cell_type=="Enterocyte")]
cells <- c(cell.1, cell.2)
X <- adult.epi@assays$RNA@data[,cells]
y <- adult.epi@meta.data[cells, "Cell_type"]
de.res <- wilcoxauc(X=X, y = y)
de.res$pct_diff <- de.res$pct_in - de.res$pct_out
sig.res <- de.res[which(de.res$logFC>0.2 & de.res$auc>0.7 & de.res$pct_diff>20 &  de.res$pct_in>25 & de.res$padj<0.01),]
g.e <- sig.res$feature[which(sig.res$group=="Enterocyte")]
g.sc <- sig.res$feature[which(sig.res$group=="Stem-cell")]
g.list <- list("Stem_cell"=g.sc, "Enterocyte"=g.e)
saveRDS(g.list, file="Res_stem_cell_vs_enterocyte_feature_genes_from_adult_duo.rds")
adult.e <- g.list$Enterocyte

## take the union set of enterocyte markers of fetal and adult duodenum as the gene set to be excluded from ISC phase-dependent genes
duo.e.markers <- union(adult.e, fetal.e)

sc.stage.marker <- readRDS("Res_HIO_fetal_and_adult_stem_cell_maturation_phase_markers.rds")
sig.res <- sc.stage.marker$sig.res
top.res <- sc.stage.marker$top.res
sig.res <- sig.res[which(!sig.res$feature%in%duo.e.markers),]
top.res <- top.res[which(!top.res$feature%in%duo.e.markers),]
res <- list("sig.res"=sig.res, "top.res"=top.res)
saveRDS(res, file="Res_HIO_fetal_and_adult_stem_cell_maturation_stage_markers.rds")
write.table(sig.res, file="Table_HIO_to_adult_ISC_state_dependent_genes.csv", sep=",", row.names = F, quote=F)

## plot HIO-to-adult cell bin information
cells <- rownames(hio.fetal.qp)
sc.prop <- hio.fetal.qp[cells,"Float_adjusted_SC_prop"]
age.vec <- hfa@meta.data[cells, "Age"]
bin.vec <- hfa@meta.data[cells, "Stem_cell_maturation_bin"]
info <- data.frame("Age"=age.vec, "Cell_bin"=bin.vec, "Stage"=stage.vec, "ISC_maturity"=sc.prop, stringsAsFactors = F)
age.order <- paste0("d", sort(as.numeric(sub("d","",unique(info$Age)))))
### make stacked bar plot to show age distribution in each cell bin
n1 <- sapply(age.order, function(age){
  sapply(sort(unique(info$Cell_bin)), function(bin){
    sum(info$Cell_bin==bin & info$Age==age)
  })
})
p1 <- t(n1/rowSums(n1))
colnames(p1) <- paste("bin", seq(nrow(n1)), sep="_")
sep <- rep(0, nrow(p1))
p2 <- rbind(cbind(p1[,1:32], sep, p1[,33:42], sep, p1[,43:57], sep, p1[,58:64], sep, p1[,65:72], matrix(0, ncol=(87-72+1), nrow=length(age.order))),
            rep(c(0,1), c(77,(87-72))))
idx <- which(colnames(p2)=="")
colnames(p2)[idx] <- c("sep", paste("bin", 72+seq(length(idx)-1), sep="_"))
rownames(p2)[nrow(p2)] <- "adult"
isc.age.cols <- setNames(c("#bdd7e7", "#b4d2e5", "#accee3","#92c1de", "#69acd5", "#2e7eba", 
                           "#fbb4b9", "#f994af", "#f775a5", "#f05d9d", "#e54d99", "#c61d8a", "#91097c", "#810378", "#7a0177",
                           "#6C3483"), 
                         c("d0", "d3", "d7", "d14", "d28", "d30", 
                           "d47", "d59", "d72", "d80", "d85", "d101", "d122", "d127", "d132", 
                           "adult"))
pdf("Plot_stacked_barplot_HIO_and_sample_age_distribution_HIO_to_adult.pdf")
barplot(p2, col=isc.age.cols[rownames(p2)], border=NA)
dev.off()

### side bar to indicate phase index of each cell bin
phase.cols <- setNames(c("#accee3","#76D7C4","#1ABC9C","#f05d9d","#91097c","#6C3483", "#ffffff"),
                       c(paste("Phase", seq(6), sep="-"), "sep"))
mat <- unique(hfa@meta.data[,c("Stem_cell_maturation_phase", "Stem_cell_maturation_bin")])
link.vec <- setNames(c(mat[,1],"sep"), c(mat[,2], "sep"))
mat[,2] <- paste("bin", mat[,2], sep="_")
vec <- link.vec[colnames(p2)]
vec.cols <- stage.cols[vec]

pdf("Plot_HIO_to_adult_ISC_cell_bin_colored_by_stage.pdf")
barplot(rep(1,ncol(p2)), col=vec.cols, border = NA)
dev.off()


