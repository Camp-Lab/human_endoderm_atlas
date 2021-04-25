library(Seurat)
library(uwot)
source("Script_functions.R")

# load fetal data
## Please get the fetal reference data from our Mendeley data repository at http://dx.doi.org/10.17632/x53tts3zfr
## load fetal atlas highly variable genes
fetal.hvg <- readLines("Res_fetal_atlas_highly_variable_genes.csv")
## load fetal atlas meta.data after decompression, which contains CSS-based UMAP embedding of the fetal atlas and organ identity
fetal.meta <- read.csv("Table_fetal_atlas_meta_data.csv")
fetal.embeddings <- fetal.meta[,c("UMAP_X", "UMAP_Y")]
ref.idx <- fetal.meta$Corrected_organ_group
## load fetal Cluster Similarity Spectrum (CSS) model
css.model <- readRDS("Res_fetal_CSS_model.rds")
## load fetal CSS-based UMAP model
fetal.umap.res <- load_uwot("Res_fetal_all_cell_type_CSS_based_UMAP_res.uwot")


# use tHIO as an example
t.hio <- readRDS("Res_H9-tHIO.rds")
## use the fetal highly variable genes that are also detected in tHIO data to calculate similarity between tHIO cells and fetal clusters per sample
## because here the fetal cluster average profile is reference, and the tHIO cells are query, the obtained similarity spectrum is called reference similarity spectrum (RSS)
shared.genes <- intersect(fetal.hvg, rownames(t.hio))
que.data <- as.matrix(t.hio@assays$RNA@data[shared.genes,])
rss.list <- lapply(seq(length(css.model$model$profiles)), function(i){
  ref <- css.model$model$profiles[[i]][shared.genes,]
  cor.mat <- cor(que.data, ref, method = "spearman")
  cor.z <- t(scale(t(cor.mat)))
  return(cor.z)
})
rss.mat <- do.call('cbind', rss.list)
t.hio[["rss"]] <- CreateDimReducObject(embeddings = rss.mat, key="RSS_", assay=DefaultAssay(t.hio))

## calculate the Euclidean distance between tHIO cells and fetal cells on the space of similarities to fetal clusters per sample, instead of on the space of expression. correlation is potentially more robust to batch variations than expression 
## after the distance calculation, get the 20-nearest fetal cells for each tHIO cell
## assign the majority of organ identity of the 20-nearest cells to be the inferred organ identity of tHIO cells 
css <- css.model$sim2profiles 
knn <- RANN::nn2(css, rss.mat, k = 20)$nn.idx
## map the organ identity of fetal cells to tHIO cells
nn.idx <- matrix(ref.idx[as.vector(knn)], nrow=nrow(knn))
trans.id <- apply(nn.idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})
t.hio@meta.data$Mapped_fetal_organ_after_correction <- trans.id

## project tHIO cells to the CSS-based UMAP embedding of fetal atlas data and visualize the projection result
thio.umap <- umap_transform(rss.mat, fetal.umap.res)
rownames(thio.umap) <- rownames(rss.mat)
t.hio[["rss_umap"]] <- CreateDimReducObject(embeddings = thio.umap, key="RSSUMAP_", assay=DefaultAssay(t.hio))
saveRDS(t.hio, file="Res_tHIO_with_CSS_and_fetal_projection.rds")

mapped.organ <- t.hio$Mapped_fetal_organ_after_correction
col.vec <- organ.cols[mapped.organ]

png("Plot_UMAP_project_tHIO_to_fetal_atlas.png", height=2000, width=2000)
plot(Embeddings(fetal, reduction = "umap_css"), pch=16, col="#e0e0e0", bty="n", xaxt="n", yaxt="n", xlab="", ylab="", cex=3)
points(Embeddings(t.hio, reduction = "rss_umap"), pch=16, col=col.vec, cex=3)
dev.off()

## plot mapped organ proportion by major cell type with stacked bar plot
ct.vec <- t.hio@meta.data$Major_cell_type
organ.vec <- t.hio@meta.data$Mapped_fetal_organ_after_correction
num.mat <- sapply(c("Epithelial", "Mesenchymal"), function(ct){
  sapply(c("Intestine", "Stomach", "Lung", "Esophagus", "Pancreas"), function(organ){
    sum(ct.vec==ct & organ.vec==organ)
  })
})
p.mat <- t(t(num.mat)/apply(num.mat, 2, sum))
pdf("Plot_barplot_tHIO_fetal_organ_projection_after_correction_for_each_major_cell_type.pdf")
par(mfrow=c(1,2))
barplot(p.mat[,"Epithelial"], col=organ.cols[rownames(p.mat)], border = NA, beside = T, main="Epithelial", las=2)
barplot(p.mat[,"Mesenchymal"], col=organ.cols[rownames(p.mat)], border = NA, beside = T, main="Mesenchymal", las=2)
dev.off()

## plot mapped organ proportion by each epithelial subtype with stacked bar plot
subtype.vec <- t.hio@meta.data$Cell_type
selected.subtype <- c("Stem cell", "Enterocyte progenitor", "Enterocyte", "Enteroendocrine", "MUC2+/MUC5AC-/MUC5B-_goblet",
                      "SPIB+/CA7+_M-cell", "MUC2-/MUC5AC+/MUC5B+_goblet", "Basal-cell-like-1", "Basal-cell-like-2",
                      "Ciliated-like", "CXCL14+", "LINC01082+", "Pro.M.", "SFRP2+", "VSM")
num.mat <- sapply(selected.subtype, function(ct){
  sapply(c("Intestine", "Stomach", "Lung", "Esophagus", "Pancreas"), function(organ){
    sum(subtype.vec==ct & organ.vec==organ)
  })
})
p.mat <- t(t(num.mat)/apply(num.mat, 2, sum))
pdf("Plot_barplot_tHIO_fetal_organ_projection_after_correction_for_each_subtype.pdf")
barplot(p.mat, col=organ.cols[rownames(p.mat)], border = NA, las=2)
dev.off()

