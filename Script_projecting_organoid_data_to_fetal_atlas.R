library(Seurat)
source("Script_functions.R")

# load fetal data
## load processed fetal atlas seurat object, which has the CSS-based UMAP embedding of the fetal atlas
fetal <- readRDS("Res_all_fetal_cells_combined_no_MT_Ribo_Sex_genes_CSS.rds")
## load fetal Cluster Similarity Spectrum (CSS) model
css.model <- readRDS("~/Work/Endoderm/include_new_SI/human_endoderm_map/remove_sex_MT_ribo_genes/Res_fetal_CSS_model.rds")
## load fetal CSS-based UMAP model
fetal.umap.res <- readRDS("Res_fetal_all_cell_type_CSS_based_UMAP_res.uwot")


# tHIO
t.hio <- readRDS("~/Work/Endoderm/used_seurat_objects/Res_H9-tHIO.rds")
## use the fetal highly variable genes that are also detected in tHIO data to calculate similarity between tHIO cells and fetal clusters per sample
## because here the fetal cluster average profile is reference, and the tHIO cells are query, the obtained similarity spectrum is called reference similarity spectrum (RSS)
fetal.hvg <- VariableFeatures(fetal)
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
css <- fetal@reductions$css@cell.embeddings
knn <- RANN::nn2(css, rss.mat, k = 20)$nn.idx
## map the organ identity of fetal cells to tHIO cells
ref.idx <- fetal@meta.data$Corrected_organ_group
nn.idx <- matrix(ref.idx[as.vector(knn)], nrow=nrow(knn))
saveRDS(nn.idx, file="Res_20_fetal_NN_organ_id_for_each_tHIO_cell_after_organ_identity_correction.rds")
trans.id <- apply(nn.idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})
t.hio@meta.data$Mapped_fetal_organ_after_correction <- trans.id
saveRDS(t.hio, file="Res_tHIO_with_CSS_and_fetal_projection.rds")

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



## CDX2 control and KO
cdx2 <- readRDS("~/Work/Endoderm/used_seurat_objects/Res_cdx2_ctrl_KO.rds")
DimPlot(cdx2, reduction = "umap_css")
rss.mat <- cdx2@reductions$rss@cell.embeddings
knn <- RANN::nn2(css, rss.mat, k = 20)$nn.idx
# map the organ identity of fetal cells to CDX2 control and KO in vitro HIO cells
ref.idx <- fetal@meta.data$Corrected_organ_group
nn.idx <- matrix(ref.idx[as.vector(knn)], nrow=nrow(knn))
saveRDS(nn.idx, file="Res_20_fetal_NN_organ_id_for_each_CDX2-control_KO_HIO_cell_after_organ_identity_correction.rds")

trans.id <- apply(nn.idx, 1, function(vec){
  freq <- table(vec)
  names(which.max(freq))
})

cdx2@meta.data$Mapped_fetal_endoderm_organ_group_after_correction <- trans.id
ct.vec <- cdx2@meta.data$Major_cell_type
condition.vec <- cdx2@meta.data$Tissue
organ.vec <- cdx2@meta.data$Mapped_fetal_endoderm_organ_group_after_correction
p.list <- lapply(c("Epithelial", "Mesenchymal"), function(ct){
  mat <- sapply(c("CDX2-WT-HIO", "CDX2-KO-HIO"), function(condition){
    sapply(c("Intestine", "Stomach", "Lung", "Esophagus", "Pancreas"), function(organ){
      sum(ct.vec==ct & condition.vec==condition & organ.vec==organ)
    })
  })
  p.mat <- t(t(mat)/apply(mat, 2, sum))
  return(p.mat)
})
names(p.list) <- c("Epithelial", "Mesenchymal")
p.mat.combined <- do.call('cbind', p.list)

# plot mapped organ proportion by major cell type by condition
pdf("Plot_barplot_cdx2_HIO_fetal_organ_projection_after_correction_for_each_major_cell_type_and_condition.pdf")
barplot(p.mat.combined, col=organ.cols[rownames(p.mat.combined)], border = NA)
dev.off()

subtype.vec <- cdx2@meta.data$Cell_type
DimPlot(cdx2, reduction = "umap_css", group.by = "Cell_type", label=T)
selected.subtype <- c("Proliferative_epi", "MUC5AC+_goblet", "FABP1+_enterocyte", "MMP7+_epi", "SPP1+_epi", "MUC16+_epi")
num.mat <- sapply(selected.subtype, function(ct){
  sapply(c("Intestine", "Stomach", "Lung", "Esophagus", "Pancreas"), function(organ){
    sum(subtype.vec==ct & organ.vec==organ)
  })
})
p.mat <- t(t(num.mat)/apply(num.mat, 2, sum))
# plot mapped organ proportion by each epithelial subtype
pdf("Plot_barplot_cdx2_HIO_fetal_organ_projection_after_correction_for_each_epi_subtype.pdf")
barplot(p.mat, col=organ.cols[rownames(p.mat)], border = NA, las=2)
dev.off()
saveRDS(cdx2, file="Res_CDX2_ctrl_KO.rds")

# plot corrected organ identity distribution
png("Plot_UMAP_CSS_fetal_atlas_organ_identity.png", height=2000, width=2000)
plotFeature2(Embeddings(fetal, reduction = "umap_css"), values = fetal$Corrected_organ_group, gCols = organ.cols, cex=2)
dev.off()


