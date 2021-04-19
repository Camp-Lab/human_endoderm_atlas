library(Seurat)
library(presto)
library(irlba)
library(dplyr)
library(gplots)
library(destiny)
library(uwot)
source("~/Work/commonScript/Script_functions.R")
library(simspec)

# set working folder
setwd("~/Work/Endoderm/include_new_SI/human_endoderm_map/remove_sex_MT_ribo_genes/")

# set colors
expr.cols <- c("#d0d1e6","#a6bddb", "#74a9cf", "#2b8cbe", "#045a8d")
age.cols <- c("#bdbdbd", "#252525")
tissue.group.cols <- setNames(c("#fcb163","#4da7b0","#5d4ea2","#9d0142","#2471A3","#a9489a","#1b7635", "#d9d9d9", "#9d0142"), c("Lung", "Liver", "Stomach","Duodenum_SI", "Other_SI_region", "Colon", "Esophagus", "UD", "SI"))
organ.cols <- setNames(c("#fcb163","#4da7b0","#5d4ea2","#9d0142","#9d0142","#9d0142","#1b7635"), c("Lung", "Liver", "Stomach","SI", "Intestine","Colon", "Esophagus"))
col.vec <- c("#F9E79F","#F4D03F", "#F5B041", "#DC7633","#1b7635","#BB8FCE","#A569BD","#6C3483","#EC7063","#EC7063","#9d0142","#2471A3","#3498DB","#02818a", "#17A589")
tissue.vec <- c("Lung-tracheal-epi","Lung-airway-trachea", "Lung-airway", "Lung-distal", "Esophagus", "Stomach","Stomach-corpus", "Stomach-antrum", "Small-intestine", "Proximal-small-intestine", "Duodenum", "Jejunum", "Ileum", "Colon", "Liver")
tissue.cols <- setNames(col.vec,tissue.vec)
cell.group.cols <- setNames(c("#8E44AD", "#DD3497", "#08519C", "#A50F15", "#006D2C","#FACB12", "#F39C12"),c("Endothelial", "Neuronal", "Immune", "Epithelial", "Mesenchymal", "Hepatocyte", "Erythroid"))

# read 10x count matrix of each sample
dir <- "/home/yuq/Work/processed_data/endoderm_atlas_fetal_data_link/"
all.files <- list.files(dir, pattern="Sample_")
seu.obj.list <- list()
for(i in seq(length(all.files))){
  cat(paste("Sample", i, "start\n"))
  sample.name <- all.files[i]
  file.path <- paste0(dir, sample.name)
  vec <- strsplit(sample.name, split="_")[[1]]
  age <- vec[3]
  individual <- vec[2]
  tissue <- vec[4]
  info.mat <- cbind(c("Age", "Tissue", "Individual", "Sample"), c(age, tissue, individual, sample.name))
  seu.obj.list[[i]] <- prepareSeuratObject2(rawData = Read10X(data.dir = file.path), namePrefix = paste0("S",i), additional.sample.info.mat=info.mat, seu.version=3, mito.cutoff=Inf, min.cell.cutoff=1, min.gene.cutoff=500, gene.num.low=500, gene.num.high=Inf)
}

# combine data of all samples and filter cells based on percent.mito and nFeature
combined <- merge(x=seu.obj.list[[1]], y=seu.obj.list[-1])
combined <- subset(combined, subset = percent.mito < 0.1 & nFeature_RNA > 1000)

# remove confounding genes and make a new combined Seurat object
dir <- "~/Work/Annotation/confound_genes/" 
all.files <- list.files(dir, pattern=".txt")
confound.genes <- lapply(seq(length(all.files)), function(i){
  file <- paste0(dir, all.files[i])
  g <- readLines(file)
  return(g)
})
names(confound.genes) <- c("Cell_cycle", "Experiment_induced", "HB", "MT", "Ribosome", "Sex_chr") 
all.confound.genes <- unique(unlist(confound.genes))
genes.to.remove <- unique(c(confound.genes[["MT"]], confound.genes[["Ribosome"]], confound.genes[["Sex_chr"]]))
counts <- combined@assays$RNA@counts
counts <- counts[setdiff(rownames(counts), genes.to.remove),]
meta <- combined@meta.data[,4:8]
fetal <- CreateSeuratObject(counts=counts, meta=meta)

# general analysis on combined object without batch effect correction
tissue.vec <- fetal@meta.data$Tissue
organ.vec <- tissue.vec
organ.vec[tissue.vec%in%c("Duodenum", "Ileum", "Jejunum", "Small-intestine", "Proximal-small-intestine")] <- "SI"
organ.vec[grep("Lung", tissue.vec)] <- "Lung"
organ.vec[grep("Stomach", tissue.vec)] <- "Stomach"
fetal@meta.data[["Organ"]] <- organ.vec
age.vec <- fetal@meta.data$Age
age.week.vec <- round(as.numeric(sub("d", "", fetal@meta.data$Age))/7)
fetal@meta.data$Age_week <- age.week.vec
fetal <- NormalizeData(object = fetal, normalization.method = "LogNormalize", scale.factor = 1e4)
fetal <- FindVariableFeatures(object = fetal, selection.method = "vst", nfeatures = 3000)
fetal <- ScaleData(object = fetal, verbose = T)
fetal <- RunPCA(fetal)
usefulPCs <- 1:20
fetal <- RunUMAP(fetal, dims = usefulPCs)
saveRDS(fetal, file=paste0("Res_all_fetal_cells_combined_no_MT_Ribo_Sex_genes_",max(usefulPCs),"PC.rds"))

# CSS (Cluster Similarity Spectrum) integration
## run CSS integration to build a CSS model
## so that other datasets could be projected to fetal data embedding or compared to fetal data on CSS space
## CSS of all fetal cells are also returned, which could be used as input for UMAP and clustering
###################################### if one wants to return the css.model, run this chunk of code #########################################
css.model <- cluster_sim_spectrum(fetal, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, 
                                  return_seuratObj = FALSE)
saveRDS(css.model, file="Res_fetal_CSS_model.rds")
css <- css.model$sim2profiles
colnames(css) <- paste("CSS", colnames(css), sep="_")
fetal[["css"]] <- CreateDimReducObject(embeddings = css, key="CSS_", assay=DefaultAssay(fetal))
##############################################################################################################################################

########################### if one only wants to get CSS, without returning the css.model, run this chunk of code ############################
#fetal <- cluster_sim_spectrum(fetal, label_tag = "orig.ident", cluster_resolution = 0.6, merge_spectrums = F, verbose = T, 
#                              return_seuratObj = TRUE)
##############################################################################################################################################

# use CSS as input for UMAP and clustering
fetal <- RunUMAP(fetal, reduction = "css", dims = 1:ncol(fetal@reductions$css@cell.embeddings), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
fetal <- FindNeighbors(object = fetal, reduction = "css", dims = 1:ncol(fetal@reductions$css@cell.embeddings), force.recalc = T, k.param = 15)
fetal <- FindClusters(object = fetal, resolution = 0.6)
saveRDS(fetal, file="Res_all_fetal_cells_combined_no_MT_Ribo_Sex_genes_CSS.rds")

# train a UMAP model using CSS of the fetal atlas data as inputs
# when Reference Similarity Spectrum (RSS) of tHIO cells to fetal atlas data is calculated,
# one could use this pre-trained UMAP model to project tHIO cells to the CSS-based UMAP embedding of the fetal atlas data
fetal.umap.res <- get_umap_models(fetal, reduction = "css", dims=1:ncol(fetal@reductions$css@cell.embeddings))
fetal.umap.res$embedding <- Embeddings(fetal, reduction = "umap_css")
save_uwot(fetal.umap.res, file="Res_fetal_all_cell_type_CSS_based_UMAP_res.uwot")


# plot general information
sc.coor=Embeddings(fetal, reduction="umap_css")
sc.coor<-fetal.umap.res$embedding
sc.cl= paste0("Cluster", fetal@meta.data$RNA_snn_res.0.6)
sc.organ=fetal@meta.data$Organ
sc.tissue=fetal@meta.data$Tissue
sc.age=as.numeric(fetal@meta.data$Age_week)
row.num=2
col.num=2
tissue.order <- c("Esophagus", "Lung-tracheal-epi","Lung-airway-trachea", "Lung-airway", "Lung-distal",  "Liver","Stomach","Stomach-corpus", "Stomach-antrum", "Small-intestine", "Proximal-small-intestine", "Duodenum", "Jejunum", "Ileum", "Colon")
organ.order <- c("Esophagus", "Lung", "Liver","Stomach", "SI", "Colon")
png("Plot_UMAP_CSS_all_fetal_RNA_LOG_basic_info-2.png", height=2000*row.num, width=2000*col.num)
par(mfrow=c(row.num,col.num))
plotFeature2(sc.coor, values=sc.organ, main="Organ", cex.main = 7, cex=2.5, gCols=organ.cols, add.legend=T, legend.cex=4, legend.pos="bottomleft", group.order=organ.order)
plotFeature2(sc.coor, values=sc.tissue, main="Tissue", cex.main = 7, cex=2.5, gCols=tissue.cols, add.legend=T, legend.cex=4, legend.pos="bottomleft", group.order=tissue.order)
plotFeature2(sc.coor, values=sc.age, main="Age",cex.main = 7,  cex=2.5, nCols=age.cols)
plotFeature2(sc.coor, values=sc.cl, main="Cluster", cex.main = 7, cex=2.5, add.label=T, label.cex = 8)
dev.off()

# plot marker gene expression in batch
plotFeature.batch(seu.obj=fetal, dr="umap_css", genes.to.plot=c("EPCAM","COL1A2","CDH5","ALB","HEMGN","PTPRC","ASCL1"), col.num = 4, plot.name = "Plot_UMAP_CSS_major_cell_type_marker_expr.png")

# classify clusters into cell classes
known.markers <- read.table("Table_major_cell_type_markers_from_literature_search.txt", head=T, sep="\t",stringsAsFactors=F)  
cell.type.order <- c("Neuronal", "Epithelial", "Mesenchymal", "Endothelial", "Immune", "Erythroid", "Hepatocyte")
pan.cm <- unique(unlist(lapply(cell.type.order, function(ct){
  intersect(known.markers$Gene[which(known.markers$Used_pan_cell_type_markers==ct)], rownames(fetal))
}))) 

fetal.expr <- getAveExpr(seu.obj=fetal, feature.to.calc = "Major_cell_type", specified.order = cell.type.order, genes=pan.cm, colname.prefix = NULL)
fetal.prop <- getExpressedProp(seu.obj = fetal, feature.to.calc = "Major_cell_type", specified.order = cell.type.order, genes=pan.cm, colname.prefix = NULL)
getDotPlot(cl.expr.mat=fetal.expr, cl.prop.mat=fetal.prop, gene.reorder=FALSE, specified.order=NULL, cl.reorder.by.hc=FALSE, genes=pan.cm, colors=c("#d9d9d9", "#252525"), point.size.factor=4.5, plot.name="Plot_dot_plot_cell_group_marker_RNA_LOG_expr.pdf", plot.height=20, plot.width=30, plot.margin=c(8,12,5,8), plot.cex=2, max.diag=FALSE, col.space.factor=0.45, row.space.factor=0.12)

meso.cl <- c(5,15,2,0,1,3,6,16,26)
epi.cl <- c(24,18,10,11,17,4,23,21)
neuron.cl <- c(13,14)
endo.cl <- 9
erythroid.cl <- 12
immune.cl <- c(19,20,7,22,8)
hepatocyte.cl <- 25
cl.mat <- cbind(c(rep("Mesenchymal", length(meso.cl)), rep("Epithelial", length(epi.cl)), rep("Neuronal", length(neuron.cl)), 
                  rep("Endothelial", length(endo.cl)), rep("Erythroid", length(erythroid.cl)),
                  rep("Immune", length(immune.cl)), rep("Hepatocyte", length(hepatocyte.cl))),
                c(meso.cl, epi.cl, neuron.cl, endo.cl, erythroid.cl, immune.cl, hepatocyte.cl))
cl.vec <- fetal@meta.data$RNA_snn_res.0.6
ct.vec <- rep(NA, length(cl.vec))
for(x in unique(cl.mat[,1])){
  cl.x <- cl.mat[cl.mat[,1]==x,2]
  ct.vec[which(cl.vec%in%cl.x)] <- x
}
fetal@meta.data$Major_cell_type <- ct.vec
DimPlot(fetal, reduction = "umap_css", group.by = "Major_cell_type")
saveRDS(fetal, file="Res_all_fetal_cells_combined_no_MT_Ribo_Sex_genes_CSS.rds")