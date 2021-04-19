combined <- readRDS("Res_HIO_tHIO_fetal_adult_duo_stem_cell_LOG_RNA_data.rds")
hio.fetal.qp <- readRDS("Res_organoid_fetal_QP_result.rds")

# visualize the stem cell maturity analysis results for organoid and fetal stem cells 
## assign blue colors to organoid cells, pink colors to fetal. Light colors indicate younger ages, while dark colors indicate older ages 
idx <- order(hio.fetal.qp[,"frxn_fetLatSC"])
sorted.sc <- hio.fetal.qp[idx,"frxn_fetLatSC"]
sorted.cl12 <- hio.fetal.qp[idx,"frxn_cl12"]
sorted.age <- as.numeric(sub("d", "", combined@meta.data[rownames(hio.fetal.qp)[idx], "Age"]))
sorted.sys <- combined@meta.data[rownames(hio.fetal.qp)[idx], "Stage"]
colorPal.hio <- grDevices::colorRampPalette(c("#bdd7e7", "#6baed6", "#3182bd", "#08519c"))
colorPal.fetal <- grDevices::colorRampPalette(c("#fbb4b9", "#f768a1", "#c51b8a", "#7a0177"))
cellColor <- rep(NA, ncol(que.expr))
idx.hio <- which(sorted.sys%in%c("thio", "hio"))
cellColor[idx.hio] <- adjustcolor(colorPal.hio(30), alpha=.8)[as.numeric(cut(sorted.age[idx.hio], breaks=30, right=F,include.lowest=T))]
idx.fetal <- which(sorted.sys=="fetal")
cellColor[idx.fetal] <- adjustcolor(colorPal.fetal(30), alpha=.8)[as.numeric(cut(sorted.age[idx.fetal], breaks=30, right=F,include.lowest=T))]

## show the increase of stem cell maturity among all examined cells
pdf("Plot_quadratic_programming_stem_cell_fraction_fetal_as_sc_ref_frame_only.pdf", width=7, height=5)
plot(c(1, length(sorted.sc)), c(0,1), type="n", pch=16,  main="", ylab="Stem cell component proportion", xlab="Sorted cell index")
dev.off()
png("Plot_quadratic_programming_stem_cell_fraction_fetal_as_sc_ref.png", width=15, height=10, unit="cm", res=500)
plot(seq(length(sorted.sc)), sorted.sc, col=cellColor, pch=16, main="", ylab="", xlab="", xaxt="n", yaxt="n",  ylim=c(0,1), cex=1.5, bty="n")
points(seq(length(sorted.sc)), sorted.cl12, col=cellColor, pch=16, cex=1.5)
dev.off()

## show the distribution of stem cell maturity score in each organoid and fetal sample  
sorted.age <- combined@meta.data[rownames(hio.fetal.qp)[idx], "Age"]
age.orders <- c("d0", "d3", "d7", "d14", "d47", "d28", "d30", "d59", "d56", "d72", "d80", "d84", "d85","d101","d127","d122", "d132")
age.cols <- sapply(age.orders,function(x){
  unique(cellColor[which(sorted.age==x)])
})
sc.prop.list <- lapply(age.orders, function(x){
  sorted.sc[which(sorted.age==x)]
})
names(sc.prop.list) <- age.orders
pdf("Plot_boxplot_stem_cell_fraction_by_age.pdf", height=5, width=8)
boxplot(sc.prop.list, las=2, col = age.cols, border=substr(age.cols, 1, 7), lwd=3)
dev.off()

# identify genes correlated with stem cell identity 
ref.vec <- hio.fetal.qp[, "frxn_fetLatSC"]
que.expr <- as.matrix(combined@assays$RNA@data[,rownames(hio.fetal.qp)])
cor.mat <- cor(ref.vec, t(que.expr), method="spearman")
saveRDS(cor.mat, file="Res_SCC_to_stem_cell_identity_prop.rds")

# use stem cell identity fraction as a pseudotime representing stem cell maturation, and identify genes showing variable expression along the pseudodtime
cells <- rownames(hio.fetal.qp)[which(hio.fetal.qp[,"frxn_fetLatSC"]>0 & hio.fetal.qp[,"frxn_fetLatSC"]<1)]
sc.prop <- hio.fetal.qp[cells,"frxn_fetLatSC"]
expr.mat <- as.matrix(combined@assays$RNA@data[,cells])

expr.by.pt.bin <- getExprByPt(pt.vec=sc.prop, expr.mat=expr.mat, cell.num.per.bin=50)
saveRDS(expr.by.pt.bin, file="Res_expr_by_stem_cell_identity_fraction_bin.rds")
age.test.p <- splineBasedAgeTest(pseudotime.vec=seq(ncol(expr.by.pt.bin)), expr.mat=expr.by.pt.bin, df=5, mode="stringent")
padj <- p.adjust(age.test.p, method="BH")
logfc <- apply(expr.by.pt.bin, 1, function(vec){
  max(vec)-min(vec)
})
age.test.res <- data.frame("Residual_comparison_P"=age.test.p, "BH_corrected_P"=padj, "Log(max-min)"=logfc, stringsAsFactors = F)
saveRDS(age.test.res, file="Res_age_test_along_stem_cell_identity_fraction_pseudotime.rds")  

sc.genes <- rownames(age.test.res)[which(age.test.res$BH_corrected_P<0.01 & age.test.res$Log.max.min.>0.5)]
length(sc.genes)
sc.expr <- expr.by.pt.bin[sc.genes,]
hc.sc <- hclust(as.dist(1-cor(t(sc.expr), method="spearman")), method="ward.D2")
plot(hc.sc, hang=-1, cex=0.1)
abline(h=5)
g=cutree(hc.sc, h=5)
pt.gene.modules <- plotClusterExprProfile(expr=sc.expr, time.vec=seq(ncol(sc.expr)), group.vec=rep("fetal", ncol(sc.expr)), cluster.vec=g, group.cols="#31a354", return.value=T, to.plot=T, plot.name="Plot_stem_cell_progession_related_gene_module_expr_profile_in_organoid_and_fetal.pdf", add.legend=F, legend.pos="topleft", cex.legend=2)
saveRDS(pt.gene.modules, file="Res_stem_cell_progession_related_gene_module_expr_profile_in_organoid_and_fetal.rds")

# Among genes belonging to cluster 2, which turns on in HIO stage and keep on in fetal, choose those showing strongest correlation with stem cell identity proportion
sc.cor.genes <- sapply(2, function(cl.x){
  g.x <- names(g)[which(g==cl.x)]
  cor.vec <- cor.mat[1,g.x]
  g.cor.idx <- order(cor.vec, decreasing=T)[1:100]
  g.cor <- names(cor.vec)[g.cor.idx]
  return(g.cor)
})
cl2.genes <- as.vector(sc.cor.genes)

# choose the top stem cell fraction correlated genes
cells <- rownames(hio.fetal.qp)[which(hio.fetal.qp[,"frxn_fetLatSC"]>0 & hio.fetal.qp[,"frxn_fetLatSC"]<1)]
sc.prop <- hio.fetal.qp[cells,"frxn_fetLatSC"]
expr.mat <- as.matrix(combined@assays$RNA@data[,cells])
cor.mat <- cor(sc.prop, t(expr.mat), method="spearman")

g0 <- colnames(cor.mat)[which(!is.na(cor.mat[1,]))]
vec0 <- cor.mat[1,g0]
idx.up <- order(vec0, decreasing = T)[1:100]
idx.down <- order(vec0)[1:100]
genes.up <- setdiff(g0[idx.up], cl2.genes)
genes.down <- setdiff(g0[idx.down], cl2.genes)
selected.genes <- c(c("LGR5","OLFM4"), genes.down, cl2.genes, genes.up)
