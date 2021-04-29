fetal.pt.genes <- readLines("Ref_data_for_ISC_score/List_fetal_pseudotime_dependent_genes_for_quadratic_programming.txt")
X <- readRDS("Ref_data_for_ISC_score/Dat_HIO_Cl12_and_d132_stem_cell_fetal_pt_gene_average_expr_as_ref_for_qp.rds")

# load query seurat object and specify cells for analysis
que.obj
que.cells
expressed.fetal.pt.genes <- intersect(fetal.pt.genes, que.obj)
## prepare expression matrix of query cells
que.expr <- as.matrix(que.objd@assays$RNA@data[expressed.fetal.pt.genes, que.cells])
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

isc.score <- identity.que[,2]
