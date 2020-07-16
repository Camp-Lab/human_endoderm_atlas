organ.cols <- setNames(c("#9D0142", "#5D4EA2", "#FCB163", "#1B7635"), c("Intestine", "Stomach", "Lung", "Esophagus"))
blue.expr.cols <- c("#FFF7FB", "#ECE7F2", "#D0D1E6", "#A6BDDB", "#74A9CF", "#3690C0", "#023858")
organ.cols <- setNames(c("#9D0142", "#5D4EA2", "#FCB163", "#1B7635", "#4da7b0"), c("Intestine","Stomach","Lung", "Esophagus","Liver"))


#' Plot feature's gradients/classes across all samples (cells) given the plotting coordinates
#'
#' This function plot the samples based on the given coordinates, coloring
#' each dot based on the value of its corresponding feature value/class.
#'
#' @param coord 	The plotting coordinates, expected to be of two columns. Each row represents one dot for one sample (cell).
#' @param value 	The values to be shown, with each value representing one sample (cell) with the same order as coord.
#' @param emphasize 	The Indices of samples (cells) to be emphasized. When it is set, the colors of all the other cells are set to be #bdbdbd30.
#' @param col 	The customized colors of dots. Either col or value must be given.
#' @param ... 	Other arguments passing to the plot function.
#' @export
plotFeature2 <- function(coor, values, label.only=FALSE, point.order= NULL, knn.pairs=NULL, emphasize=NULL, nCols=NULL, gCols=NULL, group.order=NULL, zeroAsGray=T, add.legend=F, add.label=F, label.cex=1, legend.cex=1, legend.pos="topright", add.cluster=FALSE, cl.vec=NULL, ...){
  if(label.only){
    plot(c(min(coor[,1]), max(coor[,1])), c(min(coor[,2]),max(coor[,2])),type="n", cex.main=5, xaxt="n", yaxt="n", xlab="", ylab="", bty="n", ...)
    for(x in unique(values)){
      cell.idx=which(values==x)
      if(length(cell.idx)>1){
        label.coor=apply(coor[cell.idx,], 2, median)
      }else{
        label.coor=coor[cell.idx,]
      }
      if(grepl("Cluster_", x, fixed=T)){
        text(label.coor[1], label.coor[2], label=sub("Cluster_", "", x), cex=label.cex)
      }else{
        text(label.coor[1], label.coor[2], label=sub("Cluster", "", x), cex=label.cex)
      }
    }
    return(0)
  }
  
  
  if(is.numeric(values)){
    if(is.null(nCols)){
      colorPal <- grDevices::colorRampPalette(c("darkgreen", "yellow","red"))
    }else{
      colorPal <- grDevices::colorRampPalette(nCols)
    }
    cellColor <- adjustcolor(colorPal(30), alpha=.8)[as.numeric(cut(values, breaks=30, right=F,include.lowest=T))]
    if(zeroAsGray){
      if(min(values, na.rm=T) == 0) cellColor[values == 0] <- "#bdbdbd"
    }
  }else{
    if(is.null(gCols)){
        gCols <- setNames(scales::hue_pal()(length(unique(values))), unique(values))
    }
    cellColor <- gCols[values]
  }
  
  if(is.null(knn.pairs)){
    if(is.null(emphasize)){
      if(is.null(point.order)){
        point.coor=coor
      }else if(point.order=="sorted"){
        idx <- order(values)
        cellColor <- cellColor[idx]
        point.coor <- coor[idx,]
      }else if(!is.null(group.order)){
        idx <- unlist(sapply(group.order, function(x){
          which(values==x)
        }))
        cellColor <- cellColor[idx]
        point.coor <- coor[idx,]
      }else if(point.order=="random"){
        idx <- sample(seq(length(values)))
        cellColor <- cellColor[idx]
        point.coor <- coor[idx,]
      }
      plot(point.coor, col=cellColor, pch=16, cex.main=5, xaxt="n", yaxt="n", xlab="", ylab="", bty="n", ...)
    }else{
      plot(coor, col="#bdbdbd", pch=16, cex.main=5, xaxt="n", yaxt="n", xlab="", ylab="", bty="n",  ...)
      if(is.null(point.order)){
        points(coor[emphasize,], col=cellColor[emphasize], pch=16, ...)
      }else if(point.order=="sorted"){
        idx <- order(values[emphasize])
        points(coor[emphasize,][idx,], col=cellColor[emphasize][idx], pch=16, ...)
      }else if(point.order=="random"){
        idx <- sample(seq(length(values[emphasize])))
        points(coor[emphasize,][idx,], col=cellColor[emphasize][idx], pch=16, ...)
        
      }
    }
  }else{
    plot(c(min(coor[,1]), max(coor[,1])), c(min(coor[,2]),max(coor[,2])),type="n", cex.main=5, xaxt="n", yaxt="n", xlab="", ylab="", bty="n", ...)
    for(i in seq(nrow(knn.pairs))){
      x.coor <- coor[knn.pairs[i,],1]
      y.coor <- coor[knn.pairs[i,],2]
      lines(x.coor, y.coor, col="#bdbdbd30", lwd=0.5, ...)
    }
    if(is.null(emphasize)){
      points(coor, col=cellColor, pch=16, ...)
    }else{
      points(coor, col=cellColor, col="#bdbdbd", pch=16, ...)
      points(coor[emphasize,], col=cellColor[emphasize], pch=16, ...)
    }
  }

  if(add.cluster){
    for(k in unique(cl.vec)){
      cell.idx=which(cl.vec==k)
      if(length(cell.idx)>1){
        label.coor=apply(coor[cell.idx,], 2, median)
      }else{
        label.coor=coor[cell.idx,]
      }
      if(grepl("Cluster_", k, fixed=T)){
        text(label.coor[1], label.coor[2], label=sub("Cluster_", "", k), cex=label.cex)
      }else{
        text(label.coor[1], label.coor[2], label=sub("Cluster", "", k), cex=label.cex)
      }
    }
    
  }
  
  if(add.label){
    for(x in unique(values)){
      cell.idx=which(values==x)
      if(length(cell.idx)>1){
        label.coor=apply(coor[cell.idx,], 2, median)
      }else{
        label.coor=coor[cell.idx,]
      }
      if(grepl("Cluster_", x, fixed=T)){
        text(label.coor[1], label.coor[2], label=sub("Cluster_", "", x), cex=label.cex)
      }else{
        text(label.coor[1], label.coor[2], label=sub("Cluster", "", x), cex=label.cex)
      }
    }
  }
  if(is.null(group.order)){
    idx <- intersect(names(gCols), unique(values)) 
    sorted.cols <- sort(gCols[idx])
  }else{
    sorted.cols <- gCols[group.order]
  }
  if(add.legend){
    legend(legend.pos, bty="n", text.col=sorted.cols, legend=names(sorted.cols), cex=legend.cex)
  }
  
}


get.DE.genes <- function(seu.obj=hio, g1=7, g2=12, logFC.cutoff=0.25, in.prop.cutoff=0.25, padj.cutoff=0.05){
  g1.idx <- which(hio@meta.data$RNA_snn_res.1==g1)
  g2.idx <- which(hio@meta.data$RNA_snn_res.1==g2)
  group.vec <- rep(c(g1, g2), c(length(g1.idx), length(g2.idx)))
  expr.mat <- as.matrix(cbind(hio@assays$RNA@data[, c(g1.idx, g2.idx)]))
  library("presto")
  mat <- wilcoxauc(expr.mat, group.vec)
  idx <- which(mat$logFC>logFC.cutoff & mat$pct_in>in.prop.cutoff & mat$padj<padj.cutoff)
  res <- mat[idx,]
  return(res)
}


plotGroupFeature <- function(cm.mat=top.cm, groups=NULL, gene.column="gene_name", dr="umap",group.column="cluster", add.cluster=T, cluster.col="RNA_snn_res.1", seu.obj=duo, group.prefix="Cluster", plot.name=NULL, do.par=FALSE, core.num=5, nCols=c("dark green", "yellow", "red"), point.order= "random"){
  if(is.null(groups)){
    groups <- unique(cm.mat[[group.column]])
  }
  expr.mat = seu.obj@assays$RNA@data
  if(do.par){
    library(doParallel)
    registerDoParallel(core.num)
    score.mat <- foreach(group=groups, .multicombine = T, .combine='cbind')%dopar%{
      idx <- which(cm.mat[[group.column]]==group)
      gene <- cm.mat[[gene.column]][idx]
      score <- apply(scale(t(as.matrix(expr.mat[gene,]))), 1, sum)
      return(score)
    }
  }else{
    score.mat <- sapply(groups, function(group){
      idx <- which(cm.mat[[group.column]]==group)
      gene <- cm.mat[[gene.column]][idx]
      score <- apply(scale(t(as.matrix(expr.mat[gene,]))), 1, sum)
      return(score)
    })
  }
  if(!is.null(group.prefix)){
    colnames(score.mat) <- paste(group.prefix, groups, sep="_")
  }else{
    colnames(score.mat) <- groups
  }
  
  col.num=8
  row.num=ceiling(length(groups)/col.num)
  if(is.null(plot.name)){
    plot.name <- "Plot_summarized_group_feature_pattern.png"
  }
  plot.coor=Embeddings(seu.obj, reduction=dr)
  png(plot.name, width=1000*col.num, height=1000*row.num)
  par(mfrow=c(row.num, col.num))
  for(x in colnames(score.mat)){
    plotFeature2(plot.coor, values=score.mat[,x], main=x, cex.main=5, cex=2, add.cluster=add.cluster, cl.vec=seu.obj@meta.data[[cluster.col]], label.cex=3, nCols = nCols, point.order= point.order)
  }
  dev.off()
}

plotFeature.batch <- function(seu.obj=hio, assay.type="RNA", data.type="data", add.group.info=FALSE, group.by="RNA_snn_res.1", dr="umap", genes.to.plot=NULL, plot.name=NULL, nCols=c("darkgreen", "yellow","red"), point.order="sorted",col.num=6,label.cex = 3,  cex.main=7, cex=3, per.plot.size=1000, plot.unit="px", plot.res=NA, emphasize=NULL){
  genes.to.plot <- intersect(genes.to.plot, rownames(seu.obj))
  row.num=ceiling(length(genes.to.plot)/col.num)
  plot.coor=Embeddings(seu.obj, reduction = dr)
  if(is.null(plot.name)){
    plot.name <- "Plot_UMAP_selected_gene_feature_plot.png"
  }
  if(add.group.info){
    group.vec=seu.obj@meta.data[[group.by]]
    if(is.factor(group.vec)){
      group.vec = paste0("Cluster", group.vec)
    }
  }
  png(plot.name, height=per.plot.size*row.num, width=per.plot.size*col.num, unit=plot.unit, res=plot.res)
  par(mfrow=c(row.num, col.num))
  for(g in genes.to.plot){
    #plotFeature2(plot.coor, values=slot(seu.obj@assays[[assay.type]], data.type)[g,], main=g, add.cluster=add.group.info, cl.vec=group.vec, label.cex = 3,  cex.main=5,cex=2, nCols = nCols, point.order= point.order)	
    plotFeature2(plot.coor, values=slot(seu.obj@assays[[assay.type]], data.type)[g,], main=g, add.cluster=add.group.info, cl.vec=group.vec, label.cex = label.cex,  cex.main=cex.main, cex=cex, nCols = nCols, point.order= point.order, emphasize = emphasize)	
  }
  dev.off()
}

plotClusterHeatmap <- function(gene.to.plot, cl.expr, plot.name=NULL, gene.reorder=TRUE,cl.reorder.by.hc=TRUE, cl.reorder=NULL,labCol.subset=TRUE, max.diag=TRUE,expr.cols=c("#404040", "#bababa", "#ffffff", "#f4a582", "#ca0020"), row.names=FALSE, do.scale="none",col.angle=30, plot.cexCol=1.2,plot.cexRow=1, plot.height=6, plot.width=5, plot.cl.tree=TRUE, tree.plot.name=NULL, return.tree=TRUE, return.value=FALSE){
  cm.expr <- cl.expr[gene.to.plot,]
  if(gene.reorder){
    print("Reorder genes based on hierarchical clustering")
    hc.gene <- hclust(as.dist(1-cor(t(cm.expr))))
    gene.order = hc.gene$order
  }else{
    print("No reordering for genes")
    gene.order <- rownames(cm.expr)
  }
  
  hc.cl <- hclust(as.dist(1-cor(cm.expr)))
  if(plot.cl.tree){
	if(is.null(tree.plot.name)){
		tree.plot.name <- "Plot_hc_selected_clusters.pdf"
	}
	pdf(tree.plot.name)
	plot(hc.cl, hang=-1)
	dev.off()
  }

  if(cl.reorder.by.hc){
    print("Reorder clusters based on hierarchical clustering")
    cl.order = hc.cl$order
  }else if(length(cl.reorder)==ncol(cl.expr)){
  	print("Reorder clusters based on user provided order")
	  cl.order = cl.reorder
  }else{
    print("No reordering for clusters")
    cl.order = colnames(cm.expr)
  }
  
  input <- cm.expr[gene.order, cl.order]
  if(max.diag){
    max.cl.idx <- colnames(input)[apply(input, 1, which.max)]
    idx <- unlist(lapply(colnames(input), function(i){
      which(max.cl.idx==i)
    }))
    input <- input[idx,] 
  }
  
  if(is.null(plot.name)){
    plot.name <- "Plot_heatmap_selected_genes_cluster_average_expr.pdf"
  }
  if(row.names){
  	row.names=rownames(input)
  }
  colors <- colorRampPalette(expr.cols)(50)
  library(gplots)
  if(labCol.subset==TRUE){
    colname.vec = sapply(colnames(input), function(x){strsplit(x,"_")[[1]][2]})
  }else{
    colname.vec = colnames(input)
  }
  pdf(plot.name, height=plot.height, width=plot.width)
  heatmap.2(input,trace="none",main="", density.info="density",dendrogram="none",labCol = colname.vec, scale=do.scale, labRow=row.names, col=colors,Rowv=FALSE, cexRow=plot.cexRow, cexCol = plot.cexCol, Colv=FALSE, srtCol=col.angle, margins = c(8, 5))
  dev.off()
  if(return.value){
    res <- list("hc.cl"=hc.cl, "rownames"=rownames(input), "colnames"=colnames(input))
    return(res)
  }
}

prepareTreeAndHeatmapInput <- function(expr.mat=hio.expr, cor.method="pearson", hc.method="average", genes.to.highlight=NULL){
  hc.col <- hclust(as.dist(1-cor(expr.mat, method=cor.method)), method=hc.method)
  col.orders <- hc.col$labels[hc.col$order]
  hc.row <- hclust(as.dist(1-cor(t(expr.mat), method=cor.method)), method=hc.method)
  row.orders <- hc.row$labels[hc.row$order]
  
  expr.mat <- expr.mat[row.orders, col.orders]
  max.ct <- colnames(expr.mat)[apply(expr.mat, 1, which.max)]
  idx <- unlist(lapply(colnames(expr.mat), function(x){
    which(max.ct==x)
  }))
  expr.mat <- expr.mat[idx,]
  normed.expr <- t(apply(expr.mat, 1, function(vec){(vec-min(vec))/(max(vec)-min(vec))}))
  highlight.idx <- ifelse(rownames(normed.expr)%in%genes.to.highlight, 1, NA)
  input <- cbind(normed.expr, highlight.idx)
  res <- list("heatmap_input"=input, "hc_row"=hc.row, "hc_col"=hc.col, "highlighted_genes"=intersect(rownames(normed.expr), genes.to.highlight))
  return(res)
}




getKNNEdge <- function(dis.mat, k=15){
        idx1 <- apply(dis.mat, 2, function(vec){
                order(vec)[2:(k+1)]
        })      
        knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(dis.mat)), each=k))
        return(knn.idx)
}

getDotPlot <- function(cl.expr.mat, cl.prop.mat, gene.reorder=TRUE, specified.order=NULL, cl.reorder.by.hc=TRUE, genes=NULL, colors=c("#d9d9d9", "#252525"), point.size.factor=3.5, plot.name=NULL, plot.height=20, plot.width=30, plot.margin=c(8,12,5,8), plot.cex=1.5, max.diag=TRUE, col.space.factor=1, row.space.factor=1){
  if(is.null(genes)){
    genes <- rownames(cl.expr.mat)
  }
  expr <- cl.expr.mat[genes,]
  sd.vec <- apply(expr, 1, sd)
  expr <- expr[which(sd.vec>0),]
  if(gene.reorder){
    print("Reorder genes based on hierarchical clustering")
    hc.gene <- hclust(as.dist(1-cor(t(expr))))
    gene.order = hc.gene$order
  }else{
    print("No reordering for genes")
    gene.order <- rownames(expr)
  }
  expr <- expr[gene.order,] 
  
  if(!is.null(specified.order)){
    print("Reorder clusters according to user provided order")
    expr <- expr[,specified.order]
  }else{
    if(cl.reorder.by.hc){
      print("Reorder clusters based on hierarchical clustering of provided expression matrix")
      hc <- hclust(as.dist(1-cor(expr)))
      expr <- expr[,hc$order]
    }else{
      print("Do not reorder clusters")
      expr <- expr
    }
  }

  if(max.diag){
    max.cl.idx <- colnames(expr)[apply(expr, 1, which.max)]
    idx <- unlist(lapply(colnames(expr), function(i){
      which(max.cl.idx==i)
    }))
    expr <- expr[idx,]
  }
  scale.cl.expr <- t(scale(t(expr)))
  scale.cl.expr <- rbind(scale.cl.expr, rep(c(-1,1), c(ncol(scale.cl.expr)-4,4)))
  colorPal <- grDevices::colorRampPalette(colors)
  color.break.num <- 30
  cellColor <- as.vector(t(apply(scale.cl.expr, 1, function(values){
    adjustcolor(colorPal(color.break.num), alpha=.8)[as.numeric(cut(values, breaks=color.break.num, right=F, include.lowest=T))]
  })))
  cl.num <- ncol(expr)
  prop <- cl.prop.mat[rownames(expr),colnames(expr)]
  new.vec <- rep(0, ncol(prop))
  new.vec[ncol(prop)-0:3] <- seq(from=1,to=0.25,length=4)
  prop <- rbind(prop, new.vec)
  point.size <- as.vector(prop)*point.size.factor
  
  g1 <- c(rownames(expr), "")
  if(is.null(plot.name)){
    plot.name <- "DotPlot_selected_genes.pdf"
  }
  pdf(plot.name, height=plot.height, width=plot.width)
  par(mar=plot.margin, xpd=TRUE)
  plot(c(1, length(g1)), c(1, cl.num), type="n", bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
  points(rep(seq(length(g1))*col.space.factor, cl.num), 1+rep(seq(cl.num)*row.space.factor, each=length(g1)), pch=16, col=cellColor, cex=point.size)
  mtext(g1, side=1, at=seq(length(g1))*col.space.factor, las=2, cex=plot.cex)
  mtext(sub("Cluster_", "", colnames(expr)), side=2, at=1+seq(cl.num)*row.space.factor, las=1, cex=plot.cex)
  #points(rep((length(g1)+2)*col.space.factor, 4), c(0.5*cl.num-2.9, 0.5*cl.num-1.75, 0.5*cl.num-0.5, 0.5*cl.num+0.75), pch=16, cex=c(0.25,0.5,0.75,1)*point.size.factor,xpd=TRUE)
  text(rep((length(g1)+2)*col.space.factor, 4)+0.5, cl.num-3:0, labels=c(0.25,0.5,0.75,1), cex=plot.cex, pos=4)
  dev.off()
}

getAveExpr <- function(input.type="Seurat", seu.obj, assay.type="RNA", data.type="data", feature.to.calc, X=NULL, y=NULL, method="mean", specified.order=NULL, size.cutoff=5, genes=NULL, do.par=TRUE, core.num=5, add.cluster=T, colname.prefix="Cluster", check.barcode=T){
	if(input.type=="matrix"){
	  library(doParallel)
	  registerDoParallel(core.num)
	  if(is.null(specified.order)){
	    if(is.numeric(class(y))){
	      sorted.features <- sort(as.numeric(unique(y)))
	    }else{
	      sorted.features <- sort(unique(y))
	    }
	  }else{
	    sorted.features <- specified.order
	  }
	  size.vec <- sapply(sorted.features, function(x){
	    sum(y==x)
	  })
	  sorted.features <- sorted.features[which(size.vec>size.cutoff)]
	  if(is.null(genes)){
	    genes <- rownames(X)
	  }
	  ave.expr <- foreach(k=sorted.features, .multicombine=T, .combine='cbind')%dopar%{
	    sample.idx <- which(y==k)
	    if(method=="mean"){
	      #res <- apply(X[genes, sample.idx], 1, mean)
	      res <- rowMeans(X[genes, sample.idx])
	    }else if(method=="median"){
	      res <- apply(X[genes, sample.idx], 1, median)
	    }
	    return(res)
	  }
	  stopImplicitCluster() 
	  rownames(ave.expr) <- genes
	  if(is.null(colname.prefix)){
	    colnames(ave.expr) <- sorted.features
	  }else{
	    colnames(ave.expr) <- paste(colname.prefix, sorted.features, sep="_")
	  }
	  return(ave.expr)
	}else{
    meta.data=seu.obj@meta.data
	  expr=slot(seu.obj@assays[[assay.type]], data.type)
	  if(!is.null(genes)){
	    expr=expr[genes,]
	  }
    if(check.barcode){
	  	if(sum(rownames(meta.data)!=colnames(expr))>0){
	  		return(cat("Inconsistent meta.data and expression sample names\n"))
	  	}
	  }
	  if(is.null(genes)){
	  	genes=rownames(expr)
	  }
  
	  if(is.null(specified.order)){
	  	if(is.numeric(class(meta.data[,feature.to.calc]))){
	  		sorted.features <- sort(as.numeric(unique(meta.data[,feature.to.calc])))
	  	}else{
	  		sorted.features <- sort(unique(meta.data[,feature.to.calc]))
	  	}
	  }else{
	  	sorted.features <- specified.order
	  }
	  size.vec <- sapply(sorted.features, function(x){
	  	sum(meta.data[, feature.to.calc]==x)
	  })
	  sorted.features <- sorted.features[which(size.vec>size.cutoff)]
	  if(do.par){
	    library(doParallel)
	    registerDoParallel(core.num)
	    ave.expr <- foreach(k=sorted.features, .multicombine=T, .combine='cbind')%dopar%{
	      sample.idx <- which(meta.data[,feature.to.calc]==k)
	      #apply(as.matrix(expr[genes, sample.idx]), 1, mean)
	      e <- rowMeans(expr[genes, sample.idx])
	    }
	    stopImplicitCluster()
	  }else{
	    ave.expr <- c()
	    for(k in sorted.features){
	      sample.idx <- which(meta.data[,feature.to.calc]==k)
	      #e <- apply(as.matrix(expr[genes, sample.idx]), 1, mean)
	      e <- rowMeans(expr[genes, sample.idx])
	      ave.expr <- cbind(ave.expr, e)
	    }
	  }
	  rownames(ave.expr) <- genes
	  if(is.null(colname.prefix)){
	  	colnames(ave.expr) <- sorted.features
	  }else{
	  	colnames(ave.expr) <- paste(colname.prefix, sorted.features, sep="_")
	  }
	  return(ave.expr)
  }
} 

findAllMarkers <- function(seu.obj, seu.version=2, selected.column, core.num=5, p.value.cutoff=1, bg=NULL, cluster.to.test=NULL, pos.only=T, do.par=T){
	meta.data <- seu.obj@meta.data
	if(is.factor(meta.data[,selected.column])){
	  cluster.num <- max(as.numeric(as.character(meta.data[,selected.column])))
	  if(is.null(cluster.to.test)){
	    cat("Attention: Cluster start from 0. \n")
	    cluster.to.test <- 0:cluster.num
	  }
	}else if(is.numeric(meta.data[,selected.column])){
	  if(is.null(cluster.to.test)){
	    cat(paste("Attention: Cluster start from ",min(meta.data[,selected.column]),". \n"))
	    cluster.to.test <- sort(unique(meta.data[,selected.column]))
	  }
	}else if(is.character(meta.data[,selected.column])){
	  if(is.null(cluster.to.test)){
	    cluster.to.test <- sort(unique(meta.data[,selected.column]))
	  }
	}else{
	  print("Error: Selected column is not number, factor or character. Please check.")
	  return()
	}
	
	if(is.null(bg)){
		cat("Using unique(meta.data[,selected.column]) as background.\n")
		bg=unique(meta.data[,selected.column])
	}
	
	if(seu.version==3){
		Idents(seu.obj) <- selected.column
	}else if(seu.version==2){
		original.identity <- seu.obj@ident
		seu.obj@ident <- as.factor(meta.data[,selected.column])
		names(seu.obj@ident) <- seu.obj@cell.names
	}
	
	if(do.par){
		cat("Find markers using Parallel\n")
		library(doParallel)
		registerDoParallel(core.num)
		combined.cm <- foreach(k=cluster.to.test, .multicombine=T, .combine='rbind')%dopar%{
			cat(paste("Find markers for cluster", k, "\n"))
			bg.2 = setdiff(bg, k)
			mat <- FindMarkers(seu.obj, ident.1=k, ident.2=bg.2, only.pos=pos.only, min.pct=0.25, logfc.threshold=0.25)
			dat <- data.frame(mat, "cluster"=rep(k, nrow(mat)), "gene_name"=rownames(mat), stringsAsFactors=F)
			rownames(dat) <- paste(rownames(dat), k, sep="_")
			return(dat)
		}
		stopImplicitCluster()
	}else{
		cat("Find markers in simple loop\n")
		combined.cm <- c()
		for(k in cluster.to.test){
			bg.2 = setdiff(bg, k)
			mat <- FindMarkers(seu.obj, ident.1=k, ident.2=bg.2, only.pos=pos.only, min.pct=0.25, logfc.threshold=0.25)
			dat <- data.frame(mat, "cluster"=rep(k, nrow(mat)), "gene_name"=rownames(mat), stringsAsFactors=F)
			rownames(dat) <- paste(rownames(dat), k, sep="_")
			combined.cm <- rbind(combined.cm, dat)
		}
	}
	combined.cm <- combined.cm[which(combined.cm$p_val<p.value.cutoff),]
	if(seu.version==2){
		seu.obj@ident <- original.identity
	}
	return(combined.cm)
}

highlightCells <- function(coor, cells.to.highlight, point.cex=0.8, bg.col="#d0d0d0", highlight.col="#252525"){
  plot(coor, pch=16, cex=point.cex, col=bg.col)
  points(coor[cells.to.highlight,], pch=16, cex=point.cex, col=highlight.col)
}

prepareSeuratObject <- function(rawData, namePrefix, age, tissue, sex, condition, isolation, sampleName, seu.version=2, mito.cutoff=0.1, min.cell.cutoff=3, min.gene.cutoff=1000, gene.num.low=1000, gene.num.high=20000) {

  ## Make the cell names unique due to reusing 10x barcodes across samples
  if(!is.null(namePrefix)){
	colnames(x = rawData) <- paste(namePrefix, colnames(x = rawData), sep = '_')
  }
  if(grepl("hg19", rownames(rawData)[1])){
  	rownames(x = rawData) <- sub("hg19_","",rownames(rawData))
  }
  # Read in the actual 10x raw data and create a Seurat object
  if(seu.version==2){
	dataset <- CreateSeuratObject(raw.data = rawData, min.cells = min.cell.cutoff, min.genes = min.gene.cutoff)
	# Generate a list of expressed Mitochondria genes
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = dataset@data), value = TRUE)
	# Calculate the percentage of genes in each cell that are mitochondrial and add number as metadata
	percent.mito <- Matrix::colSums(dataset@raw.data[mito.genes, ])/ Matrix::colSums(dataset@raw.data)
	dataset <- AddMetaData(object = dataset, metadata = percent.mito, col.name = "percent.mito")
	# Filter out cells with too few genes and too high a percentage of mitochondrial genes
	dataset <- FilterCells(object = dataset, subset.names = c("nGene", "percent.mito"), low.thresholds = c(gene.num.low, -Inf), high.thresholds = c(gene.num.high, mito.cutoff))
  }else if(seu.version==3){
	dataset <- CreateSeuratObject(counts=rawData, min.cells=min.cell.cutoff, min.features=min.gene.cutoff)
	# Generate a list of expressed Mitochondria genes
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = dataset), value = TRUE)
	# Calculate the percentage of genes in each cell that are mitochondrial and add number as metadata
	percent.mito <- Matrix::colSums(dataset@assays$RNA@counts[mito.genes, ])/ Matrix::colSums(dataset@assays$RNA@counts)
	dataset[["percent.mito"]] <- percent.mito
	# Filter out cells with too few genes and too high a percentage of mitochondrial genes
	dataset <- subset(dataset, cells=rownames(dataset[[]])[which(dataset[["nFeature_RNA"]]>gene.num.low & dataset[["nFeature_RNA"]]<gene.num.high & dataset[["percent.mito"]]<mito.cutoff)])
  }
 

  # Add in metadata for sample name, sample age and sample tissue type
  dataset@meta.data[, "Sample"] <- sampleName
  dataset@meta.data[, "Age"] <- age
  dataset@meta.data[, "Tissue"] <- tissue
  dataset@meta.data[, "Sex"] <- sex
  dataset@meta.data[, "Condition"] <- condition
  dataset@meta.data[, "Isolation"] <- isolation
  # Change the Seurat object to sparse format to save memory
  if(seu.version==2){
	dataset <- MakeSparse(object = dataset)
  }
  # Return the finished sparse format Seurat Object 
  return(dataset)
}

prepareSeuratObject2 <- function(rawData, namePrefix, additional.sample.info.mat, seu.version=2, mito.cutoff=0.1, min.cell.cutoff=3, min.gene.cutoff=1000, gene.num.low=1000, gene.num.high=20000) {

  ## Make the cell names unique due to reusing 10x barcodes across samples
  if(!is.null(namePrefix)){
	colnames(x = rawData) <- paste(namePrefix, colnames(x = rawData), sep = '_')
  }
  if(grepl("hg19", rownames(rawData)[1])){
  	rownames(x = rawData) <- sub("hg19_","",rownames(rawData))
  }
  # Read in the actual 10x raw data and create a Seurat object
  if(seu.version==2){
	dataset <- CreateSeuratObject(raw.data = rawData, min.cells = min.cell.cutoff, min.genes = min.gene.cutoff)
	# Generate a list of expressed Mitochondria genes
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = dataset@data), value = TRUE)
	# Calculate the percentage of genes in each cell that are mitochondrial and add number as metadata
	percent.mito <- Matrix::colSums(dataset@raw.data[mito.genes, ])/ Matrix::colSums(dataset@raw.data)
	dataset <- AddMetaData(object = dataset, metadata = percent.mito, col.name = "percent.mito")
	# Filter out cells with too few genes and too high a percentage of mitochondrial genes
	dataset <- FilterCells(object = dataset, subset.names = c("nGene", "percent.mito"), low.thresholds = c(gene.num.low, -Inf), high.thresholds = c(gene.num.high, mito.cutoff))
  }else if(seu.version==3){
	dataset <- CreateSeuratObject(counts=rawData, min.cells=min.cell.cutoff, min.features=min.gene.cutoff)
	# Generate a list of expressed Mitochondria genes
	mito.genes <- grep(pattern = "^MT-", x = rownames(x = dataset), value = TRUE)
	# Calculate the percentage of genes in each cell that are mitochondrial and add number as metadata
	percent.mito <- Matrix::colSums(dataset@assays$RNA@counts[mito.genes, ])/ Matrix::colSums(dataset@assays$RNA@counts)
	dataset[["percent.mito"]] <- percent.mito
	# Filter out cells with too few genes and too high a percentage of mitochondrial genes
	dataset <- subset(dataset, cells=rownames(dataset[[]])[which(dataset[["nFeature_RNA"]]>gene.num.low & dataset[["nFeature_RNA"]]<gene.num.high & dataset[["percent.mito"]]<mito.cutoff)])
  }

  # Add in metadata for sample name, sample age and sample tissue type
  for(i in seq(nrow(additional.sample.info.mat))){
	key=additional.sample.info.mat[i,1]
	value=additional.sample.info.mat[i,2]
	dataset@meta.data[, key] <- value
  }
  # Change the Seurat object to sparse format to save memory
  if(seu.version==2){
	dataset <- MakeSparse(object = dataset)
  }
  # Return the finished sparse format Seurat Object 
  return(dataset)
}

plotTSNEs <- function(dataset, plot.name=NULL, only.plot.cluster=cluster.only){
	if(!only.plot.cluster){
		p2 <- TSNEPlot(object=dataset, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
		p3 <- TSNEPlot(object=dataset, pt.size=1, group.by="orig.ident", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
		p4 <- TSNEPlot(object=dataset, pt.size=1, group.by="Sex", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
		p5 <- TSNEPlot(object=dataset, pt.size=1, group.by="Isolation", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
		p6 <- TSNEPlot(object=dataset, pt.size=1, group.by="Condition", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
		p7 <- TSNEPlot(object=dataset, pt.size=1, group.by="Age", do.return=T, label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
		if(is.null(plot.name)){
			plot.name="Plot_tSNE.png"
		}
		png(plot.name, width=1500, height=1000)
		print(plot_grid(p2, p3, p4, p5, p6, p7, ncol=3, nrow=2))
		dev.off()
	}else{
		p2 <- TSNEPlot(object=dataset, pt.size=1, do.return=T, do.label=T,label.size=5)+theme(axis.text=element_text(size=15), axis.ticks.length=unit(.25, "cm"), axis.title=element_text(size=18),legend.text=element_text(size=15))
		if(is.null(plot.name)){
			plot.name="Plot_tSNE.png"
		}
		png(plot.name, width=1500, height=1000)
		print(plot(p2))
		dev.off()
	}
}

preprocessSubset <- function(dataset, seu.version=2, field.to.select="orig.ident", subset.idx=idx, clean=TRUE, regressed.vars=c("percent.mito", "nUMI"), reso=0.8, pc.num=20, hvg.to.exclude=NULL, tsne.plot.name=NULL, cluster.only=FALSE, plot.known.markers=TRUE, marker.list.file=NULL, marker.plot.name=NULL, save.rds=TRUE, rds.file.name=NULL, do.return=TRUE){

	if(seu.version==2){
	eu.obj <- SubsetData(dataset, cells.use = rownames(dataset@meta.data)[which(dataset@meta.data[[field.to.select]]==subset.idx)], do.clean=clean, subset.raw=T)
	}else if(seu.version==3){
		return(cat("Error: Currently, this function is only applicable with Seurat version 2.X\n"))
	}
	if(clean){
		seu.obj <- ScaleData(object=seu.obj, vars.to.regress=regressed.vars, do.par=T)	
	}
	seu.obj <- FindVariableGenes(object=seu.obj, do.plot=FALSE, mean.function = ExpMean, dispersion.function =LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
	if(is.null(hvg.to.exclude)){
		selected.hvg <- seu.obj@var.genes
	}else{
		selected.hvg <- setdiff(seu.obj@var.genes, hvg.to.exclude)
	}
	seu.obj <- RunPCA(object=seu.obj, pc.genes=selected.hvg, do.print=FALSE,pcs.compute=pc.num)
	seu.obj <- FindClusters(object=seu.obj, reduction.type="pca", dims.use=1:pc.num, resolution=reso, force.recalc=TRUE, print.output=0, save.SNN=TRUE)
	seu.obj <- RunTSNE(object=seu.obj, dims.use=1:pc.num, do.fast=TRUE)

	if(is.null(tsne.plot.name)){
		tsne.plot.name = paste0("Plot_tSNE_res",reso,"_",pc.num,"PC.png")
	}
	plotTSNEs(dataset=seu.obj, plot.name=tsne.plot.name, only.plot.cluster=cluster.only)
	
	# plot canonical cell type marker expression on tSNE
	if(plot.known.markers){
		markers <- read.table(marker.list.file,sep="\t",stringsAsFactors=F)
		markers <- markers[markers[,2]%in%rownames(seu.obj@data),]
		g1 <- markers[,2]
		if(is.null(marker.plot.name)){
			marker.plot.name <- "Plot_marker_expression_on_tSNE_hvg.png"
		}
		png(marker.plot.name, width=1600, height=4800)
		FeaturePlot(object=seu.obj, features.plot=g1, cols.use=c("gray", "blue"), no.legend=F)
		dev.off()
	}

	if(save.rds){
		if(is.null(rds.file.name)){
			rds.file.name <- paste0("Res_",idx,".rds")
		}
		saveRDS(seu.obj, file=rds.file.name)
	}
	
	if(do.return){
		return(seu.obj)
	}
}

plotSPRING <- function(spring.coor.dir, dis.mat=NULL, genes=NULL, marker.plot.name=NULL, expr.mat=NULL, hvg.info=NULL, plot.sample.info=TRUE, sample.info.column="Sample", sample.plot.name=NULL, plot.cluster.info=TRUE, cluster.info.column="merge.cl.idx", cluster.plot.name=NULL){
	
	cat("Read spring coordinate file\n")
	spring.coor <- read.csv(paste0(spring.coor.dir,"coordinates.txt"),head=F,row.names=1)
	rownames(spring.coor) <- colnames(hvg.info)
	# get the kNN network
	knn.idx=NULL
	if(!is.null(dis.mat)){
		cat("Get k nearest neighbors, k=15\n")
		k <- 15
		idx1 <- apply(dis.mat, 2, function(vec){
			order(vec)[2:(k+1)]
		})
		knn.idx <- cbind(as.vector(idx1), rep(seq(ncol(dis.mat)), each=k))
	}
	
	# plot the marker gene expression 	
	if(!is.null(genes)){
		if(sum(!colnames(expr.mat)%in%colnames(hvg.info))>0){
			cat("Cell numbers in expression matrix is not equal to that in hvg.info\n")
			return(NA)
		}
		expr.mat <- expr.mat[,colnames(hvg.info)]
		if(min(expr.mat)<0){
			data.type="scaled data"
		}else{
			data.type="normalized data"
		}
		cat(paste0("To plot gene expression levels using ", data.type, "\n"))
		g1 <- intersect(genes, rownames(expr.mat))
		sd.vec <- apply(expr.mat[g1,],1,sd)
		g2 <- g1[which(sd.vec>0)]
		row.num <- ceiling(length(g2)/4)
		spring.file <- ifelse(is.null(marker.plot.name),  paste0(spring.coor.dir,"Plot_SPRING_marker_gene_expr.png"), marker.plot.name)
		png(spring.file, width=4*400, height=row.num*400)
		par(mfrow=c(row.num,4))
		for(x in g2){
			plotFeature2(coor=spring.coor, values=expr.mat[x,], knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=2, bty="n", xaxt="n", yaxt="n", cex=1.2)
		}
		dev.off()
	}
	
	# plot sample distribution
	if(plot.sample.info){
		cat(paste0("To plot sample distribution: ", sample.info.column, "\n"))
		spring.file <- ifelse(is.null(sample.plot.name), paste0(spring.coor.dir,"Plot_SPRING_sample_info.png"), sample.plot.name)
		sample.info <- hvg.info[sample.info.column,]
		sample.num <- length(unique(sample.info))
		row.num <- ceiling((sample.num+1)/4)
		png(spring.file, width=4*400, height=row.num*400)
		par(mfrow=c(row.num, 4))
		plotFeature2(coor=spring.coor, values=sample.info, knn.pairs=knn.idx, xlab="", ylab="", main=sample.info.column, cex.lab=1.3, cex.axis=1.3, cex.main=2, bty="n", xaxt="n", yaxt="n", cex=1.2)
		for(x in unique(sample.info)){
			plotFeature2(coor=spring.coor, values=sample.info, emphasize=which(sample.info==x), knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=2, bty="n", xaxt="n", yaxt="n", cex=1.2)
		}
		dev.off()
	}

	# plot cluster distribution
	if(plot.cluster.info){
		cat(paste0("To plot cluster distribution: ", cluster.info.column, "\n"))
		spring.file <- ifelse(is.null(cluster.plot.name), paste0(spring.coor.dir,"Plot_SPRING_cluster_info.png"), cluster.plot.name)
		cluster.info <- as.factor(as.numeric(hvg.info[cluster.info.column, ]))
		cluster.num <- length(unique(cluster.info))+1
		row.num <- ceiling(cluster.num/4)
		png(spring.file, height=row.num*400, width=4*400)
		par(mfrow=c(row.num, 4))
		plotFeature2(coor=spring.coor, values=cluster.info, knn.pairs=knn.idx, xlab="", ylab="", main=cluster.info.column, cex.lab=1.3, cex.axis=1.3, cex.main=2, bty="n", xaxt="n", yaxt="n", cex=1.2)
		for(x in seq(cluster.num-1)){
			cell.idx <- which(cluster.info==x)
			vec <- apply(spring.coor[cell.idx,], 2, median)
			text(vec[1], vec[2], label=x, cex=1.3)
		}
		for(x in sort(unique(cluster.info))){
			plotFeature2(coor=spring.coor, values=cluster.info, emphasize=which(cluster.info==x), knn.pairs=knn.idx, xlab="", ylab="", main=x, cex.lab=1.3, cex.axis=1.3, cex.main=2, bty="n", xaxt="n", yaxt="n", cex=1.2)
		}
		dev.off()
	}
}

quantifyClusterLinkage <- function(knn.idx=knn.idx, group.info=cluster.info, group.prefix="Cluster", p.adjust.method="bonferroni", enrichment.cutoff=0.05){
	groups <- sort(unique(group.info))
	p.value <- matrix(NA, nrow=length(groups), ncol=length(groups))
	or <- matrix(NA, nrow=length(groups), ncol=length(groups))
	link.num <- matrix(NA, nrow=length(groups), ncol=length(groups))
	rownames(p.value) <- paste0(group.prefix, groups)
	colnames(p.value) <- paste0(group.prefix, groups)
	rownames(or) <- paste0(group.prefix, groups)
	colnames(or) <- paste0(group.prefix, groups)
	rownames(link.num) <- paste0(group.prefix, groups)
	colnames(link.num) <- paste0(group.prefix, groups)
	for(group.i in groups){
		for(group.j in groups){
			if(group.i==group.j){
				next
			}
			cl.i.idx <- which(group.info==group.i)
			cl.j.idx <- which(group.info==group.j)
			cl.other.idx <- which(group.info!=group.i & group.info!=group.j)
			a <- sum(knn.idx[,1]%in%cl.i.idx & knn.idx[,2]%in%cl.j.idx) + sum(knn.idx[,2]%in%cl.i.idx & knn.idx[,1]%in%cl.j.idx)
			b <- sum(knn.idx[,1]%in%cl.i.idx & knn.idx[,2]%in%cl.other.idx) + sum(knn.idx[,1]%in%cl.other.idx & knn.idx[,2]%in%cl.i.idx)
			c <- length(cl.j.idx)
			d <- length(cl.other.idx)
			mod <- fisher.test(matrix(c(a,b,c,d),c(2,2)), alternative="g")
			p.value[which(groups==group.i), which(groups==group.j)] <- mod$p.value
			or[which(groups==group.i), which(groups==group.j)] <- mod$estimate	
			link.num[which(groups==group.i), which(groups==group.j)] <- a
		}
	}
	fdr <- matrix(p.adjust(p.value, method=p.adjust.method), nrow=nrow(p.value))
	rownames(fdr) <- paste0(group.prefix, groups)
	colnames(fdr) <- paste0(group.prefix, groups)
	enrichment <- fdr<enrichment.cutoff | t(fdr<enrichment.cutoff)
	cat(paste("Use",p.adjust.method,"correction for P value adjustment\n"))
	cat(paste("Use",enrichment.cutoff,"as corrected P value cutoff\n"))
	cat(paste("In total", sum(enrichment, na.rm=T)/2, "significant between group links\n"))
	res <- list("p.value"=p.value, "or"=or, "link.num"=link.num, "fdr"=fdr, "enrichment"=enrichment)
	return(res)

}
		
plotGroupLinks <- function(enrichment, or, group.info, coor, plot.file){	
	input.i <- c()
	input.j <- c()
	input.lw <- c()
	group.pairs <- c()
	groups <-  sort(unique(group.info))
	for(i in 1:(ncol(enrichment)-1)){
		for(j in (i+1):nrow(enrichment)){
			if(enrichment[j,i]){
				cl.i <- groups[i]
				cl.j <- groups[j]
				cell.idx.i <- which(group.info==cl.i)
				cell.idx.j <- which(group.info==cl.j)
				coor.i <- apply(coor[cell.idx.i, ], 2, median)
				coor.j <- apply(coor[cell.idx.j, ], 2, median)
				lw=or[i,j]
				input.i <- rbind(input.i, coor.i)
				input.j <- rbind(input.j, coor.j)
				input.lw <- c(input.lw, lw)
				group.pairs <- rbind(group.pairs, c(cl.i,cl.j))
			}

		}
	}
	if(is.factor(groups)){
		group.pairs = group.pairs-1
	}
	
	pdf(plot.file)
	plot(rbind(input.i, input.j), pch=16, xlab="", ylab="", main="", bty="n", xaxt="n", yaxt="n", type="n")
	for(k in seq(length(input.lw))){
		lines(c(input.i[k,1], input.j[k,1]), c(input.i[k,2], input.j[k,2]), lwd=sqrt(input.lw[k]), col="#bdbdbd")
	}
	for(k in seq(length(input.lw))){
		text(input.i[k,1], input.i[k,2], label=group.pairs[k,1])
		text(input.j[k,1], input.j[k,2], label=group.pairs[k,2])
	}
	dev.off()
}

quantifyAndPlotClusterLinkage <- function(data=NULL, info=NULL, dis.mat=NULL, data.dir, cell.coor.dir, distance.type, k, cluster.info.idx, p.adjust.method="bonferroni", enrichment.cutoff=0.05, return=FALSE, to.plot=TRUE, plot.name=NULL){
	
	if(is.null(data)){
		cat("Read data matrix\n")
		data.file <- paste0(data.dir, "Table_data.csv")
		data <- read.csv(data.file, header=FALSE, row.names=1)

	}else{
		cat("Data matrix provided\n")
	}
	
	if(is.null(info)){
		cat("Read sample information matrix\n")
		info.file <- paste0(data.dir, "Table_meta_data.csv")
		info <- read.csv(info.file, header=FALSE, row.names=1)
	}else{
		cat("Cell information matrix provided\n")
	}

	#calculate distance matrix 
	if(is.null(dis.mat)){
		if(distance.type=="cor"){
			cat("Calculate correlation distance\n")
			dis.mat <- 1-cor(data)
		}else if(distance.type=="eucl"){
			cat("Calculate Euclidean distance\n")
			dis.mat <- as.matrix(dist(t(data)))
		}else if(is.null(distance.type)){
			return(cat("Please specify distance calculation method. Correlation distance and Euclidean distance are supported.\n"))
		}
	}else{
		cat("Distance matrix provided\n")
	}
	
	cat(paste("Choose top", k, "nearest neighbors\n"))
	knn.idx <- getKNNEdge(dis.mat, k=15)
	cluster.info <- as.numeric(as.matrix(info[cluster.info.idx,]))
	
	res <- quantifyGroupLinkage(knn.idx=knn.idx, group.info=cluster.info,  p.adjust.method=p.adjust.method, enrichment.cutoff=enrichment.cutoff)
	p.value <- res$p.value
	or <- res$or
	link.num <- res$link.num
	fdr <- res$fdr
	enrichment $res$enrichment

	if(to.plot){
		cat("To plot\n")
		cat("Read cell coordinate file\n")
		spring.coor <- read.csv(paste0(cell.coor.dir,"coordinates.txt"),head=F,row.names=1)
		plotGroupLinks(enrichment=enrichment, group.info=cluster.info, coor=spring.coor, plot.file=ifelse(is.null(plot.name), paste0(cell.coor.dir,"Plot_cluster_linkage.pdf")))
	}

	
	if(return){
		cat("To return values\n")
		res <- list("data"=data, "cell_info"=info, "dis_mat"=dis.mat, "knn.idx"=knn.idx, "p.value"=p.value, "odds_ratio"=or, "fdr"=fdr, "enrichment"=enrichment, "link.num"=link.num)
		return(res)
	}
}

splineBasedAgeTest <- function(detected.genes=NULL, pseudotime.vec, expr.mat, df=6, mode=c("relax", "stringent"), core.num=20){
	library(splines)
  library(doParallel)
  registerDoParallel(20)
  if(is.null(detected.genes)){
	  detected.genes <- rownames(expr.mat)
	}
  if(mode=="relax"){
    cat("To compare variance explained by models")
    p.value <- foreach(gene=detected.genes, .combine=c)%dopar%{
      e <- as.vector(as.matrix(expr.mat[gene,]))
      res <- anova(lm(e~ns(pseudotime.vec, df=df)))
      p <- res$`Pr(>F)`[1]
      return(p)
    }
  }else{
    cat("To compare residuals NOT explained by models")
    p.value <- foreach(gene=detected.genes, .combine=c)%dopar%{
      e <- as.vector(as.matrix(expr.mat[gene,]))
      m0 <- lm(e~1)
      m1 <- lm(e~ns(pseudotime.vec, df=df))
      a0 <- anova(m0)
      a1 <- anova(m1)
      p <- pf(a0["Residuals", "Mean Sq"]/a1["Residuals", "Mean Sq"], df1=a0["Residuals", "Df"], df2=a1["Residuals","Df"], lower.tail = F)
      return(p)
    }
  }
	names(p.value) <- detected.genes
	stopImplicitCluster()
	return(p.value)
}

getExpressedProp <- function(seu.obj,assay.type="RNA", data.type="data", feature.to.calc, expr.cutoff=0, specified.order=NULL, size.cutoff=10, genes=NULL, core.num=5, add.cluster=T, colname.prefix="Cluster"){
      meta.data=seu.obj@meta.data
      expr=slot(seu.obj@assays[[assay.type]], data.type)
      if(sum(rownames(meta.data)!=colnames(expr))>0){
                return(cat("Inconsistent meta.data and expression sample names\n"))
      }
     if(is.null(genes)){
             genes=rownames(expr)
     }
     library(doParallel)
     registerDoParallel(core.num)
     if(is.null(specified.order)){
             if(is.numeric(class(meta.data[,feature.to.calc]))){
                     sorted.features <- sort(as.numeric(unique(meta.data[,feature.to.calc])))
             }else{
                     sorted.features <- sort(unique(meta.data[,feature.to.calc]))
             }
     }else{
             sorted.features <- specified.order
     }
     size.vec <- sapply(sorted.features, function(x){
       sum(meta.data[, feature.to.calc]==x)
     })
     sorted.features <- sorted.features[which(size.vec>size.cutoff)]
     expr.prop <- foreach(k=sorted.features, .multicombine=T, .combine='cbind')%dopar%{
             sample.idx <- which(meta.data[,feature.to.calc]==k)
             apply(expr[genes, sample.idx]>expr.cutoff, 1, sum)/length(sample.idx)
     }
     stopImplicitCluster()
     rownames(expr.prop) <- genes
     if(is.null(colname.prefix)){
       colnames(expr.prop) <- sorted.features
     }else{
       colnames(expr.prop) <- paste(colname.prefix, sorted.features, sep="_")
     }
     return(expr.prop)
}

#p_f_detected <- foreach(gene = detected_genes, .combine=c) %dopar%{
#    e <- c(as.numeric(dat_pseudo_human$expr[gene,idx_pt_human]), as.numeric(dat_pseudo_chimp$expr[gene,idx_pt_chimp]))
#    pt <- c(pt_human_aligned, pt_chimp_aligned)
#    g <- as.factor(rep(c("human","chimp"), c(length(idx_pt_human), length(idx_pt_chimp))))
#    m0 <- lm(e ~ ns(pt, df=6))
#    m1 <- lm(e ~ g * ns(pt, df=6))
#    f <- anova(m1)$"Mean Sq"[4] / anova(m0)$"Mean Sq"[2]
#    return(pf(f, m1$df, m0$df, lower.tail=T))
#}
#names(p_f_detected) <- detected_genes


ageTest <- function(expr.mat, v1.vec, v0.vec){

    library(doParallel)
    registerDoParallel(20)
    library(splines)
    pval <- foreach(i = 1:nrow(expr.mat), .combine = rbind) %dopar%{
        e <- as.numeric(expr.mat[i,])
        if (var(e) == 0)
        return(rep(NaN, 2))
             
        if(length(unique(v0.vec))>1){
            m0 <- lm(e ~ v0.vec)
            m1 <- lm(e ~ v0.vec + ns(v1.vec, df = 5))
            a0 <- anova(m0)
            a1 <- anova(m1)
            p1 <- a1$Pr[2]
            p2 <- pf(a0["Residuals","Mean Sq"] / a1["Residuals","Mean Sq"], df1 = a0["Residuals","Df"], df2 = a1["Residuals","Df"], lower.tail = F)
        }else{
            m1 <- lm(e ~ ns(v1.vec, df=5))
            a1 <- anova(m1)
            p1 <- a1$Pr[1]
            p2 <- NA
        }    
        return(c(p1, p2)) 
    }    
    rownames(pval) <- rownames(expr.mat)
    colnames(pval) <- c("ANCOVA","Resi")
    return(pval)

}

plotClusterExprProfile <- function(expr=combined.expr, time.vec=combined.time.points, group.vec=combined.group, cluster.vec=g, group.cols=c("#31a354","#3182bd"), return.value=T, to.plot=T, plot.name=NULL, add.legend=T, legend.pos="topleft", cex.legend=2, ...){
    names(group.cols) <- unique(group.vec)
    scale.expr = t(scale(t(expr)))
    g.list <- lapply(sort(as.numeric(unique(cluster.vec))), function(i){names(which(cluster.vec==i))})
    g.size = sapply(g.list, length)

    mean.list <- list()
    sd.list <- list()

    for(group in unique(group.vec)){
        print(paste(group, "start"))
        e <- scale.expr[,which(group.vec==group)]
        t <- time.vec[which(group.vec==group)]
        mean.mat <- sapply(seq(length(g.list)), function(i){
            apply(e[g.list[[i]],], 2, mean)
        })
        mean.list[[group]] <- mean.mat
        sd.mat <- sapply(seq(length(g.list)), function(i){
            apply(e[g.list[[i]],], 2, sd)
        })
        sd.list[[group]] <- sd.mat
    }
    mean.combined.mat <- do.call("rbind", mean.list)
    sd.combined.mat <- do.call("rbind", sd.list)

    res <- list("g.list"=g.list, "g.size"=g.size, "scale.expr"=scale.expr, "mean.combined.mat"=mean.combined.mat, "sd.combined.mat"=sd.combined.mat, "group.vec"=group.vec, "time.vec"=time.vec)

    if(to.plot){
        if(is.null(plot.name)){
            plot.name <- paste("Plot",paste(unique(group.vec),collapse="_"),"cluster_average_expr.pdf",sep="_")
        }
        col.num=4
        row.num=ceiling(length(g.list)/col.num)
        pdf(plot.name, height=5*row.num, width=5*col.num)
        par(mfrow=c(row.num, col.num))
        for(i in seq(length(g.list))){
            mean.vec <- mean.combined.mat[,i]
            sd.vec <- sd.combined.mat[,i]
            ymax <- max(mean.vec+sd.vec)
            ymin <- min(mean.vec-sd.vec)
            plot(c(min(time.vec), max(time.vec)), c(ymax, ymin), type="n", xlab="Pseudotime quantile",ylab="Relative expr", main=paste("Cluster",i,"(",g.size[i],")"), ...)
            for(group in unique(group.vec)){
                g.idx <- which(group.vec==group)
                g.time <- time.vec[g.idx]
                g.mean <- mean.vec[g.idx]
                g.sd <- sd.vec[g.idx]
                polygon(c(g.time, rev(g.time)), c(g.mean+g.sd, rev(g.mean-g.sd)), border=NA, col=paste0(group.cols[group],"30"))
            }
            for(group in unique(group.vec)){
                g.idx <- which(group.vec==group)
                g.time <- time.vec[g.idx]
                g.mean <- mean.vec[g.idx]
                lines(g.time, g.mean, col=group.cols[group], lwd=2)
            }
            if(add.legend) legend(legend.pos, legend=names(group.cols), text.col=group.cols, bty="n", cex=cex.legend)
        }
        dev.off()

    }

    if(return.value){
        return(res)
    }

}


do.css <- function(que.obj=duo, split.factor="orig.ident", ref.expr=NULL, pc.num=20, cluster.reso=2, ref.gene.to.remove=cc.genes, hc.cutoff=0.01, merge.css=TRUE, k=20, create.dr=TRUE, do.return=TRUE){
  flag=0
  if(is.null(ref.expr)){
    print("Reference cluster expression matrix is NOT provided")
    print("Perform de novo clustering for each sample")
    flag=1
    cl.expr.mat <- c()
    print(paste("To split query object by", split.factor))
    seu.list <- SplitObject(que.obj, split.by = split.factor)
    samples <- names(seu.list)
    for(sample in samples){
      print(paste(sample, "start"))
      seu.obj = seu.list[[sample]] 
      seu.obj = FindVariableFeatures(seu.obj, selection.method="vst", nfeatures=2000)
      seu.obj <- ScaleData(object = seu.obj, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = T)
      sample.hvg <- setdiff(VariableFeatures(seu.obj), ref.gene.to.remove)
      seu.obj <- RunPCA(object=seu.obj, features=sample.hvg, do.print=FALSE)	
      seu.obj <- FindNeighbors(object = seu.obj, dims = 1:pc.num, force.recalc = T, k.param = 15)
      seu.obj <- FindClusters(object = seu.obj, resolution = cluster.reso)
      seu.obj <- RunUMAP(object = seu.obj, dims = 1:pc.num)
      seu.list[[sample]] <- seu.obj
      ave.expr <- getAveExpr(meta.data=seu.obj@meta.data, feature.to.calc=paste0("RNA_snn_res.",cluster.reso), expr=seu.obj@assays$RNA@data, core.num=30, colname.prefix=sample)
      cl.expr.mat <- cbind(cl.expr.mat, ave.expr)
    }
    css.hvg <- SelectIntegrationFeatures(seu.list, nfeatures = 2000)
    css.hvg <- setdiff(css.hvg, ref.gene.to.remove)
    ref.expr <- cl.expr.mat[css.hvg,]
  }else{
    print("Reference cluster expression matrix is provided")
  }
  
  print("Calculate CSS")
  que.expr <- as.matrix(que.obj@assays$RNA@data[rownames(ref.expr),])
  cell2cluster.hvg <- cor(que.expr, ref.expr)
  sample.vec <- sapply(colnames(cell2cluster.hvg), function(x){strsplit(x,"_")[[1]][1]})
  print("Scale CSS within each sample")
  cell2cluster.hvg.norm <- c()
  for(sample in samples){
    mat <- cell2cluster.hvg[,which(sample.vec==sample)]
    norm.mat <- t(apply(mat,1,scale))
    cell2cluster.hvg.norm <- cbind(cell2cluster.hvg.norm, norm.mat)
  }
  colnames(cell2cluster.hvg.norm) <- colnames(cell2cluster.hvg)
  
  # merge similar cluster similarity vector to reduce the influence of dominating but similar clusters
  ## calculate cluster similarity vector using correlation distance
  print("Merge CSS with similar patterns across cells")
  hc <- hclust(as.dist(1-cor(cell2cluster.hvg.norm)))
  g <- cutree(hc, h=hc.cutoff)
  cell2cluster.merge <- sapply(unique(g), function(i){
    idx <- names(which(g==i))
    if(length(idx)==1){
      return(cell2cluster.hvg.norm[,idx])
    }else{
      vec <- apply(cell2cluster.hvg.norm[,idx], 1, mean,na.rm=T)
      return(vec)
    }
  })
  
  ## calculate cluster similarity vector using Euclidean distance
  if(merge.css){
    css.mat <- cell2cluster.merge
  }else{
    css.mat <- cell2cluster.hvg.norm
  }
  library(RANN)
  print("Construct KNN network between cells")
  knn <- RANN::nn2(css.mat, css.mat, k = k+1)$nn.idx[,-1]
  print("Get adjacent matrix between cells")
  adjacent.mat <- .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', knn, 1/15)
  adjacent.mat@Dimnames <- list(rownames(css.mat), rownames(css.mat))
  print("Identify clusters with Seurat inbuilt function")
  cl <- Seurat::FindClusters(Seurat::as.Graph(adjacent.mat), resolution = cluster.reso, verbose = TRUE)[,1]
  names(cl) <- colnames(adjacent.mat)
  print("Run UMAP on adjacent matrix")
  umap.res <- RunUMAP(as.Graph(adjacent.mat))
  umap.coor = umap.res@cell.embeddings
  
  print("Integrate CSS-based results into query Seurat object")
  que.obj[["css"]] <- CreateDimReducObject(embeddings = umap.coor, key="CSS_", assay=DefaultAssay(que.obj))
  que.obj@meta.data[["css_cluster"]] = cl
  
  if(flag){
    css.res <- list("integrated.obj"=que.obj, "seu.list"=seu.list, "ref.expr"=ref.expr, "cell_embeddings"=umap.coor, "cluster"=cl, "merged.css.mat"=cell2cluster.merge, "norm.css.mat"=cell2cluster.hvg.norm, "raw.css.mat"=cell2cluster.hvg)
  }else{
    css.res <- list("integrated.obj"=que.obj, "ref.expr"=ref.expr,"cell_embeddings"=umap.coor, "cluster"=cl, "merged.css.mat"=cell2cluster.merge, "norm.css.mat"=cell2cluster.hvg.norm, "raw.css.mat"=cell2cluster.hvg)
  }
  if(do.return){
    return(css.res)
  }
}

merge.cluster <- function(seu.obj, hc, assay.type=NULL, data.type=NULL, feature.to.merge="SCT_snn_res.2", feature.after.merge="Merged_cluster_idx", front.cl.prop=50, bg.cl.prop=10, auc.cutoff=0.6, logfc.cutoff=0.5, de.num.cutoff=10, all.confound.genes=NULL){
  if(is.null(assay.type)){
    assay.type=seu.obj@active.assay
    print(paste("Used assay:",assay.type))
  }
  if(is.null(data.type)){
    data.type="data"
    print(paste("Used data:",data.type))
  }
  hc.merge <- hc$merge
  get.merged.cl.idx <- function(row.idx, cl.before.merge, cell.cl.before.merge){
    input.pair=hc.merge[row.idx,]
    if(input.pair[1]<0 & input.pair[2]<0){
      old.cl <- abs(input.pair)-1
    }else if(input.pair[1]>0 & input.pair[2]>0){
      old.cl <- cl.before.merge[input.pair]
    }else{
      old.cl <- c(abs(min(input.pair))-1, cl.before.merge[max(input.pair)])
    }
    new.cl <- min(old.cl)
    cell.cl.before.merge[which(cell.cl.before.merge%in%old.cl)] <- new.cl
    cl.before.merge[row.idx] <- new.cl
    cell.cl.after.merge <- cell.cl.before.merge
    cl.after.merge <- cl.before.merge
    res <- list("cl.after.merge"=cl.after.merge, "cell.cl.after.merge"=cell.cl.after.merge)
    return(res)
  }
  
  cl.idx <- rep(NA, nrow(hc.merge))
  cell.cl.vec <- seu.obj@meta.data[[feature.to.merge]]
  to.check=T
  count=0
  while(to.check & count<nrow(hc.merge)){
    count <- count+1
    print(paste("Iteration", count))
    cl.vec <- hc.merge[count,]
    cl.1 <- ifelse(cl.vec[1]<0, abs(cl.vec[1])-1, cl.idx[cl.vec[1]])
    cl.2 <- ifelse(cl.vec[2]<0, abs(cl.vec[2])-1, cl.idx[cl.vec[2]])
    cell.idx <- c(which(cell.cl.vec==cl.1), which(cell.cl.vec==cl.2))
    expr <- slot(seu.obj@assays[[assay.type]], data.type)[,cell.idx]
    group.vec <- cell.cl.vec[cell.idx]
    de.res <- wilcoxauc(X=expr, y=group.vec)
    fil.idx <- which(de.res$logFC>logfc.cutoff & de.res$auc>auc.cutoff & de.res$pct_in>front.cl.prop & de.res$pct_out<bg.cl.prop & de.res$auc>0.5)
    if(length(fil.idx)<1){
      print("No DE")
      to.check=T
      updated.cl <- get.merged.cl.idx(row.idx=count, cl.before.merge=cl.idx, cell.cl.before.merge=cell.cl.vec)
      cl.idx = updated.cl$cl.after.merge
      cell.cl.vec = updated.cl$cell.cl.after.merge
      print(cl.idx)
      next
    }
    de.fil <- de.res[fil.idx,]
    g <- setdiff(de.fil$feature, all.confound.genes)
    if(length(g)<de.num.cutoff){
      print(paste("No more than", de.num.cutoff, "non-confounded DE genes"))
      to.check=T
      updated.cl <- get.merged.cl.idx(row.idx=count, cl.before.merge=cl.idx, cell.cl.before.merge=cell.cl.vec)
      cl.idx = updated.cl$cl.after.merge
      cell.cl.vec = updated.cl$cell.cl.after.merge
      print(cl.idx)
      next
    }
    de.fil <- de.fil[which(de.fil$feature%in%g),]
    cm.cl.num <- length(unique(de.fil$group))
    if(cm.cl.num<2) {
      print("Only one cluster has DE")
      to.check=T
      updated.cl <- get.merged.cl.idx(row.idx=count, cl.before.merge=cl.idx, cell.cl.before.merge=cell.cl.vec)
      cl.idx = updated.cl$cl.after.merge
      cell.cl.vec = updated.cl$cell.cl.after.merge
      print(cl.idx)
      next
    }else{
      print("No further merging")
      to.check=F    
    }
  }
  cl.num <- length(unique(cell.cl.vec))
  print(paste(cl.num, "Clusters remain"))
  
  cl.orig <- seu.obj@meta.data[[feature.to.merge]]
  cl.1 = cutree(hc, k=cl.num)
  merged.cl <- rep(NA, ncol(seu.obj))
  for(new.idx in unique(cl.1)){
    old.idx <- sub("Cluster_", "", names(which(cl.1==new.idx)))
    merged.cl[which(cl.orig%in%old.idx)] <- new.idx
  }
  seu.obj@meta.data[[feature.after.merge]] <- merged.cl
  return(seu.obj)
}

scatterplot3d_rotating <- function(coord, pt_colors = "black", pt_cex = 1, box = F, axis = F, tick.marks = F, label.tick.marks = F, grid = F,
                                   height = 1000, width = 1000,
                                   num_frames = 72, delay = 7, output_file = paste0(getwd(), "/plot.gif")){
  require(scatterplot3d)
  coord <- scale(coord)
  rotate_coord <- function(coord, angle, dims_rotate = c(1,2)){
    coord_rotated <- coord
    coord_rotated[,dims_rotate[1]] <- cos(angle) * coord[,1] - sin(angle) * coord[,2]
    coord_rotated[,dims_rotate[2]] <- cos(angle) * coord[,2] + sin(angle) * coord[,1]
    return(coord_rotated)
  }
  
  temppng <- paste0(sapply(seq(num_frames), function(x) tempfile()), ".png")
  for(i in 1:num_frames){
    png(temppng[i], height=height, width=width)
    par(cex=1)
    scatterplot3d(rotate_coord(coord, pi*2/num_frames * i), box=box, axis=axis, tick.marks = tick.marks, label.tick.marks = label.tick.marks, grid=grid, pch=16, mar=c(0,0,0,0),
                  xlim = c(-max(sqrt(rowSums(coord[,1:2]^2))), max(sqrt(rowSums(coord[,1:2]^2)))), ylim=c(-max(sqrt(rowSums(coord[,1:2]^2))), max(sqrt(rowSums(coord[,1:2]^2)))),
                  color = pt_colors, cex.symbols = pt_cex)
    dev.off()
  }
  system(paste0("~/anaconda3/bin/convert -delay ", delay, " -loop 0 ", paste(temppng, collapse = " "), " ", output_file))
  file.remove(temppng)
  return(NULL)
}

mergeCluster <- function(hc, seu.obj, assay.type=NULL, data.type="data", cluster.before.merge="SCT_snn_res.2", cluster.after.merge="Merged_cluster_idx",front.cl.prop=50, bg.cl.prop=10, padj.cutoff=0.05, logfc.cutoff=0.5, de.num.cutoff=5){
  library(presto)
  hc.merge <- hc$merge
  
  get.merged.cl.idx <- function(row.idx, cl.before.merge, cell.cl.before.merge){
    input.pair=hc.merge[row.idx,]
    if(input.pair[1]<0 & input.pair[2]<0){
      old.cl <- abs(input.pair)-1
    }else if(input.pair[1]>0 & input.pair[2]>0){
      old.cl <- cl.before.merge[input.pair]
    }else{
      old.cl <- c(abs(min(input.pair))-1, cl.before.merge[max(input.pair)])
    }
    new.cl <- min(old.cl)
    cell.cl.before.merge[which(cell.cl.before.merge%in%old.cl)] <- new.cl
    cl.before.merge[row.idx] <- new.cl
    cell.cl.after.merge <- cell.cl.before.merge
    cl.after.merge <- cl.before.merge
    res <- list("cl.after.merge"=cl.after.merge, "cell.cl.after.merge"=cell.cl.after.merge)
    return(res)
  }
  
  cl.idx <- rep(NA, nrow(hc.merge))
  cell.cl.vec <- seu.obj@meta.data[[cluster.before.merge]]
  to.check=T
  count=0
  while(to.check & count<nrow(hc.merge)){
    count <- count+1
    print(paste("Iteration", count))
    cl.vec <- hc.merge[count,]
    cl.1 <- ifelse(cl.vec[1]<0, abs(cl.vec[1])-1, cl.idx[cl.vec[1]])
    cl.2 <- ifelse(cl.vec[2]<0, abs(cl.vec[2])-1, cl.idx[cl.vec[2]])
    cell.idx <- c(which(cell.cl.vec==cl.1), which(cell.cl.vec==cl.2))
    expr <- slot(seu.obj@assays[[assay.type]], data.type)[,cell.idx]
    group.vec <- cell.cl.vec[cell.idx]
    de.res <- wilcoxauc(X=expr, y=group.vec)
    fil.idx <- which(de.res$logFC>logfc.cutoff & de.res$padj<padj.cutoff & de.res$pct_in>front.cl.prop & de.res$pct_out<bg.cl.prop & de.res$auc>0.5)
    if(length(fil.idx)<1){
      print("No DE")
      to.check=T
      updated.cl <- get.merged.cl.idx(row.idx=count, cl.before.merge=cl.idx, cell.cl.before.merge=cell.cl.vec)
      cl.idx = updated.cl$cl.after.merge
      cell.cl.vec = updated.cl$cell.cl.after.merge
      print(cl.idx)
      next
    }
    de.fil <- de.res[fil.idx,]
    g <- setdiff(de.fil$feature, all.confound.genes)
    if(length(g)<de.num.cutoff){
      print(paste("No more than", de.num.cutoff, "non-confounded DE genes"))
      to.check=T
      updated.cl <- get.merged.cl.idx(row.idx=count, cl.before.merge=cl.idx, cell.cl.before.merge=cell.cl.vec)
      cl.idx = updated.cl$cl.after.merge
      cell.cl.vec = updated.cl$cell.cl.after.merge
      print(cl.idx)
      next
    }
    de.fil <- de.fil[which(de.fil$feature%in%g),]
    cm.cl.num <- length(unique(de.fil$group))
    if(cm.cl.num<2) {
      print("Only one cluster has DE")
      to.check=T
      updated.cl <- get.merged.cl.idx(row.idx=count, cl.before.merge=cl.idx, cell.cl.before.merge=cell.cl.vec)
      cl.idx = updated.cl$cl.after.merge
      cell.cl.vec = updated.cl$cell.cl.after.merge
      print(cl.idx)
      next
    }else{
      print("No further merging")
      to.check=F    
    }
  }
  cl.num <- length(unique(cell.cl.vec))
  print(paste(cl.num, "Clusters remain"))
  
  cl.orig <- seu.obj@meta.data[[cluster.before.merge]]
  cl.1 = cutree(hc, k=cl.num)
  merged.cl <- rep(NA, ncol(meso))
  for(new.idx in unique(cl.1)){
    old.idx <- sub("Cluster_", "", names(which(cl.1==new.idx)))
    merged.cl[which(cl.orig%in%old.idx)] <- new.idx
  }
  seu.obj@meta.data[[cluster.after.merge]] <- merged.cl
  return(seu.obj)
}

get_umap_models <- function(
  object,
  dims = 1:20,
  reduction = 'pca',
  n.neighbors = 30L,
  n.components = 2L,
  metric = 'cosine',
  n.epochs = NULL,
  learning.rate = 1,
  min.dist = 0.3,
  spread = 1,
  set.op.mix.ratio = 1,
  local.connectivity = 1L,
  repulsion.strength = 1,
  negative.sample.rate = 5L,
  a = NULL,
  b = NULL,
  uwot.sgd = FALSE,
  n_threads = 1,
  verbose = TRUE
) {
  require(uwot)
  data.use <- Embeddings(object[[reduction]])[, dims]
  umap(
    X = data.use,
    n_threads = n_threads,
    n_neighbors = as.integer(x = n.neighbors),
    n_components = as.integer(x = n.components),
    metric = metric,
    n_epochs = n.epochs,
    learning_rate = learning.rate,
    min_dist = min.dist,
    spread = spread,
    set_op_mix_ratio = set.op.mix.ratio,
    local_connectivity = local.connectivity,
    repulsion_strength = repulsion.strength,
    negative_sample_rate = negative.sample.rate,
    a = a,
    b = b,
    fast_sgd = uwot.sgd,
    ret_model = TRUE,
    verbose = verbose
  )
}

getClusterMarkers <- function(seu.obj, assay.type="RNA", data.type="data", feature.to.calc=NULL, file.name=NULL, genes.to.remove=NULL, logFC.cutoff=0.25, padj.cutoff=0.05, pct_in.cutoff=25, n=50){
  library(presto)
  X=slot(seu.obj@assays[[assay.type]], data.type)
  cm <- wilcoxauc(X=X, y=seu.obj@meta.data[[feature.to.calc]])
  cm <- cm[which(cm$logFC>logFC.cutoff & cm$padj<padj.cutoff & cm$pct_in>pct_in.cutoff),]
  top.cm <- cm %>% group_by(group) %>% top_n(n, wt=logFC)
  top.cm.genes <- unique(top.cm$feature)
  top.cm.genes <- setdiff(top.cm.genes, genes.to.remove)
  cm.res <- list("all.sig.cm"=cm, "top.cm"=top.cm, "top.cm.genes.non.confound"=top.cm.genes)
  saveRDS(cm.res, file=file.name)
  return(cm.res)
}

getExprByPt <- function(pt.vec, expr.mat, cell.num.per.bin=50){
  start.idx <- seq(from=1, to=length(pt.vec), by=cell.num.per.bin)
  end.idx <- c(start.idx[-1]-1,length(pt.vec))
  idx <- cbind(start.idx, end.idx)
  pt.rank <- rank(pt.vec)
  expr.by.pt.g <- sapply(seq(nrow(idx)), function(i){
    idx.start <- idx[i,1]
    idx.end <- idx[i,2]
    idx.cell <- which(pt.rank>=idx.start & pt.rank<=idx.end)
    apply(expr.mat[,idx.cell], 1, mean)
  })
  return(expr.by.pt.g)
}
