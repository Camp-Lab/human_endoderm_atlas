library(Matrix)
library(Seurat)
build_knn_graph <- function(data, k = 20, dist_type = c("euclidean", "pearson", "raw"), use_seurat_snn = TRUE,
                            mutual = TRUE, jaccard_weighted = TRUE, jaccard_prune = 1/15,
                            return_igraph = FALSE, threads_jaccard = 1, verbose = TRUE){
  dist_type <- match.arg(dist_type)
  
  nn <- matrix(0, nrow = ncol(data), ncol = k)
  if (dist_type == "raw"){
    diag(data) <- NA
    nn <- t(apply(data, 2, function(x) order(x)[1:k]))
  } else if (dist_type %in% c("euclidean", "pearson")){
    if (dist_type == "pearson")
      data <- scale(data)
    nn <- RANN::nn2(t(data), t(data), k = k+1)$nn.idx[,-1]
  }
  if (verbose)
    cat("found nearest neighbors.\n")
  
  i <- rep(1:nrow(nn), each = ncol(nn))
  j <- as.integer(t(nn))
  adj <- Matrix::sparseMatrix(i, j, dims = c(ncol(data), ncol(data)))
  
  if(use_seurat_snn){
    library(Seurat)
    if (verbose)
      cat("revoke Seurat to compute SNN.\n")
    
    adj <- .Call('_Seurat_ComputeSNN', PACKAGE = 'Seurat', nn, jaccard_prune)
    if (verbose)
      cat("done.\n")
    
  } else{
    if (mutual){
      adj <- adj * t(adj)
    } else{
      adj <- adj + t(adj)
      adj <- Matrix::sparseMatrix(i = adj@i, p = adj@p, dims = adj@Dim, index1 = FALSE) + 0
    }
    if (verbose)
      cat("got adjacent matrix.\n")
    
    if (jaccard_prune > 0 | jaccard_weighted){
      if (verbose)
        cat("start calculating Jaccard indices...\n")
      nonzero_idx <- summary(adj)[,1:2]
      
      library(doParallel)
      registerDoParallel(threads_jaccard)
      jaccard <- foreach(i = 1:nrow(nonzero_idx), .combine = c) %dopar%{
        x <- as.integer(nonzero_idx[i,])
        n1 <- nonzero_idx[nonzero_idx[,1] == x[1],2]
        n2 <- nonzero_idx[nonzero_idx[,1] == x[2],2]
        return(length(intersect(n1, n2)) / length(union(n1, n2)))
      }
      stopImplicitCluster()
      remained <- which(jaccard >= jaccard_prune)
      adj <- Matrix::sparseMatrix(i = nonzero_idx[remained,1], j = nonzero_idx[remained,2], x = jaccard[remained], dims = c(nrow(nn), nrow(nn))) + 0
      if (verbose)
        cat("done pruning.\n")
    }
  }
  adj@Dimnames <- list(colnames(data), colnames(data))
  
  if (verbose)
    cat("done. returning result...\n")
  if (return_igraph) return(igraph::graph_from_adjacency_matrix(adj, mode = "undirected"))
  return(adj)
}

construct_pseudocells <- function(knn_adj, ratio, init_dist = 1, max_dist = 2, min_pooled_cells = 2, min_seed_num = 1){
  rownames(knn_adj) <- 1:nrow(knn_adj)
  colnames(knn_adj) <- 1:ncol(knn_adj)
  capitals <- which(runif(nrow(knn_adj)) <= ratio)
  if (length(capitals) == 0) capitals <- sample(1:nrow(knn_adj), min_seed_num)
  
  graph <- igraph::graph_from_adjacency_matrix(knn_adj, mode = "undirected")
  dist_to_capitals <- igraph::distances(graph, v = as.character(1:ncol(knn_adj)), to = as.character(capitals))
  
  selected <- dist_to_capitals <= init_dist
  selected <- t(sapply(1:nrow(selected), function(i){
    sel <- selected[i,]
    if (sum(sel) == 1) return(sel)
    if (sum(sel) > 1) return(1:length(sel) == sample(which(sel),1))
    if (sum(sel) == 0 & min(dist_to_capitals[i,]) <= max_dist) return(1:length(sel) == sample(which(dist_to_capitals[i,] <= max_dist),1))
    return(sel)
  }))
  if (ncol(dist_to_capitals) == 1) selected <- t(selected)
  
  sel_df <- data.frame(idx_raw = 1:nrow(knn_adj), idx_pseudo = -1, stringsAsFactors=F)
  sel_df$idx_pseudo[apply(selected, 1, sum) == 1] <- apply(data.frame(selected[apply(selected, 1, sum) == 1,]), 1, which)
  sel_df$idx_pseudo[sel_df$idx_pseudo %in% unique(sel_df$idx_pseudo)[sapply(unique(sel_df$idx_pseudo), function(x) sum(sel_df$idx_pseudo == x)) < min_pooled_cells]] <- -1
  return(sel_df)
}

cluster_sim_spectrum.default <- function(object, # expression matrix
                                         dr = NULL, dr_input = NULL, num_pcs_compute = 50, num_pcs_use = 20, # for dimension reduction
                                         labels, redo_pca = FALSE, k = 20, ..., # how to separate samples and whether or not to do DR separately
                                         cluster_resolution = 1, spectrum_type = c("corr_ztransform","corr_kernel","corr_raw","nnet","lasso"), # clustering and types of spectrum
                                         corr_method = c("spearman","pearson"), lambda = 1, threads = 1,  # spectrum related parameters
                                         train_on = c("raw","pseudo","rand"), downsample_ratio = 1/10, k_pseudo = 10, logscale_likelihood = F, # parameters of likelihood spectrum
                                         merge_spectrums = FALSE, merge_height_prop = 1/10, spectrum_dist_type = c("pearson", "euclidean"), spectrum_cl_method = "complete", # parameters of spectrum merging
                                         return_css_only = T, verbose = T){
  spectrum_type <- match.arg(spectrum_type)
  corr_method <- match.arg(corr_method)
  train_on <- match.arg(train_on)
  spectrum_dist_type <- match.arg(spectrum_dist_type)
  data <- object
  if (is.null(dr_input))
    dr_input <- data
  
  if (is.null(dr)){
    if (verbose)
      cat("No dimension reduction is provided. Start to do truncated PCA...\n")
    
    t_pca <- irlba::irlba(t(dr_input), nv = num_pcs_compute)
    dr <- t_pca$u %*% diag(t_pca$d)
    dr <- dr[,1:num_pcs_use]
    rownames(dr) <- colnames(data)
    
    if (verbose)
      cat("PCA finished.\n")
  }
  
  if (verbose)
    cat("Start to do clustering for each sample...\n")
  labels <- as.factor(labels)
  cl <- lapply(levels(labels), function(x){
    idx <- which(labels == x)
    dr_x <- dr[idx,]
    if (redo_pca){
      if (verbose)
        cat(paste0(">> Redoing truncated PCA on sample ", x, "...\n"))
      t_pca <- irlba::irlba(t(dr_input[,idx]), nv = num_pcs_compute)
      dr_x <- t_pca$u %*% diag(t_pca$d)
      dr_x <- dr_x[,1:num_pcs_use]
      rownames(dr_x) <- colnames(data)[idx]
    }
    knn <- build_knn_graph(t(dr_x), k = k, ..., verbose = verbose)
    rownames(knn) <- colnames(data)[idx]
    colnames(knn) <- colnames(data)[idx]
    cl <- Seurat::FindClusters(Seurat::as.Graph(knn), resolution = cluster_resolution, verbose = verbose)[,1]
    if (verbose)
      cat(paste0(">> Done clustering of sample ", x, ".\n"))
    return(cl)
  })
  names(cl) <- levels(labels)
  if (verbose)
    cat("Finished clustering.\n")
  
  if (spectrum_type %in% c("corr_ztransform","corr_kernel","corr_raw")){
    cl_profiles <- lapply(levels(labels), function(x){
      idx <- which(labels == x)
      profiles <- sapply(levels(as.factor(cl[[x]])), function(cl_x)
        apply(data[,idx[as.factor(cl[[x]])==cl_x]], 1, mean))
      return(profiles)
    })
    names(cl_profiles) <- levels(labels)
    if (verbose)
      cat("Obtained average profiles of clusters.\n")
    
    if (verbose)
      cat("Start to calculate standardized similarities to clusters...\n")
    sim2profiles <- lapply(cl_profiles, function(profiles){
      cor(as.matrix(data), profiles, method=corr_method)
    })
    
    if (spectrum_type == "corr_ztransform"){
      if (verbose)
        cat("Doing z-transformation...\n")
      sim2profiles <- lapply(sim2profiles, function(sim)
        t(scale(t(sim))))
    } else if (spectrum_type == "corr_kernel"){
      if (verbose)
        cat("Doing kernel transformation...\n")
      sim2profiles <- lapply(sim2profiles, function(sim)
        t(apply(exp(sim * lambda) * exp(-lambda), 1, function(x) x/sum(x))))
    }
  } else if (spectrum_type %in% c("nnet","lasso")){
    if (verbose)
      cat("Start to build multinomial logistic regression models...\n")
    
    models <- lapply(levels(labels), function(x){
      idx_x <- which(labels == x)
      cl_x <- cl[[x]]
      train_x <- t(data[,idx_x])
      train_y <- cl_x
      if (train_on == "rand"){
        sel_idx <- lapply(levels(cl_x), function(cl_this){
          cl_idx <- idx_x[cl_x == cl_this]
          sel_idx <- cl_idx[which(runif(length(cl_idx)) <= downsample_ratio)]
          if (length(sel_idx) == 0) sel_idx <- sample(cl_idx, 1)
          return(sel_idx)
        })
        train_x <- as.matrix(do.call(rbind, lapply(sel_idx, function(idx) t(data[,idx]))))
        train_y <- as.factor(rep(levels(cl_x), sapply(sel_idx, length)))
        if (verbose)
          cat(paste0(">> Randomly selected cells for model training: sample ", x, ".\n"))
        
      } else if (train_on == "pseudo"){
        pseudo_idx <- lapply(levels(cl_x), function(cl_this){
          cl_idx <- idx_x[cl_x == cl_this]
          knn_cl <- build_knn_graph(t(dr[cl_idx,]), k = k_pseudo,
                                    use_seurat_snn = F, mutual = F, jaccard_weighted = F, jaccard_prune = 0,
                                    ..., verbose = F)
          pseudo_idx <- construct_pseudocells(knn_cl, ratio = downsample_ratio, min_seed_num = 2)
          return(pseudo_idx)
        })
        names(pseudo_idx) <- levels(cl_x)
        train_x <- do.call(rbind, lapply(levels(cl_x), function(cl_this){
          cl_idx <- idx_x[cl_x == cl_this]
          idx <- pseudo_idx[[cl_this]]
          data_pseudo <- sapply(sort(unique(idx[idx[,2]!=-1,2])), function(ix) apply(data[,cl_idx[idx[,2] == ix]], 1, mean) )
          return(t(data_pseudo))
        }))
        train_y <- rep(levels(cl_x), sapply(pseudo_idx, function(idx) length(unique(idx[idx[,2]!=-1,2]))))
        if (verbose)
          cat(paste0(">> Constructed pseudocells for training: sample ", x, ".\n"))
      }
      idx <- which(train_y %in% unique(train_y)[sapply(unique(train_y), function(x) sum(train_y == x)) > 2])
      if (verbose & length(idx) != length(train_y))
        cat(paste0(">> Dropped clusters ",unique(train_y[-idx]), " in sample ", x, " due to too few training data.\n"))
      train_x <- train_x[idx,]
      train_y <- train_y[idx]
      
      # train model: nnet::multinom
      m <- NULL
      if (spectrum_type == "nnet"){
        require(nnet)
        train_dat <- data.frame(y = train_y, train_x)
        m <- nnet::multinom(y ~ ., data = train_dat, MaxNWts = 1E8, trace = verbose)
      } else if (spectrum_type == "lasso"){
        require(glmnet)
        if (threads > 1){
          library(doMC)
          registerDoMC(threads)
        }
        m <- glmnet::cv.glmnet(x = train_x, y = train_y, family = "multinomial", alpha = 1)
      }
      if (verbose)
        cat(paste0(">> Model trained for sample ", x, ".\n"))
      return(m)
    })
    names(models) <- levels(labels)
    if (verbose)
      cat("Models trained, start to produce likelihood spectrum...\n")
    
    sim2profiles <- lapply(models, function(m){
      if (spectrum_type == "lasso"){
        pred <- predict(m, t(data), type = "response")[,,1]
      } else if (spectrum_type == "nnet"){
        pred <- predict(m, data.frame(t(data)), "probs")
      }
      return(pred)
    })
    
    if (logscale_likelihood)
      sim2profiles <- lapply(sim2profiles, function(pred){
        pred <- t(scale(t(log(pred))))
      })
    
    if (verbose)
      cat("Done likelihood estimation.\n")
  }
  
  sim2profiles <- do.call(cbind, sim2profiles)
  sim2profiles_raw <- sim2profiles
  if (merge_spectrums){
    if (verbose)
      cat("Start to merge similar spectrums...\n")
    
    dist_css <- NULL
    if (spectrum_dist_type == "pearson"){
      dist_css <- as.dist(1 - cor(sim2profiles))
    } else if (spectrum_dist_type == "euclidean"){
      dist_css <- dist(t(sim2profiles))
    }
    cl <- hclust(dist_css, method = spectrum_cl_method)
    cl_spectrums <- cutree(cl, h = max(cl$height) * merge_height_prop)
    sim2profiles <- sapply(sort(unique(cl_spectrums)), function(x){
      if(sum(cl_spectrums==x)==1)
        return(sim2profiles[,cl_spectrums==x])
      return(apply(sim2profiles[,cl_spectrums==x], 1, mean))
    })
    #sim2profiles <- t(scale(t(sim2profiles)))
  }
  rownames(sim2profiles) <- colnames(data)
  
  if (verbose)
    cat("Done.\n")
  
  if (return_css_only)
    return(sim2profiles)
  
  model <- list(spectrum_type = spectrum_type)
  if (spectrum_type %in% c("corr_ztransform","corr_kernel")){
    model$profiles <- cl_profiles
    model$args <- c(corr_method = corr_method, lambda = lambda)
  } else if (spectrum_type %in% c("nnet","lasso")){
    model$models <- models
    model$args <- c(logscale_likelihood = logscale_likelihood)
  }
  model$merged_spectrum <- seq(1, ncol(sim2profiles_raw))
  if (merge_spectrums)
    model$merged_spectrum <- cl_spectrums
  res <- list(model = model, sim2profiles = sim2profiles)
  return(res)
}

cluster_sim_spectrum.Seurat <- function(object, var_genes = NULL, use_scale = F, use_dr = "pca", dims_use = 1:20,
                                        label_tag, redo_pca = FALSE, redo_pca_with_data = FALSE, k = 20, ...,
                                        cluster_resolution = 1, spectrum_type = c("corr_ztransform","corr_kernel","corr_raw","nnet","lasso"),
                                        corr_method = c("spearman","pearson"), lambda = 1, threads = 1,
                                        train_on = c("raw","pseudo","rand"), downsample_ratio = 1/10, k_pseudo = 10, logscale_likelihood = F,
                                        merge_spectrums = FALSE, merge_height_prop = 1/10, spectrum_dist_type = c("pearson", "euclidean"), spectrum_cl_method = "complete",
                                        return_seuratObj = T, verbose = T){
  spectrum_type <- match.arg(spectrum_type)
  corr_method <- match.arg(corr_method)
  train_on <- match.arg(train_on)
  spectrum_dist_type <- match.arg(spectrum_dist_type)
  if (is.null(var_genes))
    var_genes <- Seurat::VariableFeatures(object)
  
  data <- object[[object@active.assay]]@data[var_genes,]
  if (use_scale)
    data <- object[[object@active.assay]]@scale.data[var_genes,]
  
  dr_input <- NULL
  if (redo_pca){
    if (redo_pca_with_data){
      dr_input <- object[[object@active.assay]]@data[var_genes,]
    } else{
      dr_input <- object[[object@active.assay]]@scale.data[var_genes,]
    }
  }
  
  dr <- object@reductions[[use_dr]]@cell.embeddings[,dims_use]
  labels <- object@meta.data[,label_tag]
  css <- cluster_sim_spectrum.default(object = data, dr = dr, dr_input = dr_input, labels = labels, redo_pca = redo_pca, k = k, ..., cluster_resolution = cluster_resolution,
                                      spectrum_type = spectrum_type, corr_method = corr_method, lambda = lambda, threads = threads,
                                      train_on = train_on, downsample_ratio = downsample_ratio, k_pseudo = k_pseudo, logscale_likelihood = logscale_likelihood,
                                      merge_spectrums = merge_spectrums, merge_height_prop = merge_height_prop, spectrum_dist_type = spectrum_dist_type, spectrum_cl_method = spectrum_cl_method,
                                      return_css_only = return_seuratObj, verbose = verbose)
  if (return_seuratObj){
    object[["css"]] <- CreateDimReducObject(embeddings = css, key = "CSS_", assay = DefaultAssay(object))
    return(object)
  } else{
    return(css)
  }
}

cluster_sim_spectrum <- function(object, ...) {
  UseMethod(generic = 'cluster_sim_spectrum', object = object)
}
