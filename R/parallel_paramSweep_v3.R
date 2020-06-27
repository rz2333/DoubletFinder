.Sample <- function(real.cells, n_doublets, replace) {
  cell <- list()
  cell[[1]] <- sample(real.cells, n_doublets, replace = TRUE)
  cell[[2]] <- sample(real.cells, n_doublets, replace = TRUE)
  return cell 
  }

parallel_paramSweep_v3 <- function(n, n.real.cells, real.cells, pK, pN, data,
                                   orig.commands, PCs, sct, verbose, seed) {
  sweep.res.list <- list()
  list.ind <- 0

  ## Make merged real-artifical data
  if (verbose) {
    print(paste("Creating artificial doublets for pN = ", pN[n] * 100, "%", sep = ""))
  }
  n_doublets <- round(n.real.cells / (1 - pN[n]) - n.real.cells)
#   set.seed(seed)
#   real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
#   set.seed(seed)
#   real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  withr::with_seed(seed = seed, 
                   real.cells <- .Sample(real.cells, n_doublets, replace=TRUE))
  real.cells1 <- real.cells[[1]]
  real.cells2 <- real.cells[[2]]
  
  doublets <- (data[, real.cells1] + data[, real.cells2]) / 2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)

  ## Pre-process Seurat object
  if (sct == FALSE) {
    if (verbose) {
      print("Creating Seurat object...")
    }
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    if (verbose) {
      print("Normalizing Seurat object...")
    }
    seu_wdoublets <- NormalizeData(seu_wdoublets,
      normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
      scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
      margin = orig.commands$NormalizeData.RNA@params$margin,
      verbose = verbose
    )

    if (verbose) {
      print("Finding variable genes...")
    }
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
      selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
      loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
      clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
      mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
      dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
      num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
      binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
      nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
      mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
      dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff,
      verbose = verbose
    )

    if (verbose) {
      print("Scaling data...")
    }
    seu_wdoublets <- ScaleData(seu_wdoublets,
      features = orig.commands$ScaleData.RNA$features,
      model.use = orig.commands$ScaleData.RNA$model.use,
      do.scale = orig.commands$ScaleData.RNA$do.scale,
      do.center = orig.commands$ScaleData.RNA$do.center,
      scale.max = orig.commands$ScaleData.RNA$scale.max,
      block.size = orig.commands$ScaleData.RNA$block.size,
      min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block,
      verbose = verbose
    )

    if (verbose) {
      print("Running PCA...")
    }
    seu_wdoublets <- RunPCA(seu_wdoublets,
      features = orig.commands$ScaleData.RNA$features,
      npcs = length(PCs),
      rev.pca = orig.commands$RunPCA.RNA$rev.pca,
      weight.by.var = orig.commands$RunPCA.RNA$weight.by.var, seed.use = seed,
      verbose = FALSE
    )
  }

  if (sct == TRUE) {
    require(sctransform)
    if (verbose) {
      print("Creating Seurat object...")
    }
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)

    if (verbose) {
      print("Running SCTransform...")
    }
    seu_wdoublets <- SCTransform(seu_wdoublets)

    if (verbose) {
      print("Running PCA...")
    }
    seu_wdoublets <- RunPCA(seu_wdoublets,
      npcs = length(PCs), seed.use = seed,
      verbose = verbose
    )
  }

  ## Compute PC distance matrix
  if (verbose) {
    print("Calculating PC distance matrix...")
  }
  nCells <- nrow(seu_wdoublets@meta.data)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]
  rm(seu_wdoublets)
  gc()
  dist.mat <- fields::rdist(pca.coord)[, 1:n.real.cells]

  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  if (verbose) {
    print("Defining neighborhoods...")
  }

  for (i in 1:n.real.cells) {
    dist.mat[, i] <- order(dist.mat[, i])
  }

  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK)) + 5
  dist.mat <- dist.mat[1:ind, ]

  ## Compute pANN across pK sweep

  if (verbose) {
    print("Computing pANN across all pK...")
  }
  for (k in 1:length(pK)) {
    if (verbose) {
      print(paste("pK = ", pK[k], "...", sep = ""))
    }
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1

    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1), i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells)) / pk.temp
    }

    sweep.res.list[[list.ind]] <- pANN
  }

  return(sweep.res.list)
}
