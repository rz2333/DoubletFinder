paramSweep_v3 <- function(seu, PCs = 1:10, sct = FALSE,
                          verbose = FALSE, num.cores, seed) {
  require(Seurat)
  require(fields)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01, 0.3, by = 0.01))
  pN <- seq(0.05, 0.3, by = 0.05)
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data) / (1 - 0.05) - nrow(seu@meta.data))
  pK.test <- round(pK * min.cells)
  pK <- pK[which(pK.test >= 1)]

  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands

  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    set.seed(seed)
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data),
                                                 10000, replace = FALSE)]
    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)
  }

  if (nrow(seu@meta.data) <= 10000) {
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  ## Iterate through pN, computing pANN vectors at varying pK
  # no_cores <- detectCores()-1
  if (is.null(num.cores)) {
    num.cores <- 1
  }

  if (num.cores > 1) {
    require(parallel)
    cl <- makeCluster(num.cores)

    output2 <- withr::with_seed(
      seed,
      mclapply(as.list(1:length(pN)),
        FUN = parallel_paramSweep_v3,
        n.real.cells,
        real.cells,
        pK,
        pN,
        data,
        orig.commands,
        PCs,
        sct, mc.cores = num.cores,
        verbose = verbose,
        seed = seed,
        mc.set.seed = FALSE
      )
    )
    stopCluster(cl)
  } else {
    output2 <- lapply(as.list(1:length(pN)),
      FUN = parallel_paramSweep_v3,
      n.real.cells,
      real.cells,
      pK,
      pN,
      data,
      orig.commands,
      PCs,
      sct,
      seed = seed,
      verbose = verbose
    )
  }

  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for (i in 1:length(output2)) {
    for (j in 1:length(output2[[i]])) {
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_"))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
}
