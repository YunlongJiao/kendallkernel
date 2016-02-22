
# compute kernel matrix ---------------------------------------------------

kernmat <- function(dat, dimsize, kf = "kendall_top", procf = NULL, acc.dup = FALSE)
{
  # @param dat a matrix with observations in rows and item ranks in columns
  # @param dimsize a vector indicating multivariate positions, e.g. c(4,4) indicates 1-4 and 5-8 are bi-variate top-k rankings of totally 4 items
  # @param kf psd kernel function, examples are found in package `kernlab`
  # @param procf preprocess function applied to dat before evaluation
  # @param acc.dup accelerate computation by unique-fy duplicated instances
  # @return a symmetric square kernel matrix of size nsample x nsample of class `kernelMatrix`
  
  stopifnot(is.matrix(dat))
  if (is.character(kf)) kf <- get(kf, mode = "function")
  class(kf) <- "kernel"
  
  dm <- c(0, cumsum(dimsize))
  dm <- cut(1:max(dm), dm, include.lowest = FALSE, right = TRUE)
  dm <- lapply(levels(dm), function(u) which(dm == u))
  
  if (!is.null(procf)) {
    message("pre-process dat ...")
    if (is.character(procf)) procf <- get(procf, mode = "function")
    dat <- procf(dat = dat, dm = dm)
  }
  
  if (nrow(dat) > 1e4) {
    acc.dup <- TRUE
  }
  if (acc.dup) {
    message("accelerate kernel matrix computation for duplicated instances ...")
    datfull <- dat
    dat <- unique(dat)
    ind <- rep(0, length = nrow(datfull)) # index alignment
    invisible(lapply(seq(nrow(dat)), function(i){
      if (i %% 100 == 0) message("aligning ", i, "-th sample")
      ind[apply(datfull, 1, identical, dat[i, , drop = TRUE])] <<- i
    }))
  }
  
  message("computing kernel matrix ...")
  if (length(dimsize) == 1) {
    mat <- kernelMatrix(kernel = kf, x = dat)
  } else {
    mat <- unlist(lapply(dm, function(idx) kernelMatrix(kernel = kf, x = dat[ ,idx])))
  }
  message("done!")
  
  if (length(dm) == 1) {
    dim(mat) <- c(nrow(dat), nrow(dat))
  } else {
    dim(mat) <- c(nrow(dat), nrow(dat), length(dm))
    mat <- apply(mat, c(1, 2), mean)
  }
  if (acc.dup) {
    mat <- mat[ind, ind]
  }
  as.kernelMatrix(mat)
}

procf_Rankcluster <- function(dat, dm = list(1:ncol(dat)))
{
  # example procf where typically two arguments are taken
  # @param dat multivariate top-k datasets from package `Rankcluster` in ranking representation where small values indicate higher preference and 0 indicate missing values
  # @param dm a list of indices indicating groups of columns as multivariate
  # @return processed datasets so that larger values indicate higher preference (plus 0 replaced with NA)
  
  res <- apply(dat, 1, function(v){
    unlist(lapply(dm, function(idx){
      ind <- which(v[idx] == 0)
      vb <- length(idx) + 1 - v[idx]
      vb[ind] <- NA
      vb
    }))
  })
  t(res)
}

# compute squared-distance matrix -------------------------------------------------

distmat <- function(kmatrix)
{
  # @param kmatrix a symmetric square kernel matrix of size nsample x nsample
  # @return a squared-distance matrix of size nsample x nsample where entries are $d_{ij} = k_{ii} + k_{jj} - 2 k_{ij}$
  
  if (is.data.frame(kmatrix)) {
    kmatrix <- as.matrix(kmatrix)
  }
  dmatrix <- sweep(sweep(-2*kmatrix, 1, diag(kmatrix), "+"), 2, diag(kmatrix), "+")
  dmatrix
}

# kkmeans_multi -----------------------------------------------------------

kkmeans_multi <- function(kmatrix, klist, nrepeats, dmatrix = NULL, seed = NULL, mc.cores = NULL, verbose = FALSE, simple.output = FALSE)
{
  # @param kmatrix a kernel matrix of size nsample x nsample of class `kernelMatrix`
  # @param klist a vector of nclust
  # @param nrepeats number of random restarts for kkmeans
  # @param mc.cores Enable multicore runs. Not yet implemented.
  # @param simple.output Return simply the best repeat for each k and the best k wrt average silhouette score
  # @return a 3-layer nested list of klist -> seq(repeats) -> original result returned by kkmeans
  
  if (!is.null(seed)) set.seed(seed)
  
  kmatrix <- as.kernelMatrix(as.matrix(kmatrix))
  if (is.null(dmatrix)){
    dmatrix <- distmat(kmatrix)
  }
  
  if (simple.output) {
    # initialize for nclust
    sil.max <- -Inf
    cls <- NULL
    
    # choose best k wrt max ave sil
    for (k in klist) {
      # initialize for repeats
      dss.min <- Inf
      cl.min <- NULL
      
      # choose best repeat wrt min dss
      for (i in seq(nrepeats)) {
        if (verbose) message("k = ", k, "\t nrepeats = ", i, "\n")
        cl <- try(kkmeans(x = kmatrix, centers = k), silent = TRUE)
        if (!inherits(cl, "try-error")) {
          # distsum
          dss <- distsum(clustering = cl@.Data, kmatrix = kmatrix)
          if (dss < dss.min) {
            dss.min <- dss
            cl.min <- cl
          }
        }
      }
      
      if (!is.null(cl.min)) {
        # silhouette
        sil <- cluster::silhouette(x = cl.min@.Data, dmatrix = dmatrix)
        sil <- mean(sil[ ,"sil_width"])
        if (sil > sil.max) {
          sil.max <- sil
          cls <- cl.min
          attr(cls, "silhouette") <- sil
        }
      }
    }
  } else {
    cls <- lapply(klist, function(k){
      x <- lapply(seq(nrepeats), function(i){
        if (verbose) message("k = ", k, "\t nrepeats = ", i, "\n")
        cl <- try(kkmeans(x = kmatrix, centers = k), silent = TRUE)
        if (!inherits(cl, "try-error")) {
          # distsum
          dss <- distsum(clustering = cl@.Data, kmatrix = kmatrix)
          attr(cl, "distsum") <- dss
          # silhouette
          sil <- cluster::silhouette(x = cl@.Data, dmatrix = dmatrix)
          attr(cl, "silhouette") <- mean(sil[ ,"sil_width"])
        }
        cl
      })
      
      # remove try-error repeats
      idx <- sapply(x, inherits, what = "try-error")
      x[idx] <- NULL
      
      x
    })
    
    names(cls) <- klist
  }
  
  cls
}

distsum <- function(clustering, kmatrix)
{
  # @param clustering a cluster partition of length nsample
  # @param kmatrix a kernel matrix of size nsample x nsample
  # @return sum of squared-distances of each instance to corresp center, i.e., objective function in kkmeans
  
  stopifnot(all(dim(kmatrix) == (nsample <- length(clustering))))
  
  inds <- split(1:nsample, clustering)
  dss <- sapply(inds, function(ind){
    sum(diag(kmatrix[ind, ind, drop=FALSE])) - sum(kmatrix[ind, ind, drop=FALSE])/length(ind)
  })
  
  sum(dss)
}

# kkmeans_stats -----------------------------------------------------------

kkmeans_stats <- function(cls)
{
  # @param cls a NAMED (by nclust) list returned from kkmeans_multi
  # @return a data.frame containing statistics of clustering runs
  
  res <- list()
  for (k in names(cls)) {
    if (length(cls[[k]]) == 0) {
      next
    }
    res[[k]] <- data.frame(nclust = k, 
                           rep = seq(length(cls[[k]])), 
                           silhouette = unlist(lapply(cls[[k]], attr, which = "silhouette")), 
                           distsum = unlist(lapply(cls[[k]], attr, which = "distsum")))
  }
  res <- do.call("rbind", res)
  res$nclust <- ordered(res$nclust)
  
  res
}

kkmeans_stats_best <- function(clstats)
{
  # @param clstats a data.frame returned from kkmeans_stats
  # @return a subset of data.frame corresp to best fit of each nclust wrt min distsum
  
  res <- lapply(split(clstats, clstats$nclust), function(x){
    i <- which.min(x$distsum)
    x[i, ]
  })
  res <- do.call("rbind", res)
  
  res
}


# plot.silhouette_modified ------------------------------------------------
# modified from cluster:::plot.silhouette()

plot.silhouette_modified <-
  function(x, nmax.lab = 40, max.strlen = 5,
           main = NULL, sub = NULL,
           xlab = expression("Silhouette width " * s[i]),
           col = "gray", do.col.sort = length(col) > 1,
           border = 0, cex.names = par("cex.axis"),
           do.n.k = TRUE, do.clus.stat = TRUE, ...)
  {
    if(!is.matrix(x) || ncol(x) != 3)
      stop("No valid silhouette information (#{clusters} =? 1)")
    n <- nrow(x)
    x <- sortSilhouette(x)
    s <- rev(x[, "sil_width"])
    space <- c(0, rev(diff(cli <- x[, "cluster"])))
    space[space != 0] <- 0.5 # gap between clusters
    axisnames <- (n < nmax.lab)
    if(axisnames)
      names <- substring(rev(rownames(x)), 1, max.strlen)
#     if(is.null(main)) {
#       main <- "Silhouette plot"
#       if(!is.null(cll <- attr(x,"call"))) { # drop initial "silhouette":
#         if(!is.na(charmatch("silhouette", deparse(cll[[1]]))))
#           cll[[1]] <- as.name("FF")
#         main <-  paste(main, "of", sub("^FF","", deparse(cll)))
#       }
#     }
    smry <- summary(x)
    k <- length(nj <- smry$clus.sizes) # k clusters
#     if(is.null(sub))
#       sub <- paste("Average silhouette width : ",
#                    round(smry$avg.width, digits = 2))
    if(do.col.sort && (lc <- length(col)) > 1) {
      if(lc == k)# cluster wise coloring
        col <- col[cli]
      else ## unit wise coloring
        if(lc != n)
          col <- rep(col, length = n)
        col <- rev(col) # was rev(col[attr(x, "iOrd")])
    }
    y <- barplot(s, space = space, names = names, xlab = xlab,
                 
                 # xlim = c(min(0, min(s)), 1),
                 xlim = c(min(0, min(s)), max(s)),
                 
                 horiz = TRUE, las = 1, mgp = c(3, 1, 0),
                 col = col, border = border, cex.names = cex.names,
                 axisnames = axisnames, ...)
    title(main = main, sub = sub, adj = 0)
    if(do.n.k) {
      mtext(paste("n =", n),	adj = 0)
      mtext(substitute(k ~~ "clusters" ~~ C[j], list(k=k)), adj= 1)
    }
    if(do.clus.stat) {
      mtext(expression(paste(j," :  ", n[j]," | ", ave[i %in% Cj] ~~ s[i])),
            adj = 1.04, line = -1.2)
      y <- rev(y)
      hasCodes <- !is.null(cx <- attr(x,"codes"))
      for(j in 1:k) {
        j. <- if(hasCodes) cx[j] else j
        yj <- mean(y[cli == j.])
        text(1, yj,
             paste(j.,":  ", nj[j]," | ",
                   format(smry$clus.avg.widths[j], digits = 1, nsmall = 2)),
             xpd = NA, adj = 0.8)
      }
    }
  }

