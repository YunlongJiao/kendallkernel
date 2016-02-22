
# kmeans ------------------------------------------------------------------

kmeans_validate <- function(type, clres, keylist, klist, nrepeats, ..., 
                            func = mget(paste0("w", type), mode = "function", envir = globalenv()))
{
  # @param type a character vector indicating the type(s) of validation measure
  # @param func validation function to call is by default named paste0("w",type)
  # @param clres a 3-layer list representing nclust, key (type of kmeans), repeat respectively
  # @return validation score data.frame
  
  if (is.null(names(clres)) && max(klist) == length(clres)) {
    names(clres) <- 1:length(clres)
  }
  
  clstats <- list()
  for (k in klist) {
    clstats[[k]] <- data.frame(nclust = k, 
                               rep = rep(1:nrepeats, times = length(keylist)), 
                               method = rep(keylist, each = nrepeats))
    for (t in 1:length(type)) {
      clstats[[k]][[type[t]]] <- sapply(unlist(clres[[as.character(k)]], recursive = FALSE), func[[t]], ...)
    }
  }
  clstats <- as.data.frame(do.call("rbind", clstats), stringsAsFactors = TRUE)
  clstats$nclust <- ordered(clstats$nclust)
  clstats
}

# specific validation functions
wdistsum <- function(cl, weights = NULL, ...)
{
  # @param cl an kcca obj
  # @param weights frequency of instances
  # @return weighted sum of Kendall distance of instance to corresp center
  
  if (is.null(weights)) {
    weights <- rep(1, length(cl@cluster))
  }
  sum(cl@cldist[ ,1] * weights)
}

wdistintra <- function(cl, distmat, weights = NULL, clustering = NULL, ...)
{
  # @return weighted sum of Kendall distance between instances lying in the same cluster
  
  if (is.null(clustering)) {
    clustering <- cl@cluster
  }
  if (length(unique(clustering)) == 1) {
    return(NA)
  }
  if (is.null(weights)) {
    weights <- rep(1, length(clustering))
  }
  diag(distmat) <- 0
  dintra <- sapply(unique(clustering), function(i){
    idx <- which(clustering == i)
    imat <- upper.tri(distmat[idx,idx], diag = TRUE)
    wmat <- tcrossprod(weights[idx])
    diag(wmat) <- (diag(wmat) - weights[idx])/2
    sum((wmat * distmat[idx,idx])[imat])/sum(wmat[imat])
  })
  mean(dintra)
}

wdistinter <- function(cl, distmat, weights = NULL, clustering = NULL, ...)
{
  # @return weighted sum of Kendall distance between instances lying in different clusters
  
  if (is.null(clustering)) {
    clustering <- cl@cluster
  }
  if (length(unique(clustering)) == 1) {
    return(NA)
  }
  if (is.null(weights)) {
    weights <- rep(1, length(clustering))
  }
  diag(distmat) <- 0
  dinter <- apply(matrix(combn(unique(clustering), 2), nrow = 2), 2, function(i){
    idx1 <- which(clustering == i[1])
    idx2 <- which(clustering == i[2])
    wmat <- tcrossprod(weights[idx1], weights[idx2])
    sum(wmat * distmat[idx1,idx2])/sum(wmat)
  })
  mean(dinter)
}

wdunn <- function(cl, distmat, weights = NULL, clustering = NULL, ...)
{
  # @return weighted dunn index
  
  if (is.null(clustering)) {
    clustering <- cl@cluster
  }
  if (length(unique(clustering)) == 1) {
    return(NA)
  }
  if (is.null(weights)) {
    weights <- rep(1, length(clustering))
  }
  diag(distmat) <- 0
  sep <- apply(matrix(combn(unique(clustering), 2), nrow = 2), 2, function(i){
    idx1 <- which(clustering == i[1])
    idx2 <- which(clustering == i[2])
    min(distmat[idx1,idx2])
  })
  diam <- sapply(unique(clustering), function(i){
    idx <- which(clustering == i)
    max(distmat[idx,idx])
  })
  min(sep)/max(diam)
}

wsilhouette <- function(cl, distmat, weights = NULL, clustering = NULL, ...)
{
  # @return weighted average silhouette width
  
  if (is.null(clustering)) {
    clustering <- cl@cluster
  }
  if (length(unique(clustering)) == 1) {
    return(NA)
  }
  if (is.null(weights)) {
    weights <- rep(1, length(clustering))
  }
  diag(distmat) <- 0
  b <- sapply(1:length(clustering), function(inst){
    b_other <- sapply(setdiff(unique(clustering), clustering[inst]), function(i){
      idx <- which(clustering == i)
      sum(weights[idx] * distmat[inst,idx])/sum(weights[idx])
    })
    min(b_other)
  })
  a <- sapply(1:length(clustering), function(inst){
    idx <- which(clustering == clustering[inst])
    sum(weights[idx] * distmat[inst,idx])/(sum(weights[idx]) - 1)
  })
  sum(weights * (b-a)/pmax(a,b)) / sum(weights)
}

# mallows -----------------------------------------------------------------

mallows_validate <- function(type, mallowsres, keylist, klist, nrepeats, ..., 
                             func = mget(paste0("w", type), mode = "function", envir = globalenv()))
{
  # @param mallowsres a 3-layer list representing nclust, key (type of mallows), repeat respectively followed by result (which is a list too) returned from RMallow::Mallows
  # @return validation score data.frame
  
  if (is.null(names(mallowsres)) && max(klist) == length(mallowsres)) {
    names(mallowsres) <- 1:length(mallowsres)
  }
  
  mallowsstats <- list()
  for (k in klist) {
    mallowsstats[[k]] <- data.frame(nclust = k, 
                                    rep = rep(1:nrepeats, times = length(keylist)), 
                                    method = rep(keylist, each = nrepeats))
    for (t in 1:length(type)) {
      mallowsstats[[k]][[type[t]]] <- sapply(unlist(mallowsres[[as.character(k)]], recursive = FALSE), function(mod){
        func[[t]](mod, clustering = mod$datas$clust, ...)
      })
    }
  }
  mallowsstats <- as.data.frame(do.call("rbind", mallowsstats), stringsAsFactors = TRUE)
  mallowsstats$nclust <- ordered(mallowsstats$nclust)
  mallowsstats
}

# specific validation functions
wloglike <- function(model, ...)
{
  # @param model object returned by RMallow::Mallows, ess. having $min.like
  # @return maximized likelihood
  
  if (is.null(model)) {
    NULL
  } else {
    model[["min.like"]][max(which(model[["min.like"]] != 0))]
  }
}

wbic <- function(model, G = length(model$R), weights = NULL, 
                 df = if (isTRUE(grepl("kernel", model$key, ignore.case = TRUE))) (G-1+G+G*length(model$R[[1]])) else (3*G-1), 
                 ...)
{
  # @param model object returned by RMallow::Mallows, ess. having $min.like and $R and $datas
  # @param G number of components
  # @param nsample number of samples
  # @param df number of free parameters of mixture model
  # @return BIC penalized (maximized) likelihood
  
  if (is.null(model)) {
    NULL
  } else {
    if (is.null(weights)) {
      weights <- rep(1, nrow(model$datas))
    }
    nsample <- sum(weights)
    2 * wloglike(model) - df * log(nsample)
  }
}

wicl <- function(model, weights = NULL, ...)
{
  # @param model object returned by RMallow::Mallows, ess. having $min.like and $R and $datas$pvals
  # @return ICL penalized complete (maximized) likelihood
  
  if (is.null(model)) {
    NULL
  } else {
    if (is.null(weights)) {
      weights <- rep(1, nrow(model$datas))
    }
    bic <- wbic(model = model, weights = weights, ...)
    z <- model[["datas"]][ , grep("pvals", names(model$datas))]
    bic + 2 * sum(weights * z * log(z))
  }
}

wloglikeML <- function(model, ...)
{
  # @param model object returned by RMallow::Mallows, ess. having $min.like
  # @return mallows loglik for GMM; same as wloglike
  
  if (is.null(model[["min.like.mallows"]])) {
    wloglike(model = model, ...)
  } else {
    model[["min.like.mallows"]][max(which(model[["min.like.mallows"]] != 0))]
  }
}
