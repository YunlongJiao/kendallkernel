
# modified version of flexclust::kcca -------------------------------------
# so that weights indicate simply how many duplicates of same instances are present
# equiv to running duplicated data with option weights=NULL

kcca_modified <- function(x, k, family=kccaFamily("kmeans"), weights=NULL, 
                          group=NULL, control=NULL, simple=FALSE, save.data=FALSE)
{
    MYCALL <- match.call()
    control <- as(control, "flexclustControl")
    x <- as(x, "matrix")
    x <- family@preproc(x)
    N <- nrow(x)
    
    if(control@classify=="auto"){
        control@classify="hard"
    }

    if(!is.null(group))
    {       
        if(length(group)>N)
            warning("group vector longer than nrow(x)")
        
        group <- rep(group, length=N)
        group <- as.factor(group)
    }
    
    if(!is.null(weights))
    {
        control@classify="weighted"

        if(!family@weighted)
            stop("Centroid function of family cannot use weights")

        if(!is.null(group))
            stop("Weighted clustering for grouped data is not implemented yet.")
        ## recycle to number of observations
        weights <- rep(weights, length=N)
    }
    
    centers <- initCenters(x, k, family, control)
    cluster <- centers$cluster
    k <- centers$k
    centers <- centers$centers

    sannprob <- control@simann[1]
    if(control@classify %in% c("hard", "simann"))
    {
        for(iter in 1:control@iter.max){
            
            clustold <- cluster
            distmat <- family@dist(x, centers)
            cluster <- family@cluster(x, distmat=distmat)

            if(control@classify == "simann"){
                cluster <- perturbClusters(cluster, sannprob)
                sannprob <- sannprob * control@simann[2]
            }

            if(!is.null(group))
                cluster <- family@groupFun(cluster, group, distmat)
            
            centers <- family@allcent(x, cluster=cluster, k=k)

            ## NAs in centers are empty clusters
            centers <- centers[complete.cases(centers),,drop=FALSE]
            k <- nrow(centers)
            
            changes <- sum(cluster!=clustold)
            if(control@verbose && (iter%%control@verbose==0)){
                td <- sum(distmat[cbind(1:N, cluster)])
                printIter(iter, paste(changes, format(td), sep=" / "),
                          "Changes / Distsum")
            }

            if(changes==0) break
        }
    }
    else if(control@classify=="weighted")
    {
        td <- -1
        for(iter in 1:control@iter.max)
        {
            td.old <- td
            clustold <- cluster ### added
            distmat <- family@dist(x, centers)
            cluster <- family@cluster(distmat=distmat)
            
            changes <- sum(cluster!=clustold) ### added
            if(changes==0) break ### added
            
#             td <- sum(distmat[cbind(1:N, cluster)])
            td <- sum(distmat[cbind(1:N, cluster)] * weights) ### added
            
            if(control@verbose && (iter%%control@verbose==0))
                printIter(iter, td, "Sum of distances")
                                
            if(abs(td-td.old)/(abs(td)+0.1) < control@tolerance) break

#             ## for weight computation we need similarities
#             distmat <- mlogit(distmat)
# 
#             for(n in 1:k){
#                 w <- weights*distmat[,n]
#                 w[cluster==n] <- w[cluster==n]+control@gamma
#                 centers[n,] <- family@wcent(x, weights=w)
#             }
            
            centers <- family@allcent(x, cluster=cluster, k=k, weights = weights) ### added
        }
#         ## make sure the final centers are centroids of their partitions
#         for(n in 1:k){
#             centers[n,] <- family@cent(x[cluster==n,,drop=FALSE])
#         }
    }
    else
        stop("Unknown classification method")

    centers <- centers[complete.cases(centers),,drop=FALSE]
    
    z <- newKccaObject(x=x, family=family, centers=centers, group=group,
                       iter=iter,
                       converged=(iter<control@iter.max),
                       call=MYCALL,
                       control=control,
                       simple=simple)

    if(save.data)
        z@data <- ModelEnvMatrix(designMatrix=x)

    z
}

environment(kcca_modified) <- asNamespace("flexclust")

setClass("kcca_modified",
         contains="kccasimple",
         representation(second="integer",
                        xrange="ANY",           # range of all data
                        xcent="ANY",            # centroid of all data
                        totaldist="numeric",    # total dist data<->xcent
                        clsim="matrix"))

# modified version of flexclust::bootFlexclust ----------------------------
# add option family mainly for the family@preproc for ncol of processed x
# add weights option so that duplicates are taken into account
# (minor) plus add default option FUN=kcca_modified and pass to stepFlexclust

bootFlexclust_modified <- function(x, k, nboot=100, correct=TRUE, seed=NULL,
                                   FUN=kcca_modified, family=kccaFamily("kmeans"), weights = NULL, # added
                                   multicore=FALSE, verbose=FALSE, ...)
{
    MYCALL <- match.call()
    
    if(!is.null(seed)) set.seed(seed)
    seed <- round(2^31 * runif(nboot, -1, 1))
    
    nk <- length(k)
    nx <- nrow(x)
    
    index1 <- matrix(integer(1), nrow=nx, ncol=nboot)
    index2 <- index1
    
    ## empirical experiments show parallization does not pay for this
    ## (sample is too fast)
    for(b in 1:nboot){
        index1[,b] <- sample(1:nx, nx, replace=TRUE)
        index2[,b] <- sample(1:nx, nx, replace=TRUE)
    }
    
    BFUN <- function(b){
        if(verbose){
            if((b %% 100) == 0)
                cat("\n")
            if((b %% 10) == 0)
                cat(b, "")
        }
        
        set.seed(seed[b])
        s1 <- stepFlexclust(x[index1[,b],,drop=FALSE], k=k, verbose=FALSE,
                            FUN=FUN, family=family, weights = weights[index1[,b]], # added
                            simple=TRUE, multicore=FALSE, ...)
        s2 <- stepFlexclust(x[index2[,b],,drop=FALSE], k=k, verbose=FALSE,
                            FUN=FUN, family=family, weights = weights[index2[,b]], # added
                            simple=TRUE, multicore=FALSE, ...)
        
        clust1 <- clust2 <- matrix(integer(1), nrow=nx, ncol=nk)
        cent1 <- cent2 <- list()
        rand <- double(nk)
        
        for(l in 1:nk)
        {
            if(nk>1){
                cl1 <- getModel(s1, l)
                cl2 <- getModel(s2, l)
            }
            else{
                cl1 <- s1
                cl2 <- s2
            }
            
            clust1[,l] <- clusters(cl1, newdata=x)
            clust2[,l] <- clusters(cl2, newdata=x)
            
            cent1[[l]] <- cl1@centers
            cent2[[l]] <- cl2@centers
            
#             rand[l] <- randIndex(table(clust1[,l], clust2[,l]),
#                                  correct=correct)
            if (is.null(weights)) { # added
                rand[l] <- randIndex(table(clust1[,l], clust2[,l]), # added
                                     correct=correct) # added
            } else { # added
                stopifnot(length(weights) == nrow(x)) # added
                rand[l] <- randIndex(table(rep(clust1[,l], weights), rep(clust2[,l], weights)), # added
                                     correct=correct) # added
            } # added
        }
        list(cent1=cent1, cent2=cent2, clust1=clust1, clust2=clust2,
             rand=rand)
        
    }
    
    ## empirical experiments show parallization does not pay for the 
    ## following (element extraction from list is too fast)
    z <- MClapply(as.list(1:nboot), BFUN, multicore=multicore)
    
    clust1 <- unlist(lapply(z, function(x) x$clust1))
    clust2 <- unlist(lapply(z, function(x) x$clust2))
    dim(clust1) <- dim(clust2) <- c(nx, nk, nboot)
    
    cent1 <- cent2 <- list()
    for(l in 1:nk){
        cent1[[l]] <- unlist(lapply(z, function(x) x$cent1[[l]]))
        cent2[[l]] <- unlist(lapply(z, function(x) x$cent2[[l]]))
        # dim(cent1[[l]]) <- dim(cent2[[l]]) <- c(k[l], ncol(x), nboot)
        dim(cent1[[l]]) <- dim(cent2[[l]]) <- c(k[l], ncol(family@preproc(x)), nboot) # added
    }
    
    if(nk > 1)
        rand <- t(sapply(z, function(x) x$rand))
    else
        rand <- as.matrix(sapply(z, function(x) x$rand))
    
    colnames(rand) <- k
    
    if(verbose) cat("\n")
    
    new("bootFlexclust", k=as.integer(k), centers1=cent1, centers2=cent2,
        cluster1=clust1, cluster2=clust2, index1=index1, index2=index2,
        rand=rand, call=MYCALL)
}

environment(bootFlexclust_modified) <- asNamespace("flexclust")

setClass("bootFlexclust_modified",
         representation(k="integer",
                        centers1="list",
                        centers2="list",
                        cluster1="array",
                        cluster2="array",
                        index1="matrix",
                        index2="matrix",
                        rand="matrix",
                        call="call"))

setMethod("show", "bootFlexclust_modified",
          function(object){
              cat("An object of class", sQuote(class(object)),"\n\n")
              cat("Call:\n")
              print(object@call)
              cat("\nNumber of bootstrap pairs:", nrow(object@rand),"\n")
          })

setMethod("summary", "bootFlexclust_modified",
          function(object){
              cat("Call:\n")
              print(object@call)
              cat("\nSummary of Rand Indices:\n")
              print(summary(object@rand))
          })

setMethod("plot", signature("bootFlexclust_modified","missing"),
          function(x, y, ...){
              boxplot(x, ...)
          })

setMethod("boxplot", "bootFlexclust_modified",
          function(x, ...){
              boxplot(as.data.frame(x@rand), ...)
          })

setMethod("densityplot", "bootFlexclust_modified",
          function(x, data, ...){
              Rand <- as.vector(x@rand)
              k <- rep(colnames(x@rand), rep(nrow(x@rand), ncol(x@rand)))
              k <- factor(k, levels=colnames(x@rand))
              
              densityplot(~Rand|k, as.table=TRUE, to=1, ...)
          })

# modified function body of family@cluster --------------------------------
# tie breaking method fixed when multiple centers are within equal dist to an instance to avoid convergence issue

cluster.body <- expression({
    if (is.null(distmat)) 
        distmat <- z@dist(x, centers)
    if (n == 1) {
        return(max.col(-distmat, ties.method = "first")) ### modified
    }
    else {
        r <- t(matrix(apply(distmat, 1, rank, ties.method = "first"), ### modified
                      nrow = ncol(distmat)))
        z <- list()
        for (k in 1:n) z[[k]] <- apply(r, 1, function(x) which(x == 
                                                                   k))
    }
    return(z)
})

# modified function body of family@allcent --------------------------------
# incorporate weights option without using family@wcent

allcent.body <- expression({
    centers <- matrix(NA, nrow = k, ncol = ncol(x))
    for (n in 1:k) {
        if (sum(cluster == n, na.rm = TRUE) > 0) {
            centers[n, ] <- z@cent(x[cluster == n, , drop = FALSE], weights = weights[cluster == n]) # modified
        }
    }
    centers
})


# specific function to transfer an obj of family --------------------------

allow_weights <- function(family = NULL, name = NULL, dist = NULL, cent = NULL, preproc = NULL)
{
    if (is.null(family)) {
        family <- flexclust::kccaFamily(which = NULL, dist = dist, cent = cent, preproc = preproc)
        family@name <- name
        # body(family@cluster) <- cluster.body
    }
    
    family@weighted <- TRUE
    family@wcent <- function() stop("should not try accessing family@wcent at all! do use function ", dQuote("kcca_modified"))
    formals(family@allcent) <- c(formals(family@allcent), alist(weights = NULL))
    body(family@allcent) <- allcent.body
    family
}

# define dist, cent, and preproc functions for different versions  --------

# brute
dist_bruteKmeans <- function(x, centers)
{
  # DO NOT RUN for over-sized permutation
  
  if (ncol(x) != ncol(centers)) {
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  }
  AllKendall(r = x, seqs = centers)
}

cent_bruteKmeans <- function(x, weights = NULL)
{
  # DO NOT RUN for over-sized permutation
  
  if (is.null(weights)) {
    weights <- rep(1, nrow(x))
  }
  perm <- do.call("rbind", combinat::permn(ncol(x)))
  dist <- AllKendall(r = x, seqs = perm)
  idx <- which.min(colSums(dist * weights))
  perm[idx, ]
}

preproc_bruteKmeans <- NULL # do nothing by default

# copeland
dist_copelandKmeans <- function(x, centers)
{
  # kendall distance
  if (ncol(x) != ncol(centers)) {
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  }
  AllKendall(r = x, seqs = centers)
}

cent_copelandKmeans <- function(x, weights = NULL)
{
  if (is.null(weights)) {
    weights <- rep(1, nrow(x))
  }
  
  abils <- ncol(x)
  inds <- combinat::combn(abils, 2)
  num <- c(0, cumsum((abils-1):1))
  if (is.vector(x)) {
    x <- as.matrix(x, nrow = 1)
  } else if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  infos <- sign(x[, inds[1, ], drop = FALSE] - x[, inds[2, ], drop =  FALSE])
  cent.info <- sign(colSums(infos * weights))
  cent <- 0*(1:abils)
  for (i in 1:(abils-1)) {
    cent[i] <- cent[i] + sum(cent.info[(num[i]+1):(num[i+1])])
  }
  for (i in 2:abils) {
    cent[i] <- cent[i] - sum(cent.info[num[2:i]-(abils-i)])
  }
  rank(cent, ties.method = "first")
}

preproc_copelandKmeans <- NULL # do nothing by default

# borda
dist_bordaKmeans <- function(x, centers)
{
  # kendall distance
  if (ncol(x) != ncol(centers)) {
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  }
  AllKendall(r = x, seqs = centers)
}

cent_bordaKmeans <- function(x, weights = NULL)
{
  if (is.null(weights)) {
    weights <- rep(1, nrow(x))
  }
  cent <- apply(x * weights, 2, mean)
  rank(cent, ties.method = "first")
}

preproc_bordaKmeans <- NULL # do nothing by default

# kernel
dist_kernelKmeans <- function(x, centers)
{
  # kendall distance
  if (ncol(x) != ncol(centers)) {
    stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
  }
  kernrank:::distL2(r = x, centers = centers)
}

cent_kernelKmeans <- function(x, weights = NULL)
{
  if (is.null(weights)) {
    weights <- rep(1, nrow(x))
  }
  colMeans(x * (weights / mean(weights)))
}

preproc_kernelKmeans <- KendallInfo

