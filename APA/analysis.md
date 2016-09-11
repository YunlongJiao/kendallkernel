---
title: "kendallkernel - cluster full rankings"
author: "Yunlong Jiao"
date: "04 Jan 2016"
output: html_document
---


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, fig.width = 15, fig.height = 12, dev = "pdf", fig.keep = "high", fig.path = "figure/", cache.path = "cache/")
set.seed(61941263)

# utiles
library(caret)
library(parallel)
library(mvtnorm)
library(reshape2)
library(ggplot2)
library(gplots)

# clustering
library(kernrank)
library(flexclust)
library(Rankcluster)
source("flexclust_modified.R") # modified version of flexclust::kcca with weights
source("cluster_validation.R") # cluster validation criteria
```

## Data

`Rankcluster::APA` dataset is election data consisting of full rankings. An thorough analysis with this dataset is found in _Group Representations in Probability and Statistics_ by Persi Diaconis (1988).


```r
# load data
data(APA)
dat <- APA$frequency[ ,1:(APA$m)]
freq <- APA$frequency[ ,"frequency", drop = TRUE]

# compute distance matrix
distmat <- AllKendall(r = dat, seqs = dat)
rownames(distmat) <- colnames(distmat) <- apply(dat, 1, paste0, collapse = "")

# number of instances (duplicates allowed)
nrow(APA$data)
```

```
## [1] 5738
```

```r
# full rankings of size m
APA$m
```

```
## [1] 5
```

## Kendall kmeans

- **brute kmeans** minimizes sum of Kendall distance where centers are found at each iteration using brute-force search over all possible permutations. Note non-unique optimal center is possible and the implementation always picks the first permutation returned by `combinat::permn`.
- **copeland kmeans** minimizes sum of Kendall distance where centers are found at each iteration using Copeland rank aggreg (from an unweighted directed graph with directional edges as majority voting for pairwise wins). Note non-unique optimal center is possible and we pick in a greedy manner by _first-in first-choice (FIFC)_ of item labels to break cycles.
- **borda kmeans** minimizes sum of Kendall distance where centers are found at each iteration using Borda Count rank aggreg (a 5-approx to footrule optimal aggreg) where ties are broken by FIFC of item labels.
- **kernel kmeans** minimizes sum of Kendall distance where centers found at each iteration could be any point found in the feature space instead of being restricted to the image of a true permutation. Note `flexclust::kcca` needs to preprocess `x` so that rows represent images of permutations in the feature space of dim $n(n-1)$, specifically $\Phi(\sigma)=(sign(\sigma_i-\sigma_j))_{i<j}$.


```r
# define kcca family
keylistKM <- c("bruteKmeans", "copelandKmeans", "bordaKmeans", "kernelKmeans")
for (key in keylistKM) {
    assign(paste0("family_", key), 
           allow_weights(name = key, 
                         dist = get(paste0("dist_", key)), 
                         cent = get(paste0("cent_", key)), 
                         preproc = get(paste0("preproc_", key)))
           )
}

# testrun
for (key in keylistKM) {
    cat("\n=========================", key, "=========================\n")
    family <- get(paste0("family_", key))
    k <- 5
    ptm <- system.time(cl <- kcca_modified(x = dat, k = k, family = family, weights = freq, save.data = FALSE))
    cat("runtime = \n")
    print(ptm)
    summary(cl)
    cat("weighted sum of distance to corresp center", wdistsum(cl = cl, weights = freq))
}
```

```
## 
## ========================= bruteKmeans =========================
## runtime = 
##    user  system elapsed 
##   0.758   0.036   0.794 
## kcca object of family 'bruteKmeans' 
## 
## call:
## kcca_modified(x = dat, k = k, family = family, weights = freq, 
##     save.data = FALSE)
## 
## cluster info:
##   size  av_dist max_dist separation
## 1   23 2.130435        3          3
## 2   25 2.240000        4          3
## 3   26 2.461538        4          3
## 4   23 2.173913        4          3
## 5   23 2.217391        4          2
## 
## convergence after 8 iterations
## sum of within cluster distances: 270 
## weighted sum of distance to corresp center 11348
## ========================= copelandKmeans =========================
## runtime = 
##    user  system elapsed 
##   0.102   0.000   0.101 
## kcca object of family 'copelandKmeans' 
## 
## call:
## kcca_modified(x = dat, k = k, family = family, weights = freq, 
##     save.data = FALSE)
## 
## cluster info:
##   size  av_dist max_dist separation
## 1   26 2.346154        4          3
## 2   21 2.190476        4          2
## 3   25 2.280000        4          3
## 4   23 2.260870        4          2
## 5   25 2.320000        4          3
## 
## convergence after 3 iterations
## sum of within cluster distances: 274 
## weighted sum of distance to corresp center 12363
## ========================= bordaKmeans =========================
## runtime = 
##    user  system elapsed 
##   0.138   0.000   0.138 
## kcca object of family 'bordaKmeans' 
## 
## call:
## kcca_modified(x = dat, k = k, family = family, weights = freq, 
##     save.data = FALSE)
## 
## cluster info:
##   size  av_dist max_dist separation
## 1   23 2.130435        3          3
## 2   27 2.407407        4          3
## 3   21 2.142857        4          2
## 4   23 2.304348        4          3
## 5   26 2.307692        4          3
## 
## convergence after 8 iterations
## sum of within cluster distances: 272 
## weighted sum of distance to corresp center 12274
## ========================= kernelKmeans =========================
## runtime = 
##    user  system elapsed 
##   0.083   0.000   0.084 
## kcca object of family 'kernelKmeans' 
## 
## call:
## kcca_modified(x = dat, k = k, family = family, weights = freq, 
##     save.data = FALSE)
## 
## cluster info:
##   size  av_dist max_dist separation
## 1   25 1.555759 2.925321   2.125931
## 2   18 1.341529 2.221791   2.027778
## 3   28 1.640970 2.589069   2.024327
## 4   24 1.519737 2.389033   1.843796
## 5   25 1.531425 2.450160   1.792294
## 
## convergence after 8 iterations
## sum of within cluster distances: 183.748 
## weighted sum of distance to corresp center 7968.362
```

```r
# list of number of cluster
klist <- 2:10

# number of repeated runs to avoid local minima
nrepeats <- 50

# vary nclust and get clustering
clres <- list()
for (k in klist) {
    clres[[k]] <- lapply(keylistKM, function(key){
        family <- get(paste0("family_", key))
        cls <- lapply(1:nrepeats, function(i){
            kcca_modified(x = dat, k = k, family = family, weights = freq)
        })
        names(cls) <- 1:nrepeats
        cls
    })
    names(clres[[k]]) <- keylistKM
}
```

## Mallows mixtures

See [Murphy and Martin (2003)](http://www.researchgate.net/publication/4897368_Mixtures_of_distance-based_models_for_ranking_data) for mixture of Mallows models applied for clustering full rankings. We exploit simply the case where dispersion parameter $\lambda_g$ can be different for various clusters $g=1,\dots,G$ and cluster centers are found using different techniques below:

- **brute mallows**
- **copeland mallows**
- **borda mallows**
- **kernel mallows** where, at each iteration, centers are set as weighted average in feature space without being fully optimized and lambda that defines properly a proba distribution (namely a normalization constant) on permutation group is optimized once wrt those centers via a bifurcation method implemented by `stats::uniroot` since loglike $\ell$ is convex wrt lambda alone and the gradient equation wrt lambda is monotonic.
- **kernel mallows exhaustive** where, at each iteration, centers and lambda are alternately fully optimized implemented by `stats::nlm`. NOTE solutions are exhaustively searched but still NOT guaranteed global optima since loglik is not concave wrt the bivariate couple (center, lambda), although it is concave separately wrt either center or lambda.
- **xxx equal lambda** where xxx indicates any type above but lambda is optimized subjected to being equal across all groups


```r
# define mallows center-finding algorithm
keylistML <- c("bruteMallows", "copelandMallows", "bordaMallows", "kernelMallows", "kernelMallows_Exh")
keylistML <- c(keylistML, paste0(keylistML, "_Eqlam"))

# testrun
likes <- list()
for (key in keylistML) {
    cat("\n=========================", key, "=========================\n")
    k <- 5
    ptm <- system.time(cl <- Mallows(datas = dat, G = k, weights = freq, key = key))
    likes[[key]] <- cl$min.like
    cat("runtime = \n")
    print(ptm)
    cat("lambda estimate = \n")
    print(round(cl$lambda, 3))
    cat("cluster size (without weights) = \n")
    print(table(cl$datas$clust))
}
```

```
## 
## ========================= bruteMallows =========================
## runtime = 
##    user  system elapsed 
##   3.046   0.029   3.077 
## lambda estimate = 
## [1] 0.168 0.584 0.649 0.188 0.584
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 29  6 15 41 29 
## 
## ========================= copelandMallows =========================
## runtime = 
##    user  system elapsed 
##   1.968   0.005   1.975 
## lambda estimate = 
## [1] 0.168 0.584 0.649 0.188 0.584
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 29  6 15 41 29 
## 
## ========================= bordaMallows =========================
## runtime = 
##    user  system elapsed 
##   2.045   0.002   2.047 
## lambda estimate = 
## [1] 0.127 0.529 0.668 0.162 0.720
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 44  6 14 38 18 
## 
## ========================= kernelMallows =========================
## runtime = 
##    user  system elapsed 
##   1.564   0.013   1.579 
## lambda estimate = 
## [1] 1.694 1.647 1.189 1.153 1.981
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 14 13 43 41  9 
## 
## ========================= kernelMallows_Exh =========================
## runtime = 
##    user  system elapsed 
##  18.102   0.107  18.222 
## lambda estimate = 
## [1] 1.177 1.170 0.852 1.269 0.986
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 18 10 33 54  5 
## 
## ========================= bruteMallows_Eqlam =========================
## runtime = 
##    user  system elapsed 
##   1.176   0.007   1.186 
## lambda estimate = 
## [1] 0.447 0.447 0.447 0.447 0.447
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
##  2 11 44 12 51 
## 
## ========================= copelandMallows_Eqlam =========================
## runtime = 
##    user  system elapsed 
##   0.613   0.003   0.616 
## lambda estimate = 
## [1] 0.447 0.447 0.447 0.447 0.447
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
##  2 11 44 12 51 
## 
## ========================= bordaMallows_Eqlam =========================
## runtime = 
##    user  system elapsed 
##   0.802   0.000   0.802 
## lambda estimate = 
## [1] 0.44 0.44 0.44 0.44 0.44
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
##  2 27 38 12 41 
## 
## ========================= kernelMallows_Eqlam =========================
## runtime = 
##    user  system elapsed 
##   0.669   0.005   0.677 
## lambda estimate = 
## [1] 1.329 1.329 1.329 1.329 1.329
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 28 17 31 27 17 
## 
## ========================= kernelMallows_Exh_Eqlam =========================
## runtime = 
##    user  system elapsed 
##  17.766   0.096  17.862 
## lambda estimate = 
## [1] 0.961 0.961 0.961 0.961 0.961
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 18 10 33 54  5
```

```r
# plot observed likelihood over iteration
likes <- do.call("cbind", likes)
likes <- melt(likes, varnames = c("iter", "method"), value.name = "loglike")
likes$iter <- ordered(likes$iter)
ggplot(likes, aes(x = iter, y = loglike)) + 
  geom_line(aes(group = method, colour = method))
```

![plot of chunk mallows](figure/mallows-1.pdf)

```r
# vary nclust and get clustering
mallowsres <- mclapply(klist, function(k){
  res <- lapply(keylistML, function(key){
    cls <- lapply(1:nrepeats, function(i){
      Mallows(datas = dat, G = k, weights = freq, key = key)
    })
    names(cls) <- 1:nrepeats
    cls
  })
  names(res) <- keylistML
  res
}, mc.cores = 16)
names(mallowsres) <- klist
```

## Gaussian mixtures in feature space

- **kernel gaussian** builds a gaussian mixture and does clustering in feature space despite that it fails to form a proper proba distribution on permutation group, then normalization constant is sought based on $R$ and $lambda$ such that `loglik` corresp to true distribution over permutation group. Kernel trick may apply via [Wang et al. (2003)](http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.114.2451)


```r
# define gaussian mixture in feature space
keylistGS <- c("kernelGaussian")
keylistGS <- c(keylistGS, paste0(keylistGS, "_Eqlam"))
keylistML <- c(keylistML, keylistGS)

# test run
likes <- list()
for (key in keylistGS) {
    cat("\n=========================", key, "=========================\n")
    k <- 5
    ptm <- system.time(cl <- Mallows(datas = dat, G = k, weights = freq, key = key))
    likes[[key]] <- cl$min.like
    cat("runtime = \n")
    print(ptm)
    cat("lambda estimate = \n")
    print(round(cl$lambda, 3))
    cat("cluster size (without weights) = \n")
    print(table(cl$datas$clust))
}
```

```
## 
## ========================= kernelGaussian =========================
## runtime = 
##    user  system elapsed 
##   1.188   0.020   1.215 
## lambda estimate = 
## [1] 4.194 4.499 4.445 2.784 3.048
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 17 16 16 41 30 
## 
## ========================= kernelGaussian_Eqlam =========================
## runtime = 
##    user  system elapsed 
##   1.108   0.001   1.114 
## lambda estimate = 
## [1] 3.649 3.649 3.649 3.649 3.649
## cluster size (without weights) = 
## 
##  1  2  3  4  5 
## 25 18 25 25 27
```

```r
# plot observed likelihood over iteration
likes <- do.call("cbind", likes)
likes <- melt(likes, varnames = c("iter", "method"), value.name = "loglike")
likes$iter <- ordered(likes$iter)
ggplot(likes, aes(x = iter, y = loglike)) + 
  geom_line(aes(group = method, colour = method))
```

![plot of chunk gaussian](figure/gaussian-1.pdf)

```r
# vary nclust and get clustering
for (k in klist) {
  for (key in keylistGS) {
    cls <- lapply(1:nrepeats, function(i){
      Mallows(datas = dat, G = k, weights = freq, key = key)
    })
    names(cls) <- 1:nrepeats
    mallowsres[[as.character(k)]][[key]] <- cls
  }
}
```

## Results

Compare methods with respect to several internal validation criteria. Say $m$ distinct instances $\{x_i\}$ with frequency $\{w_i\}$, in total $n = \sum_{i=1}^m w_i$ instances allowing duplicates, are divided into $K$ groups $\{G_j\}$ with center $\{c_j\}$ (and with fuzzy proba $z_i$ for each instance). Define for **kmeans models**

- `distsum` is the objective function to minimize in kmeans as the sum of Kendall distance between instances and corresp centers, i.e., $\sum_{k=1}^K \sum_{x_i \in G_k} w_i d(x_i,c_k)$
- `distintra` is average of distance between instances lying in the same cluster, i.e., $\frac{1}{K} \sum_{k=1}^K \frac{\sum_{x_i\neq x_j \in G_k} w_i w_j d(x_i,x_j)}{\sum_{x_i\neq x_j \in G_k} w_i w_j + \sum_{x_i \in G_k} {w_i \choose 2}}$
- `distinter` is average of distance between instances lying in different clusters, i.e., $\frac{1}{{K \choose 2}} \sum_{k<l} \frac{\sum_{x_i\in G_k, x_j\in G_l} w_i w_j d(x_i,x_j)}{\sum_{x_i\in G_k, x_j\in G_l} w_i w_j}$
- `dunn` is the Dunn index measuring the ratio of inter-cluster minimum distance to intra-cluster maximum distance, i.e., $\frac{\min_{k<l} \min_{x_i\in G_k, x_j\in G_l} d(x_i,x_j)}{\max_{k} \max_{x_i\neq x_j \in G_k} d(x_i,x_j)}$
- `silhouette` is Silhouette width averaged over all instances contrasting their average distance to instances in the same cluster with their average distance to instances in the "neighbouring" (lowest of other) cluster, i.e., $\frac{\sum_{i=1}^m w_i (b_i-a_i)/\max\{a_i,b_i\}}{\sum_{i=1}^m w_i}$ where, provided that $x_i \in G_k$ (fixed), $a_i = \frac{\sum_{x_j \in G_k, x_j\neq x_i} w_j d(x_i,x_j)}{(w_i - 1) + \sum_{x_j \in G_k, x_j\neq x_i} w_j}$ and $b_i = \min_{l\neq k} \frac{\sum_{x_j \in G_l} w_j d(x_i,x_j)}{\sum_{x_j \in G_l} w_j}$

Boxplot below illustrates variance of repeated experiments from different initial cluster centers, and curves correspond to the optimal repeat with respect to min `distsum`.


```r
# for kmeans get clustering statistics and best repeat with min distsum
typeKM <- c("distsum", "distintra", "distinter", "dunn", "silhouette")
clstats <- kmeans_validate(type = typeKM, clres = clres, keylist = keylistKM, 
                           klist = klist, nrepeats = nrepeats, weights = freq, distmat = distmat)
clstats$method <- factor(clstats$method, levels = sort(levels(clstats$method)))

clstats_best <- c()
for (k in klist) {
    for (key in keylistKM) {
        x <- subset(clstats, method == key & nclust == k)
        clstats_best <- rbind(clstats_best, x[which.min(x$distsum),])
    }
}

# distsum (lower is better)
ggplot(clstats, aes(x = nclust, y = distsum)) + 
    geom_boxplot(aes(fill = method), alpha = 0.5) + 
    stat_summary(fun.y = min, geom = "line", aes(group = method, colour = method)) + 
    stat_summary(fun.y = min, geom = "point", aes(colour = method, shape = method), size = 3)
```

![plot of chunk kmeans_plot](figure/kmeans_plot-1.pdf)

```r
# distintra (lower is better)
ggplot() + 
    geom_boxplot(data = clstats, aes(x = nclust, y = distintra, fill = method), alpha = 0.5) + 
    geom_line(data = clstats_best, aes(x = nclust, y = distintra, group = method, colour = method)) + 
    geom_point(data = clstats_best, aes(x = nclust, y = distintra, colour = method, shape = method), size = 3)
```

![plot of chunk kmeans_plot](figure/kmeans_plot-2.pdf)

```r
# distinter (higher is better)
ggplot() + 
    geom_boxplot(data = clstats, aes(x = nclust, y = distinter, fill = method), alpha = 0.5) + 
    geom_line(data = clstats_best, aes(x = nclust, y = distinter, group = method, colour = method)) + 
    geom_point(data = clstats_best, aes(x = nclust, y = distinter, colour = method, shape = method), size = 3)
```

![plot of chunk kmeans_plot](figure/kmeans_plot-3.pdf)

```r
# dunn (higher is better)
ggplot() + 
    geom_boxplot(data = clstats, aes(x = nclust, y = dunn, fill = method), alpha = 0.5) + 
    geom_line(data = clstats_best, aes(x = nclust, y = dunn, group = method, colour = method)) + 
    geom_point(data = clstats_best, aes(x = nclust, y = dunn, colour = method, shape = method), size = 3)
```

![plot of chunk kmeans_plot](figure/kmeans_plot-4.pdf)

```r
# silhouette (higher is better)
ggplot() + 
    geom_boxplot(data = clstats, aes(x = nclust, y = silhouette, fill = method), alpha = 0.5) + 
    geom_line(data = clstats_best, aes(x = nclust, y = silhouette, group = method, colour = method)) + 
    geom_point(data = clstats_best, aes(x = nclust, y = silhouette, colour = method, shape = method), size = 3)
```

![plot of chunk kmeans_plot](figure/kmeans_plot-5.pdf)

Define for **mallows mixture models**

- `loglike` is the objective function to minimize in mallows mixture models or the marginal likelihood of observed variables
- `bic` is Bayesian information criterion which penalizes the model complexity, i.e., $BIC(K) = 2\ell(K) - df(K)\log(n)$
- `icl` is integrated complete likelihood which adds to `bic` an extra entropy term penalizing cases with uncertain membership assignment, i.e., $ICL(K) = BIC(K) + 2\sum_{k=1}^K \sum_{i=1}^n \hat{z}^{(i)}_k \log(\hat{z}^{(i)}_k) = BIC(K) - 2\sum_{i=1}^n Entropy(\hat{z}^{(i)})$

Boxplot below illustrates variance of repeated experiments from different initial cluster centers and curves correspond to the optimal repeat from various with respect to max `loglike`.


```r
# for mallows get clustering statistics
typeML <- c("loglike", "bic", "icl", "distintra", "distinter", "dunn", "silhouette", "loglikeML")
mallowsstats <- mallows_validate(type = typeML, mallowsres = mallowsres, keylist = keylistML, 
                                 klist = klist, nrepeats = nrepeats, weights = freq, distmat = distmat)
mallowsstats$method <- factor(mallowsstats$method, levels = sort(levels(mallowsstats$method)))

# check HOW Gaussian likelihood is related to Mallows likelihood (for same alpha/center/lambda)
for (key in keylistGS) {
  p1 <- ggplot(subset(mallowsstats, method == key), aes(x = loglike, y = loglikeML)) + 
    geom_line(aes(group = nclust, colour = nclust)) + 
    labs(title = key, x = "loglikeGaussian", y = "loglikeMallows")
  print(p1)
}
```

![plot of chunk mallows_plot](figure/mallows_plot-1.pdf)![plot of chunk mallows_plot](figure/mallows_plot-2.pdf)

```r
# get best repeat with max likelihood
mallowsstats_best <- c()
for (k in klist) {
    for (key in keylistML) {
        x <- subset(mallowsstats, method == key & nclust == k)
        mallowsstats_best <- rbind(mallowsstats_best, x[which.max(x$loglike),])
    }
}
# AND change loglike to loglikeML after best fit for GMM is selected wrt Gaussian loglike
mallowsstats$loglike <- mallowsstats$loglikeML
mallowsstats$loglikeML <- NULL
mallowsstats_best$loglike <- mallowsstats_best$loglikeML
mallowsstats_best$loglikeML <- NULL

# loglike (higher is better)
ggplot() + 
    geom_boxplot(data = mallowsstats, aes(x = nclust, y = loglike, fill = method), alpha = 0.5) + 
    geom_line(data = mallowsstats_best, aes(x = nclust, y = loglike, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-3.pdf)

```r
# bic (higher is better)
ggplot() + 
    geom_boxplot(data = mallowsstats, aes(x = nclust, y = bic, fill = method), alpha = 0.5) + 
    geom_line(data = mallowsstats_best, aes(x = nclust, y = bic, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-4.pdf)

```r
# icl (higher is better)
ggplot() + 
    geom_boxplot(data = mallowsstats, aes(x = nclust, y = icl, fill = method), alpha = 0.5) + 
    geom_line(data = mallowsstats_best, aes(x = nclust, y = icl, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-5.pdf)

```r
# loglike (higher is better; without kernelGaussian)
ggplot() + 
    geom_boxplot(data = subset(mallowsstats, !method %in% keylistGS), aes(x = nclust, y = loglike, fill = method), alpha = 0.5) + 
    geom_line(data = subset(mallowsstats_best, !method %in% keylistGS), aes(x = nclust, y = loglike, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-6.pdf)

```r
# bic (higher is better; without kernelGaussian)
ggplot() + 
    geom_boxplot(data = subset(mallowsstats, !method %in% keylistGS), aes(x = nclust, y = bic, fill = method), alpha = 0.5) + 
    geom_line(data = subset(mallowsstats_best, !method %in% keylistGS), aes(x = nclust, y = bic, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-7.pdf)

```r
# icl (higher is better; without kernelGaussian)
ggplot() + 
    geom_boxplot(data = subset(mallowsstats, !method %in% keylistGS), aes(x = nclust, y = icl, fill = method), alpha = 0.5) + 
    geom_line(data = subset(mallowsstats_best, !method %in% keylistGS), aes(x = nclust, y = icl, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-8.pdf)

```r
# distintra (lower is better)
ggplot() + 
    geom_boxplot(data = mallowsstats, aes(x = nclust, y = distintra, fill = method), alpha = 0.5) + 
    geom_line(data = mallowsstats_best, aes(x = nclust, y = distintra, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-9.pdf)

```r
# distinter (higher is better)
ggplot() + 
    geom_boxplot(data = mallowsstats, aes(x = nclust, y = distinter, fill = method), alpha = 0.5) + 
    geom_line(data = mallowsstats_best, aes(x = nclust, y = distinter, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-10.pdf)

```r
# dunn (higher is better)
ggplot() + 
    geom_boxplot(data = mallowsstats, aes(x = nclust, y = dunn, fill = method), alpha = 0.5) + 
    geom_line(data = mallowsstats_best, aes(x = nclust, y = dunn, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-11.pdf)

```r
# silhouette (higher is better)
ggplot() + 
    geom_boxplot(data = mallowsstats, aes(x = nclust, y = silhouette, fill = method), alpha = 0.5) + 
    geom_line(data = mallowsstats_best, aes(x = nclust, y = silhouette, group = method, colour = method))
```

![plot of chunk mallows_plot](figure/mallows_plot-12.pdf)

Finally, graphical illustration of clustering results (nclust selected by max of average silhouette width) with heatmaps.


```r
# for kmeans, optimal number of clusters is selected wrt average silhouette width
for (key in keylistKM){
    x <- subset(clstats_best, method == key)
    x$nclust <- as.integer(as.character(x$nclust))
    i <- which.max(x$silhouette)
    best_fit <- clres[[x$nclust[i]]][[key]][[x$rep[i]]]
    clustering <- best_fit@cluster
    o <- order(clustering)
    sep <- cumsum(table(clustering))
    heatmap.2(x = distmat[o,o], Rowv = "as-is", main = key, 
              symm = TRUE, colsep = sep, rowsep = sep, 
              dendrogram = "none", scale = "none", trace = "none")
}
```

![plot of chunk heatmap](figure/heatmap-1.pdf)![plot of chunk heatmap](figure/heatmap-2.pdf)![plot of chunk heatmap](figure/heatmap-3.pdf)![plot of chunk heatmap](figure/heatmap-4.pdf)

```r
# for mallows mixture models, optimal number of clusters is selected wrt BIC
for (key in keylistML) {
    x <- subset(mallowsstats_best, method == key)
    x$nclust <- as.integer(as.character(x$nclust))
    i <- which.max(x$bic)
    best_fit <- mallowsres[[as.character(x$nclust[i])]][[key]][[x$rep[i]]]
    clustering <- best_fit$datas$clust
    o <- order(clustering)
    sep <- cumsum(table(clustering))
    heatmap.2(x = distmat[o,o], Rowv = "as-is", main = key, 
              symm = TRUE, colsep = sep, rowsep = sep, 
              dendrogram = "none", scale = "none", trace = "none")
}
```

![plot of chunk heatmap](figure/heatmap-5.pdf)![plot of chunk heatmap](figure/heatmap-6.pdf)![plot of chunk heatmap](figure/heatmap-7.pdf)![plot of chunk heatmap](figure/heatmap-8.pdf)![plot of chunk heatmap](figure/heatmap-9.pdf)![plot of chunk heatmap](figure/heatmap-10.pdf)![plot of chunk heatmap](figure/heatmap-11.pdf)![plot of chunk heatmap](figure/heatmap-12.pdf)![plot of chunk heatmap](figure/heatmap-13.pdf)![plot of chunk heatmap](figure/heatmap-14.pdf)![plot of chunk heatmap](figure/heatmap-15.pdf)![plot of chunk heatmap](figure/heatmap-16.pdf)

## More for kmeans...

Good clustering solutions should be robust to perturbation in the input data from which centers are derived and further used to cluster the entire dataset, accouting for how similar clustering results are to bootstrap data.


```r
# number of repeated bootstrapping experiments
nboot <- 100

# get clustering from bootstrap input data
bclres <- list()
seed <- runif(length(keylistKM)) * 10^8; names(seed) <- keylistKM
for (key in keylistKM) {
    cat("\n=========================", key, "=========================\n")
    family <- get(paste0("family_", key))
    ptm <- proc.time()
    bclres[[key]] <- bootFlexclust_modified(x = dat, k = klist, nboot = nboot, multicore = FALSE, seed = seed[key], correct=TRUE, 
                                            FUN = kcca_modified, weights = freq, family = family)
    print(proc.time() - ptm)
}
```

```
## 
## ========================= bruteKmeans =========================
##     user   system  elapsed 
## 1697.756    8.464 1706.196 
## 
## ========================= copelandKmeans =========================
##    user  system elapsed 
## 248.863   0.425 249.188 
## 
## ========================= bordaKmeans =========================
##    user  system elapsed 
## 209.330   0.227 209.440 
## 
## ========================= kernelKmeans =========================
##    user  system elapsed 
## 133.326   0.363 133.615
```

Compare kmeans methods with respect to stability index. Define for **kmeans models**

- `rand` is the [adjusted rand index](https://en.wikipedia.org/wiki/Rand_index#Adjusted_Rand_index) accounting for percentage of accurate decisions made in terms of instance pairs falling in the same or in different clusters, implemented via `flexclust::randIndex`.

Boxplot below illustrates variance of repeated bootstrap experiments.


```r
# get stability statistics
bclstats <- list()
for (key in keylistKM) {
    bclstats[[key]] <- melt(bclres[[key]]@rand, varnames = c("rep", "nclust"), value.name = "rand")
    bclstats[[key]]$method <- key
}
bclstats <- as.data.frame(do.call("rbind", bclstats), stringsAsFactors = TRUE)
bclstats$nclust <- ordered(bclstats$nclust)

# rand (higher is better)
ggplot(bclstats, aes(x = nclust, y = rand, fill = method)) + 
    geom_boxplot(alpha = 0.5)
```

![plot of chunk stability_plot](figure/stability_plot-1.pdf)

## More for mallows...

Since number of free parameters is biased towards Mallows models where centers are true permutations (1 for each center) compared to centers are arbitrary in feature space ($d^2$ for each center), cross validate which bic should approx (up to an $O(1)$ term) is preferred to compare models fairly.


```r
# get test loglik based on cross-validation
mallowsresCV <- mclapply(klist, function(k){
  cls <- lapply(keylistML, function(key){
    MallowsCV(datas = dat, G = k, weights = freq, key = key)
  })
  names(cls) <- keylistML
  cls
}, mc.cores = 16)
names(mallowsresCV) <- klist
```

Define for **mallows mixture models**

- `cv.loglik` is cross-validated `loglik`, i.e., $cv.loglik(K) = \ell_{test}(K)$

Boxplot below illustrates variance of cross validated experiments and curves correspond to the average score.


```r
# get clustering statistics for each cv fold
mallowsstatsCV <- list()
i <- 0
for (k in klist) {
  for (key in keylistML) {
    i <- i + 1
    mallowsstatsCV[[i]] <- data.frame(nclust = k, method = key, rep = names(mallowsresCV[[as.character(k)]][[key]]), 
                                      cv.loglik = sapply(mallowsresCV[[as.character(k)]][[key]], function(u) u[["cv.loglik"]]))
  }
}
mallowsstatsCV <- do.call("rbind", mallowsstatsCV)
mallowsstatsCV$nclust <- ordered(mallowsstatsCV$nclust)
mallowsstatsCV$method <- factor(mallowsstatsCV$method, levels = sort(levels(mallowsstatsCV$method)))

# cv.loglik (perplexity; higher is better)
ggplot(mallowsstatsCV, aes(x = nclust, y = cv.loglik)) + 
    geom_boxplot(aes(fill = method), alpha = 0.5) + 
    stat_summary(fun.y = mean, geom = "line", aes(group = method, colour = method))
```

![plot of chunk mallowscv_plot](figure/mallowscv_plot-1.pdf)

```r
# cv.loglik (perplexity; higher is better; without kernelGaussian)
ggplot(subset(mallowsstatsCV, !method %in% keylistGS), aes(x = nclust, y = cv.loglik)) + 
    geom_boxplot(aes(fill = method), alpha = 0.5) + 
    stat_summary(fun.y = mean, geom = "line", aes(group = method, colour = method))
```

![plot of chunk mallowscv_plot](figure/mallowscv_plot-2.pdf)

## Reproduce figures from paper in section clusterAPA


```r
# set seed
set.seed(19206992)

# reduce contents to include in plots
klist_red <- seq(2, 10, 2)
keylistML_red <- c("bruteMallows", "copelandMallows", "bordaMallows", "kernelMallows_Exh", "kernelGaussian")

# kmeans time
ptmkmeans <- list()
num <- 1
for (key in keylistKM) {
  for (k in klist_red) {
    family <- get(paste0("family_", key))
    ptm <- sapply(1:nrepeats, function(i)
      system.time(cl <- kcca_modified(x = dat, k = k, family = family, weights = freq, save.data = FALSE))[1]
    )
    ptmkmeans[[num]] <- data.frame(nclust = k, method = key, time = mean(ptm))
    num <- num + 1
  }
}
ptmkmeans <- do.call("rbind", ptmkmeans)
ptmkmeans$nclust <- ordered(ptmkmeans$nclust)
ptmkmeans$method <- ordered(ptmkmeans$method)
ggplot(ptmkmeans, aes(x = nclust, y = time, colour = method)) + 
  geom_line(aes(group = method), size = 1.5) + 
  geom_point(aes(shape = method), size = 5) + 
  theme_bw() + 
  theme(legend.justification = c(0,1), legend.position = c(0,1), legend.title = element_blank(), 
        legend.key.size = unit(0.6, "in"), legend.background = element_rect(alpha("white", 0)),
        legend.text = element_text(size = 36, colour = "black", face = "bold"), 
        axis.text = element_text(size = 36, colour = "black", face = "bold"), 
        axis.title = element_text(size = 36, colour = "black", face = "bold"), 
        plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        panel.border = element_rect(colour = "black", size = 1))
```

![plot of chunk paper_plot](figure/paper_plot-1.pdf)

```r
# kmeans silhouette
ggplot(data = subset(clstats_best, nclust %in% klist_red), aes(x = nclust, y = silhouette, colour = method)) + 
  geom_line(aes(group = method), size = 1.5) + 
  geom_point(aes(shape = method), size = 5) + 
  theme_bw() + 
  theme(legend.justification = c(0,1), legend.position = c(0,1), legend.title = element_blank(), 
        legend.key.size = unit(0.6, "in"), legend.background = element_rect(alpha("white", 0)),
        legend.text = element_text(size = 36, colour = "black", face = "bold"), 
        axis.text = element_text(size = 36, colour = "black", face = "bold"), 
        axis.title = element_text(size = 36, colour = "black", face = "bold"), 
        plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        panel.border = element_rect(colour = "black", size = 1))
```

![plot of chunk paper_plot](figure/paper_plot-2.pdf)

```r
# mallows silhouette
d <- subset(mallowsstats_best, nclust %in% klist_red & method %in% keylistML_red)
stopifnot(!("kernelMallows" %in% d$method) && ("kernelMallows" %in% levels(d$method)))
d$method[d$method == "kernelMallows_Exh"] <- "kernelMallows"
ggplot(data = d, aes(x = nclust, y = silhouette, colour = method)) + 
  geom_line(aes(group = method), size = 1.5) + 
  geom_point(aes(shape = method), size = 5) + 
  theme_bw() + 
  theme(legend.justification = c(1,0.85), legend.position = c(1,0.85), legend.title = element_blank(), 
        legend.key.size = unit(0.6, "in"), legend.background = element_rect(alpha("white", 0)),
        legend.text = element_text(size = 36, colour = "black", face = "bold"), 
        axis.text = element_text(size = 36, colour = "black", face = "bold"), 
        axis.title = element_text(size = 36, colour = "black", face = "bold"), 
        plot.margin = unit(c(1, 1, 1, 1), "lines"), 
        panel.border = element_rect(colour = "black", size = 1))
```

![plot of chunk paper_plot](figure/paper_plot-3.pdf)

## Reproduce figures from supp


```r
ggplot(bclstats, aes(x = nclust, y = rand, fill = method)) + 
    geom_boxplot() + ylim(0, 0.75) + 
    theme_bw() + 
    theme(legend.justification = c(1,0), legend.position = c(1,0), legend.title = element_blank(), 
          legend.key.size = unit(0.6, "in"), legend.background = element_rect(alpha("white", 0)),
          legend.text = element_text(size = 36, colour = "black", face = "bold"), 
          axis.text = element_text(size = 36, colour = "black", face = "bold"), 
          axis.title = element_text(size = 36, colour = "black", face = "bold"), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"), 
          panel.border = element_rect(colour = "black", size = 1))
```

![plot of chunk supp_plot](figure/supp_plot-1.pdf)

## session info


```r
devtools::session_info()
```

```
## Session info --------------------------------------------------------------
```

```
##  setting  value                       
##  version  R version 3.2.3 (2015-12-10)
##  system   x86_64, linux-gnu           
##  ui       X11                         
##  language (EN)                        
##  collate  en_US.UTF-8                 
##  tz       <NA>                        
##  date     2016-09-11
```

```
## Packages ------------------------------------------------------------------
```

```
##  package      * version date       source        
##  bitops         1.0-6   2013-08-17 CRAN (R 3.2.2)
##  car            2.1-0   2015-09-03 CRAN (R 3.2.2)
##  caret        * 6.0-62  2015-11-23 CRAN (R 3.2.2)
##  caTools        1.17.1  2014-09-10 CRAN (R 3.2.2)
##  codetools      0.2-14  2015-07-15 CRAN (R 3.2.2)
##  colorspace     1.2-6   2015-03-11 CRAN (R 3.2.2)
##  combinat       0.0-8   2012-10-29 CRAN (R 3.2.2)
##  devtools       1.9.1   2015-09-11 CRAN (R 3.2.2)
##  digest         0.6.8   2014-12-31 CRAN (R 3.2.2)
##  evaluate       0.8     2015-09-18 CRAN (R 3.2.2)
##  flexclust    * 1.3-4   2013-07-02 CRAN (R 3.2.2)
##  foreach        1.4.3   2015-10-13 CRAN (R 3.2.2)
##  formatR        1.2.1   2015-09-18 CRAN (R 3.2.2)
##  gdata          2.17.0  2015-07-04 CRAN (R 3.2.2)
##  ggplot2      * 2.1.0   2016-03-01 CRAN (R 3.2.3)
##  gplots       * 2.17.0  2015-05-02 CRAN (R 3.2.2)
##  gtable         0.1.2   2012-12-05 CRAN (R 3.2.2)
##  gtools         3.5.0   2015-05-29 CRAN (R 3.2.2)
##  iterators      1.0.8   2015-10-13 CRAN (R 3.2.2)
##  kernrank     * 1.0.2   2016-04-28 local         
##  KernSmooth     2.23-15 2015-06-29 CRAN (R 3.2.2)
##  knitr          1.12.3  2016-01-22 CRAN (R 3.2.3)
##  labeling       0.3     2014-08-23 CRAN (R 3.2.2)
##  lattice      * 0.20-33 2015-07-14 CRAN (R 3.2.2)
##  lme4           1.1-10  2015-10-06 CRAN (R 3.2.2)
##  magrittr       1.5     2014-11-22 CRAN (R 3.2.2)
##  MASS           7.3-43  2015-07-16 CRAN (R 3.2.2)
##  Matrix         1.2-2   2015-07-08 CRAN (R 3.2.2)
##  MatrixModels   0.4-1   2015-08-22 CRAN (R 3.2.2)
##  memoise        0.2.1   2014-04-22 CRAN (R 3.2.2)
##  mgcv           1.8-7   2015-07-23 CRAN (R 3.2.2)
##  minqa          1.2.4   2014-10-09 CRAN (R 3.2.2)
##  modeltools   * 0.2-21  2013-09-02 CRAN (R 3.2.2)
##  munsell        0.4.2   2013-07-11 CRAN (R 3.2.2)
##  mvtnorm      * 1.0-3   2015-07-22 CRAN (R 3.2.2)
##  nlme           3.1-121 2015-06-29 CRAN (R 3.2.2)
##  nloptr         1.0.4   2014-08-04 CRAN (R 3.2.2)
##  nnet           7.3-10  2015-06-29 CRAN (R 3.2.2)
##  pbkrtest       0.4-2   2014-11-13 CRAN (R 3.2.2)
##  plyr           1.8.3   2015-06-12 CRAN (R 3.2.2)
##  quantreg       5.19    2015-08-31 CRAN (R 3.2.2)
##  Rankcluster  * 0.92.9  2014-07-25 CRAN (R 3.2.2)
##  Rcpp           0.12.2  2015-11-15 CRAN (R 3.2.2)
##  reshape2     * 1.4.1   2014-12-06 CRAN (R 3.2.2)
##  scales         0.3.0   2015-08-25 CRAN (R 3.2.2)
##  SparseM        1.7     2015-08-15 CRAN (R 3.2.2)
##  stringi        1.0-1   2015-10-22 CRAN (R 3.2.2)
##  stringr        1.0.0   2015-04-30 CRAN (R 3.2.2)
```
