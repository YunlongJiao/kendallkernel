---
title: "kendallkernel - cluster (multivariate) partial rankings"
author: "Yunlong Jiao"
date: "16 Sep 2015"
output: html_document
---


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, fig.width = 15, fig.height = 12, dev = "pdf", fig.keep = "high", fig.path = "figure/", cache.path = "cache/")
set.seed(94151402)

# utiles
library(kernrank) # fast kendall kernel
library(combinat)
library(reshape2)
library(ggplot2)
library(parallel)

# clustering
library(kernlab)
library(Rankcluster)
library(cluster)
source("kkmeans.R") # complementary functions for kernel kmeans for top-k rankings
```

## Multivariate partial rankings

Each instance $x=(x^1,\dots,x^p)$ is represented by a multivariate rank in which $x^j=(x^{j1},\dots,x^{jm_j})$ for $1\leq j \leq p$ is a (partial) rank of $m_j$ items. For example, foot index $j$ should indicate the ranking preference of the same customer over several years or on several incomparable item sets, in which case we would like to exploit information in parallel when performing clustering. The `Rankcluster::eurovision` data is a typical example of such data.

- **Mixture of multivariate isr** Implemented by `Rankcluster::rankclust` with mostly defaulted parameters (except repeat times), see [Jacques and Biernacki (2014)](https://hal.inria.fr/hal-00743384/). This mixture model is originally derived to deal with heterogeneous multivariate (partial) ranks via modeling the voting process using ISR model. NOTE parameter setting is inherited from Section 4 _Numerical experiments_ of the paper but repeated `run` times. However, we notice that the naive `set.seed` does not guarantee reproducible results here!

- **Kendall kernel kmeans** Implemented by `kernlab::kkmeans` with defaulted parameters where kernel similarity values are kendall kernel for partial rankings further averaged over multivariates present. Specificallly, this method is model-free where similarities of partial ranks are adopted by convolution Kendall kernel and multivariate information is utilized by a simple average of kernel values from variates in parallel (NOTE other multiple kernel learning approaches apply).


```r
# load data - eurovision song contest data
data(eurovision)
dat <- eurovision$data
countries <- rownames(dat)
nrow(dat) # number of countries
```

```
## [1] 34
```

```r
dimsize <- eurovision$m
dimsize # dimension of multivariate partial rankings listed in format c(m_1,m_2,...,m_p)
```

```
## [1] 8 8 8 8 8 8
```

```r
# number of clusters
klist <- 2:6
```


```r
key <- "mpr_rankcluster"

# set param
RjSE <- RjM <- rep(10, length(dimsize))
Qsem <- Ql <- 100
Bsem <- Bl <- 10
maxTry <- 20
run <- 10

# get clustering results
ptm <- proc.time()
set.seed(94121899)
assign(key, 
       rankclust(data = dat, m = dimsize, K = klist, criterion = "bic", 
                 Qsem = Qsem, Bsem = Bsem, RjSE = RjSE, RjM = RjM, Ql = Ql, Bl = Bl,
                 maxTry = maxTry, run = run))
proc.time() - ptm
```

```
##     user   system  elapsed 
## 6897.737    0.018 6898.851
```


```r
key <- "mpr_kkmeans"

# set param
nrepeats <- 100

# get clustering results
ptm <- proc.time()
kmatrix <- kernmat(dat = dat, dimsize = dimsize, kf = "kendall_top", procf = "procf_Rankcluster")
```

```
## pre-process dat ...
```

```
## computing kernel matrix ...
```

```
## done!
```

```r
dmatrix <- distmat(kmatrix)
assign(key, kkmeans_multi(kmatrix = kmatrix, klist = klist, nrepeats = nrepeats, dmatrix = dmatrix, seed = 31286963))
proc.time() - ptm
```

```
##    user  system elapsed 
##  16.679   0.037  16.799
```

## Results


```r
key <- "mpr_rankcluster"
cat("-----------------> \t", key, "\t <-----------------")

# get clustering statistics
clstats_best <- lapply(seq(length(get(key)@K)), function(i){
  data.frame(nclust = get(key)@K[i], 
             loglike = get(key)@results[[i]]@ll, 
             bic = -(get(key)@results[[i]]@bic), 
             icl = -(get(key)@results[[i]]@icl))
})
clstats_best <- do.call("rbind", clstats_best)
clstats_best$nclust <- ordered(clstats_best$nclust)

# loglike (higher is better)
ggplot(clstats_best, aes(x = nclust, y = loglike)) + 
  geom_line(aes(group = 1)) + 
  geom_point(size = 3)
```

![plot of chunk mpr_rankcluster_analysis](figure/mpr_rankcluster_analysis-1.pdf)

```r
# bic (higher is better)
ggplot(clstats_best, aes(x = nclust, y = bic)) + 
  geom_line(aes(group = 1)) + 
  geom_point(size = 3)
```

![plot of chunk mpr_rankcluster_analysis](figure/mpr_rankcluster_analysis-2.pdf)

```r
# icl (higher is better)
ggplot(clstats_best, aes(x = nclust, y = icl)) + 
  geom_line(aes(group = 1)) + 
  geom_point(size = 3)
```

![plot of chunk mpr_rankcluster_analysis](figure/mpr_rankcluster_analysis-3.pdf)

```r
# estimated clustering by rankclust (nclust corresp to largest bic)
i_best <- which.max(clstats_best$bic)
cl_best <- get(key)@results[[i_best]]

# summary of clustering
cest <- cl_best@partition
table(cest)
split(countries, cest)

# silhouette plot for best fit
sil <- cluster::silhouette(x = cest, dmatrix = dmatrix)
rownames(sil) <- countries
par(font = 2, font.axis = 2, font.lab = 2, font.main = 2, font.sub = 2, 
    cex = 1.25, cex.axis = 1.25, cex.lab = 2, cex.main = 2, cex.sub = 2, 
    mar = c(4, 17, 0, 2) + 0.1)
c6 <- c("tomato", "light blue", "forest green", "purple2", "goldenrod4", "gray20")
plot.silhouette_modified(sil, max.strlen = max(nchar(countries)), col = c6[1:max(cest)], cex.names = 1.25,
                         do.n.k = FALSE, do.clus.stat = FALSE)
```

![plot of chunk mpr_rankcluster_analysis](figure/mpr_rankcluster_analysis-4.pdf)

```
## -----------------> 	 mpr_rankcluster 	 <-----------------cest
##  1  2  3  4  5 
## 11  3  6 13  1 
## $`1`
##  [1] "Belarus"          "Bulgaria"         "Cyprus"          
##  [4] "F.Y.R. Macedonia" "Germany"          "Moldova"         
##  [7] "Romania"          "Russia"           "Serbia"          
## [10] "Turkey"           "Ukraine"         
## 
## $`2`
## [1] "France" "Greece" "Israel"
## 
## $`3`
## [1] "Albania"              "Belgium"              "Bosnia & Herzegovina"
## [4] "Finland"              "Switzerland"          "United Kingdom"      
## 
## $`4`
##  [1] "Croatia"         "Denmark"         "Estonia"        
##  [4] "Iceland"         "Ireland"         "Latvia"         
##  [7] "Lithuania"       "Malta"           "Norway"         
## [10] "Slovenia"        "Spain"           "Sweden"         
## [13] "The Netherlands"
## 
## $`5`
## [1] "Portugal"
```


```r
key <- "mpr_kkmeans"
cat("-----------------> \t", key, "\t <-----------------")

# get clustering statistics
clstats <- kkmeans_stats(get(key))
clstats_best <- kkmeans_stats_best(clstats)

# distsum (lower is better)
ggplot() + 
    geom_boxplot(data = clstats, aes(x = nclust, y = distsum), alpha = 0.5) + 
    geom_line(data = clstats_best, aes(x = nclust, y = distsum, group = 1)) + 
    geom_point(data = clstats_best, aes(x = nclust, y = distsum), size = 3)
```

![plot of chunk mpr_kkmeans_analysis](figure/mpr_kkmeans_analysis-1.pdf)

```r
# silhouette (higher is better)
ggplot() + 
    geom_boxplot(data = clstats, aes(x = nclust, y = silhouette), alpha = 0.5) + 
    geom_line(data = clstats_best, aes(x = nclust, y = silhouette, group = 1)) + 
    geom_point(data = clstats_best, aes(x = nclust, y = silhouette), size = 3)
```

![plot of chunk mpr_kkmeans_analysis](figure/mpr_kkmeans_analysis-2.pdf)

```r
# estimated clustering by kkmeans (nclust corresp to highest silhouette)
i_best <- clstats_best[which.max(clstats_best$silhouette), ]
cl_best <- get(key)[[as.character(i_best$nclust)]][[i_best$rep]]

# summary of clustering
cest <- cl_best@.Data
table(cest)
split(countries, cest)

# silhouette plot for best fit
sil <- cluster::silhouette(x = cest, dmatrix = dmatrix)
rownames(sil) <- countries
par(font = 2, font.axis = 2, font.lab = 2, font.main = 2, font.sub = 2, 
    cex = 1.25, cex.axis = 1.25, cex.lab = 2, cex.main = 2, cex.sub = 2, 
    mar = c(4, 17, 0, 2) + 0.1)
c6 <- c("tomato", "light blue", "forest green", "purple2", "goldenrod4", "gray20")
plot.silhouette_modified(sil, max.strlen = max(nchar(countries)), col = c6[1:max(cest)], cex.names = 1.25,
                         do.n.k = FALSE, do.clus.stat = FALSE)
```

![plot of chunk mpr_kkmeans_analysis](figure/mpr_kkmeans_analysis-3.pdf)

```
## -----------------> 	 mpr_kkmeans 	 <-----------------cest
##  1  2 
## 20 14 
## $`1`
##  [1] "Belarus"          "Croatia"          "Denmark"         
##  [4] "Estonia"          "F.Y.R. Macedonia" "Finland"         
##  [7] "Greece"           "Ireland"          "Israel"          
## [10] "Latvia"           "Lithuania"        "Malta"           
## [13] "Moldova"          "Norway"           "Portugal"        
## [16] "Russia"           "Slovenia"         "Spain"           
## [19] "Turkey"           "Ukraine"         
## 
## $`2`
##  [1] "Albania"              "Belgium"              "Bosnia & Herzegovina"
##  [4] "Bulgaria"             "Cyprus"               "France"              
##  [7] "Germany"              "Iceland"              "Romania"             
## [10] "Serbia"               "Sweden"               "Switzerland"         
## [13] "The Netherlands"      "United Kingdom"
```

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
##  date     2016-08-11
```

```
## Packages ------------------------------------------------------------------
```

```
##  package     * version date       source        
##  cluster     * 2.0.3   2015-07-21 CRAN (R 3.2.2)
##  colorspace    1.2-6   2015-03-11 CRAN (R 3.2.2)
##  combinat    * 0.0-8   2012-10-29 CRAN (R 3.2.2)
##  devtools      1.9.1   2015-09-11 CRAN (R 3.2.2)
##  digest        0.6.8   2014-12-31 CRAN (R 3.2.2)
##  evaluate      0.8     2015-09-18 CRAN (R 3.2.2)
##  formatR       1.2.1   2015-09-18 CRAN (R 3.2.2)
##  ggplot2     * 2.1.0   2016-03-01 CRAN (R 3.2.3)
##  gtable        0.1.2   2012-12-05 CRAN (R 3.2.2)
##  kernlab     * 0.9-22  2015-08-05 CRAN (R 3.2.2)
##  kernrank    * 1.0.2   2016-04-28 local         
##  knitr         1.12.3  2016-01-22 CRAN (R 3.2.3)
##  labeling      0.3     2014-08-23 CRAN (R 3.2.2)
##  magrittr      1.5     2014-11-22 CRAN (R 3.2.2)
##  memoise       0.2.1   2014-04-22 CRAN (R 3.2.2)
##  munsell       0.4.2   2013-07-11 CRAN (R 3.2.2)
##  plyr          1.8.3   2015-06-12 CRAN (R 3.2.2)
##  Rankcluster * 0.92.9  2014-07-25 CRAN (R 3.2.2)
##  Rcpp          0.12.2  2015-11-15 CRAN (R 3.2.2)
##  reshape2    * 1.4.1   2014-12-06 CRAN (R 3.2.2)
##  scales        0.3.0   2015-08-25 CRAN (R 3.2.2)
##  stringi       1.0-1   2015-10-22 CRAN (R 3.2.2)
##  stringr       1.0.0   2015-04-30 CRAN (R 3.2.2)
```
