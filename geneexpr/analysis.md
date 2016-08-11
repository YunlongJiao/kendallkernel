---
title: "kendallkernel - classify gene expression data"
author: "Yunlong Jiao"
date: "26 February 2016"
output: html_document
---


```r
knitr::opts_chunk$set(error = FALSE, warning = FALSE, fig.width = 8, fig.height = 8, dev = "pdf", fig.keep = "high", fig.path = "figure/", cache.path = "cache/")
set.seed(35875954)

# utiles
library(kernlab) # for kernel svm
library(pcaPP) # for fast kendall tau
library(caret) # for data split
library(parallel) # for parallel cross-validation
source("func.R")
dyn.load("src/tsp.so") # for tsp-related C codes
dyn.load("src/utiles.so") # for other C codes
```

## Datasets

All datasets are taken from publicly available sources. Briefly, for datasets with two independent parts (marked `indepval`), predictions are reported on the test set whilst 5-fold cv training is done on training set; for datasets with only one part (marked `crossval`), predictions are reported by 10 times of 5-fold cv whilst (nested) 5-fold cv is done on each training fold for parameter tuning. See [the paper](https://hal.archives-ouvertes.fr/hal-01279273) for detailed summary of the datasets.


```r
# datasets
indepvallist <- c('bc', 'lcbwh', 'pc') # two sets are available
crossvallist <- c('ct', 'ocpbsii', 'cns', 'pcout', 'transbig', 'lung_annarbor_outcome', 'Wang_Breastcancer') # only one set is available so that double-nested cv is necessary
prefixlist <- c(indepvallist, crossvallist)

# read in 10 datasets from folder data/
fnames <- list.files(path = "data/")
for (fname in fnames) load(paste0("data/", fname))

# dataset alias as they appear in the paper
cbind(dataset = (namebefore <- c("bc", "transbig", "Wang_Breastcancer", "ct", "lung_annarbor_outcome", 
                                 "lcbwh", "cns", "ocpbsii", "pc", "pcout")), 
      alias = (nameafter <- c("BC1", "BC2", "BC3", "CT", "LA1", "LC2", "MB", "OC", "PC1", "PC2")))
```

```
##       dataset                 alias
##  [1,] "bc"                    "BC1"
##  [2,] "transbig"              "BC2"
##  [3,] "Wang_Breastcancer"     "BC3"
##  [4,] "ct"                    "CT" 
##  [5,] "lung_annarbor_outcome" "LA1"
##  [6,] "lcbwh"                 "LC2"
##  [7,] "cns"                   "MB" 
##  [8,] "ocpbsii"               "OC" 
##  [9,] "pc"                    "PC1"
## [10,] "pcout"                 "PC2"
```

## Model performance comparison

Models come from 3 categories (presented in different ways for ease of coding implementation or for ease of scoring and plotting):

1. A baseline model with no tuning parameter that is all-pairs-majority-votes (or n(n-1)/2-TSP)
2. Models involving only tuning C that are SVM with linear, Gaussian RBF, (2nd-order homogeneous) polynomial and Kendall kernel where KFD (Kernel Fisher Discriminant) are penetrated in SVM codes as simple reference kernel machines
3. Models involving tuning C and k that are SVM with top-k pairs of features with aforementioned kernels


```r
# set list of C parameters for SVM-based models
Cpara_list <- 10^(-2:3)
names(Cpara_list) <- paste('C',1:length(Cpara_list),sep='')

# set list of #genepairs for corresponding models
max_nodes = 5000; npairs_out = 30;
npairs_list <- floor(exp(seq(0,1,length.out=npairs_out)*log(max_nodes)))
evenidx <- npairs_list %% 2 == 0
npairs_list[evenidx] <- npairs_list[evenidx] - 1 # keep odd numbers only
npairs_list <- unique(npairs_list)
names(npairs_list) <- paste('k',1:length(npairs_list),sep='')

# categorize models for ease of training
modelsNOpara <- c("APMV")
modelsConly <- c("SVMlinearALL", "SVMkdtALL", "SVMpolynomialALL", "SVMrbf") # plus KFD coded within each
modelsTSPrelated <- c("TSP", "kTSP", "SVMlinearTOP", "SVMkdtTOP", "SVMpolynomialTOP")
# OR reorganise for ease of plotting
models0 <- c("TSP", "APMV")
modelsConly <- c("SVMlinearALL", "SVMkdtALL", "SVMpolynomialALL", "SVMrbf") # same as before!
modelsKonly <- c("kTSP")
modelsCandK <- c("SVMlinearTOP", "SVMkdtTOP", "SVMpolynomialTOP")
modelslist <- c(models0, modelsConly, modelsKonly, modelsCandK)
# OR reorganise for feature selection plot varying K
modelsVary <- c("SVMlinearTOP", "SVMkdtTOP", "SVMpolynomialTOP", "kTSP")
modelsStatic <- c("SVMlinearALL", "SVMkdtALL", "SVMpolynomialALL", "TSP", "APMV")
```


```r
# indepval datasets
res_indepval <- mclapply(indepvallist, function(prefixname){
  xtr <- get(prefixname)$xtrain; ytr <- get(prefixname)$ytrain
  xtst <- get(prefixname)$xtest; ytst <- get(prefixname)$ytest
  if(is.null(xtst) || is.null(ytst)) stop(paste('dataset error',prefixname,sep=':'))
  
  set.seed(206)
  res <- perfClassification(NULL, prefixname, xtr, ytr, xtst, ytst,
                            Cpara_list, npairs_list, modelsConly, modelsTSPrelated, modelsNOpara,
                            nfolds = 5, nrepeats = 1, seed = 206)
  return(res)
}, mc.cores = 8)
names(res_indepval) <- indepvallist
```


```r
# crossval datasets
res_crossval <- mclapply(crossvallist, function(prefixname){
  xtr <- get(prefixname)$xtrain; ytr <- get(prefixname)$ytrain
  xtst <- get(prefixname)$xtest; ytst <- get(prefixname)$ytest
  if(!is.null(xtst) || !is.null(ytst)) stop(paste('dataset error',prefixname,sep=':'))
  
  set.seed(1226)
  outterFoldIndices <- createMultiFolds(1:nrow(xtr), k=5, times=10)
  sig <- sigest(xtr,scaled=F)['50%']
  
  res <- lapply(outterFoldIndices, function(outterFold){
    return(perfClassification(NULL, prefixname, xtr[outterFold,,drop=F], ytr[outterFold], xtr[-outterFold,,drop=F], ytr[-outterFold], 
                              Cpara_list, npairs_list, modelsConly, modelsTSPrelated, modelsNOpara, 
                              nfolds = 5, nrepeats = 1, seed = 206, sigma=sig))
  })
  return(res)
}, mc.cores = 8)
names(res_crossval) <- crossvallist
```

We report classification accuracy across different datasets and different models.


```r
modelsKFD <- sub("SVM", "KFD", modelsConly)
table_acc <- matrix(-100, 
                    nrow = length(prefixlist), ncol = length(c(modelslist,modelsKFD)),
                    dimnames = list(prefixlist, c(modelslist,modelsKFD)))

for (prefixname in prefixlist) {
  for (modelname in modelslist) {
    if (prefixname %in% indepvallist) {
      res <- res_indepval[[prefixname]]
      idx <- which.max(res[[modelname]]$cvacc)
      s <- res[[modelname]]$acc[idx]
      table_acc[prefixname,modelname] <- round(100*s,2)
      if (modelname %in% modelsConly) { # add KFD penetrated within
        s_kfd <- res[[modelname]]$acc_kfd
        table_acc[prefixname,sub("SVM", "KFD", modelname)] <- round(100*s_kfd,2)
      }
    } else if (prefixname %in% crossvallist) {
      s <- mean(sapply(res_crossval[[prefixname]], function(res){
        idx <- which.max(res[[modelname]]$cvacc)
        return(res[[modelname]]$acc[idx])
      }))
      table_acc[prefixname,modelname] <- round(100*s,2)
      if (modelname %in% modelsConly) { # add KFD penetrated within
        s_kfd <- mean(sapply(res_crossval[[prefixname]], function(res){
          return(res[[modelname]]$acc_kfd)
        }))
        table_acc[prefixname,sub("SVM", "KFD", modelname)] <- round(100*s_kfd,2)
      }
    } else {
      stop("Please add ", prefixname, " in either indepvallist or crossvallist")
    }
  }
}
rownames(table_acc) <- nameafter[match(rownames(table_acc), namebefore)] # re-name
table_acc <- table_acc[order(rownames(table_acc)), ] # re-order
table_acc <- rbind(AVERAGE = round(colMeans(table_acc), 2), table_acc) # add AVERAGE scores over all datasets
table_acc <- table_acc[ ,order(table_acc["AVERAGE",],decreasing = TRUE)] # re-order
# show score table
t(table_acc)
```

```
##                  AVERAGE   BC1   BC2   BC3    CT   LA1   LC2    MB     OC
## SVMkdtALL          79.39 78.95 71.31 67.34 85.78 70.98 97.99 63.67  99.48
## SVMlinearTOP       77.16 84.21 69.29 67.11 84.19 63.92 97.32 65.17  99.41
## SVMlinearALL       76.09 78.95 71.67 64.27 86.73 70.23 97.99 62.67  99.64
## SVMkdtTOP          75.50 52.63 70.61 65.81 85.46 67.70 97.99 58.33  99.92
## SVMpolynomialALL   74.54 68.42 71.62 63.66 78.43 70.53 98.66 61.17  99.28
## KFDkdtALL          74.33 63.16 59.41 67.22 85.46 59.08 99.33 59.33  98.73
## kTSP               74.03 57.89 58.22 64.47 87.23 61.70 97.99 56.00  99.92
## SVMpolynomialTOP   73.99 63.16 69.44 66.26 79.14 65.98 99.33 60.00  99.21
## KFDlinearALL       71.81 63.16 60.43 67.52 77.26 57.24 97.99 59.50 100.00
## KFDpolynomialALL   71.39 63.16 60.48 67.38 75.10 58.52 97.99 60.33 100.00
## TSP                69.71 68.42 49.58 57.80 85.61 58.96 95.97 52.67  99.80
## SVMrbf             69.31 63.16 71.41 65.87 81.18 70.84 93.96 63.83  98.85
## KFDrbf             66.50 63.16 60.38 66.17 84.33 58.62 97.32 60.17  98.34
## APMV               61.91 84.21 65.98 33.96 64.49 33.60 89.93 42.17  85.19
##                     PC1   PC2
## SVMkdtALL        100.00 58.40
## SVMlinearTOP      85.29 55.70
## SVMlinearALL      73.53 55.17
## SVMkdtTOP         97.06 59.47
## SVMpolynomialALL  79.41 54.23
## KFDkdtALL         97.06 54.57
## kTSP             100.00 56.83
## SVMpolynomialTOP  88.24 49.10
## KFDlinearALL      73.53 61.43
## KFDpolynomialALL  73.53 57.43
## TSP               76.47 51.83
## SVMrbf            26.47 57.50
## KFDrbf            26.47 50.00
## APMV              73.53 46.00
```

```r
# show boxplot
par(mar = c(10, 5, 1, 1) + 0.1, font.lab = 2, font.axis = 2, font.main = 2, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5)
boxplot(table_acc[-1, ]/100, las = 2, ylab = 'acc', col='royalblue2')
```

![plot of chunk perf_table](figure/perf_table-1.pdf)

```r
# wilcox test
pmatrix <- matrix(1, nrow = ncol(table_acc), ncol = ncol(table_acc), 
                  dimnames = list(colnames(table_acc),colnames(table_acc)))
for(i in 1:ncol(table_acc)){
  for(j in 1:ncol(table_acc)){
    pmatrix[i,j] = round(wilcox.test(x = table_acc[-1,i], y = table_acc[-1,j], alternative = "greater", paired = TRUE, mu = 0)$p.value, 4)
  }
}
pmatrix
```

```
##                  SVMkdtALL SVMlinearTOP SVMlinearALL SVMkdtTOP
## SVMkdtALL           1.0000       0.0654       0.0707    0.0290
## SVMlinearTOP        0.9473       1.0000       0.3477    0.7217
## SVMlinearALL        0.9463       0.6875       1.0000    0.3611
## SVMkdtTOP           0.9780       0.3125       0.6822    1.0000
## SVMpolynomialALL    0.9902       0.9033       0.9033    0.6152
## KFDkdtALL           0.9928       0.8125       0.8125    0.6880
## kTSP                0.9850       0.7539       0.6822    0.7794
## SVMpolynomialTOP    0.9971       0.8838       0.9199    0.7842
## KFDlinearALL        0.9780       0.9580       0.9293    0.7232
## KFDpolynomialALL    0.9911       0.9678       0.9463    0.8568
## TSP                 0.9980       0.9951       0.9932    0.9678
## SVMrbf              0.9951       0.8623       0.8623    0.6152
## KFDrbf              1.0000       0.9954       0.9971    0.9033
## APMV                0.9990       0.9968       0.9954    0.9814
##                  SVMpolynomialALL KFDkdtALL   kTSP SVMpolynomialTOP
## SVMkdtALL                  0.0137    0.0095 0.0212           0.0049
## SVMlinearTOP               0.1162    0.2158 0.2783           0.1377
## SVMlinearALL               0.1162    0.2158 0.3611           0.0967
## SVMkdtTOP                  0.4229    0.3631 0.2643           0.2461
## SVMpolynomialALL           1.0000    0.4609 0.3848           0.2783
## KFDkdtALL                  0.5771    1.0000 0.3994           0.4721
## kTSP                       0.6523    0.6394 1.0000           0.5771
## SVMpolynomialTOP           0.7539    0.5832 0.4609           1.0000
## KFDlinearALL               0.9033    0.7614 0.5000           0.8819
## KFDpolynomialALL           0.9473    0.6389 0.5000           0.8819
## TSP                        0.9710    0.9580 0.9756           0.9473
## SVMrbf                     0.6523    0.4528 0.3125           0.3611
## KFDrbf                     0.9580    0.9710 0.7217           0.9037
## APMV                       0.9863    0.9814 0.9814           0.9814
##                  KFDlinearALL KFDpolynomialALL    TSP SVMrbf KFDrbf   APMV
## SVMkdtALL              0.0290           0.0122 0.0029 0.0068 0.0010 0.0020
## SVMlinearTOP           0.0527           0.0420 0.0068 0.1611 0.0064 0.0046
## SVMlinearALL           0.0917           0.0707 0.0098 0.1611 0.0049 0.0064
## SVMkdtTOP              0.3178           0.1716 0.0420 0.4229 0.1162 0.0244
## SVMpolynomialALL       0.1162           0.0654 0.0378 0.3848 0.0527 0.0186
## KFDkdtALL              0.2768           0.4064 0.0527 0.5936 0.0378 0.0244
## kTSP                   0.5472           0.5472 0.0322 0.7217 0.3125 0.0244
## SVMpolynomialTOP       0.1432           0.1432 0.0654 0.6822 0.1181 0.0244
## KFDlinearALL           1.0000           0.3375 0.1875 0.5936 0.1869 0.0486
## KFDpolynomialALL       0.7353           1.0000 0.2461 0.6822 0.0691 0.0486
## TSP                    0.8389           0.7842 1.0000 0.7842 0.5000 0.0801
## SVMrbf                 0.4528           0.3611 0.2461 1.0000 0.0917 0.1377
## KFDrbf                 0.8432           0.9453 0.5391 0.9293 1.0000 0.2158
## APMV                   0.9622           0.9622 0.9346 0.8838 0.8125 1.0000
```

Among SVM-based models we compare the tuning sensitivity to C parameter that could impact on the performance in terms of classification accuracy, where KFD are set as reference baseline models.


```r
nConly <- length(modelsConly) # number of SVM models (KFD implemented within each)
for (prefixname in prefixlist) {
  mainname <- nameafter[which(prefixname == namebefore)]
  key <- (prefixname %in% indepvallist)
  
  if (key) {
    res <- res_indepval[[prefixname]]
    s <- sapply(res[modelsConly],function(u){c(u$acc_kfd,u$acc)}) # *** single value of acc_kfd followed by a vector of acc is necessary for implementation below
  } else {
    res <- res_crossval[[prefixname]]
    s <- lapply(res, function(resfold){
      sapply(resfold[modelsConly],function(u){c(u$acc_kfd,u$acc)})
    })
    dm <- c(dim(s[[1]]),length(s)); dn <- dimnames(s[[1]]) # save info of dim and dimnames
    s <- unlist(s); dim(s) <- dm # reform to an 3d array
    s <- apply(s, c(1,2), mean); dimnames(s) <- dn # average over cv-folds
  }
  
  # set y range for plot
  plotrange <- c(floor(10*min(s,na.rm=T))/10,ceiling(10*max(s,na.rm=T))/10)
  
  # plotting
  par(font.lab = 2, font.axis = 2, font.main = 2, font = 2, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5)
  plot(Cpara_list, rep(-100,length(Cpara_list)), main=mainname, 
       xlab="C parameter", ylab=ifelse(key,"acc","cvacc"), 
       ylim=plotrange,type='l',lwd=2,lty=1,col=1,log='x')
  for (col in seq(nConly)) {
    modelname <- modelsConly[col]
    lines(Cpara_list,s[-1,modelname],type='l',lwd=2,lty=1,col=col) # *** see notes above
    points(Cpara_list,rep(s[1,modelname],length(Cpara_list)),type='b',lty=5,lwd=1,pch=col,col=col,cex=1)
    if (key) { # for indepval mark cv-tuned parameter
      idx <- which.max(res[[modelname]]$cvacc)
      score <- res[[modelname]]$acc
      points(Cpara_list[idx],score[idx],lwd=2,lty=1,pch=col,col=col,cex=3)
    }
  }
  molist <- c(modelsConly,sub("SVM","KFD",modelsConly))
  molist[grep('rbf',molist)] <- paste(molist[grep('rbf',molist)],'ALL',sep='')
  legend("bottomright",legend=molist,pch=c(rep(NA,nConly),seq(nConly)),
         col=c(seq(nConly),seq(nConly)),lty=c(rep(1,nConly),rep(5,nConly)),
         lwd=c(rep(2,nConly),rep(1,nConly)),cex=1.25)
  grid(ny=16)
}
```

![plot of chunk perf_C](figure/perf_C-1.pdf)![plot of chunk perf_C](figure/perf_C-2.pdf)![plot of chunk perf_C](figure/perf_C-3.pdf)![plot of chunk perf_C](figure/perf_C-4.pdf)![plot of chunk perf_C](figure/perf_C-5.pdf)![plot of chunk perf_C](figure/perf_C-6.pdf)![plot of chunk perf_C](figure/perf_C-7.pdf)![plot of chunk perf_C](figure/perf_C-8.pdf)![plot of chunk perf_C](figure/perf_C-9.pdf)![plot of chunk perf_C](figure/perf_C-10.pdf)

Now we study the impact of feature selection onto classification accuracy. TSP-based models serve as reference models.


```r
nVary <- length(modelsVary) # number of vary-K models
nStatic <- length(modelsStatic) # number of K-independent models
for (prefixname in prefixlist) {
  mainname <- nameafter[which(prefixname == namebefore)]
  key <- (prefixname %in% indepvallist)
  
  if (key) {
    res <- res_indepval[[prefixname]]
    # acc for K-independent models
    res_static <- lapply(modelsStatic, function(modelname){
      idx <- which.max(res[[modelname]]$cvacc)
      res[[modelname]]$acc[idx]
    })
    # acc for vary-K models
    res_vary <- lapply(modelsVary, function(modelname){
      cvscore <- res[[modelname]]$cvacc
      ivscore <- res[[modelname]]$acc
      if (!modelname %in% modelsCandK) {
        cvscore <- matrix(cvscore, nrow = 1)
        ivscore <- matrix(ivscore, nrow = 1)
      }
      tl <- max.col(t(cvscore), ties.method = 'first')
      idx <- which(cvscore == max(cvscore), arr.ind = TRUE)[1,"col"] # return INDEX of smallest k
      ivscore <- sapply(1:length(tl), function(u){ivscore[tl[u],u]}) # cv-ed out regarding C
      return(list(s = ivscore, k = idx))
    })
  } else {
    res <- res_crossval[[prefixname]]
    # acc for K-independent models
    res_static <- lapply(modelsStatic, function(modelname){
      mean(sapply(res, function(resfold){
        resfold[[modelname]]$acc[which.max(resfold[[modelname]]$cvacc)]
      }))
    })
    # acc for vary-K models
    res_vary <- lapply(modelsVary, function(modelname){
      return(list(s = rowMeans(sapply(res, function(resfold){
        cvscore <- resfold[[modelname]]$cvacc
        ivscore <- resfold[[modelname]]$acc
        if (!modelname %in% modelsCandK) {
          cvscore <- matrix(cvscore, nrow = 1)
          ivscore <- matrix(ivscore, nrow = 1)
        }
        tl <- max.col(t(cvscore), ties.method = 'first')
        ivscore <- sapply(1:length(tl), function(u){ivscore[tl[u],u]}) # cv-ed out regarding C
        return(ivscore)
      })), k = NA))
    })
  }
  names(res_static) <- modelsStatic
  names(res_vary) <- modelsVary
  
  # set y range for plot
  s <- c(unlist(res_static), unlist(lapply(res_vary,function(u){u$s})))
  plotrange <- c(floor(10*min(s,na.rm=T))/10,ceiling(10*max(s,na.rm=T))/10)
  
  # plotting
  par(font.lab = 2, font.axis = 2, font.main = 2, font = 2, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5)
  plot(npairs_list, rep(-100,length(npairs_list)), main=mainname,
       xlab="Number k of top gene pairs", ylab=ifelse(key,"acc","cvacc"), 
       ylim=plotrange,type='l',lwd=2,lty=1,col=1,log='x')
  
  for(col in seq(nVary)){
    modelname <- modelsVary[col]
    score <- res_vary[[modelname]]$s
    lines(npairs_list, score, type='l',lty=1,lwd=2,col=col)
    if (key) {
      idx <- res_vary[[modelname]]$k
      points(npairs_list[idx], score[idx], lwd=2,lty=1,pch=col,col=col,cex=3)
    }
  }
  
  for(col in seq(nStatic)){
    modelname <- modelsStatic[col]
    ref <- res_static[[modelname]]
    points(npairs_list, rep(ref,length(npairs_list)), type='b',lty=5,lwd=1,pch=col,col=col,cex=1)
  }
  
  molist <- c(modelsVary,modelsStatic)
  legend("bottomleft", legend=molist,pch=c(rep(NA,nVary),seq(nStatic)),
         col=c(seq(nVary),seq(nStatic)),lty=c(rep(1,nVary),rep(5,nStatic)),
         lwd=c(rep(2,nVary),rep(1,nStatic)),cex=1.25)
  grid(ny=16)
}
```

![plot of chunk perf_fs](figure/perf_fs-1.pdf)![plot of chunk perf_fs](figure/perf_fs-2.pdf)![plot of chunk perf_fs](figure/perf_fs-3.pdf)![plot of chunk perf_fs](figure/perf_fs-4.pdf)![plot of chunk perf_fs](figure/perf_fs-5.pdf)![plot of chunk perf_fs](figure/perf_fs-6.pdf)![plot of chunk perf_fs](figure/perf_fs-7.pdf)![plot of chunk perf_fs](figure/perf_fs-8.pdf)![plot of chunk perf_fs](figure/perf_fs-9.pdf)![plot of chunk perf_fs](figure/perf_fs-10.pdf)

## Kernel approximation study

Now we turn to study the effect of stablized Kendall kernel in place of the classic one solely based on rank orders. We only focus on the case of uniform noise varying window size.


```r
prefixlist <- c("cns") # only run on one dataset for now (more models can be added to the vector)
modelname <- "SVMkdtALLquadratic" # we only study kernel jittered with uniform noise

# number of selected window sizes then specific values are chosen differently for each dataset
nEXTpara_out <- 20

# number of sampled noise perturbation
nMaxNum <- 35
MaxNum_list <- seq(nMaxNum); names(MaxNum_list) <- paste('Num', 1:length(MaxNum_list), sep='')
```


```r
for(prefixname in prefixlist){
  xdat <- rbind(get(prefixname)$xtrain,get(prefixname)$xtest)
  ydat <- factor(c(get(prefixname)$ytrain, get(prefixname)$ytest), labels = levels(get(prefixname)$ytrain))
  
  set.seed(1226)
  # split data if no separate test set available
  if(is.null(get(prefixname)$xtest) || is.null(get(prefixname)$ytest)){
    tr2tstFolds <- createMultiFolds(1:nrow(xdat), k=5, times=10)
  } else{
    tr2tstFolds <- list(tr2tst = 1:nrow(get(prefixname)$xtrain))
  }
  
  # set list of window sizes for noise perturbations
  EXTpara_list <- generateEXTpara_list(modelname, xdat, n = nEXTpara_out) # min diff to range in log scale
  
  # stabilized kernel varying window sizes
  assign(paste('res_stab',prefixname,sep='_'), 
         mclapply(EXTpara_list, function(extpara){
           # computes kernel matrix
           message('Computing kernel matrix ... prefixname = ', prefixname, ', ext_para = ', extpara)
           kf <- kdtQUADRATICdot(a = extpara)
           kernmat <- computeKernelMatrix(xdata = xdat, kf = kf)
           res <- perfSVMKernelMatrix(model = modelname, prefix = prefixname, kmat = kernmat, grp = ydat, 
                                      tr2tstFolds = tr2tstFolds, Cpara_list = Cpara_list, 
                                      extension_para = extpara)
           return(res)
         }, mc.cores = 8))
  
  # MC approx for stablized kernel varying window sizes and number of sampled noise
  assign(paste('res_approx',prefixname,sep='_'), 
         mclapply(EXTpara_list, function(extpara){
           # generate noise matrices
           set.seed(101)
           noiseMatrix_list <- replicate(n = MaxNum_list[length(MaxNum_list)],
                                         expr = matrix(runif(length(xdat),min=-extpara,max=extpara),
                                                       nrow=nrow(xdat),ncol=ncol(xdat), 
                                                       dimnames = dimnames(xdat)), 
                                         simplify = FALSE)
           
           # computes kernel matrix
           kernmat <- NULL
           ss <- lapply(MaxNum_list, function(MaxNum){
             message('Computing kernel matrix ... prefixname = ', prefixname, ', ext_para = ', extpara, ', noiseNum = ', MaxNum)
             kernmat <<- approxKernelMatrix(xdata = xdat, kf = cor.fk, num = MaxNum, 
                                            noise = noiseMatrix_list, kmold = kernmat)
             res <- perfSVMKernelMatrix(model = modelname, prefix = prefixname, kmat = kernmat, grp = ydat, 
                                        tr2tstFolds = tr2tstFolds, Cpara_list = Cpara_list, 
                                        extension_para = extpara)
             return(res)
           })
           return(ss)
         }, mc.cores = 8))
}
```


```r
for (prefixname in prefixlist) {
  mainname <- nameafter[which(prefixname == namebefore)]
  key <- (prefixname %in% indepvallist)
  
  # stab-related
  res <- get(paste('res_stab',prefixname,sep='_'))
  ext <- sapply(res, function(u){u$extension_para})
  acc <- sapply(res, function(u){u$acc})
  
  # approx-related
  res <- get(paste('res_approx',prefixname,sep='_'))
  ext_approx <- sapply(res, function(u){u[[1]]$extension_para})
  acc_approx_list <- lapply(seq(length(res[[1]])), function(num){sapply(res, function(u){u[[num]]$acc})})
  names(acc_approx_list) <- names(res[[1]])
  # now pick a window size and plot MC approx performance
  idx <- names(which.max(acc))
  acc_approxSingleWindow <- sapply(res[[idx]], function(u){u$acc})
  ext_approxSingleWindow <- ext[[idx]]
  ref_stab <- acc[which(ext_approx == ext_approxSingleWindow)]
  
  if (key) {
    res <- res_indepval[[prefixname]]
    ref <- res$SVMkdtALL$acc[which.max(res$SVMkdtALL$cvacc)]
  } else {
    res <- res_crossval[[prefixname]]
    ref <- mean(sapply(res, function(foldres){foldres$SVMkdtALL$acc[which.max(foldres$SVMkdtALL$cvacc)]}))
  }
  
  # plotting with varying number of sampled noise and varying number of window sizes
  maxnum <- min(ceiling(length(acc_approx_list)/2), 5)
  plotrange <- c(min(floor(10*acc)/10), max(ceiling(10*acc)/10))
  par(font.lab = 2, font.axis = 2, font.main = 2, font = 2, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5)
  plot(2*ext, acc, main=mainname, xlab="Noise window size a", type='o', 
       ylim=plotrange, ylab=ifelse(key,"acc","cvacc"), lwd=2,lty=1,col=1,log='x')
  points(2*ext, rep(ref,length(ext)), lwd = 1, type = 'l', lty=5, col=2)
  for(col in seq(maxnum)){
    lines(2*ext_approx, acc_approx_list[[2*col-1]], lwd = 2, lty=1, col=col)
  }
  legend("topleft", legend=c('SVMkdtALLalt--exact', paste('SVMkdtALLalt--MCapprox (D=',2*seq(maxnum)-1,')', sep=''), 'SVMkdtALL'), 
         col=c(1,seq(maxnum),2),lty=c(1, rep(1,maxnum), 5), pch=c(1, rep(NA,maxnum), NA), lwd=c(2,rep(2,maxnum), 1),
         text.width=strwidth("10000000000000000000000000000000000"), cex=1.25)
  grid(ny=16)
  
  # plotting with more varying number of sampled noise for a specified window size
  maxnum <- length(acc_approxSingleWindow)
  plotrange <- c(min(floor(10*acc_approxSingleWindow)/10), max(ceiling(10*acc_approxSingleWindow)/10))
  par(font.lab = 2, font.axis = 2, font.main = 2, font = 2, cex.axis = 1.5, cex.lab = 1.5, cex.sub = 1.5)
  plot(MaxNum_list, acc_approxSingleWindow, main=mainname, 
       ylim=plotrange, xlab='Number D of random jitter', 
       ylab=ifelse(key,"acc","cvacc"), type='l',lwd=2,lty=1,col=1)
  points(MaxNum_list, rep(ref_stab, length(MaxNum_list)), type='o', lwd = 2, lty=1, col=1)
  points(MaxNum_list, rep(ref,length(MaxNum_list)), lwd = 1, type = 'l', lty=5, col=2)
  legend("topleft", 
         legend=c(paste('SVMkdtALLalt--exact (a=',round(2*ext_approxSingleWindow,0),')', sep=''),paste('SVMkdtALLalt--MCapprox (a=',round(2*ext_approxSingleWindow,0),')', sep=''), 'SVMkdtALL'),
         col=c(1,1,2),lty=c(1,1,5),lwd=c(2,2,1),pch=c(1,NA,NA),text.width = strwidth("100000000000000000000000000000000000000"),cex=1.25)
  grid(ny=16)
}
```

![plot of chunk approx_plot](figure/approx_plot-1.pdf)![plot of chunk approx_plot](figure/approx_plot-2.pdf)

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
##  package      * version date       source        
##  car            2.1-0   2015-09-03 CRAN (R 3.2.2)
##  caret        * 6.0-62  2015-11-23 CRAN (R 3.2.2)
##  codetools      0.2-14  2015-07-15 CRAN (R 3.2.2)
##  colorspace     1.2-6   2015-03-11 CRAN (R 3.2.2)
##  devtools       1.9.1   2015-09-11 CRAN (R 3.2.2)
##  digest         0.6.8   2014-12-31 CRAN (R 3.2.2)
##  evaluate       0.8     2015-09-18 CRAN (R 3.2.2)
##  foreach        1.4.3   2015-10-13 CRAN (R 3.2.2)
##  formatR        1.2.1   2015-09-18 CRAN (R 3.2.2)
##  ggplot2      * 2.1.0   2016-03-01 CRAN (R 3.2.3)
##  gtable         0.1.2   2012-12-05 CRAN (R 3.2.2)
##  iterators      1.0.8   2015-10-13 CRAN (R 3.2.2)
##  kernlab      * 0.9-22  2015-08-05 CRAN (R 3.2.2)
##  knitr          1.12.3  2016-01-22 CRAN (R 3.2.3)
##  lattice      * 0.20-33 2015-07-14 CRAN (R 3.2.2)
##  lme4           1.1-10  2015-10-06 CRAN (R 3.2.2)
##  magrittr       1.5     2014-11-22 CRAN (R 3.2.2)
##  MASS           7.3-43  2015-07-16 CRAN (R 3.2.2)
##  Matrix         1.2-2   2015-07-08 CRAN (R 3.2.2)
##  MatrixModels   0.4-1   2015-08-22 CRAN (R 3.2.2)
##  memoise        0.2.1   2014-04-22 CRAN (R 3.2.2)
##  mgcv           1.8-7   2015-07-23 CRAN (R 3.2.2)
##  minqa          1.2.4   2014-10-09 CRAN (R 3.2.2)
##  munsell        0.4.2   2013-07-11 CRAN (R 3.2.2)
##  mvtnorm        1.0-3   2015-07-22 CRAN (R 3.2.2)
##  nlme           3.1-121 2015-06-29 CRAN (R 3.2.2)
##  nloptr         1.0.4   2014-08-04 CRAN (R 3.2.2)
##  nnet           7.3-10  2015-06-29 CRAN (R 3.2.2)
##  pbkrtest       0.4-2   2014-11-13 CRAN (R 3.2.2)
##  pcaPP        * 1.9-60  2014-10-22 CRAN (R 3.2.2)
##  plyr           1.8.3   2015-06-12 CRAN (R 3.2.2)
##  quantreg       5.19    2015-08-31 CRAN (R 3.2.2)
##  Rcpp           0.12.2  2015-11-15 CRAN (R 3.2.2)
##  reshape2       1.4.1   2014-12-06 CRAN (R 3.2.2)
##  scales         0.3.0   2015-08-25 CRAN (R 3.2.2)
##  SparseM        1.7     2015-08-15 CRAN (R 3.2.2)
##  stringi        1.0-1   2015-10-22 CRAN (R 3.2.2)
##  stringr        1.0.0   2015-04-30 CRAN (R 3.2.2)
```
