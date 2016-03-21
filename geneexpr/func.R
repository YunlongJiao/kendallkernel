###############################################################
### performance main function
###############################################################

# for each model and each dataset, run
perfClassification <- function(modelname=NULL, prefixname, xtrain, ytrain, xtest=NULL, ytest=NULL,
                               Cpara_list, npairs_list, Conly, TSPrelated, NOpara,
                               nfolds = 5, nrepeats = 10, seed = 206, sigma=NULL
                               ){
  ### modelname/prefix is of class 'character' denoting the model/dataset name, NULL means running all models there are
  ### x/y-train/test are n*p exprs data matrices or binary grp vectors of length n, x/y-test set NULL if unavailable
  ### seed is set for nfolds*nrepeats data splits
  ### Cpara_list/npairs_list give C parameter and #genepairs lists
  ### max_nodes denotes max number of top-scoring gene pairs (not necessarily disjoint!), npairs_out controls grid density (only odd integers are kept to avoid ambiguity on voting)
  ### Conly, NOpara and TSPrelated are vectors of chars denoting model names
  
  if(is.null(Conly) && is.null(TSPrelated) && is.null(NOpara))
    stop('Specify at least one of the three subtypes of models to start running!')
  
  # data split
  set.seed(seed)
  foldIndices <- createMultiFolds(1:nrow(xtrain), k=nfolds, times=nrepeats)
  
  if(!is.null(modelname)){ # run a specific model
    if(!is.null(Conly) && (modelname %in% Conly)){
      res <- perfConly(model=modelname, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices, Cpara_list=Cpara_list, sigma=sigma)
    } else if(!is.null(TSPrelated) && (modelname %in% TSPrelated)){ # models related to scoring gene pairs
      stop("Do not apply! Always call to run all TSP related models at once!")
#       res <- perfTSPrelated(model=modelname, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices, Cpara_list=Cpara_list, npairs_list=npairs_list, TSPrelated=TSPrelated)
    } else if(!is.null(NOpara) && (modelname %in% NOpara)){ # model indep of any parameters (TSP included in second condition!)
      res <- perfNOpara(model=modelname, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices)
    } else{
      stop('Indefinite model input!')
    }
  } else{ # run all models
    if(is.null(Conly)){
      resConly <- NULL
    } else{
      resConly <- lapply(Conly, function(modelname){
        message(modelname)
        return(perfConly(model=modelname, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices, Cpara_list=Cpara_list))
      })
      names(resConly) <- Conly
    }
    
    if(is.null(NOpara)){
      resNOpara <- NULL
    } else{
      resNOpara <- lapply(NOpara, function(modelname){
        message(modelname)
        return(perfNOpara(model=modelname, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices))
      })
      names(resNOpara) <- NOpara
    }
    
    if(is.null(TSPrelated)){
      resTSPrelated <- NULL
    } else{
      resTSPrelated <- perfTSPrelated(model=NULL, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices, Cpara_list=Cpara_list, npairs_list=npairs_list, TSPrelated=TSPrelated)
      # NOTE only in model==NULL means running all TSPrelated models!
    }
    
    res <- c(resConly, resNOpara, resTSPrelated)
  }
  
  return(res)
}

perfExtensions <- function(modelname=NULL, prefixname, xtrain, ytrain, xtest=NULL, ytest=NULL, 
                           Cpara_list, Extensions, ExtensionsAPPROX=NULL, 
                           MaxNum=100, MaxExt=10, EXTpara_list=NULL, 
                           kf_basic = cor.fk, nfolds = 5, nrepeats = 10, seed = 206
){
  ### modelname/prefix is of class 'character' denoting the model/dataset name, NULL means running all models there are
  ### x/y-train/test are n*p exprs data matrices or binary grp vectors of length n, x/y-test set NULL if unavailable
  ### seed is set for nfolds*nrepeats data splits
  ### Cpara_list give C parameter lists for SVM
  ### Extensions/ExtensionsAPPROX are vectors of chars denoting model names
  ### MaxNum is the number of times when data is perturbed with noise
  ### MaxExt is the number of extension parameter to generate if EXTpara_list is set different than NULL
  
  if(is.null(Extensions) && is.null(ExtensionsAPPROX))
    stop('Specify at least one of the three subtypes of models to start running!')
  
  if(!is.null(modelname))
    message('Running all models still while model specified...')
  
  # data split
  set.seed(seed)
  foldIndices <- createMultiFolds(1:nrow(xtrain), k=nfolds, times=nrepeats)
  
  # run all models
  if(is.null(Extensions)){
    resExtensions <- NULL
  } else{
    resExtensions <- lapply(Extensions, function(modelname){
      if(is.null(EXTpara_list)){
        EXTpara_list <- generateEXTpara_list(model=modelname, x=rbind(xtrain,xtest), n=MaxExt, factor=1)
      }
      
      res <- lapply(EXTpara_list, function(extpm){
        message('dataset = ', prefixname, ', model = ', modelname,', Extension parameter = ', round(extpm,6))
        return(perfConly(model=modelname, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices, Cpara_list=Cpara_list, extension_para=extpm, sigma=sigma))
      })
      names(res) <- names(EXTpara_list)
      res[['EXTpara_list']] <- EXTpara_list
      return(res)
    })
    names(resExtensions) <- Extensions
  }
  
  if(is.null(ExtensionsAPPROX)){
    resExtensionsAPPROX <- NULL
  } else{
    resExtensionsAPPROX <- lapply(ExtensionsAPPROX, function(modelname){
      if(is.null(EXTpara_list)){
        EXTpara_list <- generateEXTpara_list(model=modelname, x=rbind(xtrain,xtest), n=MaxExt, factor=1)
      }
      
      res <- lapply(EXTpara_list, function(extpm){
        message('dataset = ', prefixname, ', model = ', modelname,', Extension parameter = ', round(extpm,6))
        return(perfApprox(model=modelname, prefix=prefixname, xtrain=xtrain, ytrain=ytrain, xtest=xtest, ytest=ytest, foldIndices=foldIndices, Cpara_list=Cpara_list, extension_para=extpm, MaxNum=MaxNum, kf=kf_basic))
      })
      names(res) <- names(EXTpara_list)
      res[['EXTpara_list']] <- EXTpara_list
      return(res)
    })
    names(resExtensionsAPPROX) <- ExtensionsAPPROX
  }
  
  res <- c(resExtensions, resExtensionsAPPROX)
  
  return(res)
}




###############################################################
### performance sub-functions
###############################################################

perfConly <- function(model, prefix, xtrain, ytrain, xtest, ytest, foldIndices, Cpara_list, extension_para=NULL, degree=2, scale=1, offset=0, sigma=NULL){
  
  # mostly SVM-based models plus C-independent-but-kernel-matrix-dependent KFD method
  # KDT (kernel fisher discriminant accompanies) SVM-based models as simple reference kernel machines
  # extension_para denotes the parameter for noise-perturbed dot kernel
  
  ntr <- nrow(xtrain); num <- 1; cvacc = acc = NULL; cvacc_kfd = acc_kfd = NULL;
  
  # compute kernel matrix
  if(model == "SVMrbf" && is.null(sigma)){
    sigma <- sigest(xtrain,scaled=F)['50%'] # use median-trick (on full training data despite cross- or indep-validation) to determine sigma in case of gaussian rbf kernel
  }
  
  kf <- switch(model,
               SVMrr = vanilladot(), # further refining on data matrix is defined below!
               SVMlinearALL = vanilladot(),
               SVMrbf = rbfdot(sigma = sigma), 
               SVMkdtALL = cor.fk, 
               SVMpolynomialALL = polydot(degree = degree, scale = scale, offset = offset), # pass on user-defined arguments
               SVMkdtALLstep = kdtSTEPdot(d = extension_para),
               SVMkdtALLquadratic = kdtQUADRATICdot(a = extension_para),
               SVMkdtALLerf = kdtERRORdot(sigma = extension_para),
               stop("Indefinite model input!"))
  
  if(is.null(xtest)){
    xdata <- xtrain
  } else{
    xdata <- rbind(xtrain,xtest)
  }
  
  message('computing kernel matrix ...')
  if(model == "SVMrr"){
    xdata <- t(apply(xdata, 1, rank, ties.method = 'average')) # here is further refining!
  }
  kmat <- computeKernelMatrix(xdata = xdata, kf = kf)
  
  rm(xdata)
  if(!is.null(xtest)){
    kmatTOT <- kmat # kmatTOT denotes the kernel matrix corresponding to both training and test set
    kmat <- kmatTOT[1:ntr,1:ntr] # kmat denotes the kernel matrix corresponding to training part
  }
  
  # cross-validating
  message("Cross-validating")
  foldscore <- lapply(foldIndices, function(fold){
    message('Fold No.',num); num <<- num + 1;
    kmcs <- centerScaleKernelMatrix(kmat, fold)
    
    ### KFD
    message('+', appendLF = F)
    pred_kfd <- classifierKFD(km = kmcs, trainidx = fold, traingrp = ytrain[fold])
    s_kfd <- evaluateAcc(pred_kfd,ytrain[-fold])
    
    ### SVM
    s <- sapply(Cpara_list, function(cpm){
      message('.', appendLF = F)
      pred <- classifierSVM(km=kmcs, trainidx=fold, traingrp=ytrain[fold], cpm=cpm)
      return(evaluateAcc(pred,ytrain[-fold]))
      })
    message('DONE!')
    return(c(s_kfd,s))
  })
  cvacc <- apply(as.matrix(as.data.frame(foldscore)),1,function(u){mean(u,na.rm=TRUE)})
  cvacc_kfd <- cvacc[1]; names(cvacc_kfd) <- NULL # KFD cv score
  cvacc <- cvacc[-1]; names(cvacc) <- names(Cpara_list) # SVM cv score (varying C)
  
  # indep-validating
  if(!is.null(xtest)){
    message("Independent-validating")
    kmcs <- centerScaleKernelMatrix(kmatTOT, 1:ntr)
    
    ### KFD
    message('+', appendLF = F)
    pred_kfd <- classifierKFD(km = kmcs, trainidx = 1:ntr, traingrp = ytrain)
    acc_kfd <- evaluateAcc(pred_kfd,ytest)
    
    ### SVM
    acc <- sapply(Cpara_list, function(cpm){
      message('.', appendLF = F)
      pred <- classifierSVM(km=kmcs, trainidx=1:ntr, traingrp=ytrain, cpm=cpm)
      return(evaluateAcc(pred,ytest))
    })
    message('DONE!')
  }
  
  return(list(model=model, prefix=prefix, Cpara_list=Cpara_list, extension_para=extension_para, cvacc=cvacc, acc=acc, cvacc_kfd=cvacc_kfd, acc_kfd=acc_kfd))
}



perfSVMKernelMatrix <- function(model, prefix, kmat, grp, tr2tstFolds, Cpara_list, extension_para = NULL, nfolds = 5, nrepeats = 1, seed = 206){
  ## perfSVMKernelMatrix basically runs the same thing as perfConly while taking FULL km and tr2tstFolds as input
  ## tr2tstFolds are (simulated) separate tr-to-tst indices corresp to kmat
  ## returns the average acc over tr2tstFolds with C parameter tuned by additional inner loop (seed set)
  
  # outter loop
  tr2tstfoldscore <- lapply(tr2tstFolds, function(indepfold){
    set.seed(seed)
    foldIndices <- createMultiFolds(indepfold, k=nfolds, times=nrepeats)
    
    message("Cross-validating ")
    foldscore <- lapply(foldIndices, function(fold){
      # inner loop
      kmcs <- centerScaleKernelMatrix(kmat[indepfold, indepfold], fold)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km = kmcs, trainidx = fold, traingrp = grp[indepfold[fold]])
      s_kfd <- evaluateAcc(pred_kfd,grp[indepfold[-fold]])
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmcs, trainidx=fold, traingrp=grp[indepfold[fold]], cpm=cpm)
        return(evaluateAcc(pred,grp[indepfold[-fold]]))
      })
      
      return(c(s_kfd,s))
    })
    cvacc <- apply(as.matrix(as.data.frame(foldscore)),1,function(u){mean(u,na.rm=TRUE)})
    cvacc_kfd <- cvacc[1]; names(cvacc_kfd) <- NULL # KFD cv score
    cvacc <- cvacc[-1]; names(cvacc) <- names(Cpara_list) # SVM cv score (varying C)
    
    
    message("Independent-validating ")
    kmcs <- centerScaleKernelMatrix(kmat, indepfold)
    
    ### KFD
    message('+', appendLF = F)
    pred_kfd <- classifierKFD(km = kmcs, trainidx = indepfold, traingrp = grp[indepfold])
    acc_kfd <- evaluateAcc(pred_kfd,grp[-indepfold])
    
    ### SVM
    acc <- sapply(Cpara_list, function(cpm){
      message('.', appendLF = F)
      pred <- classifierSVM(km=kmcs, trainidx=indepfold, traingrp=grp[indepfold], cpm=cpm)
      return(evaluateAcc(pred,grp[-indepfold]))
    })
    
    message('DONE!')
    return(list(acc=acc[which.max(cvacc)], acc_kfd=acc_kfd))
  })
  acc_kfd <- mean(sapply(tr2tstfoldscore, function(u){u$acc_kfd}))
  acc <- mean(sapply(tr2tstfoldscore, function(u){u$acc}))
  
  return(list(model=model, prefix=prefix, Cpara_list=Cpara_list, extension_para=extension_para, acc=acc, acc_kfd=acc_kfd))
}



perfApprox <- function(model, prefix, xtrain, ytrain, xtest, ytest, foldIndices, Cpara_list, extension_para, MaxNum, kf){
  # kernel machines with approx kernel matrices rather than expected version
  # MaxNum denotes max number of sampled data for kernel matrix computation
  # extension_para denotes the parameter for noise that perturbs Kendall kernel!!
  
  ntr <- nrow(xtrain); cvacc = acc = NULL; cvacc_kfd = acc_kfd = NULL;
  
  if(is.null(xtest)){
    xdata <- xtrain
  } else{
    xdata <- rbind(xtrain,xtest)
  }
  kmat <- matrix(0, nrow = nrow(xdata), ncol = nrow(xdata))
  
  numscore <- lapply(1:MaxNum, function(num){
    message('Noise sampling No. ',num)
    kmat <<- kmat + computeKernelMatrix(kf = kf,
                                        xdata = switch(model,
                                                       SVMkdtALLquadraticAPPROX = xdata + matrix(runif(n = length(xdata), min = -extension_para, max = extension_para), nrow = nrow(xdata)),
                                                       SVMkdtALLerfAPPROX = xdata + matrix(rnorm(n = length(xdata), mean = 0, sd = extension_para), nrow = nrow(xdata)),
                                                       stop('Indefinite model input!')))
    
    # cross-validating
    message("Cross-validating ")
    foldscore <- lapply(foldIndices, function(fold){
      kmcs <- centerScaleKernelMatrix(kmat[1:ntr,1:ntr], fold)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km = kmcs, trainidx = fold, traingrp = ytrain[fold])
      s_kfd <- evaluateAcc(pred_kfd,ytrain[-fold])
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmcs, trainidx=fold, traingrp=ytrain[fold], cpm=cpm)
        return(evaluateAcc(pred,ytrain[-fold]))
      })
      
      return(c(s_kfd,s))
    })
    cvacc <- apply(as.matrix(as.data.frame(foldscore)),1,function(u){mean(u,na.rm=TRUE)})
    cvacc_kfd <- cvacc[1]; names(cvacc_kfd) <- NULL # KFD cv score
    cvacc <- cvacc[-1]; names(cvacc) <- names(Cpara_list) # SVM cv score (varying C)
    
    # indep-validating
    if(!is.null(xtest)){
      message("Independent-validating ")
      kmcs <- centerScaleKernelMatrix(kmat, 1:ntr)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km = kmcs, trainidx = 1:ntr, traingrp = ytrain)
      acc_kfd <- evaluateAcc(pred_kfd,ytest)
      
      ### SVM
      acc <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmcs, trainidx=1:ntr, traingrp=ytrain, cpm=cpm)
        return(evaluateAcc(pred,ytest))
      })
    }
    
    return(list(model=model, prefix=prefix, Cpara_list=Cpara_list, extension_para=extension_para, cvacc=cvacc, acc=acc, cvacc_kfd=cvacc_kfd, acc_kfd=acc_kfd))
  })
  
  cvacc <- as.matrix(as.data.frame(lapply(numscore,function(u){u$cvacc})));
  dimnames(cvacc) <- list(names(Cpara_list),1:MaxNum)
  cvacc_kfd <- unlist(lapply(numscore,function(u){u$cvacc_kfd}));
  names(cvacc_kfd) <- 1:MaxNum
  
  if(!is.null(xtest)){
    acc <- as.matrix(as.data.frame(lapply(numscore,function(u){u$acc})));
    dimnames(acc) <- list(names(Cpara_list),1:MaxNum)
    acc_kfd <- unlist(lapply(numscore,function(u){u$acc_kfd}));
    names(acc_kfd) <- 1:MaxNum
  }
  
  return(list(model=model, prefix=prefix, Cpara_list=Cpara_list, noiseApproxMaxNum=MaxNum, extension_para=extension_para, cvacc=cvacc, acc=acc, cvacc_kfd=cvacc_kfd, acc_kfd=acc_kfd))
}


perfNOpara <- function(model, prefix, xtrain, ytrain, xtest, ytest, foldIndices){
  
  ntr <- nrow(xtrain); nts <- nrow(xtest); num <- 1; cvacc = acc = NULL;
  
  # NO cross-validating
  cvacc <- c(-100)
#   message("Cross-validating")
#   foldscore <- lapply(foldIndices, function(fold){
#     message('Fold No.',num); num <<- num + 1;
#     pred <- switch(model,
#                    APMV = classifierAPMV(xtrain=xtrain[fold,,drop=F],ytrain=ytrain[fold],xtest=xtrain[-fold,,drop=F]),
#                    stop('Indefinite model input!'))
#     message('DONE!')
#     return(evaluateAcc(pred,ytrain[-fold]))
#   })
#   cvacc <- mean(unlist(foldscore),na.rm=TRUE)
  
  # indep-validating
  if(!is.null(xtest)){
    message("Independent-validating")
    pred <- switch(model,
                   APMV = classifierAPMV(xtrain=xtrain,ytrain=ytrain,xtest=xtest),
                   stop('Indefinite model input!'))
    message('DONE!')
    acc <- evaluateAcc(pred,ytest)
  }
  
  return(list(model=model, prefix=prefix, cvacc=cvacc, acc=acc))
}


perfTSPrelated <- function(model=NULL, prefix, xtrain, ytrain, xtest, ytest, foldIndices, Cpara_list, npairs_list, TSPrelated){
  
  # This script is temporarily running all models and returning their prediction accuracies all at once!
  if(!all(TSPrelated==c("TSP", "kTSP", "SVMlinearTOP", "SVMkdtTOP", "SVMpolynomialTOP")))
    stop('Scripts of perfTSPrelated() need modification for adapting new models!')
  
  # C&k-dependent models are mostly SVM-based models plus C-independent-but-kernel-matrix-dependent KFD method
  # KDT (kernel fisher discriminant accompanies) SVM-based models as simple reference kernel machines
  
  for(tempname in TSPrelated){
    assign(paste('cvacc', tempname, sep='_'), NULL)
    assign(paste('acc', tempname, sep='_'), NULL)
    assign(paste('cvacc_kfd', tempname, sep='_'), NULL)
    assign(paste('acc_kfd', tempname, sep='_'), NULL)
  }
  
  ntr <- nrow(xtrain); nts <- nrow(xtest); num <- 1; 
  
  xtrain_rank <- t(apply(xtrain, 1, rank))
  
  # cross-validating
  message("Cross-validating")
  foldscore <- lapply(foldIndices, function(fold){
    message('Fold No.',num); num <<- num + 1;
    
    message('Calculating gene pairs ranking ... ')
    res_rankGenePairs <- scoreRankGenePairsLimitedNodes(xtrain_rank[fold,,drop=F],ytrain[fold],disjoint=F,max_nodes=npairs_list[length(npairs_list)])
    print(dim(res_rankGenePairs))
    message('DONE!')
    
    # TSP
    message('TSP')
    pred_TSP <- classifierkTSP(classes=levels(ytrain), xtest=xtrain[-fold,,drop=F], res_rankGenePairs=res_rankGenePairs[1,,drop=F])
    score_TSP <- evaluateAcc(pred_TSP, ytrain[-fold])
    message('DONE!')
    
    # kTSP
    message('kTSP')
    score_kTSP <- sapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      pred <- classifierkTSP(classes=levels(ytrain), xtest=xtrain[-fold,,drop=F], res_rankGenePairs=res_rankGenePairs[1:max_pairs,,drop=F])
      return(evaluateAcc(pred, ytrain[-fold]))
    })
    message('DONE!')
    
    # SVMkdtTOP
    message('SVMkdtTOP')
    score_SVMkdtTOP <- lapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      kmatrix <- computeKernelMatrix(xdata=xtrain, kf=kdtTOPdot(res_rankGenePairs[1:max_pairs,,drop=F]))
      kmatrix <- centerScaleKernelMatrix(kmat=kmatrix, trainidx=fold)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km=kmatrix, trainidx=fold, traingrp=ytrain[fold])
      s_kfd <- evaluateAcc(pred_kfd,ytrain[-fold])
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmatrix, trainidx=fold, traingrp=ytrain[fold], cpm=cpm)
        return(evaluateAcc(pred,ytrain[-fold]))
      })
      return(c(s_kfd,s))
    })
    score_kfd_SVMkdtTOP <- as.matrix(as.data.frame(score_SVMkdtTOP))[1,]; names(score_kfd_SVMkdtTOP) <- names(npairs_list)
    score_SVMkdtTOP <- as.matrix(as.data.frame(score_SVMkdtTOP))[-1,]; dimnames(score_SVMkdtTOP) <- list(names(Cpara_list),names(npairs_list))
    message('DONE!')
    
    # SVMlinearTOP
    message('SVMlinearTOP')
    score_SVMlinearTOP <- lapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      idxlist <- res_rankGenePairs[1:max_pairs,,drop=F]
      idxlist <- unique(as.vector(idxlist))
      kmatrix <- computeKernelMatrix(xdata=xtrain[,idxlist,drop=F], kf=vanilladot())
      kmatrix <- centerScaleKernelMatrix(kmat=kmatrix, trainidx=fold)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km=kmatrix, trainidx=fold, traingrp=ytrain[fold])
      s_kfd <- evaluateAcc(pred_kfd,ytrain[-fold])
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmatrix, trainidx=fold, traingrp=ytrain[fold], cpm=cpm)
        return(evaluateAcc(pred,ytrain[-fold]))
      })
      return(c(s_kfd,s))
    })
    score_kfd_SVMlinearTOP <- as.matrix(as.data.frame(score_SVMlinearTOP))[1,]; names(score_kfd_SVMlinearTOP) <- names(npairs_list)
    score_SVMlinearTOP <- as.matrix(as.data.frame(score_SVMlinearTOP))[-1,]; dimnames(score_SVMlinearTOP) <- list(names(Cpara_list),names(npairs_list))
    message('DONE!')
    
    # SVMpolynomialTOP
    message('SVMpolynomialTOP')
    score_SVMpolynomialTOP <- lapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      idxlist <- res_rankGenePairs[1:max_pairs,,drop=F]
      idxlist <- unique(as.vector(idxlist))
      kmatrix <- computeKernelMatrix(xdata=xtrain[,idxlist,drop=F], kf=polydot(degree=2, scale=1, offset=0))
      kmatrix <- centerScaleKernelMatrix(kmat=kmatrix, trainidx=fold)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km=kmatrix, trainidx=fold, traingrp=ytrain[fold])
      s_kfd <- evaluateAcc(pred_kfd,ytrain[-fold])
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmatrix, trainidx=fold, traingrp=ytrain[fold], cpm=cpm)
        return(evaluateAcc(pred,ytrain[-fold]))
      })
      return(c(s_kfd,s))
    })
    score_kfd_SVMpolynomialTOP <- as.matrix(as.data.frame(score_SVMpolynomialTOP))[1,]; names(score_kfd_SVMpolynomialTOP) <- names(npairs_list)
    score_SVMpolynomialTOP <- as.matrix(as.data.frame(score_SVMpolynomialTOP))[-1,]; dimnames(score_SVMpolynomialTOP) <- list(names(Cpara_list),names(npairs_list))
    message('DONE!')
    
    return(list(TSP=score_TSP,
                kTSP=score_kTSP,
                SVMkdtTOP=score_SVMkdtTOP,
                kfd_SVMkdtTOP=score_kfd_SVMkdtTOP,
                SVMlinearTOP=score_SVMlinearTOP,
                kfd_SVMlinearTOP=score_kfd_SVMlinearTOP,
                SVMpolynomialTOP=score_SVMpolynomialTOP,
                kfd_SVMpolynomialTOP=score_kfd_SVMpolynomialTOP))
  })
  
  # average fold scores
  cvacc_TSP <- averageFoldscore(fs = foldscore, nm = "TSP", Clist = NULL, nlist = NULL)
  cvacc_kTSP <- averageFoldscore(fs = foldscore, nm = "kTSP", Clist = NULL, nlist = npairs_list)
  cvacc_SVMlinearTOP <- averageFoldscore(fs = foldscore, nm = "SVMlinearTOP", Clist = Cpara_list, nlist = npairs_list)
  cvacc_kfd_SVMlinearTOP <- averageFoldscore(fs = foldscore, nm = "kfd_SVMlinearTOP", Clist = NULL, nlist = npairs_list)
  cvacc_SVMkdtTOP <- averageFoldscore(fs = foldscore, nm = "SVMkdtTOP", Clist = Cpara_list, nlist = npairs_list)
  cvacc_kfd_SVMkdtTOP <- averageFoldscore(fs = foldscore, nm = "kfd_SVMkdtTOP", Clist = NULL, nlist = npairs_list)
  cvacc_SVMpolynomialTOP <- averageFoldscore(fs = foldscore, nm = "SVMpolynomialTOP", Clist = Cpara_list, nlist = npairs_list)
  cvacc_kfd_SVMpolynomialTOP <- averageFoldscore(fs = foldscore, nm = "kfd_SVMpolynomialTOP", Clist = NULL, nlist = npairs_list)
  
  # indep-validating
  if(!is.null(xtest)){
    message("Independent-validating")
    
    ntr <- nrow(xtrain); nts <- nrow(xtest)
    
    message('Calculating gene pairs ranking ... ')
    res_rankGenePairs <- scoreRankGenePairsLimitedNodes(xtrain_rank,ytrain,disjoint=F,max_nodes=npairs_list[length(npairs_list)])
    message('DONE!')
    
    # TSP
    message('TSP')
    pred_TSP <- classifierkTSP(classes=levels(ytrain), xtest=xtest, res_rankGenePairs=res_rankGenePairs[1,,drop=F])
    acc_TSP <- evaluateAcc(pred_TSP, ytest)
    message('DONE!')
    
    # kTSP
    message('kTSP')
    acc_kTSP <- sapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      pred <- classifierkTSP(classes=levels(ytrain), xtest=xtest, res_rankGenePairs=res_rankGenePairs[1:max_pairs,,drop=F])
      return(evaluateAcc(pred, ytest))
    })
    message('DONE!')
    
    # SVMSVMkdtTOP
    message('SVMkdtTOP')
    acc_SVMkdtTOP <- lapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      kmatrix <- computeKernelMatrix(xdata=rbind(xtrain,xtest), kf=kdtTOPdot(res_rankGenePairs[1:max_pairs,,drop=F]))
      kmatrix <- centerScaleKernelMatrix(kmat=kmatrix, trainidx=1:ntr)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km=kmatrix, trainidx=1:ntr, traingrp=ytrain)
      s_kfd <- evaluateAcc(pred_kfd,ytest)
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmatrix, trainidx=1:ntr, traingrp=ytrain, cpm=cpm)
        return(evaluateAcc(pred,ytest))
      })
      return(c(s_kfd,s))
    })
    acc_kfd_SVMkdtTOP <- as.matrix(as.data.frame(acc_SVMkdtTOP))[1,]; names(acc_kfd_SVMkdtTOP) <- names(npairs_list)
    acc_SVMkdtTOP <- as.matrix(as.data.frame(acc_SVMkdtTOP))[-1,]; dimnames(acc_SVMkdtTOP) <- list(names(Cpara_list),names(npairs_list))
    message('DONE!')
    
    # SVMlinearTOP
    message('SVMlinearTOP')
    acc_SVMlinearTOP <- lapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      idxlist <- res_rankGenePairs[1:max_pairs,,drop=F]
      idxlist <- unique(as.vector(idxlist))
      kmatrix <- computeKernelMatrix(xdata=rbind(xtrain,xtest)[,idxlist,drop=F], kf=vanilladot())
      kmatrix <- centerScaleKernelMatrix(kmat=kmatrix, trainidx=1:ntr)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km=kmatrix, trainidx=1:ntr, traingrp=ytrain)
      s_kfd <- evaluateAcc(pred_kfd,ytest)
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmatrix, trainidx=1:ntr, traingrp=ytrain, cpm=cpm)
        return(evaluateAcc(pred,ytest))
      })
      return(c(s_kfd,s))
    })
    acc_kfd_SVMlinearTOP <- as.matrix(as.data.frame(acc_SVMlinearTOP))[1,]; names(acc_kfd_SVMlinearTOP) <- names(npairs_list)
    acc_SVMlinearTOP <- as.matrix(as.data.frame(acc_SVMlinearTOP))[-1,]; dimnames(acc_SVMlinearTOP) <- list(names(Cpara_list),names(npairs_list))
    message('DONE!')
    
    # SVMpolynomialTOP
    message('SVMpolynomialTOP')
    acc_SVMpolynomialTOP <- lapply(npairs_list, function(max_pairs){
      message(max_pairs,' ', appendLF = F)
      idxlist <- res_rankGenePairs[1:max_pairs,,drop=F]
      idxlist <- unique(as.vector(idxlist))
      kmatrix <- computeKernelMatrix(xdata=rbind(xtrain,xtest)[,idxlist,drop=F], kf=polydot(degree=2, scale=1, offset=0))
      kmatrix <- centerScaleKernelMatrix(kmat=kmatrix, trainidx=1:ntr)
      
      ### KFD
      message('+', appendLF = F)
      pred_kfd <- classifierKFD(km=kmatrix, trainidx=1:ntr, traingrp=ytrain)
      s_kfd <- evaluateAcc(pred_kfd,ytest)
      
      ### SVM
      s <- sapply(Cpara_list, function(cpm){
        message('.', appendLF = F)
        pred <- classifierSVM(km=kmatrix, trainidx=1:ntr, traingrp=ytrain, cpm=cpm)
        return(evaluateAcc(pred,ytest))
      })
      return(c(s_kfd,s))
    })
    acc_kfd_SVMpolynomialTOP <- as.matrix(as.data.frame(acc_SVMpolynomialTOP))[1,]; names(acc_kfd_SVMpolynomialTOP) <- names(npairs_list)
    acc_SVMpolynomialTOP <- as.matrix(as.data.frame(acc_SVMpolynomialTOP))[-1,]; dimnames(acc_SVMpolynomialTOP) <- list(names(Cpara_list),names(npairs_list))
    message('DONE!')
  }
  
  # concatenating results
  modelResults <- list()
  for(tempname in TSPrelated){
    modelResults[[tempname]] <- switch(tempname,
                                       TSP = list(model=tempname, prefix=prefix, cvacc=get(paste('cvacc',tempname,sep='_')), acc=get(paste('acc',tempname,sep='_'))),
                                       kTSP = list(model=tempname, prefix=prefix, npairs_list=npairs_list, cvacc=get(paste('cvacc',tempname,sep='_')), acc=get(paste('acc',tempname,sep='_'))),
                                       SVMlinearTOP = , 
                                       SVMkdtTOP = , 
                                       SVMpolynomialTOP = list(model=tempname, prefix=prefix, Cpara_list=Cpara_list, npairs_list=npairs_list, cvacc=get(paste('cvacc',tempname,sep='_')), acc=get(paste('acc',tempname,sep='_')), cvacc_kfd=get(paste('cvacc_kfd',tempname,sep='_')), acc_kfd=get(paste('acc_kfd',tempname,sep='_')))
                                       )
  }
  
  return(modelResults)
}

###############################################################
### useful functions
###############################################################

### Kernel matrix calculating and centering/scaling
computeKernelMatrix <- function(xdata, kf){
  ### xdata n*p data matrix
  n <- nrow(xdata); samplenames <- rownames(xdata)
  
  class(kf) <- 'kernel'
  kmat <- kernelMatrix(kf, as.matrix(xdata))
  dimnames(kmat) <- list(samplenames,samplenames)
  
  return(as.kernelMatrix(kmat))
}

centerScaleKernelMatrix <- function(kmat, trainidx){
  n <- nrow(kmat); ntr <- length(trainidx); nts <- n - ntr;
  
  #centering
  ed <- matrix(0,nrow=n,ncol=n)
  ed[trainidx,] <- 1
  kmcs <- kmat - t(ed)%*%kmat/ntr - kmat%*%ed/ntr + t(ed)%*%kmat%*%ed/(ntr^2)
  #scaling
  dsr <- sqrt(diag(kmcs))
  kmcs <- sweep(
    sweep(kmcs, 1, dsr, FUN='/'),
    2, dsr, FUN='/')
  rownames(kmcs) <- rownames(kmat)
  colnames(kmcs) <- colnames(kmat)
  
  return(as.kernelMatrix(kmcs))
}

approxKernelMatrix <- function(xdata, kf, num, noise, kmold = NULL){
  ### kmold is the one with (num-1) realizations and this function updates from there, default NULL indicates to start de novo
  ### therefore num needs to increase consecutively!!!
  
  ### xdata n*p data matrix, noise = {noise^1, ..., noise^MaxNum} is a list of n*p noise matrices
  ### num denotes the first of how many noise matrices to use, kf is the basic dot function
  ### kmold is the old matrix based on which we add more sampled realizations to the estimate
  ### [approxKernelMatrix_{ij}] returns [(1/num^2) \sum_{1<=k<=num, 1<=l<=num} kf(x_i+noise^k_i, x_j+noise^l_j)]
  
  if(is(xdata,"vector"))
    xdata <- as.matrix(xdata)
  if(!is(xdata,"matrix"))
    stop("data must be a matrix")
  if(!is(noise,"list"))
    stop("noise must be given as list object")
  
  n <- nrow(xdata)
  if(is.null(kmold)){
    kmold <- matrix(0, ncol = n, nrow = n, dimnames = list(rownames(xdata),rownames(xdata)))
  }
  res1 <- matrix(0, ncol = n, nrow = n, dimnames = list(rownames(xdata),rownames(xdata)))
  
  for(i in 1:n) {
    for(j in i:n) {
      if(num > 1){
        updi <- sum(sapply(1:(num-1), function(u){ kf(xdata[i,]+noise[[num]][i,], xdata[j,]+noise[[u]][j,]) }))
        updj <- sum(sapply(1:(num-1), function(u){ kf(xdata[i,]+noise[[u]][i,], xdata[j,]+noise[[num]][j,]) }))
      } else{
        updi = updj = 0
      }
      updij <- kf(xdata[i,]+noise[[num]][i,], xdata[j,]+noise[[num]][j,])
      res1[i,j]  <- ((num-1)*(num-1)*kmold[i,j] + updi + updj + updij) / (num*num)
    }
  }
  res1 <- res1 + t(res1)
  diag(res1) <- diag(res1)/2
  
  return(as.kernelMatrix(res1))
}


### generate noise distribution parameters for extension models
generateEXTpara_list <- function(model, x, n = 10, factor = 1){
  
  ### excerpt (slightly different though) of jitter() from base package, originally designed for uniform noise
  
  if (model %in% c('SVMkdtALLerf', 'SVMkdtALLerfAPPROX')) {
    factor <- factor/3 # 3-sigma principle or empirical 68–95–99.7 rule
  }
  
  x <- as.vector(x)
  
  z <- diff(r <- range(x[is.finite(x)]))
  if (z == 0) 
    z <- abs(r[1L])
  if (z == 0) 
    z <- 1
  
  # d <- diff(xx <- unique(sort.int(round(x, 3 - floor(log10(z))))))
  d <- diff(xx <- sort(unique(x))) # no rounding!
  
  d <- if (length(d)) {
    min(d)
  } else if (xx != 0) {
    xx/10
  } else {
    z/10
  }
  
  ### end of excerpt
  
  amount_min <- factor * abs(d) # minimum difference in all values
  amount_max <- factor * (z) # range in all values
  
  # if(amount_min >= amount_max){
  #   warning('RARE: amount_min larger than amount_max thus reversed!')
  #   amount_min <- factor * (z/50)
  #   amount_max <- factor/5 * abs(d)
  # }
  
  # if(amount_max <= abs(d)){
  #   warning('Extension parameter is NOT large enough to differ from kdt!')
  # }
  
  # Case I: in equidistant-scale
  # para <- seq(amount_min, amount_max, length.out = n)
  # Case II: in log-scale
  para <- amount_min*exp(seq(0,1,length.out=n)*log(amount_max/amount_min))
  # Case III: user-defined coef
  # if(n > 2)
  #   para <- c(factor/5 * abs(d), factor*abs(d)*(1:(n-2)), factor * (z/50))
  # else
  #   para <-  c(factor/5 * abs(d), factor * (z/50))[1:n]
  
  names(para) <- paste('EXT',1:length(para),sep='')
  return(para)
}

### model assessment
evaluateAcc <- function(predictions,observations){
  if(!is(predictions,"factor")) stop('Predictions are not factors!')
  if(all(levels(predictions)==levels(observations)))
    return(sum(predictions==observations,na.rm=T)/length(predictions))
  else
    stop('Predictions and observations have different levels!')
}

### TSP scoring function
scoreRankGenePairsLimitedNodes <- function(xtrain_rank, ytrain, max_pairs=0, disjoint=FALSE, max_nodes=0, score_tie=TRUE){
  
  ### xtrain is ntr*p matrix of expression profiles
  ### ytrain is ntr-vector of binary labels
  ### max_nodes == 0 by default stands for keeping track of all pairs
  ### disjoint denotes whether select only disjoint pairs
  ### scoreRankGenePairsLimitedNodes() integrates the two func above and does everything in C with linked list structure,
  ### thus returning a 2-col matrix of top scoring gene pairs, which ALWAYS BREAK TIES with secondary scores!
  
  # convert xtrain expression data to rank data assumed done!!!
  # xtrain_rank <- t(apply(xtrain, 1, rank))
  
  p <- ncol(xtrain_rank)
  
  if(max_pairs == 0){
    if(disjoint == T){
      max_pairs <- floor(p/2)
    } else {
      max_pairs <- p*(p-1)/2
    }
  }
  
  if(max_nodes == 0){
    max_nodes <- p*(p-1)/2
  }
  
  # split the rank data.frame according to factor levels
  xtrain_split <- split(as.data.frame(xtrain_rank),ytrain)
  
  # call C functions to compute scores
  result <- .Call('scoreRankGenePairsLimitedNodesC2R',as.matrix(xtrain_split[[1]]),as.matrix(xtrain_split[[2]]),as.integer(max_pairs),as.logical(disjoint), as.integer(max_nodes))
  
  max_pairs <- nrow(result)
  colnames(result) <- c('gene1idx','gene2idx')
  rownames(result) <- 1:max_pairs
  
  if(disjoint == T){
    message('Hey! You have just found ',max_pairs, ' disjoint top gene pairs!');
  } else{
    message('Hey! You have just found ',max_pairs, ' (overlapping) top gene pairs!');
  }
  
  return(result)
}

mytspR <- function(xtrain_rank, ytrain, max_nodes=0){
  
  ################################################
  ### AVL tree proved SLOWER THAN linked list !!
  ################################################
  
  ### xtrain is ntr*p RANK matrix of expression profiles
  ### ytrain is ntr-vector of binary labels
  ### max_nodes == 0 by default stands for keeping track of all pairs
  ### mytspR() does everything in C using AVL tree structure
  ### and returns a 2-col matrix of top scoring gene pairs
  
  # convert xtrain expression data to rank data assumed done!!!
  # xtrain_rank <- t(apply(xtrain, 1, rank))
  
  p <- ncol(xtrain_rank)
  
  if(max_nodes == 0){
    max_nodes <- p*(p-1)/2
  }
  
  # split the rank data.frame according to factor levels
  xtrain_split <- split(as.data.frame(xtrain_rank),ytrain)
  
  # call C functions to compute scores
  result <- .Call('mytspC2R',as.matrix(xtrain_split[[1]]),as.matrix(xtrain_split[[2]]), as.integer(max_nodes))
  colnames(result) <- c('gene1idx','gene2idx')
  rownames(result) <- 1:max_nodes
  
  message('Hey! You have just found ',max_nodes, ' (non-disjoint) top gene pairs!');
  
  return(result)
}

### kdtTOP dot
kdtTOPdot <- function(res_rankGenePairs){
  
  # res_rankGenePairs is the returned obj from calling scoreRankGenePairsLimitedNodes()
  # tspdot() returns the kernel function that computes natural product of two augmented profiles
  
  kernf <- function(x, y) {
    x_aug <- (x[res_rankGenePairs[,1]] < x[res_rankGenePairs[,2]])
    y_aug <- (y[res_rankGenePairs[,1]] < y[res_rankGenePairs[,2]])
    res <- (x_aug == y_aug)
    res <- sum(res)/nrow(res_rankGenePairs)
    return(res)
  }
  
  return(kernf)
}

### kdtExtension1 - piece-wise STEP dot
kdtSTEPdot <- function(d){
  
  # kdtSTEPdot() returns the kernel function that computes naive kdt extension
  
  kernf <- function(x, y){
    if(length(x) != length(y))
      stop('Input vectors should be of the same size!')
    else
      p <- length(x)
    
    res <- .Call('kdtSTEPdotC2R', as.double(x),as.double(y),as.integer(p),as.double(d))
    
    return(res)
  }
  
  return(kernf)
}

### kdtExtension2 - piece-wise QUADRATIC dot
kdtQUADRATICdot <- function(a){
  
  # kdtQUADRATICdot() returns the kernel function that computes uniform-noise-stablized kdt extension
  
  kernf <- function(x, y){
    if(length(x) != length(y))
      stop('Input vectors should be of the same size!')
    else
      p <- length(x)
    
    res <- .Call('kdtQUADRATICdotC2R', as.double(x),as.double(y),as.integer(p),as.double(a))
    
    return(res)
  }
  
  return(kernf)
}

### kdtExtension3 - error function dot
kdtERRORdot <- function(sigma){
  
  # kdtERRORdot() returns the kernel function that computes Gaussian-noise-stablized kdt extension
  
  kernf <- function(x, y){
    if(length(x) != length(y))
      stop('Input vectors should be of the same size!')
    else
      p <- length(x)
    
    res <- .Call('kdtERRORdotC2R', as.double(x),as.double(y),as.integer(p),as.double(sigma))
    
    return(res)
  }
  
  return(kernf)
}


### averaging foldscores (not a generic function)
averageFoldscore <- function(fs, nm, Clist = NULL, nlist = NULL){
  # Clist/nlist given non-null indicates the dimension of each foldscore, namely 1x1 nx1 or Cxk
  
  fs_ave <- apply(as.data.frame(lapply(fs,function(u){as.vector(u[[nm]])})),1,function(v){mean(v,na.rm=TRUE)});
  
  if(is.null(Clist) && !is.null(nlist)){
    names(fs_ave) <- names(nlist)
  }
  if(!is.null(Clist) && !is.null(nlist)){
    fs_ave <- matrix(fs_ave, nrow = length(Clist), ncol = length(nlist), dimnames = list(names(Clist),names(nlist)))
  }
  
  return(fs_ave)
}

###############################################################
### classifiers
###############################################################

classifierSVM <- function(km, trainidx, traingrp, cpm){
  
  # classifierSVM returns a ntst-vector of predicted binary classes
  # NOTE those indices not in trainidx is ready for test!!
  # NOTE that (row/sample) NAMES of kernel matrix is necessary for naming predicted vector
  
  # training
  mo <- ksvm(as.kernelMatrix(km[trainidx,trainidx]), y=traingrp, C=cpm, type="C-svc")
  # predicting
  pred <- predict(mo, as.kernelMatrix(km[-trainidx,trainidx,drop=F][,SVindex(mo),drop=F]), type = "response")
  names(pred) <- rownames(km)[-trainidx]
  
  return(pred)
}

classifierAPMV <- function(xtrain, ytrain, xtest){
  
  # xtrain/xtest expr data matrix of n*p
  # classifierAPMV() predicts classes using majority votes based on all pairs
  
  ntr <- nrow(xtrain); nts <- nrow(xtest);
  trainsamples <- rownames(xtrain); testsamples <- rownames(xtest)
  classes <- levels(ytrain)
  
  pred_response <- .Call('allPairsMajorityVotesC2R', as.numeric(as.matrix(xtrain)),as.integer(ytrain),as.numeric(as.matrix(xtest)))
  pred_response <- factor(pred_response,levels=c(1,2),labels=classes)
  names(pred_response) <- testsamples
  
  return(pred_response)
}

classifierkTSP <- function(classes, xtest, res_rankGenePairs){
  # classes denotes the labels of groups in training set
  # xtest denotes gene expression profiles of dimension nts*p
  # classifierkTSP() returns a list of class predictions
  
  nts <- nrow(xtest); testsamples <- rownames(xtest)
  
  pred <- xtest[,res_rankGenePairs[,1]] < xtest[,res_rankGenePairs[,2]]
  
  pred_response <- apply(matrix(pred,nrow=nts), 1, function(u){
    tab <- table(u)
    return(names(tab)[which.max(tab)])
  })
  pred_response[pred_response=='TRUE'] <- classes[1]
  pred_response[pred_response=='FALSE'] <- classes[2]
  pred_response <- factor(pred_response, levels=classes)
  names(pred_response) <- testsamples
  
  return(pred_response)
}

classifierKFD <- function(km, trainidx, traingrp, mu = 1e-3){
  
  # mu is a small number added to diag(N) to avoid singularity
  # classifierKFD returns a ntst-vector of predicted binary classes
  # NOTE those indices not in trainidx is ready for test!!
  # NOTE that (row/sample) NAMES of kernel matrix is necessary for naming predicted vector
  
  ntot <- nrow(km); ntrn <- length(trainidx); ntst <- ntot - ntrn
  classes <- levels(traingrp)
  if(length(classes) != 2)
    stop('Training response should be (potentially) two classes!')
  
  tab <- table(traingrp)
  if(sum(tab) == 0)
    stop('Wrong input!')
  else if(tab[1] == 0){
    pred_response <- factor(rep(classes[2],ntst),levels=classes)
    return(pred_response)
  } else if(tab[2] == 0){
    pred_response <- factor(rep(classes[1],ntst),levels=classes)
    return(pred_response)
  }
  
  idx1 <- trainidx[traingrp == classes[1]]
  idx2 <- trainidx[traingrp == classes[2]]
  idxtst <- setdiff(1:ntot, trainidx); testsamples <- rownames(km)[idxtst]
  
  M1 <- apply(km[trainidx,idx1,drop=F],1,mean)
  M2 <- apply(km[trainidx,idx2,drop=F],1,mean)
  
  L1 <- matrix(1/length(idx1),nrow=length(idx1),ncol=length(idx1))
  L2 <- matrix(1/length(idx2),nrow=length(idx2),ncol=length(idx2))
  N <- km[trainidx,trainidx,drop=F] %*% km[trainidx,trainidx,drop=F] - km[trainidx,idx1,drop=F] %*% L1 %*% km[idx1,trainidx,drop=F] - km[trainidx,idx2,drop=F] %*% L2 %*% km[idx2,trainidx,drop=F]
  
  # case I: being  adjusted by within-class variance i.e. KFD
  alpha <- solve(N + diag(mu, nrow = nrow(N), ncol = ncol(N))) %*% matrix(M2 - M1, ncol=1)
  # case II: without being  adjusted by within-class variance
  # alpha <- matrix(M2 - M1, ncol=1)
  
  pred <- as.vector(km[idxtst,trainidx,drop=F] %*% matrix(alpha, ncol=1))
  pred_response <- as.vector(as.numeric(pred >= 0))
  pred_response <- factor(pred_response, levels=c(0,1), labels=classes)
  names(pred_response) <- testsamples
  
  return(pred_response)
}

