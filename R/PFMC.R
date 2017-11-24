calcvcox <-
  function(Time,
           Delta,
           X,
           ncore = 2,
           K = 2,
           nfolds = 10,
           nopenaltyresult = NULL,
           alpha = 1,
           scad = FALSE,
           adpcoef = NULL,
           nlambda = 100,
           seed=1){
    if(is.null(nopenaltyresult))
      stop("Fitting of the full finite-mixture Cox PH model is null!")
    U = nopenaltyresult$U
    pi = nopenaltyresult$pi

    n = length(Time)
    X = as.data.frame(X)
    Xm = as.matrix(X)
    Time = matrix(Time, ncol = 1)
    Delta = matrix(Delta, ncol = 1)
    p = dim(Xm)[2]
    surv = cbind(time = Time, status = Delta)
    colnames(surv) = c("time", "status")

    ####check adpcoef
    if(is.null(adpcoef)|scad)
      adpcoef = lapply(1:K, function(k)rep(1,times = p))
    if(!is.list(adpcoef))
      stop("adpcoef must be a list!")
    if(length(adpcoef) != K)
      stop("K and length of adpcoef do not match!")
    if(any(sapply(adpcoef, length) != p))
      stop("length in adpcoef and the dimension of covarites do not match!")
    adpcoef = lapply(adpcoef, function(x)(x/sum(x))*p)
    if(scad)
      alpha = 1


    corenum = detectCores()
    if(ncore > corenum){
      warning("there are(is) only ",corenum, " core(s) on this machine!")
    }
    ncore = min(ncore, corenum)
    cl = makeCluster(rep("localhost", times = ncore), type = "SOCK")
    registerDoSNOW(cl)

    k = 1
    allresult = foreach(k = 1:K)%dopar%{
      set.seed(seed)
      if(!scad){
        idx = U[,k]/pi[k] >= 1e-6
        weights = U[idx, k] / sum(U[idx, k]) * sum(idx)
        ret = cv.glmnet(x = Xm[idx,], y = surv[idx,], weights = weights, family="cox", grouped = TRUE,type.measure = "deviance", nfolds = nfolds,nlambda=nlambda,alpha = alpha, parallel = FALSE, penalty.factor = adpcoef[[k]])
        lambda = c(ret$lambda.min, ret$lambda.1se)
      }else{
        idx = U[,k]/pi[k] >= 1e-6
        weights = U[idx, k] / sum(U[idx, k]) * sum(idx)
        fit = glmnet::cv.glmnet(x = Xm[idx,], y = surv[idx,], weights = weights, family="cox", grouped = TRUE,type.measure = "deviance", nfolds = nfolds, nlambda = nlambda,thresh = 1,alpha = alpha, parallel = FALSE,maxit =100000)
        lambda1 = fit$lambda
        ret = cv.ncvsurv(X = Xm[idx,], y = surv[idx,], penalty = "SCAD", gamma = 3.7, alpha = 1,lambda = lambda1,max.iter=100000, nfolds = nfolds, penalty.factor = adpcoef[[k]],weights = weights)
        lambda = c(ret$lambda.min, ret$lambda.1se)
      }
      names(lambda) = c("lambda.min", "lambda.1se")
      return(lambda)
    }
    names(allresult) = paste0("comp", 1:K)
    stopCluster(cl)
    out <- cbind(sapply(allresult,function(x) x))
    rownames(out) <- c("lambda.min","lambda.1se")
    return(out)
  }

pmixcox <-
function(Time,
         Delta,
         X,
         K = 3,
         iter.max = 1000,
         u.init = NULL,
         tparm = 0,
         alpha = 1,
         scad = FALSE,
         adpcoef = NULL,
         abstol = 1e-2,
         reltol = 1e-4,
         seed = 1
  ){
    set.seed(seed)
    n = length(Time)
    X = as.data.frame(X)
    Xm = as.matrix(X)
    Time = matrix(Time, ncol = 1)
    Delta = matrix(Delta, ncol = 1)
    p = dim(Xm)[2]
    surv = cbind(time = Time, status = Delta)
    colnames(surv) = c("time", "status")

    ####check tparm
    if(length(tparm) == 1)
      tparm = rep(tparm, times = K)
    else
      if(length(tparm) != K)
        stop("length of tparm must be 1 or K!")

    ####check u.init and initial U matrix
    if(is.matrix(u.init)){
      if(all(dim(u.init) == c(n,K)))
        U = u.init
      else
        stop("parameter u.init wrong!")
    }else if(is.vector(u.init)){
      if(length(u.init) != n | any(!(u.init %in% (1:K))))
        stop("parameter u.init wrong!")
      else{
          U = matrix(0.001, K, K)
          diag(U) = 1 - (K-1)*0.001
          U = U[u.init, ]
          if(K==1)
            U = matrix(U, ncol=1)
        }
    }else if(is.null(u.init)){
      u.init = sample(1:K, size = n, replace = TRUE)
      U = matrix(0.001, K, K)
      diag(U) = 1 - (K-1)*0.001
      U = U[u.init, ]
      if(K==1)
        U = matrix(U,ncol=1)
    }else{
      stop("parameter u.init wrong!")
    }

    ####check adpcoef
    if(is.null(adpcoef)|scad)
      adpcoef = lapply(1:K, function(k)rep(1,times = p))
    if(!is.list(adpcoef))
      stop("adpcoef must be a list!")
    if(length(adpcoef) != K)
      stop("K and length of adpcoef do not match!")
    if(any(sapply(adpcoef, length) != p))
      stop("length in adpcoef and the dimension of covarites do not match!")
    adpcoef = lapply(adpcoef,function(x)(x/sum(x))*p)

    if(scad)
      alpha = 1

    ####initial pi
    pi = apply(U, 2, mean)

    ####EM iteration
    ploglik = NULL
    mixloglik = NULL
    ploglikpenalty = NULL
    mixloglikpenalty = NULL
    lastploglikpenalty = -Inf
    lastmixloglikpenalty = -Inf
    convergence = FALSE
    iter = 0

    while(!convergence & iter<iter.max){
      ##Mstep
      if(!scad){
        fit = lapply(1:K, function(k){
          idx = U[,k]/pi[k] >= 1e-6
          weights = U[idx, k] / sum(U[idx, k]) * sum(idx)
          result = glmnet(x = Xm[idx,], y = surv[idx,], family = c("cox"), weights = weights,offset = NULL, alpha = alpha, nlambda = 1, penalty.factor = adpcoef[[k]],lambda = tparm[k],standardize = TRUE, thresh = 1e-07)
          return(as.numeric(result$beta))
        })
      }else{
        fit = lapply(1:K, function(k){
          idx = U[,k]/pi[k] >= 1e-6
          weights = U[idx, k] / sum(U[idx, k]) * sum(idx)
          result = ncvsurv(X = Xm[idx,], y = surv[idx,], penalty = "SCAD", gamma = 3.7,alpha = 1, weights = weights,max.iter =100000, penalty.factor = adpcoef[[k]],lambda = tparm[k])
          return(as.numeric(result$beta))
        })
      }

      ##E-step
      h0 = NULL
      H0 = NULL
      EE = NULL
      PLL = NULL
      allUEXB = NULL
      for (k in 1:K) {
        Uk = U[,k]
        XB = Xm%*%fit[[k]]
        EXB = exp(XB)
        EE = cbind(EE, EXB)

        UEXB = Uk*EXB
        allUEXB = cbind(allUEXB, UEXB)
        aa = sapply(Time, function(x)sum(UEXB[Time>=x]))
        bb = 1/aa
        bb[bb == Inf] = 0

        h00 = Uk*bb*Delta
        H00 = sapply(Time, function(x)sum(h00[Time<=x]))
        h0 = cbind(h0, h00)
        H0 = cbind(H0, H00)

        PLL = c(PLL, sum((Uk*(XB-log(aa)))[Delta==1&aa>0]))
      }

      CC = (h0 * EE)^(Delta %*% t(rep(1, K)))
      CC = (CC * exp(-H0 * EE)) %*% diag(as.vector(pi), K)
      nowploglik = sum(PLL) + sum((U%*%diag(log(pi), K))[U != 0])
      nowmixloglik = sum(log(apply(CC, 1, sum)))

      if(!scad){
        penalty  = sum(n*pi*tparm * sapply(1:K, function(k){
          sum(adpcoef[[k]] * (alpha * abs(fit[[k]]) + (1 - alpha)/2*(fit[[k]]^2)))
        }))
      }else{
        penalty  = sum(pi * sapply(1:K, function(k){
          sum(sapply(fit[[k]],function(x){
            absx = abs(x)
            return(ifelse(absx <= tparm[k], tparm[k]*absx,
                          ifelse(absx <= 3.7*tparm[k],(tparm[k]*absx*3.7 - 0.5*(absx^2 + tparm[k]^2))/2.7,tparm[k]^2 * 4.7/2)))
          }) * adpcoef[[k]])
        }))
      }

      nowploglikpenalty = nowploglik - penalty
      nowmixloglikpenalty = nowmixloglik - penalty
      ploglikpenalty = c(ploglikpenalty, nowploglikpenalty)
      mixloglikpenalty = c(mixloglikpenalty, nowmixloglikpenalty)
      ploglik = c(ploglik, nowploglik)
      mixloglik = c(mixloglik, nowmixloglik)

      if(!convergence){
        U = do.call(rbind, lapply(1:n,function(i){
          x = CC[i,]
          if(any(is.na(x)))
            return(pi)
          else
            if(all(x == 0)){
              restart = rep(0, times = K)
              restart[sample(1:K, size = 1, replace = FALSE, prob = pi)] = 1
              return(restart)
            }
          return(x / sum(x))
        }))
        if(K == 1)
          U = matrix(U, ncol=1)
        pi = apply(U, 2, mean)
      }

      if(is.finite(nowploglikpenalty) & is.finite(lastploglikpenalty) &
           (abs(nowploglikpenalty-lastploglikpenalty) < abstol |
              abs(nowploglikpenalty-lastploglikpenalty) <
              abs(reltol*(nowploglikpenalty + reltol))))
        convergence = TRUE
      lastploglikpenalty = nowploglikpenalty
      lastmixloglikpenalty = nowmixloglikpenalty

      iter = iter + 1
    }
    if(iter==iter.max)
      warning("EM iterations reach iter.max!")
    class = apply(U, 1, function(x){
      which.max(x)[1]
    })
    finalmixloglik = mixloglik[length(mixloglik)]

    ret = list(U = U,
               fit = fit,
               pi = pi,
               class = class,
               ploglik = ploglik,
               mixloglik = mixloglik,
               iter = iter,
               convergence = iter < iter.max)
    return(ret)
  }



convexMin <- function(b, X, penalty, gamma, l2, penalty.factor, a, Delta=NULL, weights) {
  n <- nrow(X)
  p <- ncol(X)
  l <- ncol(b)

  if (penalty=="MCP") {
    k <- 1/gamma
  } else if (penalty=="SCAD") {
    k <- 1/(gamma-1)
  } else if (penalty=="lasso") {
    return(NULL)
  }
  if (l==0) return(NULL)

  val <- NULL
  for (i in 1:l) {
    A1 <- if (i==1) rep(1,p) else b[,i]==0
    if (i==l) {
      L2 <- l2[i]
      U <- A1
    } else {
      A2 <- b[,i+1]==0
      U <- A1&A2
      L2 <- l2[i+1]
    }
    if (sum(!U)==0) next
    Xu <- X[,!U]
    p.. <- k*(penalty.factor[!U]!=0) - L2*penalty.factor[!U]

    eta <- if (i==l) X%*%b[,i] else X%*%b[,i+1]
    haz <- drop(exp(eta))
    rsk <- rev(cumsum(rev(weights*haz)))
    h <- weights*haz*cumsum(Delta*weights/rsk)
    xwxn <- crossprod(sqrt(h) * Xu)/n
    eigen.min <- min(eigen(xwxn-diag(diag(xwxn)*p.., nrow(xwxn), ncol(xwxn)))$values)


    if (eigen.min < 0) {
      val <- i
      break
    }
  }
  val
}



cv.ncvsurv <- function(X, y, ..., cluster, nfolds=10, seed, returnY=FALSE, trace=FALSE) {

  # Complete data fit
  fit.args <- list(...)
  fit.args$X <- X
  fit.args$y <- y
  fit.args$returnX <- TRUE
  fit <- do.call("ncvsurv", fit.args)
  if(is.null(fit.args$weights)) weights = rep(1, nrow(X)) else weights <- fit.args$weights
  if(is.null(fit.args$alpha)) alpha = 1 else alpha <- fit.args$alpha
  if(is.null(fit.args$penalty.factor)) penalty.factor = rep(1,ncol(X)) else penalty.factor <- fit.args$penalty.factor

  lambda <- fit$lambda
  #cat(weights,"\n")


  # Get standardized X, y
  #X <- fit$X
  #y <- cbind(fit$time, fit$fail)
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX) fit$X <- NULL

  # Set up folds
  n <- nrow(X)
  if (!missing(seed)) set.seed(seed)
  cv.ind <- ceiling(sample(1:n)/n*nfolds)
  Y <- matrix(NA, nrow=n, ncol=length(fit$lambda))

  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$warn <- FALSE
  cv.args$convex <- FALSE
  cv.args$penalty.factor <- fit$penalty.factor
  if (!missing(cluster)) {
    if (!("cluster" %in% class(cluster))) stop("cluster is not of class 'cluster'; see ?makeCluster")
    clusterExport(cluster, c("cv.ind","fit","X", "y", "cv.args"), envir=environment())
#    clusterCall(cluster, function() require(ncvreg))
    fold.results <- parLapply(cl=cluster, X=1:nfolds, fun=cvf.surv, XX=X, y=y, cv.ind=cv.ind, cv.args=cv.args,weights=weights)
  }
  cvraw = matrix(NA, nfolds, length(fit$lambda))

  beta = list()

    for (i in seq(nfolds)) {
      which = cv.ind == i
      if (!missing(cluster)) {
        res <- fold.results[[i]]
        beta[[i]] = res$beta
      } else {
        if (trace) cat("Starting CV fold #",i,sep="","\n")
        fit.args2 <- list(...)
        fit.args2$X <- X[!which,]
        fit.args2$y <- y[!which,]
        fit.args2$returnX <- TRUE
        fit.args2$weights <- weights[!which]
        fit.args2$lambda <- lambda
        fit.args2$alpha <- alpha
        fit.args2$penalty.factor <- penalty.factor
        res <- do.call("ncvsurv", fit.args2)
        beta[[i]] = res$beta

      }
      coefmat = beta[[i]]

        plminusk = coxnet.deviance(x = X[!which, ], y = y[!which,],  weights = weights[!which],
                                   beta = coefmat)
        plfull = coxnet.deviance(x = X, y = y,
                                 weights = weights, beta = coefmat)
        cvraw[i, seq(along = plfull)] = plfull - plminusk

    }
    status = y[, "status"]
    N = nfolds - apply(is.na(cvraw), 2, sum)
    err.ind = which(N>1)
    N = N[err.ind]
    weights = as.vector(tapply(weights * status, cv.ind, sum))
    cvraw = (cvraw/weights)[,err.ind]
    cvm = apply(cvraw, 2, weighted.mean, w = weights, na.rm = TRUE)
    cvsd = sqrt(apply(scale(cvraw, cvm, FALSE)^2, 2, weighted.mean,
                      w = weights, na.rm = TRUE)/(N - 1))
    lamin=ncvgetmin(fit$lambda[err.ind],cvm,cvsd)
    val = list(cvm = cvm, cvsd = cvsd,lambda=fit$lambda[err.ind], fit=fit,
               lambda.min=lamin[[1]],lambda.1se=lamin[[2]])
    if (returnY) val$Y <- Y
    structure(val)
    #Y[cv.ind==i, 1:res$nl] <- res$yhat
    #E[cv.ind==i, 1:res$nl] <- res$loss
  }



cvf.surv <- function(i, XX, y, cv.ind, cv.args,weights=rep(1,nrow(XX))) {

  mark=which(cv.ind!=i)
  cv.args$X <- XX[mark, , drop=FALSE]
  cv.args$y <- y[mark,]
  n=nrow(cv.args$X)

  cv.args$weights=weights[mark]
  fit.i <- do.call("ncvsurv", cv.args)
  beta <- fit.i$beta


  list( nl=length(fit.i$lambda), beta = beta)
}

ncvgetmin=function(lambda,cvm,cvsd){
  cvmin=min(cvm,na.rm=TRUE)
  idmin=cvm<=cvmin
  lambda.min=max(lambda[idmin],na.rm=TRUE)
  idmin=match(lambda.min,lambda)
  semin=(cvm+cvsd)[idmin]
  idmin=cvm<=semin
  lambda.1se=max(lambda[idmin],na.rm=TRUE)
  list(lambda.min=lambda.min,lambda.1se=lambda.1se)
}

lamNames <- function(l) {
  if (length(l) > 1) {
    d <- ceiling(-log10(-max(diff(l))))
    d <- min(max(d,4), 10)
  } else {
    d <- 4
  }
  formatC(l, format="f", digits=d)
}


ncvsurv <- function(X, y, penalty=c("MCP", "SCAD", "lasso"), gamma=switch(penalty, SCAD=3.7, 3),
                    alpha=1, lambda.min=ifelse(n>p,.001,.05), nlambda=100, lambda, eps=1e-4, max.iter=10000,
                    convex=TRUE, dfmax=p, penalty.factor=rep(1, ncol(X)),weights=rep(1,length(y[,1])), warn=TRUE, returnX=FALSE, ...) {

  # Coersion
  penalty <- match.arg(penalty)
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (class(y) != "matrix") {
    tmp <- try(y <- as.matrix(y), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("y must be a matrix or able to be coerced to a matrix")
    if (ncol(y)!=2) stop("y must have two columns for survival data: time-on-study and a censoring indicator")
  }
  if (storage.mode(y)=="integer") storage.mode(y) <- "double"
  if (storage.mode(penalty.factor) != "double") storage.mode(penalty.factor) <- "double"

  ## Error checking
  if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
  if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
  if (nlambda < 2) stop("nlambda must be at least 2")
  if (alpha <= 0) stop("alpha must be greater than 0; choose a small positive number instead")
  if (length(penalty.factor)!=ncol(X)) stop("penalty.factor does not match up with X")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y to ncvreg")

  ## Set up XX, yy, lambda
  tOrder <- order(y[,1])
  yy <- as.numeric(y[tOrder,1])
  Delta <- y[tOrder,2]
  weights <- weights[tOrder]
  weights <- weights*length(weights)/sum(weights)
  n <- length(yy)
  XX <- std(X[tOrder,,drop=FALSE],weights)
  ns <- attr(XX, "nonsingular")
  penalty.factor <- penalty.factor[ns]
  p <- ncol(XX)
  if (missing(lambda)) {
    lambda <- setupLambdaCox(XX, yy, Delta, alpha, lambda.min, nlambda, penalty.factor, weights = weights)
    user.lambda <- FALSE
  } else {
    nlambda <- length(lambda)
    user.lambda <- TRUE
  }

  ## Fit
  res <- .Call("cdfit_cox_dh", XX, Delta, penalty, lambda, eps, as.integer(max.iter), as.double(gamma), penalty.factor,
               alpha, as.integer(dfmax), as.integer(user.lambda | any(penalty.factor==0)), as.integer(warn), weights)
  b <- matrix(res[[1]], p, nlambda)
  loss <- -1*res[[2]]
  iter <- res[[3]]
  Eta <- matrix(res[[4]], n, nlambda)

  ## Eliminate saturated lambda values, if any
  ind <- !is.na(iter)
  b <- b[, ind, drop=FALSE]
  iter <- iter[ind]
  lambda <- lambda[ind]
  loss <- loss[ind]
  Eta <- Eta[,ind,drop=FALSE]
  if (warn & sum(iter)==max.iter) warning("Algorithm failed to converge for some values of lambda")

  ## Local convexity?
  convex.min <- if (convex) convexMin(b, XX, penalty, gamma, lambda*(1-alpha), penalty.factor, Delta=Delta,weights=weights) else NULL

  ## Unstandardize
  beta <- matrix(0, nrow=ncol(X), ncol=length(lambda))
  bb <- b/attr(XX, "scale")[ns]
  beta[ns,] <- bb
  offset <- -crossprod(attr(XX, "center")[ns], bb)

  ## Names
  varnames <- if (is.null(colnames(X))) paste("V",1:ncol(X),sep="") else colnames(X)
  dimnames(beta) <- list(varnames, lamNames(lambda))

  ## Output
  val <- structure(list(beta = beta,
                        iter = iter,
                        lambda = lambda,
                        penalty = penalty,
                        gamma = gamma,
                        alpha = alpha,
                        convex.min = convex.min,
                        loss = loss,
                        penalty.factor = penalty.factor,
                        n = n,
                        time = yy,
                        fail = Delta,
                        order = tOrder))
  val$Eta <- sweep(Eta, 2, offset, "-")
  if (returnX) {
    val$X <- XX
    val$y <- yy
  }
  val
}


setupLambdaCox <- function(X, y, Delta, alpha, lambda.min, nlambda, penalty.factor, weights=rep(1,length(Delta))) {
  n <- nrow(X)
  p <- ncol(X)

  # Fit to unpenalized covariates
  ind <- which(penalty.factor!=0)
  if (length(ind)!=p) {
    nullFit <- coxph(Surv(y, Delta) ~ X[, -ind, drop=FALSE])
    eta <- nullFit$linear.predictors
    rsk <- rev(cumsum(rev(weights*exp(eta))))
    s <- Delta*weights - weights*exp(eta)*cumsum(Delta*weights/rsk)
  } else {
    w <- 1/(n-(1:n)+1)
    s <- Delta*weights - weights*cumsum(Delta*weights*w)
  }

  # Determine lambda.max
  zmax <- .Call("maxprod", X, s, ind, penalty.factor) / n
  lambda.max <- zmax/alpha

  if (lambda.min==0) lambda <- c(exp(seq(log(lambda.max),log(.001*lambda.max),len=nlambda-1)),0)
  else lambda <- exp(seq(log(lambda.max),log(lambda.min*lambda.max),len=nlambda))
  lambda
}


std <- function(X, weights) {
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  STD <- .Call("standardize", X, weights)
  dimnames(STD[[1]]) <- dimnames(X)
  ns <- which(STD[[3]] > 1e-6)
  if (length(ns) == ncol(X)) {
    val <- STD[[1]]
  } else {
    val <- STD[[1]][, ns, drop=FALSE]
  }
  attr(val, "center") <- STD[[2]]
  attr(val, "scale") <- STD[[3]]
  attr(val, "nonsingular") <- ns
  val
}

