### to do ###
# 1. aabc: select the model with only aabc > 0, and max

#' Function for AUC when input is X and Y.
#' @title AUC
#' @description Empirical AUC estimate
#' @param outcome binary outcome (1: desired outcome; 0: otherwise)
#' @param predict prediction score
#' @return a numeric value of empirical estimation of area under the ROC curves
#' @examples
#' # no run

AUC <- function(outcome, predict){
  if(length(predict)!=length(outcome)) stop("predict and outcome should vectors with the same length")
  Y <- predict[outcome==1]
  X <- predict[outcome==0]
  if(stats::sd(predict)==0){
    return(0.5)
  }else{
    return(mean(outer(Y,X, FUN=function(y,x) 1*(y>x))))
  }
}

#' Function for c-index when input is X and Y.
#' @title C.index
#' @description Empirical c-index estimate
#' @param yvar column name for observed time
#' @param score column name for marker value
#' @param censorvar column name for censor (1 is event, 0 is censored)
#' @param data input data matrix
#' @return a numeric value of empirical estimation of c-index
#' @examples
#' # no run

C.index <- function(yvar, score, censorvar, data){
  Cs <- data[, censorvar, drop=FALSE]
  Z <- outer(data[, score], data[, score], '-')
  Z <- Z[lower.tri(Z)]
  Ts <- outer(data[, yvar], data[, yvar], '-')
  Ts <- Ts[lower.tri(Ts)]
  Cs.inx <- outer(rep(1, nrow(data)), data[, censorvar])
  dels <- Cs.inx[lower.tri(Cs.inx)]
  sum(1*(Z<0)*(Ts>0)*dels)/(sum(1*(Ts>0)*dels))
}

#' Function for generate CV fold index
#' @title cvfolds0
#' @description internal function for generating CV fold index
#' @param X marker matrix for non-responders
#' @param Y marker matrix for responders
#' @param idx m*n by 2 matrix for row index of marker matrix, first column is row index in X; second column is for Y
#' @param nfolds the cross-validation folds
#' @return a vector containing CV fold index for each row in Z
#' @examples
#' # no run

cvfolds0 <- function(X, Y, idx, nfolds=5){
  fold.X <- sample.int(nfolds, nrow(X), replace = TRUE)
  fold.Y <- sample.int(nfolds, nrow(Y), replace = TRUE)
  fold.Z <- ifelse(fold.X[idx[,1]] == fold.Y[idx[,2]],
                   fold.X[idx[,1]], NA)
  fold.Z
}

#' Function for generate Z matrix for binary endpoint
#' @title pair.diff
#' @description internal function for generating Z matrix (binary endpoint)
#' @param X marker matrix for non-responders
#' @param Y marker matrix for responders
#' @param A Treatment arm indicator (1 is treatment, 0 is control)
#' @return A list of prepared data input for ROCSI
#' @examples
#' # no run

pair.diff <- function(X, Y, A){
  X.trt <- X[A==1 & Y==0,, drop=FALSE]
  Y.trt <- X[A==1 & Y==1,, drop=FALSE]
  X.ctl <- X[A==0 & Y==0,, drop=FALSE]
  Y.ctl <- X[A==0 & Y==1,, drop=FALSE]
  idx.trt <- data.frame(xidx=rep(1:nrow(X.trt), nrow(Y.trt)), yidx=rep(1:nrow(Y.trt), each=nrow(X.trt)))
  idx.ctl <- data.frame(xidx=rep(1:nrow(X.ctl), nrow(Y.ctl)), yidx=rep(1:nrow(Y.ctl), each=nrow(X.ctl)))
  Z.trt <- X.trt[idx.trt[,1],] - Y.trt[idx.trt[,2], ]
  Z.trt <- cbind(1,Z.trt)
  W.trt <- rep(1/(nrow(X.trt)*nrow(Y.trt)), nrow(Z.trt))
  Z.ctl <- X.ctl[idx.ctl[,1],] - Y.ctl[idx.ctl[,2], ]
  Z.ctl <- cbind(0,Z.ctl)
  W.ctl <- rep(1/(nrow(X.ctl)*nrow(Y.ctl)), nrow(Z.ctl))
  colnames(Z.trt)[1] <- colnames(Z.ctl)[1] <- "trt"
  Z <- rbind(Z.trt, Z.ctl)
  Z <- as.matrix(Z)
  W <- c(W.trt, W.ctl)
  list(Z=Z, W=W, Z.trt=Z.trt, Z.ctl=Z.ctl, X.trt=X.trt, Y.trt=Y.trt, X.ctl=X.ctl, Y.ctl=Y.ctl, idx.trt=idx.trt, idx.ctl=idx.ctl)
}

#' Function for generate Z matrix for time-to-event endpoint
#' @title pair.diff.surv
#' @description internal function for generating Z matrix (time-to-event endpoint)
#' @param X marker matrix
#' @param Y a vector for observed time
#' @param A a vector for Treatment arm indicator (1 is treatment, 0 is control)
#' @param C a vector for censor (1 is event, 0 is censored)
#' @return A list of prepared data input for ROCSI
#' @examples
#' # no run

pair.diff.surv <- function(X, Y, A, C){
  n.trt <- nrow(X[A==1, , drop=FALSE])
  idx.trt <- outer(c(1:n.trt), c(1:n.trt), 'paste', sep=",")
  idx.trt <- idx.trt[lower.tri(idx.trt)]
  idx.trt <- scan(text = idx.trt, what = integer(), sep = ',', quiet = TRUE)
  idx.trt <- matrix(idx.trt, length(idx.trt)/2, byrow = TRUE)

  Zi.trt <- Zj.trt <- X[A==1, , drop=FALSE]
  Ti.trt <- Tj.trt <- Y[A==1]
  Cs.trt <- C[A==1]
  Z.trt<-Zj.trt[idx.trt[,1], ,drop=F] - Zi.trt[idx.trt[,2], ,drop=F]
  Ts.trt <- Ti.trt[idx.trt[,1]] - Tj.trt[idx.trt[,2]]
  Cs.inx.trt <- outer(rep(1, length(Cs.trt)), Cs.trt)
  dels.trt <- Cs.inx.trt[lower.tri(Cs.inx.trt)]
  weights.trt <- 1*(Ts.trt>0)*dels.trt
  Z.trt <- -Z.trt[weights.trt==1,] # need to check the sign of Z
  W.trt <- rep(1/sum(weights.trt), nrow(Z.trt))

  n.ctl <- nrow(X[A==0, , drop=FALSE])
  idx.ctl <- outer(c(1:n.ctl), c(1:n.ctl), 'paste', sep=",")
  idx.ctl <- idx.ctl[lower.tri(idx.ctl)]
  idx.ctl <- scan(text = idx.ctl, what = integer(), sep = ',', quiet = TRUE)
  idx.ctl <- matrix(idx.ctl, length(idx.ctl)/2, byrow = TRUE)

  Zi.ctl <- Zj.ctl <- X[A==0, , drop=FALSE]
  Ti.ctl <- Tj.ctl <- Y[A==0]
  Cs.ctl <- C[A==0]
  Z.ctl<-Zj.ctl[idx.ctl[,1], ,drop=F] - Zi.ctl[idx.ctl[,2], ,drop=F]
  Ts.ctl <- Ti.ctl[idx.ctl[,1]] - Tj.ctl[idx.ctl[,2]]
  Cs.inx.ctl <- outer(rep(1, length(Cs.ctl)), Cs.ctl)
  dels.ctl <- Cs.inx.ctl[lower.tri(Cs.inx.ctl)]
  weights.ctl <- 1*(Ts.ctl>0)*dels.ctl
  Z.ctl <- -Z.ctl[weights.ctl==1,] # need to check the sign of Z
  W.ctl <- rep(1/sum(weights.ctl), nrow(Z.ctl))

  Z.trt <- cbind(1,Z.trt)
  Z.ctl <- cbind(0,Z.ctl)
  colnames(Z.trt)[1] <- colnames(Z.ctl)[1] <- "trt"
  Z <- rbind(Z.trt, Z.ctl)
  Z <- as.matrix(Z)
  W <- c(W.trt, W.ctl)
  list(Z=Z, W=W, Z.trt=Z.trt, Z.ctl=Z.ctl, Zi.trt=Zi.trt,  Zj.trt=Zj.trt, Zi.ctl=Zi.ctl, Zj.ctl=Zj.ctl, idx.trt=idx.trt[weights.trt==1,], idx.ctl=idx.ctl[weights.ctl==1,])
}

#' Function for HIC calculation
#' @title HIC
#' @description function for HIC calculation
#' @param beta  estimates of coefficient beta
#' @param Z matrix prepared for ROCSI
#' @param index m*n by 2 matrix  for the subindex for the pair difference in Z
#' @param w a vector of weights Z (can be used for inverse probability weighting for missing data, default is 1)
#' @return A numeric value with corresponding HIC
#' @examples
#' # no run

HIC <- function(beta, Z, index, w=1){
  #  inv<-try(ginv(hessAUC(beta, Z)),silent=T)
  inv<-try(solve(hessAUC(beta, Z, w=w)),silent=T)
  if(inherits(inv, "try-error"))
  {
    return(NA)
  } else {
    return(sum(diag(inv%*%gradsqr(beta, Z, index, w=w))))
  }
}

#' Internal function for HIC calculation
#' @title gradsqr
#' @description Internal function for HIC calculation
#' @param beta  estimates of coefficient beta
#' @param Z0 (m x n) x p Z matrix as prepared for ROCSI
#' @param index m*n by 2 matrix  for the subindex for the pair difference in Z
#' @param w a vector of weights Z (can be used for inverse probability weighting for missing data, default is 1)
#' @return gradient square for the GCV.
#' @examples
#' # no run
gradsqr <- function(beta, Z0, index, w=1){
  n.Z0 <- nrow(Z0)
  k.Z0 <- ncol(Z0)
  tmp <- w*apply(Z0, MARGIN=1, FUN=grad.sub, beta)
  tmp <- matrix(tmp, nrow=k.Z0)
  index.x <- index[,1]
  an0 <- cbind(t(tmp), index.x)
  an1 <- do.call("rbind", by(an0[,1:k.Z0, drop=F], an0[,"index.x"], FUN=colSums, simplify = F))
  an2 <- apply(an1, MARGIN=1, FUN=function(x) x%*%t(x))
  an2 <- matrix(an2, nrow=k.Z0^2)
  a1 <- matrix(rowSums(an2), k.Z0, k.Z0)
  index.y <- index[,2]
  am0 <- cbind(t(tmp), index.y)
  am1 <- do.call("rbind", by(am0[,1:k.Z0, drop=F], am0[,"index.y"], FUN=colSums, simplify = F))
  am2 <- apply(am1, MARGIN=1, FUN=function(x) x%*%t(x))
  am2 <- matrix(am2, nrow=k.Z0^2)
  a2 <- matrix(rowSums(am2), k.Z0, k.Z0)
  tmp.a <- apply(tmp, MARGIN=2, FUN=function(x) x%*%t(x))
  tmp.a <- matrix(tmp.a, nrow=k.Z0^2)
  a <- matrix(rowSums(tmp.a), k.Z0, k.Z0)
  a0 <- a1+a2-a
  a0
}

#' Internal function of grad_square in the GCV
#' @title grad.sub
#' @description Internal function of grad_square in the GCV
#' @param z  (m x n) x p data matrix as prepared for ROCSI
#' @param beta estimates of coefficient beta
#' @return grad_square in the GCV
#' @examples
#' # no run
grad.sub <- function(z, beta){
  (z)*(-1/(1+exp(beta*(z))))
}

#' Internal function for hessAUC
#' @title hessAUC.sub
#' @description Internal function for hessAUC
#' @param z  (m x n) x p data matrix as prepared for ROCSI
#' @param beta estimates of coefficient beta
#' @return Hessian matrix components.
#' @examples
#' # no run
hessAUC.sub <- function(z, beta){
  z%*%t(z)*c((exp(t(beta)%*%(z))/(1+exp(t(beta)%*%(z)))^2))
}

#' function for Hessian matrix of AUC
#' @title hessAUC
#' @description function for Hessian matrix of AUC
#' @param beta estimates of coefficient beta
#' @param Z  (m x n) x p data matrix as prepared for ROCSI
#' @param w a vector of weights Z (can be used for inverse probability weighting for missing data, default is 1)
#' @return Hessian matrix of AUC.
#' @examples
#' # no run

hessAUC <- function(beta, Z, w=1){
  tmp <- w*apply(Z, MARGIN=1, FUN=hessAUC.sub, beta)
  tmp <- matrix(tmp, nrow=ncol(Z)^2)
  matrix(rowSums(tmp), ncol(Z), ncol(Z))
}


#' function for ROCSI
#' @title ROCSI
#' @description function for ROCSI
#' @param Dtrain data matrix for training dataset
#' @param Dtest optional data matrix for testing dataset
#' @param yvar column name for outcome
#' @param xvars a string vector of column names for input markers
#' @param trtvar column name for treatment (the column should contain binary code with 1 being treatment and 0 being control)
#' @param cvar column name for censor (the column should contain binary code with 1 being event and 0 being censored)
#' @param nfolds n fold CV used for cv.glmnet
#' @param type outcome type ("binary" for binary outcome and "survival" for time-to-event outcome)
#' @return A list with ROCSI output
#' \describe{
#'   \item{beta.aABC}{final beta estimated from ROCSI based on \eqn{ABC^{(acv)}}}
#'   \item{beta.1se}{final beta estimated from lambda.1se based on nfold CV}
#'   \item{lambda.aABC}{optimal lambda selected by optimizing \eqn{ABC^{(acv)}}}
#'   \item{fit.cv}{fitted cv.glmnet model}
#'   \item{log}{log matrix of all lambdas and ABCs}
#'   \item{abc.test}{ABC in testing dataset based on optimal beta}
#'   \item{abc.test1se}{ABC in testing dataset based on 1se beta}
#'   \item{predScore}{a data.frame of testing data and its predictive signature scores (based on beta.aABC) for each subjects}
#'   \item{predScore.1se}{a data.frame of testing data and its predictive signature scores (based on beta.1se) for each subjects}
#' }
#' @examples
#' n <- 100
#' k <- 5
#' prevalence <- sqrt(0.5)
#' rho<-0.2
#' sig2 <- 2
#' rhos.bt.real <- c(0, rep(0.1, (k-3)))*sig2
#' y.sig2 <- 1
#' yvar="y.binary"
#' xvars=paste("x", c(1:k), sep="")
#' trtvar="treatment"
#' prog.eff <- 0.5
#' effect.size <- 1
#' a.constent <- effect.size/(2*(1-prevalence))
#' ObsData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                     sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                     rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#' TestData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                      sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                      rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#'bst.aabc <- ROCSI(Dtrain=ObsData$data, Dtest = TestData$data, yvar=yvar,
#'xvars=xvars, trtvar=trtvar, cvar=NULL, nfolds=5, type="binary")
#'bst.aabc$beta.aABC
#'bst.aabc$log
#'bst.aabc$abc.test
#'bst.aabc$beta.1se
#'bst.aabc$abc.test1se
#' @export
ROCSI <- function(Dtrain, Dtest=NULL, yvar, xvars, trtvar, cvar=NULL, nfolds=5, type="binary"){
  coef.cv.glmnet <- utils::getFromNamespace("coef.cv.glmnet", "glmnet")
  predict.cv.glmnet <- utils::getFromNamespace("predict.cv.glmnet", "glmnet")
  if(type=="binary"){
    Xtmp <- pair.diff(X=Dtrain[,xvars, drop=FALSE], Y=Dtrain[,yvar], A=Dtrain[,trtvar])
    X.trt <- Xtmp$X.trt
    Y.trt <- Xtmp$Y.trt
    X.ctl <- Xtmp$X.ctl
    Y.ctl <- Xtmp$Y.ctl
    Z.trt <- Xtmp$Z.trt
    Z.ctl <- Xtmp$Z.ctl
    W <- Xtmp$W
    Z <- rbind(Z.trt, -Z.ctl)
    idx.trt <- Xtmp$idx.trt
    idx.ctl <- Xtmp$idx.ctl
    fold.Z.trt <- cvfolds0(X=X.trt, Y=Y.trt, idx=idx.trt, nfolds=nfolds)
    fold.Z.ctl <- cvfolds0(X=X.ctl, Y=Y.ctl, idx=idx.ctl, nfolds=nfolds)
    fold.Z <- c(fold.Z.trt, fold.Z.ctl)
  }else if(type=="survival"){
    Xtmp <- pair.diff.surv(X=Dtrain[,xvars, drop=FALSE], Y=Dtrain[,yvar], A=Dtrain[,trtvar], C=Dtrain[,cvar])
    Z.trt <- Xtmp$Z.trt
    Z.ctl <- Xtmp$Z.ctl
    Z <- rbind(Z.trt, -Z.ctl)
    Zi.trt <- Xtmp$Zi.trt
    Zj.trt <- Xtmp$Zj.trt
    Zi.ctl <- Xtmp$Zi.ctl
    Zj.ctl <- Xtmp$Zj.ctl
    idx.trt <- Xtmp$idx.trt
    idx.ctl <- Xtmp$idx.ctl
    W <- Xtmp$W
    fold.Z.trt <- cvfolds0(X=Zi.trt, Y=Zj.trt, idx=idx.trt, nfolds=nfolds)
    fold.Z.ctl <- cvfolds0(X=Zi.ctl, Y=Zj.ctl, idx=idx.ctl, nfolds=nfolds)
    fold.Z <- c(fold.Z.trt, fold.Z.ctl)
  }else{
    stop("Only binary or survival types allowed")
  }
  yb.sim<-stats::rbinom(nrow(Z),1,0.3)
  data.new<-cbind(yb.sim,Z)
  Y.new<-data.new[yb.sim==1,]
  W1 <- W[yb.sim==1]
  X.new<-data.new[yb.sim==0,]
  W0 <- W[yb.sim==0]
  X.new[,-1]<--X.new[,-1]
  data.fin<-rbind(Y.new,X.new)
  W.fin <- c(W1, W0)
  data.fin <- as.data.frame(data.fin)
  y=data.fin[,"yb.sim"]
  x=as.matrix(data.fin[,xvars])
  cv.idx <- !is.na(fold.Z)
  fit.cv <- glmnet::cv.glmnet(x[cv.idx,], y[cv.idx], weights=W.fin[cv.idx], family="binomial", foldid=fold.Z[cv.idx], intercept = FALSE, type.measure="auc", alpha=0.5, parallel = TRUE)
  #fit.cv$lambda.min
  #plot(fit.cv)
  #fit.cv$glmnet.fit$beta
  #beta.1se <- as.matrix(coef(fit.cv, s = "lambda.1se")[-1,])

  ## get HIC for each glmnet solution path ###
  fit.cv.nzero.idx<-which(!duplicated(fit.cv$nzero))[-1]
  lambda.sel <- c(fit.cv$lambda[fit.cv.nzero.idx])
  fit.cv.nzero.idx <- c(1, fit.cv.nzero.idx)
  lambda.sel <- c(fit.cv$lambda[1], lambda.sel)
  fit.cv.nzero.mat<-as.matrix(fit.cv$glmnet.fit$beta[,fit.cv.nzero.idx])
  betahat.mat <- fit.cv.nzero.mat

  aABC.mat <- data.frame(lambda=lambda.sel, abc.train=NA, bias=NA)
  for(i in 1:length(lambda.sel)){
    beta.idx <- which(betahat.mat[,i]!=0)
    beta.hat <- betahat.mat[beta.idx,i]
    var.sel <- row.names(betahat.mat)[beta.idx]
    if(type=="binary"){
      pred0 <- predict.cv.glmnet(fit.cv, newx = as.matrix(Dtrain[,xvars, drop=FALSE]), s=lambda.sel[i])
      abc.train <- AUC(Dtrain[,yvar][Dtrain[,trtvar]==1], pred0[Dtrain[,trtvar]==1]) -
        AUC(Dtrain[,yvar][Dtrain[,trtvar]==0], pred0[Dtrain[,trtvar]==0])
    }else if(type=="survival"){
      score <- c(predict.cv.glmnet(fit.cv, newx = as.matrix(Dtrain[,xvars, drop=FALSE]), s=lambda.sel[i]))
      pred0 <- data.frame(yvar=Dtrain[,yvar], censorvar=Dtrain[,cvar], score=score, trt=Dtrain[,trtvar])
      abc.train <- C.index("yvar", "score", "censorvar", data=pred0[pred0[,"trt"]==1,]) -
        C.index("yvar", "score", "censorvar", data=pred0[pred0[,"trt"]==0,])
    }else{
      stop("Only binary or survival types allowed")
    }
    # if(length(beta.idx)==1){
    #   correct <- 0
    # }else{
    #   correct <- try(abs(HIC(beta=beta.hat, Z=Z.trt[,var.sel, drop=F], index=idx.trt)) +
    #                  abs(HIC(beta=beta.hat, Z=Z.ctl[,var.sel, drop=F], index=idx.ctl)),silent=T)/4 # divided by 4 only correct for 1:1 treatment ratio
    # }
    #    correct <- try(abs(HIC(beta=beta.hat, Z=Z.trt[,var.sel, drop=F], index=idx.trt)) +
    #                     abs(HIC(beta=beta.hat, Z=Z.ctl[,var.sel, drop=F], index=idx.ctl)),silent=T) # divided by 4 only correct for 1:1 treatment ratio
    if(length(beta.hat)>0){
      beta.hat0 <- ifelse(length(beta.hat)==1, beta.hat, theta2beta(beta2theta(beta.hat)))
      correct <- try(abs(HIC(beta=beta.hat0, Z=rbind(Z.trt[,var.sel, drop=F], Z.ctl[,var.sel, drop=F]), index=rbind(idx.trt, idx.ctl),w=W))) # try weight options 1 or W
    }else(
      correct <- 0
    )
    aABC.mat$abc.train[i] <- abs(abc.train)
    aABC.mat$bias[i] <- correct
  }

 # if(!all(aABC.mat$bias == cummax(aABC.mat$bias))){ aABC.mat$bias <- cummax(aABC.mat$bias) } # make sure bias is monotone increasing

  aABC.mat$aABC <- aABC.mat$abc.train - aABC.mat$bias

  lambda.aABC <- aABC.mat$lambda[which.max(aABC.mat$aABC)]
  beta.aABC <- betahat.mat[,which.max(aABC.mat$aABC)]


  # perfrmance on the test data
  if(!is.null(Dtest)){
    if(type=="binary"){

      pred1 <- predict.cv.glmnet(fit.cv, newx = as.matrix(Dtest[,xvars, drop=FALSE]), s=lambda.aABC)
      abc.test <- AUC(Dtest[,yvar][Dtest[,trtvar]==1], pred1[Dtest[,trtvar]==1]) -
        AUC(Dtest[,yvar][Dtest[,trtvar]==0], pred1[Dtest[,trtvar]==0])

      pred1se <- predict.cv.glmnet(fit.cv, newx = as.matrix(Dtest[,xvars, drop=FALSE]), s='lambda.1se')
      abc.test1se <- AUC(Dtest[,yvar][Dtest[,trtvar]==1], pred1se[Dtest[,trtvar]==1]) -
        AUC(Dtest[,yvar][Dtest[,trtvar]==0], pred1se[Dtest[,trtvar]==0])
      beta1se <- as.matrix(coef.cv.glmnet(fit.cv, s = "lambda.1se")[-1,])
    }else if(type=="survival"){
      score1 = c(predict.cv.glmnet(fit.cv, newx=as.matrix(Dtest[,xvars, drop=FALSE]), s=lambda.aABC))
      pred1 <- data.frame(yvar=Dtest[,yvar], censorvar=Dtest[,cvar], score=score1, trt=Dtest[,trtvar])
      abc.test <- C.index("yvar", "score", "censorvar", data=pred1[pred1[,"trt"]==1,]) -
        C.index("yvar", "score", "censorvar", data=pred1[pred1[,"trt"]==0,])

      score1se <- c(predict.cv.glmnet(fit.cv, newx=as.matrix(Dtest[,xvars, drop=FALSE]), s='lambda.1se'))
      pred1se <- data.frame(yvar=Dtest[,yvar], censorvar=Dtest[,cvar], score=score1se, trt=Dtest[,trtvar])
      abc.test1se <- C.index("yvar", "score", "censorvar", data=pred1se[pred1se[,"trt"]==1,]) -
        C.index("yvar", "score", "censorvar", data=pred1se[pred1se[,"trt"]==0,])
      beta1se <- as.matrix(coef.cv.glmnet(fit.cv, s = "lambda.1se")[-1,])

    }else{
      stop("Only binary or survival types allowed")
    }
  }else{
    abc.test<-abc.test1se<-0
    pred1 <- beta1se <- pred1se <- NULL
  }

  list(beta.aABC=beta.aABC, beta.1se=beta1se, lambda.aABC=lambda.aABC, fit.cv=fit.cv, log=aABC.mat, abc.test=abs(abc.test), abc.test1se=abs(abc.test1se), predScore=pred1, predScore.1se=pred1se)
}


#' Function to translate beta into theta, the n-sphere constrain
#' @title beta2theta
#' @description Function to translate beta into theta, the n-sphere constrain
#' @param beta estimates of coefficient beta
#' @return a numeric vector for theta (dimension-1)
#' @examples
#' # no run

beta2theta <- function(beta){
  dim <- length(beta)
  theta<-rep(0,dim-1)
  for (i in 1:dim-1){
    theta[i]<-atan(sqrt(sum(beta[c((i+1):dim)]^2))/beta[i])
  }
  theta
}

### Function to translate theta to beta, the n-sphere constrain ###
#' Function to translate beta into theta, the n-sphere constrain
#' @title theta2beta
#' @description Function to translate theta into beta
#' @param theta n-sphere coordination
#' @return a numeric vector for beta (dimension+1)
#' @examples
#' # no run

theta2beta <- function(theta){
  dim<-length(theta)+1
  beta<-matrix(1,dim,1)
  ### assign the n-sphere coordination #######
  for (k in 1: dim){
    if(k==1){
      beta[1,1]<-cos(theta[1])
    }else if (k== dim){
      for(ii in 1:(k-1)){
        beta[k,1]<-beta[k,1]*sin(theta[ii])}
    }else {
      for(ii in 1:(k-1)){
        beta[k,1]<-beta[k,1]*sin(theta[ii])}
      beta[k,1]<-beta[k,1]*cos(theta[k])}
  }
  beta
}


###############################################################
### function for modified covariate methods based on glmnet ###
#' function for ROCSI
#' @title MClogit
#' @description function for modified covariate methods based on glmnet
#' @param dataset data matrix for training dataset
#' @param yvar column name for outcome
#' @param xvars a string vector of column names for input markers
#' @param trtvar column name for treatment (the column should contain binary code with 1 being treatment and 0 being control)
#' @param cvar column name for censor (the column should contain binary code with 1 being event and 0 being censored)
#' @param nfolds n fold CV used for cv.glmnet
#' @param type outcome type ("binary" for binary outcome and "survival" for time-to-event outcome)
#' @param newx data matrix for testing dataset X
#' @param bestsub criteria for best lambda, used by glmnet
#' @param type.measure type of measure used by glmnet
#' @return A list with ROCSI output
#' \describe{
#'   \item{x.logit}{final beta estimated from MClogit}
#'   \item{predScore}{a data.frame of testing data and its predictive signature scores (based on beta.aABC) for each subjects}
#'   \item{abc}{ABC in testing dataset based on optimal beta}
#'   \item{fit.cv}{the fitted glmnet object}
#' }
#' @examples
#' n <- 100
#' k <- 5
#' prevalence <- sqrt(0.5)
#' rho<-0.2
#' sig2 <- 2
#' rhos.bt.real <- c(0, rep(0.1, (k-3)))*sig2
#' y.sig2 <- 1
#' yvar="y.binary"
#' xvars=paste("x", c(1:k), sep="")
#' trtvar="treatment"
#' prog.eff <- 0.5
#' effect.size <- 1
#' a.constent <- effect.size/(2*(1-prevalence))
#' ObsData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                     sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                     rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#' TestData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                      sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                      rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#' bst.mod <- MClogit(dataset=ObsData$data, yvar=yvar, xvars=xvars,
#' trtvar=trtvar, nfolds = 5, newx=TestData$data,
#' type="binary", bestsub="lambda.1se")
#' bst.mod$abc
#' bst.mod$x.logit[-1,1]
#' @export
MClogit <- function(dataset, yvar, xvars, trtvar, cvar=NULL, nfolds = 5, type="binary", newx=NULL, bestsub="lambda.1se", type.measure = "auc"){
  coef.cv.glmnet <- utils::getFromNamespace("coef.cv.glmnet", "glmnet")
  predict.cv.glmnet <- utils::getFromNamespace("predict.cv.glmnet", "glmnet")
  t.mod <- (dataset[,trtvar]==1)*1 - (dataset[,trtvar]==0)*1
  data.bin <- data.frame(responder=dataset[,yvar], t.mod*dataset[,xvars])

  if(type=="binary"){
    cvglmnet.fit <- glmnet::cv.glmnet(x=as.matrix(data.bin[,xvars]), y=data.bin$responder,
                              family = "binomial", intercept = FALSE, type.measure = type.measure)
    x.logit <- coef.cv.glmnet(cvglmnet.fit, s = bestsub)
    if(!is.null(newx)){
      pred1 <- predict.cv.glmnet(cvglmnet.fit, newx = as.matrix(newx[,xvars]), s=bestsub)
      abc <- AUC(newx[,yvar][newx[,trtvar]==1], pred1[newx[,trtvar]==1]) -
        AUC(newx[,yvar][newx[,trtvar]==0], pred1[newx[,trtvar]==0])
    }else{
      abc <- NULL
    }

  }else if(type=="survival"){
    y <- cbind(time=dataset[,yvar],status=dataset[,cvar])
    cvglmnet.fit <- glmnet::cv.glmnet(x=as.matrix(data.bin[,xvars]), y=y, family = "cox", type.measure = type.measure)
    x.logit <- coef.cv.glmnet(cvglmnet.fit, s = bestsub)
    if(!is.null(newx)){
      score = predict.cv.glmnet(cvglmnet.fit, newx = as.matrix(newx[,xvars]), s=bestsub)
      pred1 <- data.frame(yvar=newx[,yvar], censorvar=newx[,cvar], score=c(score), trt=newx[,trtvar])
      abc <- C.index("yvar", "score", "censorvar", data=pred1[pred1[,"trt"]==1,]) -
        C.index("yvar", "score", "censorvar", data=pred1[pred1[,"trt"]==0,])
    }else{
      abc <- NULL
    }
  }else{
    stop("Only binary or survival types allowed")
  }
  list(x.logit=x.logit, predScore=pred1, abc=abs(abc), fit.cv=cvglmnet.fit)
}


#' Function for simulated data generation
#' @title data.gen
#' @description Function for simulated data generation
#' @param n Total sample size
#' @param k Number of markers
#' @param prevalence prevalence of predictive biomarkers with values above the cutoff
#' @param prog.eff  effect size \eqn{beta} for prognostic biomarker
#' @param sig2 standard deviation of each marker
#' @param y.sig2 Standard Deviation of the error term in the linear component
#' @param rho rho*sig2 is the entries for covariance matrix between pairs of different k markers
#' @param rhos.bt.real correlation between each prognostic and predictive markers
#' @param a.constent a constant is set such that there is no overall treatment effect
#' @return A list of simulated clinical trial data with heterogeneous prognostic and predictive biomarkers
#' @examples
#' n <- 500
#' k <- 10
#' prevalence <- sqrt(0.5)
#' rho<-0.2
#' sig2 <- 2
#' rhos.bt.real <- c(0, rep(0.1, (k-3)))*sig2
#' y.sig2 <- 1
#' prog.eff <- 0.5
#' effect.size <- 1
#' a.constent <- effect.size/(2*(1-prevalence))
#' ObsData <- data.gen(n=n, k=k, prevalence=prevalence, prog.eff=prog.eff,
#'                     sig2=sig2, y.sig2=y.sig2, rho=rho,
#'                     rhos.bt.real=rhos.bt.real, a.constent=a.constent)
#' @export
data.gen <- function(n, k, prevalence=sqrt(0.5), prog.eff=1, sig2, y.sig2, rho, rhos.bt.real, a.constent){
  covm <- matrix(rho*sig2,k,k)
  diag(covm) <- sig2
  covm[1,3:k] <- rhos.bt.real
  covm[2,3:k] <- rhos.bt.real
  covm[3:k,1] <- rhos.bt.real
  covm[3:k,2] <- rhos.bt.real
  dich.cutoff <- stats::qnorm(prevalence, sd = sqrt(sig2))
  x <- MASS::mvrnorm(n, rep(0,k), covm)
  x.dich <- 1*(x < dich.cutoff)
  w <- x.dich[,3:5]
  trt <- stats::rbinom(n, 1, prob=0.5)
  #trt <- rbinom(n, size = 1, prob = plogis(-1+0.05*w[,1] + 0.25*w[,2] + 0.6*w[,3] + 0.4*w[,1]*w[,3])) # confounding variables

  # predictive effect: x1, x2, prognostic effect x3
  prog.part <- prog.eff*w[,1] + stats::rnorm(n, sd=sqrt(y.sig2))
  pred.part <- a.constent*(-2*prevalence + x.dich[,1] + x.dich[,2]) # simulation based on binary biomarker
  #prog.part <- x[,3] + rnorm(n, sd=sqrt(y.sig2))
  #pred.part <- a.constent*(-2*prevalence + x[,1] + x[,2]) # simulation based on continous biomarker
  y.0 <- pred.part*0 + prog.part
  y.1 <- pred.part*1 + prog.part

  # continuous outcome
  y <- y.1*trt + y.0*(1-trt)

  # create binary outcome
  y.binary.0 <- 1*(stats::plogis(y.0)>prevalence)
  y.binary.1 <- 1*(stats::plogis(y.1)>prevalence)
  y.binary <- y.binary.1*trt + y.binary.0*(1-trt)

  # create time-to-event outcome
  surv.time.0 <- exp(y.0)
  surv.time.1 <- exp(y.1)
  cens.time <- exp(stats::rnorm(n, sd = 3))
  y.time.to.event.0 <- pmin(surv.time.0, cens.time)
  y.time.to.event.1 <- pmin(surv.time.1, cens.time)
  y.time.to.event <- y.time.to.event.1*trt + y.time.to.event.0*(1-trt)
  status.0 <- 1*(surv.time.0 <= cens.time)
  status.1 <- 1*(surv.time.1 <= cens.time)
  status <- status.1*trt + status.0*(1-trt)

  data <- cbind(y, y.binary, y.time.to.event, status,
                y.0, y.1, y.binary.0, y.binary.1,
                y.time.to.event.0, y.time.to.event.1,
                status.0, status.1, trt, x)
  colnames(data) <- c("y", "y.binary", "y.time.to.event", "status",
                      "y.0", "y.1", "y.binary.0", "y.binary.1",
                      "y.time.to.event.0", "y.time.to.event.1",
                      "status.0", "status.1","treatment", sapply(c(1:k), FUN=function(x) paste("x", x, sep="")))
  colnames(x) <- sapply(c(1:k), FUN=function(x) paste("x", x, sep=""))
  colnames(w) <- sapply(c(1:3), FUN=function(x) paste("w", x, sep=""))
  list(data=data, y=y, y.binary = y.binary,
       y.time.to.event = y.time.to.event, status = status,
       y.0=y.0, y.1=y.1, y.binary.0=y.binary.0, y.binary.1=y.binary.1,
       y.time.to.event.0=y.time.to.event.0, y.time.to.event.1=y.time.to.event.1,
       status.0=status.0, status.1=status.1, trt=trt, x=x.dich, x.continues=x,
       w=w, sigpos=(x.dich[,1]==1&x.dich[,2]==1))
}
