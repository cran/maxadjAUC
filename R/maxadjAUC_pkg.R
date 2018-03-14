##########
### maxadjAUC
### Code to maximize adjusted AUC
##########

############ Define helper functions
globalVariables('ind')

saauc <- function(betavec, Xmat, yvec, hc){
  wts <- prop.table(table(yvec, Xmat[,1]),1)[2,]
  costval <- rep(NA,length(wts))

  betadiffs <- list();
  for(i in 1:length(wts)){
    XDc <- Xmat[(Xmat[,1]==i & yvec==1),]
    XnDc <- Xmat[(Xmat[,1]==i & yvec==0),]
    diffs <- matrix(NA, nrow=nrow(XDc)*nrow(XnDc), ncol=ncol(Xmat)-1)
    for(j in 1:ncol(diffs)){
      diffs[,j] <- as.vector(outer(XDc[,j+1],XnDc[,j+1],"-"))
    }
    diffsB <- diffs %*% betavec
    costval[i] = -sum(pnorm(diffsB/hc[i]))/nrow(diffs)
  }

  totcost = (costval %*% wts)
  as.numeric(totcost)
}

"%ni%" <- Negate("%in%")

constr <- function(betavec, Xmat, yvec, hc){
  norm(matrix(betavec), type="F")
}
####################

maxadjAUC <- function(outcome, predictors, covariate, initialval="rGLM", approxh = 1/3, conditional=FALSE, tolval = 1e-6, stepsz = 1e-5){

  ###### Give error if any observations missing
  if(max(is.na(outcome))==1 | max(is.na(predictors))==1 | max(is.na(covariate))==1){
    stop("missing values not allowed")
  }
  ###### Give error if wrong variable types
  if(!is.matrix(predictors)){
    stop("predictors must be a matrix")
  }
  if(!is.numeric(covariate) | !is.numeric(outcome) | !is.numeric(predictors)){
    stop("outcome, covariate, and predictors must be numeric")
  }
  ###### Give error if outcome is not 0/1
  outcomevalue <- unique(outcome)
  if((length(outcomevalue) != 2) | (min(outcomevalue) != 0) | (max(outcomevalue) != 1)){
    stop("outcome must be 0/1 variable")
  }
  ###### Give error if any inputs have different lengths
  if(length(outcome) != length(covariate) | length(outcome) != nrow(predictors)){
    stop("outcome, covariate, and predictors are not the same length")
  }

  ### Re-level covariate so it's 1,...,M
  covval <- as.factor(covariate)
  levels(covval) <- rank(as.numeric(levels(covval)))
  cord <- as.numeric(covval)

  ###### Remove covariate values that are concordant
  numcases <- table(outcome, cord)[2,]
  numobs <- table(cord)
  dropcovar <- which(pmax(numcases,numobs-numcases)==numobs)
  onew <- outcome[which(cord %ni% dropcovar)]
  pnew <- predictors[which(cord %ni% dropcovar), ]
  cnew <- cord[which(cord %ni% dropcovar)]
  numcov <- length(table(cnew))

  ### Re-level covariate so it's 1,...,numcov
  covvalnew <- as.factor(cnew)
  levels(covvalnew) <- rank(as.numeric(levels(covvalnew)))
  cnew <- as.numeric(covvalnew)

  ###### Get initial estimates
  data <- data.frame(onew, cnew, pnew)
  names(data) <- c("outcome","covariate",paste("V",c(1:ncol(pnew)),sep=""))
  varnames <- names(data)
  prednames <- varnames[-c(1:2)]
  ### Logistic regression
  if(conditional==FALSE){
    glmcoef <- glm(as.formula(paste(varnames[1], " ~ ", paste(prednames, collapse=" + "),
                                    " + factor(",varnames[2],")", sep="")),
                   family="binomial", data=data)$coef[2:(length(prednames)+1)]
  }else{
    glmcoef <- survival::clogit(as.formula(paste(varnames[1], " ~ ", paste(prednames, collapse=" + "),
                                                 " + strata(",varnames[2],")", sep="")), data=data)$coef[1:length(prednames)]
  }
  ### Robust logistic regression
  robustmod <- aucm::rlogit(as.formula(paste(varnames[1], " ~ ", paste(prednames, collapse=" + "),
                                             " + factor(",varnames[2],")", sep="")), dat=data)
  if(robustmod$convergence==TRUE){
    rglmcoef <- robustmod$coef[2:(length(prednames)+1)]
  }else{
    rglmcoef <- glmcoef
  }

  if(initialval=="rGLM"){
    beta0 = rglmcoef/norm(matrix(rglmcoef),type="F")
  }else{
    beta0 = glmcoef/norm(matrix(glmcoef),type="F")
  }
  normglm = glmcoef/norm(matrix(glmcoef),type="F")
  normrglm = rglmcoef/norm(matrix(rglmcoef),type="F")
  hvalc <- rep(NA,numcov)
  for(i in 1:numcov){
    hvalc[i] = sd(as.matrix(pnew[which(cnew==i),]) %*% matrix(beta0,ncol=1))/((table(cnew)[i])^approxh)
  }

  #### Results
  rslt <- Rsolnp::solnp(beta0, saauc, eqfun=constr, eqB=1, Xmat=data[,-1], yvec=data[,1], hc=hvalc,
                                control=list("outer.iter"=10^3, "inner.iter"=10^4, "delta"=stepsz, "tol"=tolval,"trace"=0))

  #### Performance in training data
  AUCcTRrglm <- rep(NA, numcov)
  AUCcTRglm <- rep(NA, numcov)
  AUCcTRsuppl <- rep(NA, numcov)
  wtvals <-  prop.table(table(onew, cnew),1)[2,]

  for(i in 1:numcov){
    AUCcTRrglm[i] = Hmisc::somers2(pnew[which(cnew==i),] %*% normrglm, onew[which(cnew==i)])[1]
    AUCcTRglm[i] = Hmisc::somers2(pnew[which(cnew==i),] %*% normglm, onew[which(cnew==i)])[1]
    AUCcTRsuppl[i] = Hmisc::somers2(pnew[which(cnew==i),] %*% rslt$pars, onew[which(cnew==i)])[1]
  }
  aAUCTRrglm <- AUCcTRrglm %*% wtvals
  varTRrglm <-(wtvals %*% (AUCcTRrglm - AUCcTRrglm %*% wtvals)^2)

  aAUCTRsuppl <- AUCcTRsuppl %*% wtvals
  varTRsuppl <- (wtvals %*% (AUCcTRsuppl - AUCcTRsuppl %*% wtvals)^2)

  aAUCTRglm <- AUCcTRglm %*% wtvals
  varTRglm <- (wtvals %*% (AUCcTRglm - AUCcTRglm %*% wtvals)^2)

  if(rslt$convergence > 0){
    warning("SaAUC algorithm failed to converge")
  }
  if(!(robustmod$convergence)){
    warning("rGLM algorithm failed to converge; GLM used instead")
  }


  return(list(NumCov=numcov, FittedCombs=list(InitialVal=beta0, NormGLM=normglm, NormrGLM=normrglm,
                                                MaxSaAUC=rslt$pars),
                aAUCTR=c("aAUCTRrGLM"=aAUCTRrglm, "aAUCTRGLM"=aAUCTRglm, "aAUCTRsuppl"=aAUCTRsuppl),
                varTR=c("varTRrGLM"=varTRrglm,"varTRGLM"=varTRglm,"varTRsuppl"=varTRsuppl)))
}
