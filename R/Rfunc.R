#' @import Rcpp
#' @import glmnet
#' @import stats
#' @import mvtnorm
#' @import doParallel
#' @import MASS
#' @import ggplot2
#' @import foreach
#' @importFrom caret createMultiFolds
#' @importFrom matrixcalc is.singular.matrix
#' @useDynLib SA24204174
#' @title Calculate GIC
#' @description Calculate GIC of estimator
#' @param beta the estimator
#' @param n the number of samples
#' @param an the function using in gic
#' @param u the U-statistics
#' @return gic
#' @examples
#' \dontrun{
#' gic(c(1,2),5,log,c(1,2))
#'
#' }
#' @export
gic <- function(beta,n,an,u){
  return(u%*%beta-an(n)/n*sum(beta!=0))
}

#' @title weight_cos
#' @description Computes a weighted coefficient estimate.
#' @param prior A vector or matrix representing prior information.
#' @param sigmab A vector or matrix of scale factors for the prior.
#' @param u A matrix used in the computation, often representing features or data transformations.
#' @param betas A vector or matrix of coefficients to combine with prior information.
#' @param lambda A scalar weighting parameter to balance the prior and the data-driven component.
#' @return A vector or matrix of weighted coefficient estimates.
#' @examples
#' \dontrun{
#' # Example usage
#' prior <- c(1, 2, 3)
#' sigmab <- c(0.5, 0.5, 0.5)
#' u <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow = 3)
#' betas <- c(0.2, 0.3, 0.4)
#' lambda <- 0.1
#' result <- weight_cos(prior, sigmab, u, betas, lambda)
#' print(result)
#' }
#' @export
weight_cos <- function(prior,sigmab,u,betas,lambda){
  un_temp <- u%*%betas+lambda*prior%*%sigmab
  return(betahat(sigmab,un_temp))
}

#' @title betascale
#' @description Scales a coefficient vector `beta` using the quadratic form induced by `sigma`.
#' @param beta A numeric vector of coefficients.
#' @param sigma A covariance matrix used for scaling.
#' @return A numeric vector representing the scaled version of `beta`.
#' @examples
#' # Example usage
#' beta <- c(1, 2, 3)
#' sigma <- matrix(c(2, 1, 1, 1, 2, 1, 1, 1, 2), nrow = 3)
#' scaled_beta <- betascale(beta, sigma)
#' print(scaled_beta)
#'
#' # Verify normalization
#' scaled_beta %*% sigma %*% scaled_beta  # Should equal 1
#' @export
betascale<- function(beta,sigma){
  return(as.numeric(beta/as.numeric(sqrt(beta%*%sigma%*%beta))))
}

#' @title lambda_cross
#' @description
#' This function performs cross-validation to select the optimal regularization parameter (`lambda`)
#' for a given dataset and prior. It uses a custom scoring method to evaluate `lambda` values across
#' multiple folds and repetitions.
#' @param prior A numeric vector representing the prior information for weighting.
#' @param X A numeric matrix representing the feature data.
#' @param y A numeric vector representing the response variable.
#' @param betas A matrix where each column is a candidate coefficient vector.
#' @param lambda A numeric vector of candidate `lambda` values to evaluate.
#' @param ntimes An integer indicating the number of cross-validation repetitions.
#' @param k An integer specifying the number of folds in each cross-validation repetition.
#' @return The optimal `lambda` value that maximizes the cross-validation score.
#' @examples
#' # Generate synthetic data
#' set.seed(123)
#' prior <- c(1, 0.5, 0.2)
#' X <- matrix(rnorm(300), nrow = 100, ncol = 3)
#' betas <- t(matrix(rnorm(9), nrow = 3, ncol = 3))
#' y <- X%*%betas[,1]+rnorm(100)
#' lambda <- seq(0.01, 1, by = 0.1)
#' ntimes <- 5
#' k <- 5
#'
#' # Run the function
#' best_lambda <- lambda_cross(prior, X, y, betas, lambda, ntimes, k)
#' print(best_lambda)
#'
#' @export
lambda_cross <- function(prior,X,y,betas,lambda,ntimes,k){
  beta_num <- dim(betas)[2]
  dimx <- dim(X)[2]
  prior <- as.numeric(prior)
  if(beta_num != length(prior)){
    stop("Illegal prior length!")
  }
  if(ntimes < 1){
    stop("Illegal ntimes!")
  }
  if(is.null(lambda)){
    stop('NULL lambda!')
  }
  lambda <- as.numeric(lambda)
  ntimes <- floor(ntimes)
  folds <-createMultiFolds(y=X[,1],k=k,times=ntimes)
  lambda_score <- numeric(length(lambda))
  temp_score <- 0
  for(i in 1:ntimes){
    for(j in 1:length(lambda)){
      timess <- 0
      temp_score <- 0
      for(z in 1:k){
        X_temp <- X[folds[[z+k*(i-1)]],]
        y_temp <- y[folds[[z+k*(i-1)]]]
        u_temp <- un(X_temp,y_temp)
        sigmab_temp <- as.matrix(t(betas)%*%cov(X_temp)%*%betas)
        if(is.singular.matrix(sigmab_temp)){
          cat('2noni','\n')
          next
        }else{
          timess <- timess + 1
        }
        weight_temp <- weight_cos(prior,sigmab_temp,u_temp,betas,lambda[j])
        kuso <- as.numeric(betas%*%weight_temp)
        temp_score <- temp_score + as.numeric(un(X[-folds[[z+k*(i-1)]],],y[-folds[[z+k*(i-1)]]])%*%kuso)
      }
      temp_score <- temp_score/timess
      lambda_score[j] <- lambda_score[j] + temp_score
    }
    #cat('times:',timess,'\n')
  }
  cat('lambda_score:',order(-abs(lambda_score)),'\n')
  index <- order(-lambda_score)[1]
  return(lambda[index])
}


#' @title GIC to select
#' @description Use GIC (Generalized Information Criterion) to select the best regularization parameter
#' @param lambda a vector of candidate regularization parameters
#' @param u the U-statistics vector
#' @param covx the covariance matrix of the predictors
#' @param n the number of samples
#' @param p the number of predictors
#' @param an the penalty adjustment function used in GIC
#' @param start the starting point for the optimization
#' @param ori the initial estimator using in debias penalty term if theres no need for debias set it 0
#' @param weight the weights applied in the lasso optimization
#' @param s sign of the first term or u can use sign of U-statistics
#' @param L Lipschitz constant used in the optimization algorithm
#' @param alpha the learning rate for the optimization
#' @param iter the maximum number of iterations for the optimization
#' @param tol the tolerance level for convergence
#' @return the best lambda value that minimizes the GIC
#' @examples
#' \dontrun{
#' gic_select(lambda = c(0.1, 0.2, 0.3), u = c(1, 2),
#'            covx = matrix(c(1, 0.5, 0.5, 1), 2, 2),
#'            n = 100, p = 2, an = log, start = c(1, 1),
#'            ori = c(1, 1), weight = c(1, 1))
#' }
#' @export
gic_select <- function(lambda,u,covx,n,p,an,start,ori,weight,s=1,L=1.2,alpha=0.11,iter=2e2,tol=1e-6){
  lambda_score <- numeric(length(lambda))
  for(i in 1:length(lambda)){
    bbet <- lmrc_lasso(n,p,covx,L,s,u,start,ori,weight,lambda[i],alpha,iter,tol)
    bbet <- as.numeric(bbet)
    bbet <- bbet/as.numeric(sqrt(t(bbet)%*%covx%*%bbet))
    lambda_score[i] =  gic(bbet,n,an,u)
  }
  return(lambda[order(-lambda_score)[1]])
}

#' @title GIC select with debias term
#' @description Use GIC (Generalized Information Criterion) to select the best parameters with a debias term
#' @param parameters a matrix where each row contains a combination of gamma, lambda1, and lambda2 to evaluate with at most 3 gammas.
#' @param n the number of samples
#' @param an the penalty adjustment function used in GIC
#' @param u the U-statistics vector
#' @param covx the covariance matrix of the predictors
#' @param p the number of predictors
#' @param start the starting point for the optimization
#' @param ori the original data or parameter used in the optimization process
#' @param weight the weights applied in the lasso optimization
#' @param s scaling factor for the solution
#' @param L Lipschitz constant used in the optimization algorithm
#' @param alpha the learning rate for the optimization
#' @param iter the maximum number of iterations for the optimization
#' @param tol the tolerance level for convergence
#' @return the debiased coefficients selected by the best combination of gamma, lambda1, and lambda2
#' @examples
#' \dontrun{
#' parameters <- matrix(c(0.5, 0.1, 0.01, 1, 0.2, 0.02), ncol = 3, byrow = TRUE)
#' debias_gic_select(parameters, n = 100, an = log,
#'                   u = c(1, 2), covx = matrix(c(1, 0.5, 0.5, 1), 2, 2),
#'                   p = 2, start = c(1, 1), ori = c(1, 1),
#'                   weight = c(1, 1), s = 1)
#' }
#' @export
debias_gic_select <- function(parameters,n,an,u,covx,p,start,ori,weight,s,L=1.2,alpha=0.11,iter=2e2,tol=1e-6){
  gic_max <- -Inf
  max_gamma <- 0
  max_lambda1 <- 0
  max_lambda2 <- 0
  if(dim(parameters)[2]!=3){
    stop('Illegal parameters dim!')
  }
  for(i in 1:dim(parameters)[1]){
    temp_ori <- abs(ori)^as.numeric(parameters[i,1])
    kuso <- lmrc_lasso_debias(n,p,covx,L,s,u,start,ori ,weight,temp_ori ,as.numeric(parameters[i,2]),alpha,iter,tol,as.numeric(parameters[i,3]))
    if(gic(kuso/L,n,an,u) > gic_max){
      gic_max <- gic(kuso/L,n,an,u)
      max_gamma <- as.numeric(parameters[i,1])
      max_lambda1 <- as.numeric(parameters[i,2])
      max_lambda2 <- as.numeric(parameters[i,3])
    }
  }
  temp_ori <- abs(ori)^max_gamma
  kuso <- lmrc_lasso_debias(length(u),length(u),covx,L,s,u,start,ori ,weight,temp_ori ,max_lambda1,alpha,iter,tol,max_lambda2)
  cat('last_gamma:',max_gamma,'last_lambda1:',max_lambda1,'last_lambda2:',max_lambda2,'gic_max:',gic_max,'\n')
  return(kuso)
}


#' @title FP
#' @description Calculate the false positive rate of an estimator
#' @param beta the estimated coefficients
#' @param beta0 the true coefficients
#' @return the false positive rate, defined as the proportion of nonzero entries in `beta` corresponding to zero entries in `beta0`
#' @examples
#' \dontrun{
#' beta <- c(0, 1, 0, 1)
#' beta0 <- c(0, 0, 0, 1)
#' fp(beta, beta0)  # Returns 0.5
#' }
#' @export
fp <- function(beta,beta0){
  support <- (beta0==0)
  return(sum(beta[support]!=0)/sum(support))
}

#' @title FN
#' @description Calculate the false negative rate of an estimator
#' @param beta the estimated coefficients
#' @param beta0 the true coefficients
#' @return the false negative rate, defined as the proportion of zero entries in `beta` corresponding to nonzero entries in `beta0`
#' @examples
#' \dontrun{
#' beta <- c(0, 1, 0, 1)
#' beta0 <- c(1, 1, 0, 1)
#' fn(beta, beta0)  # Returns 0.3333333
#' }
#' @export
fn <- function(beta,beta0){
  noise <- (beta0!=0)
  return(sum(beta[noise]==0)/sum(noise))
}

#' @title CM
#' @description Calculate the correctness metric for an estimator
#' @param beta the estimated coefficients
#' @param beta0 the true coefficients
#' @return a binary value: 1 if the support of `beta` matches exactly with the support of `beta0`, and 0 otherwise
#' @examples
#' \dontrun{
#' beta <- c(0, 1, 0, 1)
#' beta0 <- c(0, 1, 0, 1)
#' cm(beta, beta0)  # Returns 1
#'
#' beta <- c(0, 1, 0, 1)
#' beta0 <- c(0, 1, 1, 1)
#' cm(beta, beta0)  # Returns 0
#' }
#' @export
cm <- function(beta,beta0){
  return(as.numeric(sum((beta!=0)==(beta0!=0))==length(beta0)))
}


## This is a function written by Lisai in Trans-Lasso.
#' @title Q Aggregation
#' @description Aggregation function to estimate coefficients using either model selection or Q-aggregation
#' @param B a matrix of estimated coefficients (rows represent predictors, columns represent models)
#' @param X.test the test predictor matrix
#' @param y.test the test response vector
#' @param total.step the maximum number of iterations for the Q-aggregation process (default is 10)
#' @param selection a boolean flag; if `TRUE`, selects the coefficients with the smallest prediction error.
#' If `FALSE`, performs Q-aggregation (default is `FALSE`)
#' @return a list containing:
#' \item{theta}{the final weight vector for the models}
#' \item{beta}{the aggregated coefficients}
#' \item{beta.ew}{the coefficients obtained at the end of the Q-aggregation process (if `selection = FALSE`)}
#' @examples
#' \dontrun{
#' B <- matrix(rnorm(20), nrow = 4)
#' X.test <- matrix(rnorm(40), nrow = 10)
#' y.test <- rnorm(10)
#' result <- agg.fun(B, X.test, y.test, total.step = 10, selection = FALSE)
#' print(result)
#' }
#' @export
agg.fun<- function(B, X.test,y.test, total.step=10, selection=F){
  if(sum(B==0)==ncol(B)*nrow(B)){
    return(rep(0,nrow(B)))
  }
  p<-nrow(B)
  K<-ncol(B)
  colnames(B)<-NULL
  if(selection){#select beta.hat with smallest prediction error
    khat<-which.min(colSums((y.test-X.test%*%B)^2))
    theta.hat<-rep(0, ncol(B))
    theta.hat[khat] <- 1
    beta=B[,khat]
    beta.ew=NULL
  }else{#Q-aggregation
    theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2)
    theta.hat=theta.hat/sum(theta.hat)
    theta.old=theta.hat
    beta<-as.numeric(B%*%theta.hat)
    beta.ew<-beta
    # theta.old=theta.hat
    for(ss in 1:total.step){
      theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2+colSums((as.vector(X.test%*%beta)-X.test%*%B)^2)/8)
      theta.hat<-theta.hat/sum(theta.hat)
      beta<- as.numeric(B%*%theta.hat*1/4+3/4*beta)
      if(sum(abs(theta.hat-theta.old))<10^(-3)){break}
      theta.old=theta.hat
    }
  }
  list(theta=theta.hat, beta=beta, beta.ew=beta.ew)
}

## This is a function written by Lisai in Trans-Lasso.
###oracle Trans-Lasso
#' @title Oracle Trans-LASSO
#' @description Oracle Trans-LASSO
#' @param X the design matrix
#' @param y the response vector
#' @param A0 a set of indices for the initial active set
#' @param n.vec a vector defining sample splitting indices; `n.vec[1]` defines the size of the first part
#' @param lam.const an optional constant for the penalty term; if `NULL`, it is determined through cross-validation
#' @param l1 a logical flag indicating whether to use lasso regression in the first stage (default is `TRUE`)
#' @return a list containing:
#' \item{beta.kA}{the estimated coefficients after the two-stage procedure}
#' \item{w.kA}{the coefficients from the first stage of lasso regression (if applicable)}
#' \item{lam.const}{the constant for the penalty term}
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(100 * 10), 100, 10)
#' y <- rnorm(100)
#' A0 <- c(1, 2, 3)
#' n.vec <- c(80, 20)
#' result <- las.kA(X, y, A0, n.vec)
#' print(result)
#' }
#' @export
las.kA<-function(X, y, A0, n.vec, lam.const=NULL, l1=T){
  p<-ncol(X)
  size.A0<- length(A0)
  if(size.A0 > 0){
    ind.kA<- ind.set(n.vec, c(1, A0+1))
    ind.1<-1:n.vec[1]
    if(l1){
      y.A<-y[ind.kA]
    }
    if(is.null(lam.const)){
      cv.init<-cv.glmnet(X[ind.kA,], y.A, nfolds=8, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.kA)))
      lam.const <- cv.init$lambda.min/sqrt(2*log(p)/length(ind.kA))
    }
    w.kA <- as.numeric(glmnet(X[ind.kA,], y.A, lambda=lam.const*sqrt(2*log(p)/length(ind.kA)))$beta)
    w.kA<-w.kA*(abs(w.kA)>=lam.const*sqrt(2*log(p)/length(ind.kA)))
    # cv.delta<-cv.glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.1)))
    #delta.kA<-predict(cv.delta, s='lambda.min', type='coefficients')[-1]
    delta.kA <- as.numeric(glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=lam.const*sqrt(2*log(p)/length(ind.1)))$beta)
    delta.kA<-delta.kA*(abs(delta.kA)>=lam.const*sqrt(2*log(p)/length(ind.1)))
    beta.kA <- w.kA + delta.kA
    lam.const=NA
  }else{
    cv.init<-cv.glmnet(X[1:n.vec[1],], y[1:n.vec[1]], nfolds=8, lambda=seq(1,0.1,length.out=20)*sqrt(2*log(p)/n.vec[1]))
    lam.const<-cv.init$lambda.min/sqrt(2*log(p)/n.vec[1])
    beta.kA <- predict(cv.init, s='lambda.min', type='coefficients')[-1]
    w.kA<-NA
  }
  list(beta.kA=as.numeric(beta.kA),w.kA=w.kA, lam.const=lam.const)

}

## This is a function written by Lisai in Trans-Lasso.
#Trans Lasso method
#' @title Trans-LASSO
#' @description Perform Transfer Learning with Lasso regression using aggregation and selection
#' @param X the design matrix
#' @param y the response vector
#' @param n.vec a vector defining the sample splitting indices; `n.vec[k]` defines the size of the `k`th split
#' @param I.til a vector of indices for the transfer dataset to be used for aggregation
#' @param l1 a boolean flag indicating whether to use lasso regression in the selection steps (default is `TRUE`)
#' @return a list containing:
#' \item{beta.hat}{the aggregated coefficients using Q-aggregation}
#' \item{theta.hat}{the weights assigned to models during Q-aggregation}
#' \item{rank.pi}{the ranks of transfer datasets based on the computed Rhat values}
#' \item{beta.pool}{the aggregated coefficients using an alternative pool method}
#' \item{theta.pool}{the weights assigned to models during the alternative pool aggregation}
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(1000), 100, 10)
#' y <- rnorm(100)
#' n.vec <- c(50, 30, 20)
#' I.til <- sample(1:50, 10)
#' result <- Trans.lasso(X, y, n.vec, I.til)
#' print(result)
#' }
#' @export
Trans.lasso <- function(X, y, n.vec, I.til, l1=T){
  M= length(n.vec)-1
  #step 1
  X0.til<-X[I.til,] #used for aggregation
  y0.til<-y[I.til]
  X<- X[-I.til,]
  y<-y[-I.til]
  #step 2
  Rhat <- rep(0, M+1)
  p<- ncol(X)
  n.vec[1]<- n.vec[1]-length(I.til)
  ind.1<-ind.set(n.vec,1)
  for(k in 2: (M+1)){
    ind.k<-ind.set(n.vec,k)
    Xty.k <- t(X[ind.k,])%*%y[ind.k]/n.vec[k] - t(X[ind.1,])%*%y[ind.1]/ n.vec[1]
    margin.T<-sort(abs(Xty.k),decreasing=T)[1:round(n.vec[1]/3)]
    Rhat[k] <-  sum(margin.T^2)
  }
  Tset<- list()
  k0=0
  kk.list<-unique(rank(Rhat[-1]))
  #cat(rank(Rhat[-1]),'\n')
  for(kk in 1:length(kk.list)){#use Rhat as the selection rule
    Tset[[k0+kk]]<- which(rank(Rhat[-1]) <= kk.list[kk])
  }
  k0=length(Tset)
  Tset<- unique(Tset)
  #cat(length(Tset),'\n')

  beta.T<-list()
  init.re<-las.kA(X=X, y=y, A0=NULL, n.vec=n.vec, l1=l1)
  beta.T[[1]] <- init.re$beta.kA
  beta.pool.T<-beta.T ##another method for comparison
  for(kk in 1:length(Tset)){#use pi.hat as selection rule
    T.k <- Tset[[kk]]
    re.k<- las.kA(X=X, y=y, A0=T.k, n.vec=n.vec, l1=l1, lam.const=init.re$lam.const)
    beta.T[[kk+1]] <-re.k$beta.kA
    beta.pool.T[[kk+1]]<-re.k$w.kA
  }
  beta.T<-beta.T[!duplicated((beta.T))]
  beta.T<- as.matrix(as.data.frame(beta.T))
  agg.re1 <- agg.fun(B=beta.T, X.test=X0.til, y.test=y0.til)
  beta.pool.T<-beta.pool.T[!duplicated((beta.pool.T))]
  beta.pool.T<- as.matrix(as.data.frame(beta.pool.T))
  agg.re2<-agg.fun(B=beta.pool.T, X.test=X0.til, y.test=y0.til)

  return(list(beta.hat=agg.re1$beta, theta.hat=agg.re1$theta, rank.pi=rank(Rhat[-1]),
              beta.pool=agg.re2$beta, theta.pool=agg.re2$theta))
}

## This is a function written by Lisai in Trans-Lasso.
#A method for comparison: Trans-Lasso(l1). It has the same pipeline of the Trans-Lasso
###but with sparsity index R_k=\|w^{(k)}-\beta\|_1 and a naive aggregation (empirical risk minimization)
#' @title Trans.lasso.sp
#' @description Perform Transfer Learning with a sparsity constraint using Lasso regression
#' @param X the design matrix
#' @param y the response vector
#' @param n.vec a vector defining the sample splitting indices; `n.vec[k]` defines the size of the `k`th split
#' @param I.til a vector of indices for the transfer dataset to be used for aggregation
#' @param l1 a boolean flag indicating whether to use lasso regression in the selection steps (default is `TRUE`)
#' @return a list containing:
#' \item{beta.sp}{the selected coefficients using model selection aggregation}
#' \item{theta.sp}{the weights assigned to models during aggregation}
#' \item{rank.pi}{the ranks of transfer datasets based on the computed Rhat values}
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(1000), 100, 10)
#' y <- rnorm(100)
#' n.vec <- c(50, 30, 20)
#' I.til <- sample(1:50, 10)
#' result <- Trans.lasso.sp(X, y, n.vec, I.til)
#' print(result)
#' }
#' @export
Trans.lasso.sp <- function(X, y, n.vec, I.til, l1=T){
  M= length(n.vec)-1
  #step 1
  X0.til<-X[I.til,] #used for aggregation
  y0.til<-y[I.til]
  X<- X[-I.til,]
  y<-y[-I.til]
  #step 2
  Rhat <- rep(0, M+1)
  p<- ncol(X)
  n.vec[1]<- n.vec[1]-length(I.til)
  ind.1<-ind.set(n.vec,1)
  init.re<-las.kA(X=X, y=y, A0=NULL, n.vec=n.vec, l1=l1)
  for(k in 2: (M+1)){
    ind.k <- ind.set(n.vec,k)
    w.init.k<- as.numeric(glmnet(X[ind.k,], y[ind.k], lambda=init.re$lam.const*sqrt(2*log(p)/length(ind.k)))$beta)
    Rhat[k] <-  sum(abs(w.init.k-init.re$beta.kA)) ##\|w^{(k)}-\beta\|_1
  }
  Tset<- list()
  k0=0
  kk.list<-unique(rank(Rhat[-1]))
  #cat(rank(Rhat[-1]),'\n')
  for(kk in 1:length(kk.list)){#use pi.hat as selection rule
    Tset[[k0+kk]]<- which(rank(Rhat[-1]) <= kk.list[kk])
  }
  k0=length(Tset)
  Tset<- unique(Tset)
  # cat(length(Tset),'\n')

  beta.T<-list()
  beta.T[[1]] <- init.re$beta.kA
  for(kk in 1:length(Tset)){#use pi.hat as selection rule
    T.k <- Tset[[kk]]
    beta.T[[kk+1]] <-las.kA(X=X, y=y, A0=T.k, lam.const=init.re$lam.const,n.vec=n.vec, l1=l1)$beta.kA
  }
  beta.T<-beta.T[!duplicated((beta.T))]
  beta.T<- as.matrix(as.data.frame(beta.T))
  agg.re <- agg.fun(B=beta.T, X.test=X0.til, y.test=y0.til, selection =T)
  return(list(beta.sp=agg.re$beta, theta.sp=agg.re$theta, rank.pi=rank(Rhat[-1])))
}

#' @title mse
#' @description Compute estimation error and prediction error (optional)
#' @param beta the true coefficients
#' @param est the estimated coefficients
#' @param X.test (optional) the test design matrix
#' @return a list containing:
#' \item{est.err}{the estimation error (squared error between `beta` and `est`)}
#' \item{pred.err}{the prediction error (mean squared error of predictions), or `NA` if `X.test` is not provided}
#' @examples
#' \dontrun{
#' beta <- c(1, 0, 0.5)
#' est <- c(1.1, 0, 0.4)
#' X.test <- matrix(c(1, 0.5, -0.2, 1.5, 0.3, 0.7), 3, 2)
#' mse.fun(beta, est, X.test)
#' }
#' @export
mse.fun<- function(beta,est, X.test=NULL){
  pred.err<-NA
  est.err<- sum((beta-est)^2)

  if(!is.null(X.test)){
    pred.err<-  mean((X.test%*%(beta-est))^2)
  }
  return(list(est.err=est.err, pred.err= pred.err))
}

## This is a function written by Lisai in Trans-Lasso.
#' @title index of set
#' @description Generate indices for specified groups based on sample splits
#' @param n.vec a vector specifying the size of each group; `n.vec[k]` represents the size of the `k`-th group
#' @param k.vec a vector of group indices for which the function generates indices
#' @return A vector of indices corresponding to the specified groups
#' @examples
#' \dontrun{
#' # Example with 3 groups of sizes 50, 30, and 20
#' n.vec <- c(50, 30, 20)
#' # Generate indices for groups 1 and 2
#' ind.set(n.vec, c(1, 2))
#' # Output: indices 1 to 80
#' }
#' @export
ind.set<- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}

#' @title rep columns
#' @description Replicate a vector as columns in a matrix
#' @param x a numeric vector to be replicated
#' @param n the number of columns to replicate
#' @return A matrix with `n` columns, each containing the replicated values of `x`
#' @examples
#' \dontrun{
#' # Replicate the vector c(1, 2, 3) into 4 columns
#' repcol(c(1, 2, 3), 4)
#' # Output:
#' #      [,1] [,2] [,3] [,4]
#' # [1,]    1    1    1    1
#' # [2,]    2    2    2    2
#' # [3,]    3    3    3    3
#' }
#' @export
repcol<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}


## This is a function written by Lisai in Trans-Lasso.
#' @title generate coeffients
#' @description Generate coefficient matrices and true coefficients for simulations
#' @param s the number of nonzero elements in the true coefficient vector (`beta0`)
#' @param h the number of modified elements for selected prior groups
#' @param q the number of modified elements for non-selected prior groups (default is 30)
#' @param size.A0 the number of prior groups with selected modifications
#' @param M the number of prior estimates to generate
#' @param sig.beta the signal magnitude of nonzero elements in the true coefficients
#' @param sig.delta1 the magnitude of modifications for selected prior groups
#' @param sig.delta2 the magnitude of modifications for non-selected prior groups
#' @param p the total number of predictors
#' @param exact a logical value; if `TRUE`, modifications are deterministic, otherwise they are random
#' @return A list containing:
#' \item{W}{A matrix of prior estimates (size `p x M`), with modifications applied}
#' \item{beta0}{The true coefficient vector of size `p`}
#' @examples
#' \dontrun{
#' # Generate coefficients for a simulation
#' result <- Coef.gen(s = 10, h = 5, q = 15, size.A0 = 3, M = 10,
#'                    sig.beta = 1, sig.delta1 = 0.5, sig.delta2 = 0.2, p = 100)
#' W <- result$W   # Prior estimates matrix
#' beta0 <- result$beta0 # True coefficient vector
#' }
#' @export
Coef.gen<- function(s, h,q=30, size.A0, M, sig.beta,sig.delta1, sig.delta2, p, exact=T){
  beta0 <- c(rep(sig.beta,s), rep(0, p - s))
  W <- repcol(beta0,  M)# ten prior estimates
  W[1,]<-W[1,]-2*sig.beta
  for(k in 1:M){
    if(k <= size.A0){
      if(exact){
        samp0<- sample(1:p, h, replace=F)
        W[samp0,k] <-W[samp0,k] + rep(-sig.delta1, h)
      }else{
        W[1:100,k] <-W[1:100,k] + rnorm(100, 0, h/100)
      }
    }else{
      if(exact){
        samp1 <- sample(1:p, q, replace = F)
        W[samp1,k] <- W[samp1,k] + rep(-sig.delta2,q)
      }else{
        W[1:100,k] <-W[1:100,k] + rnorm(100, 0, q/100)
      }
    }
  }

  return(list(W=W, beta0=beta0))
}
