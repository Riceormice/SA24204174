## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(SA24204174)
library(MASS)
library(Rcpp)
library(glmnet)
library(caret)#做交叉验证用
library(mvtnorm)
library(parallel)
library("foreach")
library("doParallel")
library(matrixcalc)

## ----dgp----------------------------------------------------------------------
nums <- 0
sims <- 0
set.seed(20231027+sims)
priors <- numeric(72)
for(i in 1:36){
  priors[(2*i-1):(2*i)] <- sort(runif(2,1,2))
}
# 0.5 去掉两个，偏大了！！！！！
sse_fit <- NULL
sse_target_only <- NULL
fp_lasso <- NULL
fn_lasso <- NULL
mse1 <- NULL
fp_oracle <- NULL
fn_oracle <- NULL

mse2 <- NULL
fp_trans <- NULL
fn_trans <- NULL

mse3 <- NULL
mse4 <- NULL
mse5 <- NULL
mse6 <- NULL
#mse7 <- NULL
#mse8 <- NULL
mse_single_trans1 <- NULL
mse_single_trans2 <- NULL
fp_target <- NULL
fn_target <- NULL
fp_trans1 <- NULL
fn_trans1 <- NULL
fp_trans2 <- NULL
fn_trans2 <- NULL

p_ini <- NULL
p_1 <- NULL
p_2 <- NULL
p_3 <- NULL
p_4 <- NULL
p_5 <- NULL
p_6 <- NULL
p_last <- NULL

fp_ini <- NULL
fp_1 <- NULL
fp_2 <- NULL
fp_3 <- NULL
fp_4 <- NULL
fp_5 <- NULL
fp_6 <- NULL
fp_last <- NULL

fn_ini <- NULL
fn_1 <- NULL
fn_2 <- NULL
fn_3 <- NULL
fn_4 <- NULL
fn_5 <- NULL
fn_6 <- NULL
fn_last <- NULL

iter_times <- NULL

mymse_l2_ini <- NULL
mymse_l2_1 <- NULL
mymse_l2_2 <- NULL
mymse_l2_3 <- NULL
mymse_l2_4 <- NULL
mymse_l2_5 <- NULL
mymse_l2_6 <- NULL
mymse_l2_last <- NULL
  trans_sse <- NULL
  ps <- NULL
  n0s <- NULL
  sizea0s <- NULL
  sigbetas <- NULL
  Ms <- NULL
  hs <- NULL
  exacts <- NULL
  ss <- NULL
  
  cm_oracle <- NULL
  cm_lasso <- NULL
  cm_trans <- NULL
  cm_my <- NULL
  
  # source data size
  size.A0 = 0
  sig.beta = 0.3
  p = 500
  M = 20
  n0 = 150
  sig.z <- 1
  n_s <- 100
  n.vec <- c(n0, rep(n_s, M))
Sig.X <- diag(1, p)
Niter = 200
l1=T
A0 = 1:size.A0
ind.1 <- 1:n0
#这里h变为12了，6.6
h=6
exact = TRUE
s=16
set.seed(20231027+ sims + nums)
                  #source("TransLasso-functions.R")
                  start <- (proc.time())[3][[1]]
                  #if(size.A0<1){next}
                  nums<-nums+1
                  
                  #if(nums !=11) next
                  print(nums)
                  beta0<-
                  coef.all <-Coef.gen( s, h = h, q = 2*s, size.A0 = size.A0,  M = M,   sig.beta = sig.beta,
                                         sig.delta1 = sig.beta, sig.delta2 = sig.beta+0.2, p = p, exact)
                  B <- cbind(coef.all$beta0, coef.all$W)
                  beta0 <- coef.all$beta0
                  
                  ###generate the data
                  X <- NULL
                  y <- NULL
                  for (k in 1:(M + 1)) {
                    X <- rbind(X, rmvnorm(n.vec[k], rep(0, p), Sig.X))
                    ind.k <- ind.set(n.vec, k)
                    y <- c(y, X[ind.k, ] %*% B[, k] + rnorm (n.vec[k], 0, 1))
                  }

## ----translasso---------------------------------------------------------------
                  ###compute init beta ####
                  mse.vec<-rep(NA,8)
                  beta.init <-
                    as.numeric(glmnet(X[1:n.vec[1], ], y[1:n.vec[1]], lambda = 0.52*sqrt(2 * log(p) / n.vec[1]))$beta)
                  beta.init <- betascale(beta.init,Sig.X)*1.2
                  mse.vec[1] = mse.fun(as.numeric(beta.init),beta0,X[1:n0,])$est.err
                  fp_lasso <- c(fp_lasso,fp(beta.init,beta0))
                  fn_lasso <- c(fn_lasso,fn(beta.init,beta0))
                  #return(list(X=X,y=y))
                  ######Oracle Trans-Lasso#######
                  if (size.A0 == 0) {
                    beta.kA <- beta.init
                    mse.vec[2] = mse.fun(as.numeric(beta.kA),beta0,X[1:n0,])$est.err
                    #mse.vec[7] =  mse.vec[2]
                    #r1 <- my_difftuning(X,y,n0,beta0,beta.init,0.04)$beta
                    #mse.vec[7] =  mse.fun(as.numeric(r1),y[1:n0], beta0)$est.err
                  } else{
                    sr <- las.kA(X, y, A0 = 1:size.A0, n.vec = n.vec, l1=l1)
                    beta.kA <- sr$beta.kA
                    w.kA <- sr$w.kA
                    #cat('oracle_lambda:',sr$lam.const,'\n')
                    beta.kA <- betascale(beta.kA,Sig.X)*1.2
                    mse.vec[2] = mse.fun(as.numeric(beta.kA),beta0,X[1:n0,])$est.err
                    #mse.vec[7] =  mse.fun(as.numeric(w.kA),y[1:n0], beta0)$est.err
                    #r1 <- my_difftuning(X,y,n0,beta0,sr$beta.kA,sr$lam.const)$beta
                    #cat('oracle with my tuning:',r1[1:16],'\n')
                    #mse.vec[7] =  mse.fun(r1,y[1:n0], beta0)$est.err
                  }
                  fp_oracle <- c(fp_oracle,fp(beta.kA,beta0))
                  fn_oracle <- c(fn_oracle,fn(beta.kA,beta0))
                  ###########Trans-Lasso#############
                  prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:50, l1 = l1)
                  prop.re2 <- Trans.lasso(X, y, n.vec, I.til = 51:100, l1=l1)
                  prop.re3 <- Trans.lasso(X, y, n.vec, I.til = 101:150, l1=l1)
                  if(size.A0 > 0 & size.A0< M){ #Rank.re characterizes the performance of the sparsity index Rk
                    Rank.re<- (sum(prop.re1$rank.pi[1:size.A0]<=size.A0) +
                                 sum(prop.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
                  }else{ Rank.re <- 1 }
                  beta.prop <- (prop.re1$beta.hat  + prop.re3$beta.hat) / 2
                  beta.prop <- beta.prop*(abs(beta.prop)>0.002)
                  prop.re1$beta.hat <- prop.re1$beta.hat*(abs(prop.re1$beta.hat)>0.002)
                  prop.re3$beta.hat <- prop.re1$beta.hat*(abs(prop.re3$beta.hat)>0.002)
                  mse_single_trans1 <- c(  mse_single_trans1,mse.fun(prop.re1$beta.hat ,beta0,X[1:n0,])$est.err)
                  mse_single_trans2 <- c(  mse_single_trans2,mse.fun(prop.re3$beta.hat ,beta0,X[1:n0,])$est.err)
                  cat('beta.kA:',beta.kA,'\n')
                  cat('beta.prop:',beta.prop,'\n')
                  beta.prop <- betascale(beta.prop,Sig.X)*1.2
                  cat('trans:',mse.fun(beta.prop,beta0,X[1:n0,])$est.err,'\n')
                  trans_sse <- c(trans_sse,mse.fun(beta.prop,beta0,X[1:n0,])$est.err)
                  fp_trans <- c(fp_trans,fp(beta.prop,beta0))
                  fn_trans <- c(fn_trans,fn(beta.prop,beta0))
                  fp_trans1 <- c(fp_trans1,fp(prop.re1$beta.hat,beta0))
                  fn_trans1 <- c(fn_trans1,fn(prop.re1$beta.hat,beta0))
                  fp_trans2 <- c(fp_trans2,fp(prop.re3$beta.hat,beta0))
                  fn_trans2 <- c(fn_trans2,fn(prop.re3$beta.hat,beta0))
                  mse.vec[3] = mse.fun(beta.prop,beta0,X[1:n0,])$est.err
                  
                  ######A method for comparison: it has the same pipeline of the Trans-Lasso
                  ###but with sparsity index R_k=\|w^{(k)}-\beta\|_1 and a naive aggregation (empirical risk minimization)
                  prop.sp.re1 <- Trans.lasso.sp(X, y, n.vec, I.til = 1:50, l1 = l1)
                  prop.sp.re2 <- Trans.lasso.sp(X, y, n.vec, I.til = 101:n.vec[1], l1=l1)
                  if(size.A0 > 0 & size.A0< M){
                    Rank.re.sp <- (sum(prop.sp.re1$rank.pi[1:size.A0]<=size.A0) +
                                     sum(prop.sp.re2$rank.pi[1:size.A0]<=size.A0))/2/size.A0
                  }else{ Rank.re.sp <-1 }
                  beta.sp <- (prop.sp.re1$beta.sp + prop.sp.re2$beta.sp) / 2
                  mse.vec[4] = mse.fun(beta.sp,beta0,X[1:n0,])$est.err
                  
                  ######another method for comparison: it is the same as Trans-Lasso except
                  ##that the bias correction step (step 2 of Oracle Trans-Lasso) is omitted
                  beta.pool<-(prop.re1$beta.pool+prop.re2$beta.pool)/2
                  mse.vec[5] = mse.fun(beta.pool,beta0,X[1:n0,])$est.err
                  #####naive translasso: simply assumes A0=1:K
                  beta.all <- las.kA(X, y, A0 = 1:M, n.vec = n.vec, l1=l1)$beta.kA#naive transLasso
                  mse.vec[6] = mse.fun(as.numeric(beta.all),beta0,X[1:n0,])$est.err
                  cat('Trans-LASSO:',mse.vec,'\n')

## ----SA24204174---------------------------------------------------------------
                  L <- as.numeric(sqrt(t(beta0)%*%Sig.X%*%beta0))
                  X1<-X
                  y1<-y
                  #if(nums < 6) next
                  u0 <- un(X[1:n.vec[1],],y[1:n.vec[1]])
                  Us <- matrix(rep(0,p*(M+1)),p,M+1)
                  for(i in 1:M){
                    inde <- (n.vec[1]+1+(i-1)*(n.vec[i+1])):(n.vec[1]+i*n.vec[i+1])
                    Us[, i] <- un(X[inde,],y[inde])
                  }
                  Us[,M+1] <- u0
                  u_source1=rep(1.2,p)
                  u_source11 <- u_source1
                  u_source11[1]=-1
                  #next
                  #if(size.A0 ==4) return(list(X=X,y=y))
                  dimx <- p
                  record<-rep(1,p)
                  beta_source_sum <- rep(0,p)
                  beta_source_sum1 <- rep(0,p)
                  beta_source_sum2 <- rep(0,p)
                  beta_sum_all <- matrix(rep(0,p*(M+1)),p,M+1)
                  dimx <- length(beta_source_sum)
                  lambda_temp <- 0.001
                  
                  
                  p_ini <- c(p_ini,p)
                  wudi_num <- 0.0
                  
                  cat('第一次循环：','\n')
                  #cat('beta:',length(beta0),'\n')
                  system.time({
                    cl<- makeCluster(1)
                    registerDoParallel(cl)       #进行进程注册
                    beta_sum_all <- foreach(
                      i=1:M,          #输入等待请求的参数
                      .combine=cbind,  #返回结果的整合
                      .packages = c("SA24204174")
                      #多个进程共享的系统环境
                    ) %dopar% {
                      an <- function(n){
                        log(log(n))/sqrt(n)
                      }
                      inde <- (n0 + 1 + (i - 1) * n_s):(n0 + i * n_s)
                      lambda_temp <- SA24204174::gic_select(seq(0.065, 0.076, 0.001), Us[, i], cov(X[inde, ]), n_s,dimx,an, u_source1, rep(0, dimx), 1/sqrt(abs(u_source1)), s=0,tol = 1e-6)
                      
                      kuso <- SA24204174::lmrc_lasso(dimx, dimx, cov(X[inde, ]), L, 0, Us[, i], u_source1, rep(0, dimx), 1/sqrt(abs(u_source11)), lambda_temp, 0.11, 2e2, 1e-6)
                      #kuso = kuso * (abs(kuso) > wudi_num)
                      kuso <- kuso / as.numeric(sqrt(t(kuso) %*% cov(X[1:n0, ]) %*% kuso))
                      return(kuso)
                    }
                  })
                  
                  #cat(beta_sum_all,'\n')
                  an <- function(n){
                    log(log(n))/sqrt(n)
                  }
                  inde <- 1:n0
                  lambda_temp <- SA24204174::gic_select(seq(0.055, 0.088, 0.001), Us[, i], cov(X[inde, ]), n_s,dimx,an, u_source1, rep(0, dimx), 1/sqrt(abs(u_source1)), s=1)
                  kuso1 <- SA24204174::lmrc_lasso(dimx, dimx, cov(X[inde, ]), L, 1, u0, u_source1, rep(0, dimx), 1/sqrt(abs(u_source1)), lambda_temp, 0.11, 2e2, 1e-6)
                  beta_fit <- kuso1*(as.numeric(kuso1%*%t(X[1:n0,])%*%y[1:n0])/as.numeric(kuso1%*%t(X[1:n0,])%*%X[1:n0,]%*%kuso1))
                  cat('target_only:',sum((beta_fit-beta0)^2),'\n')
                  sse_target_only <- c(sse_target_only,sum((beta_fit-beta0)^2))
                  beta_target <- beta_fit
                  iter_times_num <- 1
                  fp_target <- c(fp_target,fp(beta_fit,beta0))
                  fn_target <- c(fn_target,fn(beta_fit,beta0))
                  
                  
                  
                  
                  
                  
                  
                  ## ini weighting
                  
                  prior <- rep(1,M)
                  sigmab <- as.matrix(t(beta_sum_all)%*%cov(X[1:150,])%*%beta_sum_all)
                  prior <- betascale(prior,sigmab)
                  #priors <- rep(1,36)
                  #prior<-priors[(2*nums-1):(2*nums)]
                  #cat('prior:',prior,'\n')
                  #prior <- c(rep(prior[2],size.A0),rep(prior[1],M-size.A0))
                  beta_cal<-numeric(p)
                  if(is.singular.matrix(sigmab)){
                    cat('weixian,no inverse!!!!','\n')
                  }
                  
                  if(size.A0 <=4){
                    gred_lambda <- lambda_cross(prior,X[1:150,],y[1:150],beta_sum_all,seq(1.3,2,0.1),5,3)
                  }else{
                    gred_lambda <- lambda_cross(prior,X[1:150,],y[1:150],beta_sum_all,seq(1.1,2.0,0.1),5,3)
                  }
                  weight_temp <- weight_cos(prior,sigmab,u0,beta_sum_all,gred_lambda)
                  kuso1<- as.numeric(beta_sum_all%*%weight_temp)
                  kuso1 <- as.numeric(kuso1)
                  cat('weight_ini:',weight_temp,'\n')
                  cat('gred_lambda:',gred_lambda,'\n','ini_kuso1:',kuso1,'\n')
                  weight_beifen <- weight_temp
                  
                  an <- function(n){
                    log(log(n))
                  }
                  #if(size.A0==0){
                  #kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.008,0.013,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.001,min(0.005,0.002+size.A0/8*0.001),0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00001,0.00005+size.A0/8*0.00001,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  #}else{
                  #  if(size.A0==4){
                  #    kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.008,0.01,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.001,0.002,0.0005)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00001,0.00004,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  #  }else{
                  #  kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.008,0.013,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.001,min(0.005,0.002+size.A0/8*0.001),0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00001,0.00004+size.A0/8*0.00001,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  #  }
                  #}
                  if(size.A0==0){
                    kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.005,0.012,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.0005,0.0025,0.0005)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00002,0.00005,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  }else if(size.A0==4){
                    kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.005,0.012,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.0005,0.0025,0.0005)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00002,0.00005,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  }else if(size.A0==8){
                    kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.005,0.012,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.0005,0.0025,0.0005)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00002,0.00005,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  }else if(size.A0==12){
                    kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.012,0.015,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.0015,0.0026,0.0005)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00002,0.00005,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  }else if(size.A0>=16){
                    kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.005,0.012,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.0005,0.0025,0.0005)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00002,0.00005,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  }else{
                    kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.005,0.012,0.001)),expand.grid(1,seq(0.05,0.05,0.01),seq(0.0005,0.0025,0.0005)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00002,0.00005,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                  }
                  #if(nums==19) return(list(X,y,prior))
                  beta_cal <- rep(0,p)
                  beta_cal[as.logical(record)] = kiu
                  
                  cat('mse.ini:',sum((beta_cal-beta0)^2),'\n')
                  
                  cat('ini_kiu:',kiu,'\n')
                  
                  mymse_l2_ini <- c(  mymse_l2_ini , sum((beta_cal-beta0)^2))
                  
                  
                  kuso1 <- as.numeric(as.logical(kiu))*0.15
                  cat('kuso1:',length(kuso1),'\n')
                  fp_ini <- c(fp_ini,fp(kuso1,beta0))
                  fn_ini <- c(fn_ini,fn(kuso1,beta0))
                  
                  cat('ini_fn:',fn(kuso1,beta0),'\n')
                  
                  
                  l1_signal <- 0
                  
                  
                  dimx_temp <- dimx + 1
                  
                  
                  
                  iter_pp <- 500
                  
                  
                  iter_l1_true <- 0
                  iter_sqrt <- 0
                  iter_l1 <- 0
                  for(zz in 1:iter_pp){
                    cat('sizea0:',size.A0,'\n')
                    kuso1<-as.numeric(as.logical(kiu))*0.15
                    cat('iter:',zz,'\n')
                    beta_source_sum<-kuso1*(M+1)
                    beta_cal <- rep(0,p)
                    record[which(record==1)[which(abs(beta_source_sum/(M+1))<=wudi_num)]] = 0
                    cat('sump:',sum(record),'\n')
                    beta0_temp = beta0[which(abs(beta_source_sum/(M+1))>wudi_num)]
                    cat('dimx',dim(X),'\n')
                    cat(length(abs(beta_source_sum/(M+1))>wudi_num),'\n')
                    X<-X[,abs(beta_source_sum/(M+1))>wudi_num]
                    beta_source_sum <- beta_source_sum[abs(beta_source_sum/(M+1))>wudi_num]
                    dimx <- length(beta_source_sum)
                    
                    
                    
                    iter_sqrt <- iter_sqrt + 1
                    lambda_temp <- 0.001
                    u0 <- un(X[1:n0,],y[1:n0])
                    Us <- matrix(rep(0,dimx*(M+1)),dimx,(M+1))
                    for(i in 1:(M)){
                      inde <- (n0+1+(i-1)*(n_s)):(n0+i*n_s)
                      Us[, i] <- un(X[inde,],y[inde])
                    }
                    Us[,M+1] <- u0
                    beta_source_sum1 <- rep(0,dimx)
                    beta_sum_all <- matrix(rep(0,dimx*(M+1)),dimx,(M+1))
                    jiujiu <- 1
                    #if(iter_sqrt>2){
                    #  jiujiu <- 1.25
                    #}
                    if(zz < 60){
                      system.time({
                        
                        beta_sum_all <- foreach(
                          i=1:M,          #输入等待请求的参数
                          .combine=cbind,  #返回结果的整合
                          .packages = c("lmrc")
                          #多个进程共享的系统环境
                        ) %dopar% {
                          an <- function(n){
                            log(log(n))/sqrt(n)
                          }
                          inde <- (n0 + 1 + (i - 1) * n_s):(n0 + i * n_s)
                          lambda_temp <- SA24204174::gic_select(seq(0.065,0.135,0.01)/3*sqrt(p/sum(record))*1.00074^(sum(record)-p)*jiujiu,Us[,i],cov(X[inde,]),n_s,dimx,an,rep(1,dimx),rep(0,dimx),1/sqrt(abs(beta_source_sum)),-1,alpha=0.21)
                          kuso <- SA24204174::lmrc_lasso(dimx,dimx,cov(X[inde,]),L,-1,Us[,i],rep(1,dimx),rep(0,dimx),1/sqrt(abs(beta_source_sum)),lambda_temp,0.21,1e2,1e-6)
                          kuso <- kuso / as.numeric(sqrt(t(kuso) %*% cov(X[1:n0, ]) %*% kuso))
                          return(kuso)
                        }
                      })
                      
                      
                    }else{
                      stopCluster(cl)
                      for(i in 1:M){
                        inde <- (n0+1+(i-1)*(n_s)):(n0+i*n_s)
                        
                        lambda_temp <- gic_select(c(0.09,0.1,0.115,0.11,0.115,0.12)/3*sqrt(p/sum(record))*1.00075^(sum(record)-p)*jiujiu,Us[,i],cov(X[inde,]),dimx,rep(1,dimx),rep(0,dimx),1/sqrt(abs(beta_source_sum)),-1,alpha=0.21)
                        #cat(lambda_temp,' ')
                        kuso<-lmrc_lasso(dimx,dimx,cov(X[inde,]),L,-1,Us[,i],rep(1,dimx),rep(0,dimx),1/sqrt(abs(beta_source_sum)),lambda_temp,0.21,1e2,1e-6)
                        # kuso<-lmrc_lasso(dimx,dimx,cov(X[1:n0,]),L,0,un(X[1:n0,],y[1:n0]),rep(1,dimx),kuso,rep(1,dimx),0.5,0.11,2e2,1e-6)
                        #kuso<-lmrc_lasso(dimx,dimx,cov(X[inde,]),L,1,u0,rep(1,dimx),kuso,rep(1,dimx),0.10,0.11,1e2,1e-6)
                        kuso <- kuso/as.numeric(sqrt(t(kuso)%*%cov(X[1:n0,])%*%kuso))
                        # print(beta_source_sum+kuso)
                        #beta_source_sum <- (beta_source_sum+kuso)
                        beta_sum_all[,i]=kuso
                        beta_source_sum1 <- (beta_source_sum1+kuso)
                      }
                    }
                    
                    beta_cal<-numeric(p)
                    #if(size.A0 <=4){
                    #  gred_lambda <- SA24204174::greedy_cross(1000,as.numeric(prior),X[1:n0,],y[1:n0],as.matrix(beta_sum_all[,1:M]),seq(0.9,1.25,0.01),5,3)
                    #}
                    prior <- rep(1,M)
                    sigmab <- as.matrix(t(beta_sum_all)%*%cov(X[1:150,])%*%beta_sum_all)
                    if(is.singular.matrix(sigmab)){
                      weight_temp <- weight_beifen 
                      kuso1<- as.numeric(beta_sum_all%*%weight_temp)
                      kuso1 <- as.numeric(kuso1)
                      kuso1<-betascale(kuso1,cov(X[1:150,]))
                      cat('weight:',weight_temp,'\n')
                      cat('kuso:',kuso1,'\n')
                      cat('mse_l2:',sum((beta_cal-beta0)^2),'\n')
                    }else{
                      prior <- betascale(prior,sigmab)
                      lambdass <- seq(0.3,0.5,0.02)*sqrt(p)/log(sum(record))*(zz^1.62)
                      if(size.A0 <=4){
                        gred_lambda <- lambda_cross(prior,X[1:150,],y[1:150],beta_sum_all,lambdass,5,3)
                      }else{
                        gred_lambda <- lambda_cross(prior,X[1:150,],y[1:150],beta_sum_all,lambdass,5,3)
                      }
                      
                      weight_temp <- weight_cos(prior,sigmab,u0,beta_sum_all,gred_lambda)
                      weight_beifen <- weight_temp
                      kuso1<- as.numeric(beta_sum_all%*%weight_temp)
                      kuso1 <- as.numeric(kuso1)
                      cat('weight:',weight_temp,'\n')
                      cat('拉姆慕达：',lambdass,'\n')
                      cat('gred_lambda:',gred_lambda,'\n','weight_kuso1:',kuso1,'\n')
                      beta_cal[as.logical(record)] = kuso1
                      cat('mse_l2:',sum((beta_cal*L-beta0)^2),'\n')
                      
                    }

                    beta_cal<-numeric(p)
                    if(size.A0==0 & dimx >16){
                      if(zz==1) kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.015,0.023,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.001,0.004,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0002,0.0003,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 2) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.001,0.035,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.001,0.011,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0005,0.0009,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 3) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.001,0.045,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.001,0.014,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0005,0.0014,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.001,0.051,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.001,0.017,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0005,0.0017,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz > 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.0045,0.055,0.001)+(zz-4)*0.00),expand.grid(1,seq(0.03,0.05,0.01),seq(0.001,0.019,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0001,0.002,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      
                    }
                    if(size.A0==4& dimx >16){
                      if(zz==1) kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.015,0.022,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.001,0.004,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0001,0.0003,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 2) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.025,0.035,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.003,0.008,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0005,0.0009,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 3) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.031,0.045,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.005,0.012,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0009,0.0013,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.035,0.051,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.01,0.015,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.001,0.0015,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz > 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.045,0.055,0.001)+(zz-5)*0.00),expand.grid(1,seq(0.03,0.05,0.01),seq(0.015,0.019,0.001)+(zz-5)*0.00),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0015,0.002,0.0001)+(zz-5)*0.000)),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                    }
                    if(size.A0==8& dimx >16){
                      if(zz==1) kiu <- debias_gic_select(rbind(expand.grid(3,seq(0.02,0.05,0.01),seq(0.00002,0.000025,0.000001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.003,0.006,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0005,0.0006,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 2) kiu<-debias_gic_select(rbind(expand.grid(3,seq(0.02,0.05,0.01),seq(0.000035,0.00004,0.000001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.006,0.008,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0007,0.001,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 3) kiu<-debias_gic_select(rbind(expand.grid(3,seq(0.02,0.05,0.01),seq(0.000031,0.000045,0.000001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.008,0.013,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0009,0.0013,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 4) kiu<-debias_gic_select(rbind(expand.grid(3,seq(0.02,0.05,0.01),seq(0.000041,0.000052,0.000001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.012,0.016,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0014,0.0018,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz > 4) kiu<-debias_gic_select(rbind(expand.grid(3,seq(0.02,0.05,0.01),seq(0.000045,0.000055,0.000001)+(zz-5)*0.000001),expand.grid(1,seq(0.03,0.05,0.01),seq(0.014,0.019,0.001)+(zz-5)*0.001),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0015,0.002,0.0001)+(zz-5)*0.0001)),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                    }
                    if(size.A0==12& dimx >16){
                      if(zz==1) kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.021,0.026,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.003,0.006,0.001)),expand.grid(2,seq(0.04,0.05,0.01),seq(0.0003,0.0005,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 2) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.032,0.037,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.008,0.011,0.001)),expand.grid(2,seq(0.04,0.05,0.01),seq(0.0008,0.0012,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 3) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.04,0.045,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.011,0.015,0.001)),expand.grid(2,seq(0.04,0.06,0.01),seq(0.0012,0.0015,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.043,0.051,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.012,0.016,0.001)),expand.grid(2,seq(0.05,0.07,0.01),seq(0.0014,0.0018,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz > 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.045,0.055,0.001)+(zz-5)*0.001),expand.grid(1,seq(0.03,0.05,0.01),seq(0.013,0.019,0.001)+(zz-5)*0.001),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0015,0.002,0.0001)+(zz-5)*0.0001)),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                    }
                    
                    if(size.A0>12& dimx >16){
                      if(zz==1) kiu <- debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.021,0.027,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.003,0.006,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0003,0.0005,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 2) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.03,0.035,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.007,0.010,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0008,0.0012,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 3) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.035,0.041,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.012,0.016,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0012,0.0016,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz == 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.045,0.055,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.014,0.018,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0015,0.0019,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                      if(zz > 4) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.045,0.055,0.001)+(zz-5)*0.001),expand.grid(1,seq(0.03,0.05,0.01),seq(0.013,0.019,0.001)+(zz-5)*0.001),expand.grid(2,seq(0.02,0.05,0.01),seq(0.0015,0.002,0.0001)+(zz-5)*0.0001)),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                    }
                    if(dimx <=16) kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.001,0.01,0.001)),expand.grid(1,seq(0.03,0.05,0.01),seq(0.001,0.01,0.001)),expand.grid(2,seq(0.02,0.05,0.01),seq(0.00001,0.0001,0.00001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                    cat('kiu:',kiu,'\n')
                    if(sum(which(kiu==0))==0){
                      if(sum(which(kiu==0))==0){
                        if(dimx>16){
                          kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.035,0.04,0.001)),expand.grid(1,seq(0.045,0.05,0.01),seq(0.015,0.017,0.001)),expand.grid(2,seq(0.03,0.05,0.01),seq(0.0016,0.0019,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                          
                          if(sum(which(kiu==0))==0){
                            kiu<-debias_gic_select(rbind(expand.grid(0.5,seq(0.05,0.1,0.01),seq(0.035,0.04,0.001)),expand.grid(1,seq(0.045,0.05,0.01),seq(0.015,0.017,0.001)),expand.grid(2,seq(0.03,0.05,0.01),seq(0.002,0.0022,0.0001))),n0,an,u0,cov(X[1:n0,]),dimx,rep(1,dimx),kuso1,rep(1,dimx),1,1.2,0.21,2e2,1e-6)
                          }
                        }
                        if(sum(which(kiu==0))==0){
                          stopCluster(cl)
                          kuso1<-kiu
                          break
                        }
                      }
                    }
                    
                    
                    kuso1 <- kiu
                    kuso1 <- as.numeric(kuso1)
                    beta_cal=rep(0,p)
                    beta_cal[as.logical(record)] = kuso1
                    beifenkuso<- beta_cal
                    cat('mse_l2:',sum((beta_cal-beta0)^2),'\n')
            
                    if(iter_times_num == 1){
                      fp_1 <- c(fp_1,fp(beta_cal,beta0))
                      fn_1 <- c(fn_1,fn(beta_cal,beta0))
                      p_1 <- c(p_1,dimx)
                      mymse_l2_1 <- c(mymse_l2_1,sum((beta_cal-beta0)^2))
                      cat('fp_1:',fp(beta_cal,beta0),' ','fn_1:',fn(beta_cal,beta0),'\n')
                    }
                    if(iter_times_num == 2){
                      fp_2 <- c(fp_2,fp(beta_cal,beta0))
                      fn_2 <- c(fn_2,fn(beta_cal,beta0))
                      p_2 <- c(p_2,dimx)
                      mymse_l2_2 <- c(mymse_l2_2,sum((beta_cal-beta0)^2))
                      cat('fp_2:',fp(beta_cal,beta0),' ','fn_2:',fn(beta_cal,beta0),'\n')
                    }
                    if(iter_times_num == 3){
                      fp_3 <- c(fp_3,fp(beta_cal,beta0))
                      fn_3 <- c(fn_3,fn(beta_cal,beta0))
                      p_3 <- c(p_3,dimx)
                      mymse_l2_3 <- c(mymse_l2_3,sum((beta_cal-beta0)^2))
                      cat('fp_3:',fp(beta_cal,beta0),' ','fn_3:',fn(beta_cal,beta0),'\n')
                    }
                    if(iter_times_num == 4){
                      fp_4 <- c(fp_4,fp(beta_cal,beta0))
                      fn_4 <- c(fn_4,fn(beta_cal,beta0))
                      p_4 <- c(p_4,dimx)
                      mymse_l2_4 <- c(mymse_l2_4,sum((beta_cal-beta0)^2))
                      cat('fp_4:',fp(beta_cal,beta0),' ','fn_4:',fn(beta_cal,beta0),'\n')
                    }
                    if(iter_times_num == 5){
                      fp_5 <- c(fp_5,fp(beta_cal,beta0))
                      fn_5 <- c(fn_5,fn(beta_cal,beta0))
                      p_5 <- c(p_5,dimx)
                      mymse_l2_5 <- c(mymse_l2_5,sum((beta_cal-beta0)^2))
                      cat('fp_5:',fp(beta_cal,beta0),' ','fn_5:',fn(beta_cal,beta0),'\n')
                    }
                    if(iter_times_num == 6){
                      fp_6 <- c(fp_6,fp(beta_cal,beta0))
                      fn_6 <- c(fn_6,fn(beta_cal,beta0))
                      p_6 <- c(p_6,dimx)
                      mymse_l2_6 <- c(mymse_l2_6,sum((beta_cal-beta0)^2))
                    }
                    iter_times_num <- iter_times_num + 1
                    #if (min(ttm[1:16]) < 0.06 && max(ttm[1:16]) < 0.34) {
                    #对整体乘以1.2
                    # ttm <- ttm * 1.2
                    #}
                    
                  }
                  

                  beta_cal=rep(0,p)
                  beta_cal[as.logical(record)] = kuso1
                  #beta_cal <- beifenkuso#这句话用来直接调用备份保存的kuso，因为不可逆
                  beta_cal <- betascale(beta_cal,Sig.X)
                  mymse_l2_last <- c(mymse_l2_last,sum((beta_cal*L-beta0)^2))
                  p_last <- c(p_last,dimx)
                  fp_last <- c(fp_last,fp(beta_cal,beta0))
                  fn_last <- c(fn_last,fn(beta_cal,beta0))
                  iter_times <- c(iter_times,iter_times_num)


                  cat('kuso1:',kuso1,'\n')
                  cat('zhi:',(as.numeric(beta_cal%*%t(X1[1:n0,])%*%y1[1:n0])/as.numeric(beta_cal%*%t(X1[1:n0,])%*%X1[1:n0,]%*%beta_cal)),'\n')
                  beta_fit <- beta_cal*(as.numeric(beta_cal%*%t(X1[1:n0,])%*%y1[1:n0])/as.numeric(beta_cal%*%t(X1[1:n0,])%*%X1[1:n0,]%*%beta_cal))
                  cat('beta_fit:',beta_fit,'\n')

                  

                  sse_fit <- c(sse_fit,sum((beta_fit -beta0)^2))
                  cat('His methods:',mse.vec,'\n')
                  cat('sse_last:',sum((beta_cal*L-beta0)^2),'\n')
                  cat('sse_fit:',sum((beta_fit -beta0)^2),'\n')
                  cat('fp_oracle',fp(beta.kA,beta0),' ')
                  cat('fn_oracle',fn(beta.kA,beta0),'\n')
                  
                  cat('fp_trans',fp(beta.prop,beta0),' ')
                  cat('fn_trans',fn(beta.prop,beta0),'\n')
                  
                  cat('fp_last',fp(beta_cal,beta0),' ')
                  cat('fn_last',fn(beta_cal,beta0),'\n')
                  
                  if(fp(beta_cal,beta0)==0 & fn(beta_cal,beta0)==0){
                    cm_my <- c(cm_my,1)
                  }else{
                    cm_my <- c(cm_my,0)
                  }
                  
                  if(fp(beta.kA,beta0)==0 & fn(beta.kA,beta0)==0){
                    cm_oracle <- c(cm_oracle,1)
                  }else{
                    cm_oracle <- c(cm_oracle,0)
                  }
                  
                  if(fp(beta.prop,beta0)==0 & fn(beta.prop,beta0)==0){
                    cm_trans <- c(cm_trans,1)
                  }else{
                    cm_trans <- c(cm_trans,0)
                  }
    
                  mse1 <- c(mse1,mse.vec[1])
                  mse2 <- c(mse2,mse.vec[2])
                  mse3 <- c(mse3,mse.vec[3])
                  mse4 <- c(mse4,mse.vec[4])
                  mse5 <- c(mse5,mse.vec[5])
                  mse6 <- c(mse6,mse.vec[6])

                  
                  
                  
                  ps <- c(ps,p)
                  n0s <- c(n0s,n0)
                  sizea0s <- c(sizea0s,size.A0)
                  sigbetas <- c(sigbetas,sig.beta)
                  Ms <- c(Ms,M)
                  hs <- c(hs,h)
                  exacts <- c(exacts,exact)
                  ss <- c(ss,s)
                  
                  
                  
                  #print(betahat)
                  
                  end <- (proc.time())[3][[1]]
                  print(paste('time = ', (end-start), 's', sep=''))

## ----savedata-----------------------------------------------------------------
  df <- data.frame(
    ps = ps,
    n0s = n0s,
    sizea0s = sizea0s,
    sigbetas = sigbetas,
    Ms = Ms,
    hs = hs,
    exacts = exacts,
    ss = ss,
    mse1 = mse1,
    mse2 = mse2,
    mse3 = mse3,
    mse4 = mse4,
    mse5 = mse5,
    mse6 = mse6,
    trans_sse,
    mse_single_trans2,
    mse_single_trans1,
    fp_lasso,
    fn_lasso,
    fp_oracle,
    fn_oracle,
    fp_trans,
    fn_trans,
    fp_trans1,
    fn_trans1,
    fp_trans2,
    fn_trans2,
    fp_target,
    fn_target,
    cm_my,
    cm_oracle,
    cm_trans,
    iter_times,
    p_ini,p_last,
    fp_ini,fp_last,
    fn_ini,fn_last,
    sse_target_only,
    mymse_l2_ini,mymse_l2_last,sse_fit
    
  )
  data_name <-paste('newtranslearningsim_1027_ID',as.character(sims),'_seed',as.character(20231027+sims),'.csv',sep='')
  write.csv(df, file = data_name)

