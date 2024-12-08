// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

// not user level
Eigen::VectorXd theta_gradient(int n, int p, Eigen::VectorXd beta, Eigen::MatrixXd sigma, double s, Eigen::VectorXd u) {
  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(sigma.block(1, 0, p - 1, 1).data(), sigma.block(1, 0, p - 1, 1).size());
  double a = sigma(0, 0);
  Eigen::MatrixXd c = sigma.block(1, 1, p - 1, p - 1);
  double t0 = (b.dot(beta.tail(p - 1)));
  Eigen::VectorXd t1 = c*beta.tail(p - 1);
  double t2 = 0.5*a;
  double t3 = 1 / pow((pow(t0, 2)+a*(1 - (beta.tail(p - 1).dot(t1)))), 0.5);
  Eigen::VectorXd theta_grad = Eigen::VectorXd::Zero(p - 1);
  Eigen::VectorXd xt(p);
  return -(s * u(0) / a * (t0 * t3 * b - (t2 /sqrt(t0 * t0 + a * (1 - (beta.tail(p - 1).dot(t1)))) * t1 + t2 * t3 * t1)) - u(0) / a * b + u.tail(p - 1));
}


//' @title Pairwise U-statistic computation
//' @description Computes a U-statistic based on pairwise comparisons of rows in a predictor matrix `x` and response vector `y`.
//' @param x A numeric matrix where each row represents an observation and each column represents a feature.
//' @param y A numeric vector representing the response variable, with length equal to the number of rows in `x`.
//' @return A numeric vector of length equal to the number of columns in `x`, representing the computed U-statistic.
//' @examples
//' \dontrun{
//' # Example data
//' x <- matrix(rnorm(1000), ncol = 5) # 200 rows, 5 columns
//' y <- rnorm(200)
//' # Compute U-statistic
//' u_stat <- un(x, y)
//' print(u_stat)
//' }
//' @export
// [[Rcpp::export]]
Eigen::VectorXd un(Eigen::MatrixXd x, Eigen::VectorXd y) {
  int n = y.size();
  int p = x.cols();
  Eigen::VectorXd ustat = Eigen::VectorXd::Zero(p);


  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - i - 1; j++) {
      if (y[i] > y[j + i + 1]) {
        ustat +=  x.row(i)- x.row(j + i + 1);
      }
      else if (y[i] < y[j + i + 1]) {
        ustat -=  x.row(i) - x.row(j + i + 1);
      }
    }
  }
  //Rcout << "The value of n : " << n << "\n";
  //Rcout << "The value of p : " << p << "\n";

  return ustat/n/(n-1);
}

//' @title ReLU Activation Function
//' @description Computes the Rectified Linear Unit (ReLU) activation for a given numeric vector.
//' @param vec A numeric vector (Eigen::VectorXd) for which the ReLU activation is computed.
//' @return A numeric vector (Eigen::VectorXd) with the same dimensions as `vec`, where each element is replaced by `max(0, element)`.
//' @examples
//' \dontrun{
//' # Example vector
//' input_vec <- c(-2, -1, 0, 1, 2)
//' # Compute ReLU
//' output_vec <- relu(input_vec)
//' print(output_vec)  # Should return c(0, 0, 0, 1, 2)
//' }
//' @export
// [[Rcpp::export]]
Eigen::VectorXd relu(Eigen::VectorXd vec){
  return (vec.array() > 0).select(vec, 0);
}



// 退休版本
//Eigen::VectorXd relu1(Eigen::VectorXd vec){
//  Eigen::VectorXd result(vec.size());
//  for (int i = 0; i < vec.size(); i++) {
//    if (vec(i) > 0) {
//      result(i) = vec(i) ;
//    } else {
//      result(i) = 0 ;
//    }
//  }
//  return result;
//}

// not for user // [[Rcpp::export]]
double posi_nega(Eigen::VectorXd theta, Eigen::MatrixXd sigma, double s,int p){
  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(sigma.block(1, 0, p - 1, 1).data(), sigma.block(1, 0, p - 1, 1).size());
  double a = sigma(0, 0);
  Eigen::MatrixXd c = sigma.block(1, 1, p - 1, p - 1);
  return((pow(b.dot(theta),2) + a*(1-theta.transpose()*c*theta)));
}

// not for user // [[Rcpp::export]]
Eigen::VectorXd theta_to_beta(Eigen::VectorXd theta, Eigen::MatrixXd sigma, double s, int p){
  Eigen::VectorXd beta_n(p);
  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(sigma.block(1, 0, p - 1, 1).data(), sigma.block(1, 0, p - 1, 1).size());
  double a = sigma(0, 0);
  Eigen::MatrixXd c = sigma.block(1, 1, p - 1, p - 1);
  beta_n(0) = -b.dot(theta) / a +
    s / a * sqrt(pow(b.dot(theta), 2) +
    a * (1 - theta.transpose() * c * theta));
  beta_n.segment(1, p - 1) = theta;
  return (beta_n);
}

// not for user // [[Rcpp::export]]
Eigen::VectorXd sign(Eigen::VectorXd vec) {
  return (vec.array() > 0).select(Eigen::VectorXd::Ones(vec.size()),
          (vec.array() < 0).select(-Eigen::VectorXd::Ones(vec.size()),
           Eigen::VectorXd::Zero(vec.size())));
}


// not for user // [[Rcpp::export]]
Eigen::VectorXd proximal(Eigen::VectorXd beta, double lambda, double alpha, Eigen::VectorXd beta_original, Eigen::VectorXd w){
  return(beta_original+sign(beta-beta_original).cwiseProduct(relu(((beta-beta_original).cwiseAbs()-alpha*lambda*w))));
}

//' @title Compute Beta-Hat Vector
//' @description Computes the beta-hat vector using the inverse of the covariance matrix and an input vector.
//' @param sigma A covariance matrix (Eigen::MatrixXd) used in the computation.
//' @param u A numeric vector (Eigen::VectorXd) to be transformed.
//' @return A numeric vector (Eigen::VectorXd) representing the beta-hat, scaled to have a unit norm.
//' @examples
//' \dontrun{
//' n <- 100
//' p <- 5
//' beta_eg <- rep(0.5,p)
//' Sig.X <- diag(rep(1,p))
//' X <- rmvnorm(n.vec[k], rep(0, p), Sig.X)
//' y <- X%*%beta_eg + rnorm (n, 0, 1))
//' u <- un(X,y)
//' betahat(cov(X),u)
//' }
//' @export
// [[Rcpp::export]]
Eigen::VectorXd betahat(Eigen::MatrixXd sigma, Eigen::VectorXd u){
  Eigen::VectorXd k = sigma.inverse() * u ;
  return(k / sqrt(u.dot(k)));
}

// not for user
// // [[Rcpp::export]]
double target_value(Eigen::VectorXd beta_new, Eigen::VectorXd u, double lambda){
  return(-u.dot(beta_new) + lambda*(beta_new.array().cwiseAbs().sum()));
}


// 20021009 means inf


// not for user just for debug
// // [[Rcpp::export]]
Eigen::VectorXd lmrc_lasso_s(int n, int p, Eigen::MatrixXd sigma, double L, double s,Eigen::VectorXd u, Eigen::VectorXd beta_ini,
                             Eigen::VectorXd beta_ori, Eigen::VectorXd w, double lambda=0.1,
                             double alpha=0.15, int max_iter=1e2, double tol=5e-6 ){
  Eigen::MatrixXd betas(p, max_iter + 1);
  double dif = 20021009;
  double ratio = alpha/1.645 ;
  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(sigma.block(1, 0, p - 1, 1).data(), sigma.block(1, 0, p - 1, 1).size());
  double a = sigma(0, 0);
  Eigen::MatrixXd c = sigma.block(1, 1, p - 1, p - 1);
  if(w.size() == 0){
    Eigen::VectorXd w = Eigen::VectorXd::Ones(p);
  }
  if(beta_ori.size() == 0){
    Eigen::VectorXd beta_ori = Eigen::VectorXd::Zero(p);
  }
  if(beta_ini.size() == 0){
    if(n >= 1.1*p){
      // 这步有可能会有问题，因为不一定可逆，即使我用了1.1
      Eigen::VectorXd beta_ini = betahat(sigma,u);
    }
    else{
      Eigen::VectorXd beta_ini = Eigen::VectorXd::Zero(p);
      beta_ini(0) = 1 ;
    }
  }
  beta_ini = beta_ini / sqrt(beta_ini.transpose() * sigma * beta_ini);
  betas.col(0) = beta_ini ;
  //Rcpp::Rcout << "The value of   betas.col(0) : " <<   betas.col(0).transpose() << "\n";
  double alpha_now = alpha ;
  double target_min = 20021009 ;
  int target_min_index = 0 ;
  Eigen::VectorXd beta_n(p) ;
  Eigen::VectorXd beta_new(p) ;
  Eigen::VectorXd theta_n(p - 1) ;
  Eigen::VectorXd theta_n1(p - 1) ;
  Eigen::VectorXd theta_grad(p - 1) ;
  Eigen::VectorXd alpha_list(1) ;

  //NEW8.19
  Eigen::VectorXd alpha_nn(1) ;
  alpha_nn(0) = 1;   //alpha_nn(1) = 0.85;   alpha_nn(2) = 0.5;   alpha_nn(3) = 0.15;   alpha_nn(4) = 0.03;

  double alpha_target = 0 ;
  double j = 0.995 ; //投影的比例0-1
  double coef = 0 ;
  double modify_coef = 1 ;
  double temp_target = 0 ;
  int alpha_target_min_index = 0 ;
  for (int i = 0; i < max_iter; i++) {
    if(i == 0){lambda = lambda*1.2; alpha = alpha*1.2;}
    if(i == 5){lambda = lambda/1.2; alpha = alpha/1.2;}
    //Rcout << "The value of n : " << n << "\n";
    //Rcpp::Rcout << "The value of i : " << i << "\n";
    beta_n = betas.col(i) ;
    theta_grad = theta_gradient(n, p, beta_n, sigma, s, u) ;
    //Rcpp::Rcout << "The value of theta_grad  : " << theta_grad  << "\n";
    alpha_now = alpha ;
    ratio = alpha_now/1.645 ;
    alpha_target = 20021009 ;
    alpha_target_min_index = 0;
    modify_coef = 1 ;
    for(int z = 0; z < 1; z++){
      beta_n = betas.col(i) ;
      //Rcpp::Rcout << "The value of alpha_now : " << alpha_now << " ";
      alpha_now = alpha * alpha_nn(z);
      // alpha_now = alpha_now - (z > 0).select(1.0 / (z * z) * ratio, 0.0) ;
      alpha_list[z] = alpha_now ;
      theta_n = beta_n.tail(p - 1) - alpha_now*theta_grad ;

      //Rcpp::Rcout << "The value of posi_nega(theta_n,sigma,s,p)w : " << posi_nega(theta_n,sigma,s,p)<< "\n";
      if(posi_nega(theta_n,sigma,s,p) < 0){
        coef = a/(a*(theta_n).transpose() * c * theta_n-pow(b.dot(theta_n),2)) ;
        //Rcpp::Rcout << "The value of coef : " << coef << "\n";
        theta_n1 = theta_n ;
        theta_n1 = j * sqrt(coef) * theta_n ;
        beta_n = theta_to_beta(theta_n1,sigma,s,p) ;
        beta_new = proximal(beta_n,lambda,alpha_now,beta_ori,w) ;
        //Rcpp::Rcout << "The value of beta_new : " << beta_new << "\n";
        beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
        //modify_coef = j ;
        //temp_target = -u.dot(beta_new) + lambda*(beta_new.array().cwiseAbs().sum()) ;
        //if(alpha_target > temp_target ) {
        //alpha_target_min_index = z ;
        //alpha_target = temp_target ;
        modify_coef = j ;
        //}
      }
      else{
        beta_n = theta_to_beta(theta_n,sigma,s,p) ;
        beta_new = proximal(beta_n,lambda,alpha_now,beta_ori,w) ;
        beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
        //temp_target = -u.dot(beta_new) + lambda*(beta_new.array().cwiseAbs().sum()) ;
        // Rcpp::Rcout << "The value of temp_target : " <<  temp_target << "\n";
        //if(alpha_target > temp_target ) {
        //alpha_target_min_index = z ;
        //alpha_target = temp_target ;
        modify_coef = 1 ;
        //}
      }
    }
    alpha = alpha/1.05 ;
    //alpha_now = alpha_list(alpha_target_min_index) ;
    //Rcpp::Rcout << "The value of alpha_target_min_index : " << alpha_target_min_index << "\n";
    //if(modify_coef == 1){
    //  beta_n = betas.col(i) ;
    //  theta_n = beta_n.tail(p - 1) - alpha_now * theta_grad ;
    //  //Rcpp::Rcout << "The value of beta_n : " << beta_n.transpose() << "\n";
    // beta_n = theta_to_beta(theta_n,sigma,s,p) ;
    //  beta_new = proximal(beta_n,lambda,alpha_now,beta_ori,w) ;
    //  beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
    //}else{
    //  beta_n = betas.col(i) ;
    //  theta_n = beta_n.tail(p - 1) - alpha_now * theta_grad ;
    //  coef = a/(a*(theta_n).transpose() * c * theta_n-pow(b.dot(theta_n),2)) ;
    //  theta_n = modify_coef * sqrt(coef) * theta_n ;
    //  beta_n = theta_to_beta(theta_n,sigma,s,p) ;
    //  beta_new = proximal(beta_n,lambda,alpha_now,beta_ori,w) ;
    //  beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
    //}
    betas.col(i + 1) = beta_new ;
    dif = (betas.col(i)-betas.col(i + 1)).array().square().sum() ;
    temp_target = -u.dot(beta_new) + lambda*(((beta_ori.array()-beta_new.array()).cwiseAbs().sum())) ;
    //Rcpp::Rcout << "The value of temp_target  : " << temp_target  << "\n";
    //Rcpp::Rcout << "The value of beta_new.array()  : " << beta_new.array()  << "\n";

    if(temp_target < target_min){
      target_min = temp_target ;
      target_min_index = i ; // 这里进行了减一，意思是原本target__min_index应该是i+1这里自动减去1了
    }
    if(dif < tol) break;
  }
  //return(beta_new.array()*L);

  return(betas.col(target_min_index)*L) ;
}


//' @title LMRC Lasso Estimator
//' @description Computes the Lasso estimator using the LMRC method for a given set of parameters.
//' @param n Integer. The number of samples.
//' @param p Integer. The number of predictors (features).
//' @param sigma Eigen::MatrixXd. The covariance matrix of the predictors.
//' @param L Double. The Lipschitz constant for the gradient.
//' @param s Double. Indicator for positive (1) or negative (-1) direction. If 0, computes both and chooses the better one.
//' @param u Eigen::VectorXd. The gradient vector.
//' @param beta_ini Eigen::VectorXd. Initial beta estimate.
//' @param beta_ori Eigen::VectorXd. Original beta estimate.
//' @param w Eigen::VectorXd. A weighting vector.
//' @param lambda Double. Regularization parameter. Default is 0.1.
//' @param alpha Double. Step size parameter. Default is 0.15.
//' @param max_iter Integer. Maximum number of iterations. Default is 100.
//' @param tol Double. Convergence tolerance. Default is 5e-6.
//' @return Eigen::VectorXd. The Lasso beta estimate.
//' @examples
//' \dontrun{
//' n <- 100
//' p <- 50
//' beta_eg <- c(rep(0.5,10),p-10)
//' Sig.X <- diag(rep(1,p))
//' X <- rmvnorm(n.vec[k], rep(0, p), Sig.X)
//' y <- X%*%beta_eg + rnorm (n, 0, 1))
//' L <- 1.0
//' s <- 0
//' u <- un(X,y)
//' beta_ini <- rep(1, p)
//' beta_ori <- rep(0, p)
//' w <- rep(1, p)
//' lambda <- 0.1
//' alpha <- 0.15
//' max_iter <- 100
//' tol <- 5e-6
//' beta_lasso <- lmrc_lasso(n, p, cov(X), L, s, u, beta_ini, beta_ori, w, lambda, alpha, max_iter, tol)
//' print(beta_lasso)
//' }
//' @export
// [[Rcpp::export]]
Eigen::VectorXd lmrc_lasso(int n, int p, Eigen::MatrixXd sigma, double L, double s,Eigen::VectorXd u, Eigen::VectorXd beta_ini,
                           Eigen::VectorXd beta_ori, Eigen::VectorXd w, double lambda = 0.1,
                           double alpha = 0.15, int max_iter = 1e2, double tol = 5e-6 ){
  if(s == 0){
    Eigen::VectorXd beta_posi1 = lmrc_lasso_s(n, p, sigma, L, 1, u, beta_ini, beta_ori , w, lambda,alpha, max_iter, tol);
    Eigen::VectorXd beta_nega1 = lmrc_lasso_s(n, p, sigma, L, -1, u, beta_ini, beta_ori , w, lambda,alpha, max_iter, tol);
    if(-u.dot(beta_nega1) + lambda*((beta_ori.array()-beta_nega1.array()).cwiseAbs().sum()) < -u.dot(beta_posi1) + lambda*((beta_ori.array()-beta_posi1.array()).cwiseAbs().sum())
    ){
      return beta_nega1 ;
    }else{
      return beta_posi1 ;
    }
  }else{
    return(lmrc_lasso_s(n, p, sigma, L, s, u, beta_ini, beta_ori , w, lambda,alpha, max_iter, tol)) ;
  }
}



// not for user
// // [[Rcpp::export]]
double trans_target(Eigen::VectorXd beta, Eigen::VectorXd Q, Eigen::VectorXd u0, Eigen::MatrixXd Us, Eigen::VectorXd weight, Eigen::VectorXd prior,
                    double v, double lambda1, double lambda2){
  Eigen::VectorXd bbbb = (weight.cwiseQuotient(prior).log()) ;
  double result = -((v*u0 + (Us * weight).cwiseProduct(Q)).dot(beta)) + lambda1 *(beta.array().cwiseAbs().sum())
    + lambda2 * (weight.dot(bbbb)) ;
  return (result) ;
}


// not for user
Eigen::VectorXd proximal_debias(Eigen::VectorXd beta, double lambda, double alpha, Eigen::VectorXd beta_original, Eigen::VectorXd w, double l1_lambda,Eigen::VectorXd denominator){
  Eigen::VectorXd grad(beta.size()) ;
  for (int i = 0; i < beta.size(); i++) {

    if(beta_original(i) >= 0){
      if(beta(i) >= beta_original(i) + alpha*lambda*w(i) + alpha*l1_lambda){
        grad(i) = beta(i) - alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) ;
      }
      if( beta_original(i) - alpha*lambda*w(i) + alpha*l1_lambda/denominator(i)  <=  beta(i) && beta(i) <= beta_original(i) + alpha*lambda*w(i)+ alpha*l1_lambda/denominator(i)){
        grad(i) = beta_original(i) ;
      }
      if(beta(i) <= beta_original(i) - alpha*lambda*w(i) + alpha*l1_lambda/denominator(i) && -alpha*lambda*w(i) + alpha*l1_lambda/denominator(i) <= beta(i)){
        grad(i) = beta(i) + alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) ;
      }
      if(beta(i) <= - alpha*lambda*w(i) + alpha*l1_lambda/denominator(i) && -alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) <= beta(i)){
        grad(i) = 0 ;
      }
      if( beta(i) <= -alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) ){
        grad(i) = beta(i) + alpha*lambda*w(i) + alpha*l1_lambda/denominator(i) ;
      }
    }
    if(beta_original(i) < 0){
      if(beta(i) >= alpha*lambda*w(i) + alpha*l1_lambda/denominator(i)){
        grad(i) = beta(i) - alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) ;
      }
      if( alpha*lambda*w(i) - alpha*l1_lambda/denominator(i)  <=  beta(i) && beta(i) <=  alpha*lambda*w(i) + alpha*l1_lambda/denominator(i)){
        grad(i) = 0 ;
      }
      if(beta(i) <=   alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) && beta_original(i) + alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) <= beta(i)){
        grad(i) = beta(i) - alpha*lambda*w(i) + alpha*l1_lambda/denominator(i) ;
      }
      if(beta(i) <= beta_original(i) + alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) && beta_original(i) - alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) <= beta(i)){
        grad(i) = beta_original(i)  ;
      }
      if( beta(i) <= beta_original(i) - alpha*lambda*w(i) - alpha*l1_lambda/denominator(i) ){
        grad(i) = beta(i) + alpha*lambda*w(i) + alpha*l1_lambda/denominator(i) ;
      }
    }


  }
  return(grad) ;
}

// not for user
// 所以可以这么说，v只是一个初步值，如果我用双向贪心的话那么可以做到调整v，因为一个情况就是Q的每个q都变小了，那么实际上就相当于调整了v。
// // [[Rcpp::export]]
Eigen::VectorXd lmrc_lasso_s_debias(int n, int p, Eigen::MatrixXd sigma, double L, double s,Eigen::VectorXd u, Eigen::VectorXd beta_ini,
                                    Eigen::VectorXd beta_ori, Eigen::VectorXd w, Eigen::VectorXd denominator,double lambda=0.1,
                                    double alpha=0.15, int max_iter=1e2, double tol=5e-6, double l1_lambda= 0.01 ){
  Eigen::MatrixXd betas(p, max_iter + 1);
  double dif = 20021009;
  double ratio = alpha/1.645 ;
  Eigen::VectorXd b = Eigen::Map<Eigen::VectorXd>(sigma.block(1, 0, p - 1, 1).data(), sigma.block(1, 0, p - 1, 1).size());
  double a = sigma(0, 0);
  Eigen::MatrixXd c = sigma.block(1, 1, p - 1, p - 1);
  if(w.size() == 0){
    Eigen::VectorXd w = Eigen::VectorXd::Ones(p);
  }
  if(beta_ori.size() == 0){
    Eigen::VectorXd beta_ori = Eigen::VectorXd::Zero(p);
  }
  if(beta_ini.size() == 0){
    if(n >= 1.1*p){
      // 这步有可能会有问题，因为不一定可逆，即使我用了1.1
      Eigen::VectorXd beta_ini = betahat(sigma,u);
    }
    else{
      Eigen::VectorXd beta_ini = Eigen::VectorXd::Zero(p);
      beta_ini(0) = 1 ;
    }
  }
  beta_ini = beta_ini / sqrt(beta_ini.transpose() * sigma * beta_ini);
  betas.col(0) = beta_ini ;
  //Rcpp::Rcout << "The value of   betas.col(0) : " <<   betas.col(0).transpose() << "\n";
  double alpha_now = alpha ;
  double target_min = 20021009 ;
  int target_min_index = 0 ;
  Eigen::VectorXd beta_n(p) ;
  Eigen::VectorXd beta_new(p) ;
  Eigen::VectorXd theta_n(p - 1) ;
  Eigen::VectorXd theta_n1(p - 1) ;
  Eigen::VectorXd theta_grad(p - 1) ;
  Eigen::VectorXd alpha_list(1) ;

  //NEW8.19
  Eigen::VectorXd alpha_nn(1) ;
  alpha_nn(0) = 1;   //alpha_nn(1) = 0.85;   alpha_nn(2) = 0.5;   alpha_nn(3) = 0.15;   alpha_nn(4) = 0.03;

  double alpha_target = 0 ;
  double j = 0.995 ; //投影的比例0-1
  double coef = 0 ;
  double modify_coef = 1 ;
  double temp_target = 0 ;
  int alpha_target_min_index = 0 ;
  for (int i = 0; i < max_iter; i++) {
    if(i == 0){lambda = lambda*1.2; alpha = alpha*1.2;}
    if(i == 5){lambda = lambda/1.2; alpha = alpha/1.2;}
    //Rcout << "The value of n : " << n << "\n";
    //Rcpp::Rcout << "The value of i : " << i << "\n";
    beta_n = betas.col(i) ;
    theta_grad = theta_gradient(n, p, beta_n, sigma, s, u) ;
    //Rcpp::Rcout << "The value of theta_grad  : " << theta_grad  << "\n";
    alpha_now = alpha ;
    ratio = alpha_now/1.645 ;
    alpha_target = 20021009 ;
    alpha_target_min_index = 0;
    modify_coef = 1 ;
    for(int z = 0; z < 1; z++){
      beta_n = betas.col(i) ;
      //Rcpp::Rcout << "The value of alpha_now : " << alpha_now << " ";
      alpha_now = alpha * alpha_nn(z);
      // alpha_now = alpha_now - (z > 0).select(1.0 / (z * z) * ratio, 0.0) ;
      alpha_list[z] = alpha_now ;
      theta_n = beta_n.tail(p - 1) - alpha_now*theta_grad ;

      //Rcpp::Rcout << "The value of posi_nega(theta_n,sigma,s,p)w : " << posi_nega(theta_n,sigma,s,p)<< "\n";
      if(posi_nega(theta_n,sigma,s,p) < 0){
        coef = a/(a*(theta_n).transpose() * c * theta_n-pow(b.dot(theta_n),2)) ;
        //Rcpp::Rcout << "The value of coef : " << coef << "\n";
        theta_n1 = theta_n ;
        theta_n1 = j * sqrt(coef) * theta_n ;
        beta_n = theta_to_beta(theta_n1,sigma,s,p) ;
        beta_new = proximal_debias(beta_n,lambda,alpha_now,beta_ori,w,l1_lambda,denominator) ;
        //Rcpp::Rcout << "The value of beta_new : " << beta_new << "\n";
        beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
        //modify_coef = j ;
        //temp_target = -u.dot(beta_new) + lambda*(beta_new.array().cwiseAbs().sum()) ;
        //if(alpha_target > temp_target ) {
        //alpha_target_min_index = z ;
        //alpha_target = temp_target ;
        modify_coef = j ;
        //}
      }
      else{
        beta_n = theta_to_beta(theta_n,sigma,s,p) ;
        beta_new = proximal_debias(beta_n,lambda,alpha_now,beta_ori,w,l1_lambda,denominator) ;
        beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
        //temp_target = -u.dot(beta_new) + lambda*(beta_new.array().cwiseAbs().sum()) ;
        // Rcpp::Rcout << "The value of temp_target : " <<  temp_target << "\n";
        //if(alpha_target > temp_target ) {
        //alpha_target_min_index = z ;
        //alpha_target = temp_target ;
        modify_coef = 1 ;
        //}
      }
    }
    alpha = alpha/1.05 ;
    //alpha_now = alpha_list(alpha_target_min_index) ;
    //Rcpp::Rcout << "The value of alpha_target_min_index : " << alpha_target_min_index << "\n";
    //if(modify_coef == 1){
    //  beta_n = betas.col(i) ;
    //  theta_n = beta_n.tail(p - 1) - alpha_now * theta_grad ;
    //  //Rcpp::Rcout << "The value of beta_n : " << beta_n.transpose() << "\n";
    // beta_n = theta_to_beta(theta_n,sigma,s,p) ;
    //  beta_new = proximal(beta_n,lambda,alpha_now,beta_ori,w) ;
    //  beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
    //}else{
    //  beta_n = betas.col(i) ;
    //  theta_n = beta_n.tail(p - 1) - alpha_now * theta_grad ;
    //  coef = a/(a*(theta_n).transpose() * c * theta_n-pow(b.dot(theta_n),2)) ;
    //  theta_n = modify_coef * sqrt(coef) * theta_n ;
    //  beta_n = theta_to_beta(theta_n,sigma,s,p) ;
    //  beta_new = proximal(beta_n,lambda,alpha_now,beta_ori,w) ;
    //  beta_new = beta_new/sqrt(beta_new.transpose() * sigma * beta_new) ;
    //}
    betas.col(i + 1) = beta_new ;
    dif = (betas.col(i)-betas.col(i + 1)).array().square().sum() ;
    temp_target = -u.dot(beta_new)  +lambda* ((beta_ori.array()-beta_new.array()).cwiseAbs().sum()) + l1_lambda*(beta_new.array().cwiseAbs().sum()) ;
    //Rcpp::Rcout << "The value of temp_target  : " << temp_target  << "\n";
    //Rcpp::Rcout << "The value of beta_new.array()  : " << beta_new.array()  << "\n";

    if(temp_target < target_min){
      target_min = temp_target ;
      target_min_index = i ; // 这里进行了减一，意思是原本target__min_index应该是i+1这里自动减去1了
    }
    if(dif < tol) break;
  }
  //return(beta_new.array()*L);
  // 这里可能需要调整
  return(beta_new*L) ;
}




//' @title LMRC Lasso Debiased Estimator
//' @description Computes a debiased estimator using LMRC Lasso for a given set of parameters.
//' @param n Integer. The number of samples.
//' @param p Integer. The number of predictors (features).
//' @param sigma Eigen::MatrixXd. The covariance matrix of the predictors.
//' @param L Double. The Lipschitz constant for the gradient.
//' @param s Double. Indicator for positive (1) or negative (-1) direction. If 0, computes both and chooses the better one.
//' @param u Eigen::VectorXd. The gradient vector.
//' @param beta_ini Eigen::VectorXd. Initial beta estimate.
//' @param beta_ori Eigen::VectorXd. Original beta estimate.
//' @param w Eigen::VectorXd. A weighting vector.
//' @param denominator Eigen::VectorXd. A denominator vector for scaling.
//' @param lambda Double. Regularization parameter for debiasing. Default is 0.1.
//' @param alpha Double. Step size parameter. Default is 0.15.
//' @param max_iter Integer. Maximum number of iterations. Default is 100.
//' @param tol Double. Convergence tolerance. Default is 5e-6.
//' @param l1_lambda Double. Regularization parameter for L1 penalty. Default is 0.01.
//' @return Eigen::VectorXd. The debiased beta estimate.
//' @examples
//' \dontrun{
//' n <- 100
//' p <- 50
//' beta_eg <- c(rep(0.5,10),p-10)
//' Sig.X <- diag(rep(1,p))
//' X <- rmvnorm(n.vec[k], rep(0, p), Sig.X)
//' y <- X%*%beta_eg + rnorm (n, 0, 1))
//' L <- 1.0
//' u <- un(X,y)
//' ini <- rep(1, p)
//' ori <- rep(0, p)
//' w <- rep(1, p)
//' weight <- rep(1, p)
//' lam1 <- 0.1
//' alpha <- 0.15
//' iter <- 100
//' tol <- 5e-6
//' lam2 <- 0.01
//' s = sign(u[1])
//' bt <- lmrc_lasso_debias(n, p, cov(X), L, s, u, ini, ori, w, weight, lam1, alpha, iter, tol, lam2)
//' cat(bt)
//' }
//' @export
// [[Rcpp::export]]
Eigen::VectorXd lmrc_lasso_debias(int n, int p, Eigen::MatrixXd sigma, double L, double s,Eigen::VectorXd u, Eigen::VectorXd beta_ini,
                                  Eigen::VectorXd beta_ori, Eigen::VectorXd w,Eigen::VectorXd denominator, double lambda = 0.1,
                                  double alpha = 0.15, int max_iter = 1e2, double tol = 5e-6 ,double l1_lambda= 0.01){
  if(s == 0){
    Eigen::VectorXd beta_posi1 = lmrc_lasso_s_debias(n, p, sigma, L, 1, u, beta_ini, beta_ori , w,denominator, lambda,alpha, max_iter, tol, l1_lambda);
    Eigen::VectorXd beta_nega1 = lmrc_lasso_s_debias(n, p, sigma, L, -1, u, beta_ini, beta_ori , w,denominator, lambda,alpha, max_iter, tol, l1_lambda);
    if(-u.dot(beta_nega1) + lambda*((beta_ori.array()-beta_nega1.array()).cwiseAbs().sum()) + l1_lambda*(beta_nega1.array().cwiseAbs().sum()) < -u.dot(beta_posi1)
         + lambda*((beta_ori.array()-beta_posi1.array()).cwiseAbs().sum()) + l1_lambda*(beta_posi1.array().cwiseAbs().sum())
    ){
      return (beta_nega1) ;
    }else{
      return (beta_posi1) ;
    }
  }else{
    return(lmrc_lasso_s_debias(n, p, sigma, L, s, u, beta_ini, beta_ori , w,denominator, lambda,alpha, max_iter, tol, l1_lambda)) ;
  }
}


//' @title Kullback-Leibler (KL) Divergence
//' @description Computes the Kullback-Leibler (KL) divergence between two probability distributions represented by vectors.
//' @param weight A numeric vector (Eigen::VectorXd) representing the posterior distribution (weights).
//' @param prior A numeric vector (Eigen::VectorXd) representing the prior distribution.
//' @return A numeric scalar representing the KL divergence.
//' @examples
//' \dontrun{
//' # Example usage
//' weight <- c(0.4, 0.6)
//' prior <- c(0.5, 0.5)
//' kl_div <- kl(weight, prior)
//' print(kl_div)
//' }
//' @export
// [[Rcpp::export]]
double kl(Eigen::VectorXd weight, Eigen::VectorXd prior) {
  int n = weight.size();
  Eigen::VectorXd result(n);
  double sumkl = 0;
  for (int i = 0; i < n; i++) {
    if (weight(i) <  0.00001) {
      result(i) = 0;
    } else {
      result(i) = weight(i) * log(weight(i) / prior(i));
    }
    sumkl = sumkl + result(i);
  }

  return sumkl;
}


int argmax(Eigen::VectorXd vec) {
  int n = vec.size();
  if (n == 0) {
    return -1;  // Return -1 for empty vector
  }

  double max_val = vec[0];
  int max_index = 0;

  for (int i = 1; i < n; i++) {
    if (vec[i] > max_val) {
      max_val = vec[i];
      max_index = i;
    }
  }

  return max_index;
}


















