#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
// using namespace arma;

//' @title A Conjugate Gradient method using Rcpp
//' @description Solve a system of linear equations Ax=b
//' @param A the coefficient matrix
//' @param b the augmented matrix
//' @return x0 the solution
//' @import Rcpp
//' @import RcppArmadillo
//' @useDynLib StatComp22067
//' @examples
//' \dontrun{
//' A = matrix(c(4,-2,-1,-2,4,-2,-1,-2,3),3)
//' b = c(0,-2,3)
//' conjugate(A,b)
//' }
//' @export
// [[Rcpp::export]]
arma::mat conjugate(arma::mat A, arma::vec b) {
  int k = 1;
  int n = b.n_elem; 
  arma::mat v0 = arma::zeros(n,1);
  arma::mat x0 = v0;
  arma::mat x1 = v0;
  arma::mat p1 = v0;
  arma::mat r0 = v0;
  arma::mat r1 = v0;
  arma::mat p0 = v0;
  r0.col(0) = -b;
  p0.col(0) = b;
  while(k <= n+1 && arma::norm(r0,"inf") > 0.05){
    arma::mat alpha = r0.t() * r0 / (p0.t() * A * p0);
    x1 = x0 + alpha(0,0) * p0;
    r1 = r0 + alpha(0,0) * A * p0;
    arma::mat beta = r1.t() * r1 / (r0.t() * r0);
    p1 = -r1 + beta(0,0) * p0;
    p0 = p1; 
    r0 = r1;
    x0 = x1;
    k += 1;
  }
  return(x0);
}

//' @title Compute first-order gradient
//' @description Compute first-order gradient of the objective function
//' @param x a positive definite symmetric matrix
//' @param lambda a regularization parameter
//' @param s the sample covariance
//' @return d the first-order gradient
//' @import Rcpp
//' @import RcppArmadillo
//' @useDynLib StatComp22067
//' @examples
//' \dontrun{
//' f_gra1c(diag(3),1,diag(3))
//' }
//' @export
// [[Rcpp::export]]
arma::mat f_gra1c(arma::mat x, double lambda, arma::mat s){
  arma::mat d = s - arma::inv(x) - lambda * x;
  return(d);
}

//' @title Compute second-order gradient
//' @description Compute second-order gradient of the objective function
//' @param x a positive definite symmetric matrix
//' @param lambda a regularization parameter
//' @param s the sample covariance
//' @return d the second-order gradient
//' @import Rcpp
//' @import RcppArmadillo
//' @useDynLib StatComp22067
//' @examples
//' \dontrun{
//' f_gra2c(diag(3),1,diag(3))
//' }
//' @export
// [[Rcpp::export]]
arma::mat f_gra2c(arma::mat x, double lambda, arma::mat s){
  int len = arma::size(x)(0);
  arma::mat iden = arma::eye(len*len,len*len);
  arma::mat d = arma::kron(arma::inv(x),arma::inv(x)) + lambda * iden;
  return(d);
}

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the number of samples
//' @param mu1 mean of the first random variable
//' @param mu2 mean of the second random variable
//' @param sigma1 standard deviation of the first random variable
//' @param sigma2 standard deviation of the second random variable
//' @param rho correlation of the two random variables
//' @return X random sample of size \code{n}
//' @import Rcpp
//' @useDynLib StatComp22067
//' @examples
//' \dontrun{
//' gc <- GibbsC(100)
//' par(mfrow=c(2,1));
//' plot(gc[,1],type='l')
//' plot(gc[,2],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix GibbsC(int N,double mu1=0, double mu2=0, double sigma1=1, double sigma2=1, double rho=0.9) {
  NumericMatrix X(N, 2);
  double x = 0, y = 0;
  X(0,0) = x;
  X(0,1) = y;
  for(int i = 1; i < N; i++) {
    y = X(i-1,1);
    X(i, 0) = rnorm(1, mu1 + rho * (y - mu2) * sigma1/sigma2, sqrt((1 - rho * rho) * sigma1 * sigma1))[0];
    x = X(i, 0);
    X(i, 1) = rnorm(1, mu2 + rho * (x - mu1) * sigma2/sigma1, sqrt((1 - rho * rho) * sigma2 * sigma2))[0];
  }
  return(X);
}
