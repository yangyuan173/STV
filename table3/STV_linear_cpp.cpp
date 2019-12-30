//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#define pi 3.14159265
// [[Rcpp::export]]
arma::vec h_beta_t0_Rcpp(double tau, arma::vec beta_t, double alpha){
  arma::vec a1=beta_t-alpha;
  arma::vec b1=a1/tau;
  arma::vec a2=beta_t+alpha;
  arma::vec b2=a2/tau;
  arma::vec bt=0.5*(1.0+2.0/pi*atan(b1))%a1+0.5*(1.0-2.0/pi*atan(b2))%a2;
  return bt;
}

// [[Rcpp::export]]
double llk_linear_Rcpp(arma::vec par, arma::mat z, arma::vec y, int q_n, arma::mat B, 
                       arma::vec alpha, double tau){
  int nc = par.n_elem/q_n;
  int n = y.n_elem;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  double llk;
  arma::vec thetaj;
  arma::vec beta_t;
  for(int j = 0; j < nc; j++){
    int first = j*q_n;
    int last = first + q_n-1;
    thetaj = par.subvec(first,last);
    arma::vec beta_t = B*thetaj;
    arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j));
    sum_zh = sum_zh + z.col(j)%h_t;
  }
  llk = mean(square(y-sum_zh));
  return llk;
}

// [[Rcpp::export]]
double llk_linear_penal_Rcpp(arma::vec par, arma::mat z, arma::vec y, int q_n, arma::mat B, 
                       arma::vec alpha,double rho, double tau=0.01){
  int nc = par.n_elem/q_n;
  int n = y.n_elem;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  double penalty = 0;
  double llk;
  arma::vec thetaj;
  arma::vec beta_t;
  for(int j = 0; j < nc; j++){
    int first = j*q_n;
    int last = first + q_n-1;
    thetaj = par.subvec(first,last);
    arma::vec beta_t = B*thetaj;
    arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j));
    sum_zh = sum_zh + z.col(j)%h_t;
    penalty = penalty + mean(beta_t%beta_t);
  }
  llk = mean(square(y-sum_zh))+rho*penalty;
  return llk;
}
