//#include <Rcpp.h>
#include <RcppArmadillo.h>
#include <math.h>

using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]
#define pi 3.14159265
// [[Rcpp::export]]
arma::vec zeta_t_Rcpp(double tau, arma::vec beta_t, double alpha){
  arma::vec a1=beta_t-alpha;
  arma::vec a2=beta_t+alpha;
  arma::vec bt=a1%arma::conv_to<arma::vec>::from(beta_t>alpha)+a2%arma::conv_to<arma::vec>::from(beta_t< -alpha);
  return bt;
}

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

// [[Rcpp::export]]
double llk_linear_penal_coor_Rcpp(arma::vec par, int j_star, arma::mat z, arma::vec y, arma::vec zh, int q_n, arma::mat B, 
                             arma::vec alpha,double rho, double tau=0.01){
  int n = y.n_elem;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  sum_zh = sum_zh + zh;
  double penalty = 0;
  double llk;
  arma::vec thetaj;
  arma::vec beta_t;
  thetaj = par;
  beta_t = B*thetaj;
  arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j_star));
  sum_zh = sum_zh + z.col(j_star)%h_t;
  penalty = penalty + mean(beta_t%beta_t);
  llk = mean(square(y-sum_zh))+rho*penalty;
  return llk;
}

// [[Rcpp::export]]
double llk_linear_penal_coor_Rcpp2(arma::vec par, int j_star, arma::mat z, arma::vec y, arma::vec zh, int q_n, arma::mat B, 
                                  arma::vec alpha,double rho, double tau=0.01){
  int n = y.n_elem;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  sum_zh = sum_zh + zh;
  double penalty = 0;
  double llk;
  arma::vec thetaj;
  arma::vec beta_t;
  thetaj = par;
  beta_t = B*thetaj;
  arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j_star));
  sum_zh = sum_zh + z.col(j_star)%h_t;
  penalty = penalty + mean(thetaj%thetaj);
  llk = mean(square(y-sum_zh))+rho*penalty;
  return llk;
}


// [[Rcpp::export]]
double llk_linear_penal_coor_Rcpp3(arma::vec par, int j_star, arma::mat z, arma::vec y, double penals, arma::vec zh, int q_n, arma::mat B, 
                                   arma::vec alpha,double rho, double tau=0.01){
  int n = y.n_elem;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  sum_zh = sum_zh + zh;
  double penalty = 0;
  double llk;
  arma::vec thetaj;
  arma::vec beta_t;
  thetaj = par;
  beta_t = B*thetaj;
  arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j_star));
  sum_zh = sum_zh + z.col(j_star)%h_t;
  penalty = penals + mean(thetaj%thetaj);
  llk = mean(square(y-sum_zh))+rho*penalty;
  return llk;
}


// [[Rcpp::export]]
List sum_zh_Rcpp2(arma::vec par_all, int j_star, arma::mat z, arma::mat B, int q_n, arma::vec alpha, double tau=0.01){
  int nc = z.n_cols;
  int n = z.n_rows;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  arma::vec thetaj;
  arma::vec beta_t;
  double penals=0;
  for(int j = 0; j < nc; j++){
    if(j != j_star){
      int first = j*q_n;
      int last = first + q_n-1;
      thetaj = par_all.subvec(first,last);
      arma::vec beta_t = B*thetaj;
      arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j));
      sum_zh = sum_zh + z.col(j)%h_t;
      penals = penals + mean(beta_t%beta_t);
    }
  }
  return Rcpp::List::create(Rcpp::Named("sum_zh") = sum_zh,
                            Rcpp::Named("penals") = penals);
}
// [[Rcpp::export]]
arma::vec sum_zh_Rcpp(arma::vec par_all, int j_star, arma::mat z, arma::mat B, int q_n, arma::vec alpha, double tau=0.01){
  int nc = z.n_cols;
  int n = z.n_rows;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  arma::vec thetaj;
  arma::vec beta_t;
  for(int j = 0; j < nc; j++){
    if(j != j_star){
      int first = j*q_n;
      int last = first + q_n-1;
      thetaj = par_all.subvec(first,last);
      arma::vec beta_t = B*thetaj;
      arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j));
      sum_zh = sum_zh + z.col(j)%h_t;
    }
  }
  return sum_zh;
}

// [[Rcpp::export]]
arma::vec sum_zbeta_Rcpp(arma::vec par_all, arma::mat z, arma::mat B, arma::vec alpha, double tau=0.01){
  int nc = z.n_cols;
  int n = z.n_rows;
  int q_n = B.n_cols;
  arma::vec sum_zh = arma::zeros<arma::vec>(n);
  arma::vec thetaj;
  arma::vec beta_t;
  for(int j = 0; j < nc; j++){
      int first = j*q_n;
      int last = first + q_n-1;
      thetaj = par_all.subvec(first,last);
      arma::vec beta_t = B*thetaj;
      arma::vec h_t = h_beta_t0_Rcpp(tau, beta_t, alpha(j));
      sum_zh = sum_zh + z.col(j)%h_t;
  }
  return sum_zh;
}
// [[Rcpp::export]]
arma::vec norm_beta_Rcpp(arma::vec par_all,  arma::mat B,  arma::vec alpha, double tau=0.01){
  int q_n =B.n_cols;
  int nc  = par_all.n_elem/q_n;
  arma::vec norms(nc);
  arma::vec thetaj;
  arma::vec beta_t;
  for(int j = 0; j < nc; j++){
    int first = j*q_n;
    int last = first + q_n-1;
    thetaj = par_all.subvec(first,last);
    arma::vec beta_t = B*thetaj;
    arma::vec zeta_t = zeta_t_Rcpp(tau, beta_t, alpha(j));
    norms(j) = sqrt(mean(square(zeta_t)));
  }
  return norms;
}

// [[Rcpp::export]]
arma::vec norm_beta_Rcpp_other(arma::vec par_all,  arma::mat B, int q_n){
  int nc  = par_all.n_elem/q_n;
  arma::vec norms(nc);
  arma::vec thetaj;
  arma::vec beta_t;
  for(int j = 0; j < nc; j++){
    int first = j*q_n;
    int last = first + q_n-1;
    thetaj = par_all.subvec(first,last);
    arma::vec beta_t = B*thetaj;;
    norms(j) = sqrt(mean(square(beta_t)));
  }
  return norms;
}

// [[Rcpp::export]]
arma::mat create_ZB(arma::mat z, arma::mat B){
  arma::mat ZB(z.n_rows,z.n_cols*B.n_cols);
  for(int i=0; i < z.n_cols; i++){
    int k = 0;
    for(int j = i*B.n_cols; j< B.n_cols*(i+1); j++){
      ZB.col(j) = z.col(i)%B.col(k);
      k++;
    }
  }
  return ZB;
}
  
  
  

