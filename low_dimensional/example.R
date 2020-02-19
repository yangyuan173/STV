## This file is to demonstrate how our method works. 
## Below is a toy example with three varying coefficients to be estimated.

rm(list = ls())
source('STV_linear_functions.R')
source('loc_linear.R')
Rcpp::sourceCpp('STV_linear_cpp.cpp')

n = 500 #sample size
p = 3 #number of covariates
max_time = 3 #maximum of index variable [0,max_time]
err_sd_ratio =  0.1 # noise to signal ratio
z_var = 0.1 # covariance matrix parameter
f1<-function(x){
  (-x^2+3)*(x<=sqrt(3))
}
f1 =  Vectorize(f1)

f2<-function(x){
  2*log(x+0.01)*(x>=1)
}
f2 = Vectorize(f2)

f3<-function(x){
  2*(-3/(x+1)+1)*(x<=2)
}
f3 = Vectorize(f3)

f_list = list(f1,f2,f3)
f_pos = 1:length(f_list)#sample(1:p, size = length(f_list))
tau = 0.01 #paramter in the smooth approximation of threshold operator.


n_grid = 100 #number of grid points
w_grid = seq(0,3,length.out = n_grid)


sim_data = STV_linear_simulation(n, p, max_time = 3, err_sd_ratio = 0.2, f_list, f_pos, seed = 12, cov_type = "AR1",z_rho = 0.5, sigma=1)
# max_time: index variable maximum value
# err_sd_ratio: error to signal ratio
# f_list: list of true betas
# f_pos: positions of true betas
# seed: random seed
# z_rho: z covariance matrix parameter
# sigma: z covariance matrix parameter
z = sim_data$z
y = sim_data$y
w = sim_data$w


stv.fit = SoftThreshVarying(z,y,w,rho=NULL,n_folds=5,w_grid=NULL)
# z: covariates
# w: index variable
# y: linear outcome
# rho: thresholding approximation parameter, defaul is n^{1/2}
# n_folds: cross-validation to choose number of bases
# w_grid: grid points for output estimation. w_grid has to be within the range of w. Default is w. 

j = 1
betaj = stv.fit$beta[,j]
sdj = stv.fit$sds[,j]
plot_CI0(w,betaj,sdj)


j = 2
betaj = stv.fit$beta[,j]
sdj = stv.fit$sds[,j]
plot_CI0(w,betaj,sdj)

j = 3
betaj = stv.fit$beta[,j]
sdj = stv.fit$sds[,j]
plot_CI0(w,betaj,sdj)
