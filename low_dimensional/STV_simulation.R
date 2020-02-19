## This file is to generate simulation results for low dimensional case
## The appeoximate run time is 4 hours.

rm(list=ls())
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
rho = 1/n^2 #penalty coefficient

n_grid = 100 #number of grid points
w_grid = seq(0,3,length.out = n_grid) #grid points
L = 4 #number of candidate bandwith h in the local polynomial method
true_m = cbind(f1(w_grid),f2(w_grid),f3(w_grid)) #true coefficient functions
for(seed in 1:200){
  cat("loop: ", seed,"\n")
  sim_data = STV_linear_simulation(n, p, max_time = 3, err_sd_ratio = 0.2, f_list, f_pos, seed = seed, cov_type = "AR1",z_rho = z_var, sigma=1)
  z = sim_data$z
  y = sim_data$y
  w = sim_data$w
  
  #### STV ####
  t0 = proc.time()[3]
  q_n = STV_linear_cv(z,y,w,rho, n_folds = 5,llk_fun = llk_linear_penal_Rcpp) 
  B_grid = gen_basis_linear(w_grid, q_n = q_n, w_bounds = c(0,3))
  B = gen_basis_linear(w, q_n = q_n, w_bounds = c(0,3))
  gma_alpha = ini_gma(z,y,q_n)
  gma_ini = gma_alpha$gma_ini
  alpha = gma_alpha$alpha
  STV_fit = optim(par = gma_ini, llk_linear_penal_Rcpp, z=z,y=y,q_n=q_n,B = B,alpha= alpha,rho=rho, method ="BFGS",control = list(maxit=200) )
  STV_llk = STV_fit$value
  STV_gma = matrix(STV_fit$par,nrow = q_n, nc = p)
  STV_beta_grid = stv_beta_t_thre(STV_gma,q_n=q_n,alpha,B = B_grid)
  STV_theta_grid = B_grid%*%STV_gma
  STV_sd_grid = STV_sd(y, z, B, STV_gma, alpha, B_grid)
  STV_ISE = apply((STV_beta_grid  - true_m)^2,2,mean)
  STV_AISE = mean(STV_ISE)
  STV_cp_zero = STV_cp_all(STV_theta_grid, STV_sd_grid, true_m, alpha)
  STV_cp = STV_cp_zero$cp_mat
  STV_zero = STV_cp_zero$zero_mat
  STV_zero1 = STV_TFPN_all(STV_beta_grid,true_m)
  STV_zero2 = STV_TFPN_all(STV_zero,true_m)
  write.table(q_n,file = "STV_q_n.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  
  #### B-spline method ####
  t1 = proc.time()[3]
  regTV_fit = optim(par = gma_ini, llk_linear_reg, z=z,y=y,q_n=q_n,B = B,method ="BFGS",control = list(maxit=200) )
  regTV_gma = matrix(regTV_fit$par, nrow = q_n, nc = p)
  regTV_beta_grid = B_grid%*%regTV_gma
  regTV_sd_grid = regTV_sd_f(y, z, B, regTV_gma, B_grid)
  regTV_llk = regTV_fit$value
  regTV_ISE = apply((regTV_beta_grid - true_m)^2,2,mean)
  regTV_AISE = mean(regTV_ISE)
  regTV_cp = ((regTV_beta_grid-1.96*regTV_sd_grid)<= true_m) &  ((regTV_beta_grid+1.96*regTV_sd_grid)>= true_m)
  regTV_zero = !(((regTV_beta_grid-1.96*regTV_sd_grid)<= 0) &  ((regTV_beta_grid+1.96*regTV_sd_grid)>= 0))
  regTV_zero2 = STV_TFPN_all(regTV_zero,true_m)
  
  
  #### local polynomial method ####
  t2 = proc.time()[3]
  loc_fit = loc_linear_fit(z,y,w,w_grid)
  loc_beta_grid = loc_fit$beta_grid
  loc_sd_grid = sigma_loc_poly(z,w,y,h = loc_fit$h_opt, w_grid)
  loc_ISE = apply((loc_beta_grid - true_m)^2,2,mean)
  loc_AISE = mean(loc_ISE)
  loc_cp = ((loc_beta_grid-1.96*loc_sd_grid)<= true_m) &  ((loc_beta_grid+1.96*loc_sd_grid)>= true_m)
  loc_zero = !(((loc_beta_grid-1.96*loc_sd_grid)<= 0) &  ((loc_beta_grid+1.96*loc_sd_grid)>= 0))
  loc_zero2 = STV_TFPN_all(loc_zero,true_m)

  t3 = proc.time()[3]
  write.table(STV_beta_grid,file = "STV_beta_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(regTV_beta_grid,file = "regTV_beta_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(loc_beta_grid,file = "loc_beta_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_sd_grid,file = "STV_sd_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(regTV_sd_grid,file = "regTV_sd_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(loc_sd_grid,file = "loc_sd_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_cp,file = "STV_cp.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(regTV_cp,file = "regTV_cp.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(loc_cp,file = "loc_cp.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  
  write.table(STV_zero1,file = "STV_zero1.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_zero2,file = "STV_zero2.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(regTV_zero2,file = "regTV_zero2.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(loc_zero2,file = "loc_zero2.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  
  write.table(t(c(STV_ISE,STV_AISE)),file = "STV_IMSE.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(t(c(regTV_ISE,regTV_AISE)),file = "regTV_IMSE.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(t(c(loc_ISE,loc_AISE)),file = "loc_IMSE.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(t(c(t1-t0,t2-t1,t3-t2)),file = "run_time.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
}


############################
## Results Summary 
############################
IMSE_mall = NULL
IMSE_sdall = NULL
pdfname = paste("n",n,"p",p,"err",err_sd_ratio,"zcor",z_var,sep="")

prefix =  "STV"
beta_all = read.table(paste(prefix,"_beta_grid.txt",sep=""))
cp_all = read.table(paste(prefix,"_cp.txt",sep=""))
sd_all = read.table(paste(prefix,"_sd_grid.txt",sep=""))
STV_cp1 =  apply(matrix(cp_all[,1], nr = n_grid),1, mean)
STV_cp2 =  apply(matrix(cp_all[,2], nr = n_grid),1, mean)
STV_cp3 =  apply(matrix(cp_all[,3], nr = n_grid),1, mean)
STV1 = beta_all[1:n_grid,] # one simulation
STVsd1 = sd_all[1:n_grid,]
STV_m1 = apply(matrix(beta_all[,1], nr = n_grid),1, mean)
STV_m2 = apply(matrix(beta_all[,2], nr = n_grid),1, mean)
STV_m3 = apply(matrix(beta_all[,3], nr = n_grid),1, mean)
STV_m10 = apply(matrix(beta_all[,1], nr = n_grid) == 0,1, mean)
STV_m20 = apply(matrix(beta_all[,2], nr = n_grid) == 0,1, mean)
STV_m30 = apply(matrix(beta_all[,3], nr = n_grid) == 0,1, mean)
STV_zero1_all = read.table(paste(prefix,"_zero1.txt",sep="")) 
STV_zero1_1 =  apply(matrix(STV_zero1_all[,1], nr = 4),1, mean)
STV_zero1_2 =  apply(matrix(STV_zero1_all[,2], nr = 4),1, mean)
STV_zero1_3 =  apply(matrix(STV_zero1_all[,3], nr = 4),1, mean)
STV_zero1_1_sd =  apply(matrix(STV_zero1_all[,1], nr = 4),1, sd)
STV_zero1_2_sd =  apply(matrix(STV_zero1_all[,2], nr = 4),1, sd)
STV_zero1_3_sd =  apply(matrix(STV_zero1_all[,3], nr = 4),1, sd)
STV_zero2_all = read.table(paste(prefix,"_zero2.txt",sep="")) 
STV_zero2_1 =  apply(matrix(STV_zero2_all[,1], nr = 4),1, mean)
STV_zero2_2 =  apply(matrix(STV_zero2_all[,2], nr = 4),1, mean)
STV_zero2_3 =  apply(matrix(STV_zero2_all[,3], nr = 4),1, mean)
STV_zero2_1_sd =  apply(matrix(STV_zero2_all[,1], nr = 4),1, sd)
STV_zero2_2_sd =  apply(matrix(STV_zero2_all[,2], nr = 4),1, sd)
STV_zero2_3_sd =  apply(matrix(STV_zero2_all[,3], nr = 4),1, sd)

prefix =  "regTV"
beta_all = read.table(paste(prefix,"_beta_grid.txt",sep=""))
cp_all = read.table(paste(prefix,"_cp.txt",sep=""))
sd_all = read.table(paste(prefix,"_sd_grid.txt",sep=""))
regTV_cp1 =  apply(matrix(cp_all[,1], nr = n_grid),1, mean)
regTV_cp2 =  apply(matrix(cp_all[,2], nr = n_grid),1, mean)
regTV_cp3 =  apply(matrix(cp_all[,3], nr = n_grid),1, mean)
regTV1 = beta_all[1:n_grid,] # one simulation
regTVsd1 = sd_all[1:n_grid,]
regTV_m1 = apply(matrix(beta_all[,1], nr = n_grid),1, mean)
regTV_m2 = apply(matrix(beta_all[,2], nr = n_grid),1, mean)
regTV_m3 = apply(matrix(beta_all[,3], nr = n_grid),1, mean)
regTV_m10 = apply(matrix(beta_all[,1], nr = n_grid) == 0,1, mean)
regTV_m20 = apply(matrix(beta_all[,2], nr = n_grid) == 0,1, mean)
regTV_m30 = apply(matrix(beta_all[,3], nr = n_grid) == 0,1, mean)
regTV_zero2_all = read.table(paste(prefix,"_zero2.txt",sep="")) 
regTV_zero2_1 =  apply(matrix(regTV_zero2_all[,1], nr = 4),1, mean)
regTV_zero2_2 =  apply(matrix(regTV_zero2_all[,2], nr = 4),1, mean)
regTV_zero2_3 =  apply(matrix(regTV_zero2_all[,3], nr = 4),1, mean)
regTV_zero2_1_sd =  apply(matrix(regTV_zero2_all[,1], nr = 4),1, sd)
regTV_zero2_2_sd =  apply(matrix(regTV_zero2_all[,2], nr = 4),1, sd)
regTV_zero2_3_sd =  apply(matrix(regTV_zero2_all[,3], nr = 4),1, sd)

prefix =  "loc"
beta_all = read.table(paste(prefix,"_beta_grid.txt",sep=""))
cp_all = read.table(paste(prefix,"_cp.txt",sep=""))
sd_all = read.table(paste(prefix,"_sd_grid.txt",sep=""))
loc_cp1 =  apply(matrix(cp_all[,1], nr = n_grid),1, mean)
loc_cp2 =  apply(matrix(cp_all[,2], nr = n_grid),1, mean)
loc_cp3 =  apply(matrix(cp_all[,3], nr = n_grid),1, mean)
loc1 = beta_all[1:n_grid,] # one simulation
locsd1 = sd_all[1:n_grid,]
loc_m1 = apply(matrix(beta_all[,1], nr = n_grid),1, mean)
loc_m2 = apply(matrix(beta_all[,2], nr = n_grid),1, mean)
loc_m3 = apply(matrix(beta_all[,3], nr = n_grid),1, mean)
loc_m10 = apply(matrix(beta_all[,1], nr = n_grid) == 0,1, mean)
loc_m20 = apply(matrix(beta_all[,2], nr = n_grid) == 0,1, mean)
loc_m30 = apply(matrix(beta_all[,3], nr = n_grid) == 0,1, mean)
loc_zero2_all = read.table(paste(prefix,"_zero2.txt",sep="")) 
loc_zero2_1 =  apply(matrix(loc_zero2_all[,1], nr = 4),1, mean)
loc_zero2_2 =  apply(matrix(loc_zero2_all[,2], nr = 4),1, mean)
loc_zero2_3 =  apply(matrix(loc_zero2_all[,3], nr = 4),1, mean)
loc_zero2_1_sd =  apply(matrix(loc_zero2_all[,1], nr = 4),1, sd)
loc_zero2_2_sd =  apply(matrix(loc_zero2_all[,2], nr = 4),1, sd)
loc_zero2_3_sd =  apply(matrix(loc_zero2_all[,3], nr = 4),1, sd)

zero_all_1 = cbind(STV_zero1_1,STV_zero2_1,regTV_zero2_1,loc_zero2_1)  
rownames(zero_all_1) = c("TP","TN","FP","FN")
zero_all_2 = cbind(STV_zero1_2,STV_zero2_2,regTV_zero2_2,loc_zero2_2)  
rownames(zero_all_2) = c("TP","TN","FP","FN")
zero_all_3 = cbind(STV_zero1_3,STV_zero2_3,regTV_zero2_3,loc_zero2_3)  
rownames(zero_all_3) = c("TP","TN","FP","FN")
zeros_all = rbind(zero_all_1,zero_all_2,zero_all_3)

zero_all_1_sd = cbind(STV_zero1_1_sd,STV_zero2_1_sd,regTV_zero2_1_sd,loc_zero2_1_sd)  
rownames(zero_all_1_sd) = c("TP","TN","FP","FN")
zero_all_2_sd = cbind(STV_zero1_2_sd,STV_zero2_2_sd,regTV_zero2_2_sd,loc_zero2_2_sd)  
rownames(zero_all_2_sd) = c("TP","TN","FP","FN")
zero_all_3_sd = cbind(STV_zero1_3_sd,STV_zero2_3_sd,regTV_zero2_3_sd,loc_zero2_3_sd)  
rownames(zero_all_3) = c("TP","TN","FP","FN")
zeros_all_sd = rbind(zero_all_1_sd,zero_all_2_sd,zero_all_3_sd)

zeros_table = rbind(zero_all_1[1,],zero_all_1_sd[1,],
                    zero_all_1[2,],zero_all_1_sd[2,],
                    zero_all_2[1,],zero_all_2_sd[1,],
                    zero_all_2[2,],zero_all_2_sd[2,],
                    zero_all_3[1,],zero_all_3_sd[1,],
                    zero_all_3[2,],zero_all_3_sd[2,])
write.csv(zeros_table,file = "zeros_table.csv")

############# Table 2 ###################
zero_all_1 = cbind(STV_zero1_1,STV_zero2_1,regTV_zero2_1,loc_zero2_1)  
rownames(zero_all_1) = c("TP","TN","FP","FN")
zero_all_2 = cbind(STV_zero1_2,STV_zero2_2,regTV_zero2_2,loc_zero2_2)  
rownames(zero_all_2) = c("TP","TN","FP","FN")
zero_all_3 = cbind(STV_zero1_3,STV_zero2_3,regTV_zero2_3,loc_zero2_3)  
rownames(zero_all_3) = c("TP","TN","FP","FN")
zeros_all = rbind(zero_all_1,zero_all_2,zero_all_3)
xtable::xtable(zeros_all)
Hmisc::latex(round(zeros_all,digits = 3),file = "")

############# Table 1 ###################
IMSE_all = NULL
for(prefix in c("STV","regTV","loc")){  
  IMSE = read.table(paste(prefix,"_IMSE.txt",sep=""))
  IMSE_m = apply(IMSE,2,mean)
  IMSE_sd = apply(IMSE,2,sd)
  
  IMSE_all = rbind(IMSE_all,IMSE_m,IMSE_sd)
}
colnames(IMSE_all) = c("ISE1","ISE2","ISE3","AISE")
write.csv(IMSE_all,file = "IMSE_all.csv")

run_time = read.table("run_time.txt")
run_time_m = apply(run_time,2,mean)
run_time_m# 7/14: 21.136415  0.446155 26.970645 

