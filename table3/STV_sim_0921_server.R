seed = 1 ## seed: 1~100; repeat 100 times
source('STV_linear_functions.R')
Rcpp::sourceCpp('STV_linear_cpp.cpp')
#pdf(file=paste("STV_linear_NR_la1.pdf",sep=""),width=6,height=6)
n = 2000
p = 3
max_time = 3
err_sd_ratio =  0.1
z_var = 0.5

f_list = list(f1,f2,f3)
f_pos = 1:length(f_list)#sample(1:p, size = length(f_list))

tau = 0.01
rho = 0.0001

n_grid = 100
w_grid = seq(0,3,length.out = n_grid)
L = 4 #number of candidate h.
p.adjust.m = "BH" #c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
true_m = cbind(f1(w_grid),f2(w_grid),f3(w_grid))

q_n = 8
  cat("loop: ", seed,"\n")
  sim_data = STV_linear_simulation(n, p, max_time = 3, err_sd_ratio =  err_sd_ratio, f_list, f_pos, seed = seed, z_var = z_var)
  z = sim_data$z
  y = sim_data$y
  w = sim_data$w
  #plot_true(w,n,p,f_pos,f_list)
  t0 = proc.time()[3]
  B_grid = gen_basis_linear(w_grid, q_n = q_n, w_bounds = c(0,3))
  B = gen_basis_linear(w, q_n = q_n, w_bounds = c(0,3))
  gma_alpha = ini_gma(z,y,q_n)
  gma_ini = gma_alpha$gma_ini
  alpha = gma_alpha$alpha
  t1 = proc.time()[3]
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
  STV_p = apply(abs(STV_beta_grid/STV_sd_grid),2, function(x) p.adjust(2*pnorm(x,lower.tail = FALSE),method = p.adjust.m))
  STV_zero_c = STV_p<0.05
  STV_zero2 = STV_TFPN_all(STV_zero_c,true_m)
  t3 = proc.time()[3]-t1
  cat(t3,"\n")
  write.table(STV_beta_grid,file = "STV_beta_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_beta_grid,file = "STV_theta_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_sd_grid,file = "STV_sd_grid.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_cp,file = "STV_cp.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_zero1,file = "STV_zero1.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(STV_zero2,file = "STV_zero2.txt",col.names = FALSE,row.names = FALSE, append = TRUE)
  write.table(t(c(STV_ISE,STV_AISE)),file = "STV_IMSE.txt",col.names = FALSE,row.names = FALSE, append = TRUE)


