library(grpreg)
source("STV_linear_functions.R")
Rcpp::sourceCpp("STV_linear_cpp.cpp")
n = 1000
p = 1500
q_n = 12 
alpha.list = 2
p0 = 6
max_time = 3
err_sd_ratio = 0.1
cov_type = "AR1" #"AR1","CS","Ind"
z_rho = 0.5
sgma = 1
f_list = list(f4,f5,f6,f7,f8,f9)
f_pos = 1:6
true_set = rep(0,p)
true_set[f_pos] = 1
tau = 0.01
rho = min(1/n,0.0001)
n_grid = 100
w_grid = seq(0,3,length.out = n_grid)
true_m = cbind(f1(w_grid),f2(w_grid),f3(w_grid))
seed = 1
nsim_all = 100
for(seed in 1:nsim_all){
  cat("loop: ", seed,"\n")
  sim_data = STV_linear_simulation(n, p, max_time = 3, err_sd_ratio =  err_sd_ratio, f_list, f_pos, seed = seed, cov_type = cov_type,z_rho = z_rho, sigma = sgma)
  z = sim_data$z
  y = sim_data$y
  w = sim_data$w
  B = gen_basis_linear(w, q_n = q_n, w_bounds = c(0,3))
  sim_data2 = STV_linear_simulation(n, p, max_time = 3, err_sd_ratio =  err_sd_ratio, f_list, f_pos, seed = seed+nsim_all*2, cov_type = cov_type,z_rho = z_rho, sigma = sgma)
  z2 = sim_data2$z
  y2 = sim_data2$y
  w2 = sim_data2$w
  B2 = gen_basis_linear(w2, q_n = q_n, w_bounds = c(0,3))
  
  true_llk = llk_linear_stv_penalized_true(w,z,y,f_pos,f_list,rho=0)
  zB = create_ZB(z,B)
  group = rep(1:p,each = q_n)
  
  #group = rep(1:25,each=10)
  t0 = proc.time()[3]
  cvfit <- cv.grpreg(zB, y, group,penalty="grLasso",seed = 666)
  t_Lasso = proc.time()[3]-t0
  gma_grLasso = coef(cvfit)[-1]
  llk_grLasso = llk_linear_SCAD(gma_grLasso,z,y-coef(cvfit)[1],q_n,B)
  PMSE_grLasso = llk_linear_SCAD(gma_grLasso,z2,y2-coef(cvfit)[1],q_n,B2)
  norms_grLasso = norm_beta_Rcpp_other(coef(cvfit)[-1],B,q_n)
  select_grLasso = norms_grLasso>0
  TMSE_grLasso = TMSE_other(f_pos,norms_grLasso,gma = coef(cvfit)[-1],B,f_list,w )+coef(cvfit)[1]^2
  FPss_grLasso =  FPFNSeSpLik(TrueBeta=true_set,beta=select_grLasso)
  cou_grLasso = COU_model(select_grLasso,true_set)
  rslt_grLasso = c(FPss_grLasso,cou_grLasso, TMSE_grLasso,llk_grLasso-true_llk,PMSE_grLasso,t_Lasso)
  write.table(t(select_grLasso),file = "select_grLasso.txt",row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(t(rslt_grLasso),file = "rslt_grLasso.txt",row.names = FALSE,col.names = FALSE,append = TRUE)
  
  t0 = proc.time()[3]
  cvfit2 <- cv.grpreg(zB, y, group,penalty="grSCAD",seed = 666)
  t_SCAD = proc.time()[3]-t0
  gma_grSCAD = coef(cvfit2)[-1]
  llk_grSCAD = llk_linear_SCAD(gma_grSCAD,z,y-coef(cvfit2)[1],q_n,B)
  PMSE_grSCAD = llk_linear_SCAD(gma_grSCAD,z2,y2-coef(cvfit2)[1],q_n,B2)
  norms_grSCAD = norm_beta_Rcpp_other(coef(cvfit2)[-1],B,q_n)
  select_grSCAD = norms_grSCAD>0
  FPss_grSCAD =  FPFNSeSpLik(TrueBeta=true_set,beta=select_grSCAD)
  TMSE_grSCAD = TMSE_other(f_pos,norms_grSCAD,gma = gma_grSCAD,B,f_list,w )
  cou_grSCAD = COU_model(select_grSCAD,true_set)
  rslt_grSCAD = c(FPss_grSCAD,cou_grSCAD,TMSE_grSCAD,llk_grSCAD-true_llk,PMSE_grSCAD,t_SCAD)
  write.table(t(select_grSCAD),file = "select_grSCAD.txt",row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(t(rslt_grSCAD),file = "rslt_grSCAD.txt",row.names = FALSE,col.names = FALSE,append = TRUE)
  
  ######## construct z(*)B
  t0 = proc.time()[3]
  stv.fit = cv.stv(z,y,B,nfolds=10,seed=123,alpha.list = alpha.list)
  t_STV = proc.time()[3]-t0
  norms = norm_beta_Rcpp(par_all = stv.fit$gma, B, alpha = stv.fit$alpha, tau=0.01)
  select = norms > 0
  TMSE_STV = TMSE_stv(f_pos,norms,gma = stv.fit$gma,B,f_list,w,stv.fit$alpha)
  cou_STV = COU_model(select,true_set)
  llk_STV = llk_linear_stv_penalized(stv.fit$gma,z,y,q_n,B,stv.fit$alpha,rho=0)
  PMSE_STV = llk_linear_stv_penalized(stv.fit$gma,z2,y2,q_n,B2,stv.fit$alpha,rho=0)
  FPss = FPFNSeSpLik(TrueBeta=true_set,beta=select)
  rslt_STV = c(FPss,cou_STV,TMSE_STV,llk_STV-true_llk,PMSE_STV,t_STV)
  write.table(t(select),file = "select_STV.txt",row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(stv.fit$alpha[1],file = "alpha_STV.txt",row.names = FALSE,col.names = FALSE,append = TRUE)
  write.table(t(rslt_STV),file = "rslt_STV.txt",row.names = FALSE,col.names = FALSE,append = TRUE)
}

rslt_STV = read.table("rslt_STV.txt")
rslt_grSCAD = read.table("rslt_grSCAD.txt")
rslt_grLasso = read.table("rslt_grLasso.txt")
mean_all = rbind(colMeans(rslt_STV),colMeans(rslt_grSCAD),colMeans(rslt_grLasso))
sd_all = rbind(apply(rslt_STV,2,sd),apply(rslt_grSCAD,2,sd),apply(rslt_grLasso,2,sd))
colnames(mean_all) = c("FP", "FN", "TNR", "TPR","C","O","U","TISE","llk_diff","PMSE","Runtime")
rownames(mean_all) = c("STV","grSCAD","grLASSO")

#### Table 4 ####
tab = mean_sd_table(mean_all,sd_all,digits=2,multiply = 1)
tab_final = cbind(rep(n,3),rep(p,3),rep(cov_type,3),rep(z_rho,3),mean_all[,5:7],tab[,c(1:2,8,10,11)])
library(xtable)
print(xtable(tab_final),file = paste("table_",cov_type,"_n",n,"p",p,".txt",sep=""))
write.csv(tab,file = paste("taball_summary_n",n,"p",p,".csv",sep=""))

print(xtable(tab_final),file = paste("summary_1128.txt",sep=""),append = TRUE)


  
