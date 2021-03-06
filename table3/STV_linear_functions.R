pr0 <- function(theta,sds,alpha){
  a = (alpha-theta)/sds
  b = (-alpha-theta)/sds
  p0 = pnorm(a)-pnorm(b)
  return(p0)
}

plot_p0 <-function(w_grid,f0,fa,p0,main = NULL){
  ylim = range(f0,fa)+c(-0.25,0.25)
  par(mar = c(5,5,2,5))
  plot(w_grid,f0,ylim = ylim,type = "l",lwd = 2, lty = 1,ylab = expression(hat(beta)),xlab = "w",col="grey",main = main)
  lines(w_grid,fa,type = "l",lwd = 2, lty = 1)
  par(new = T)
  plot(w_grid, p0, type="l",lwd = 2, col="grey", axes=F, xlab=NA, ylab=NA, cex=1.2,ylim = c(0,1))
  axis(side = 4)
  mtext(side = 4, line = 3, expression(Pr(hat(beta)==0)))
}
plot_CP <-function(w_grid,truth,STV_cp1,regTV_cp1,loc_cp1,ylab2){
  par(mar = c(5,5,2,5))
  plot(w_grid,STV_cp1,ylim = c(0,1),type = "l",lwd = 2, lty = 1,ylab = "Coverage Probability",xlab = "w")
  abline(h = 0.95,lty=4)
  lines(w_grid,regTV_cp1,type = "l",lwd = 2, lty = 2)
  lines(w_grid,loc_cp1,type = "l",lwd = 2, lty = 3)
  par(new = T)
  plot(w_grid, truth, type="l",lwd = 2, col="grey", axes=F, xlab=NA, ylab=NA, cex=1.2)
  axis(side = 4)
  mtext(side = 4, line = 3, ylab2)
}


lm_beta_CI<-function(w,j,lm_fit2,intcpt = c(1,3,26)){
  corr2 = vcov(lm_fit2)
  coef2 = summary(lm_fit2)$coefficients
  others = matrix((1:nrow(coef2))[-intcpt],ncol = 3)
  coef.indx = rbind(intcpt,others)
  w_w2 = cbind(rep(1,length(w)),w,w^2)
  indx = coef.indx[j,]
  beta = w_w2%*%coef2[indx ,1]
  beta.sd  = sapply(1:length(w), function(x) sqrt(w_w2[x,]%*%corr2[indx,indx]%*%w_w2[x,]))
  return(cbind(beta,beta-1.96*beta.sd,beta+1.96*beta.sd))
}

STV_real <- function(data, y_name, w_name, z_names, q_n = 8,test_indx=NULL){
  data = data[order(data[,w_name]),]
  if(is.null(test_indx)){
    train_indx = 1:nrow(data)
    w_test = data[,w_name]
  }else{
    train_indx = (1:nrow(data))[-test_indx]
    w_test = data[test_indx,w_name]
  }
  y = data[train_indx,y_name]
  w = data[train_indx,w_name]
  z = as.matrix(data[train_indx,z_names])
  n = length(y)
  rho = 1/n^2
  z = cbind(rep(1,n),z)
  p = ncol(z)
  w_all = data[,w_name]
  w_bounds = range(w_all)
  knot_set = quantile(w_all,(1:(q_n-4))/(q_n-3))
  extra = min(diff(knot_set))
  boundary = c(w_bounds[1]-extra,w_bounds[2]+extra)
  B0 = splines::bs(w, knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = boundary)
  B = matrix(B0, nrow=length(w))
  B0test = splines::bs(w_test, knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = boundary)
  B_test = matrix(B0test, nrow=length(w_test))
  gma_alpha = ini_gma(z,y,q_n)
  gma_ini = gma_alpha$gma_ini
  alpha = gma_alpha$alpha
  STV_fit = optim(par = gma_ini, llk_linear_penal_Rcpp, z=as.matrix(z),y=y,q_n=q_n,B = B,alpha= alpha,rho = rho, method ="BFGS",control = list(maxit=200) )
  STV_gma = matrix(STV_fit$par,nrow = q_n, nc = p)
  STV_theta = B%*%STV_gma
  STV_sds = STV_sd(y, z, B, STV_gma, alpha, B,rho=rho)
  #test
  STV_beta = stv_beta_t_thre(STV_gma,q_n=q_n,alpha,B = B)
  STV_theta_test = B_test%*%STV_gma
  STV_sds_test = STV_sd(y, z, B, STV_gma, alpha, B_test,rho=rho)
  STV_beta_test = stv_beta_t_thre(STV_gma,q_n=q_n,alpha,B = B_test)
  return(list(STV_beta = STV_beta,STV_sds =  STV_sds, STV_theta=STV_theta,alpha = alpha, B = B, STV_fit = STV_fit, STV_gma= STV_gma,STV_beta_test = STV_beta_test,STV_sds_test =  STV_sds_test, STV_theta_test=STV_theta_test, B_test = B_test ))
}

plot_3R <-function(w_grid,f0,fa,fb,fc,main = NULL,ylab = expression(beta[1]),xlab = "w"){
  ylim = range(f0,fa,fb,fc)+c(-0.25,0.25)
  plot(w_grid,f0,main = main,ylim = ylim,type = "l",ylab = ylab,xlab = xlab,lty = 1,lwd=2,col="grey") 
  lines(w_grid,fa,lty = 2,lwd=2)
  lines(w_grid,fb,lty = 2,lwd=2)
  lines(w_grid,fc,lty = 2,lwd=2)
}


plot_CI <-function(w_grid,f0,fa,CI,main = NULL,ylab = expression(beta),xlab = "w",vw = NULL){
  fl = CI[,1]
  fu = CI[,2]
  ylim = range(f0,fa,fl,fu)+c(-0.25,0.25)
  plot(w_grid,f0,main = main,ylim = ylim,type = "l",ylab = ylab,xlab = xlab,lty = 1,lwd=2,col="grey") 
  lines(w_grid,fa,lty = 1,lwd=1)
  lines(w_grid,fl,lty = 3,lwd=1)
  lines(w_grid,fu,lty = 3,lwd=1)
  if(!is.null(vw)){
    for(v in vw){
      abline(v = v, col = "grey")
    }
  }

}

plot_3M <-function(w_grid,f0,fa,fb,fc,main = NULL,ylab = expression(beta[1]),xlab = "w"){
  ylim = range(f0,fa,fb,fc)+c(-0.25,0.25)
  plot(w_grid,f0,main = main,ylim = ylim,type = "l",ylab = ylab,xlab = xlab,lty = 1,lwd=2,col="grey") 
  lines(w_grid,fa,lty = 2,lwd=2)
  lines(w_grid,fb,lty = 3,lwd=2)
  lines(w_grid,fc,lty = 4,lwd=2)
}


STV_linear_cv <- function(z,y,w,rho, n_folds = 5, llk_fun){
  N = length(y)
  folds_i = sample(rep(1:n_folds, length.out = N))
  h = c(6,8,10,12,15)
  p = ncol(z)
  L = length(h)
  mse = rep(0,L)
  cv_tmp = matrix(NA,nrow = n_folds, ncol = length(h))
  for(fold in 1:n_folds){
    test_i = which(folds_i == fold)
    z_train = z[-test_i,]
    z_test = z[test_i,]
    y_train = y[-test_i]
    y_test = y[test_i]
    w_train = w[-test_i]
    w_test = w[test_i]
    K = length(w_test)
    for(l in 1:L){
      q_n = h[l]
      B = gen_basis_linear(w, q_n = q_n, w_bounds = c(0,3))
      gma_alpha = ini_gma(z,y,q_n)
      gma_ini = gma_alpha$gma_ini
      alpha = gma_alpha$alpha
      B_train = B[-test_i,]
      B_test = B[test_i,]
      STV_fit = optim(par = gma_ini, llk_fun, z=z_train,y=y_train,q_n=h[l],B = B_train,alpha= alpha,rho=rho, method ="BFGS",control = list(maxit=200) )
      gma =  matrix(STV_fit$par,nrow = q_n, nc = p)
      theta = B_test%*%gma
      beta = stv_beta_t_thre(gma,q_n=q_n,alpha,B = B_test)
      y_hat = rowSums(sapply(1:p,function(x) z_test[,x]*beta[,x]))
      cv_tmp[fold,l] = mean((y_test-y_hat)^2)
    }
  }
  cv <- colMeans(cv_tmp)
  h_opt = h[cv == min(cv)]
  return(h_opt)
}




STV_TFPN_all<-function(beta,beta0){
  rslt = matrix(0,ncol = ncol(beta),nrow = 4)
  for(j in 1:ncol(beta)){
    rslt[,j] = TFPN(beta[,j],beta0[,j])
  }
  return(rslt)
}

TFPN <-function(estimate,truth){
  n.zero = sum(truth == 0)
  n.nonzero = sum(truth !=0)
  p.TP = sum(truth !=0 & estimate != 0)/n.nonzero
  p.TN = sum(truth ==0 & estimate == 0)/n.zero
  p.FP = sum(truth ==0 & estimate != 0)/n.zero
  p.FN = sum(truth !=0 & estimate == 0)/n.nonzero
  return(c(p.TP,p.TN,p.FP,p.FN))
}

# CI for one variable
STV_CI <- function(theta_grid, sd_grid, alpha){
  L = theta_grid-1.96*sd_grid
  R = theta_grid+1.96*sd_grid
  L1 = L
  R1 = R
  for(i in 1:length(theta_grid)){
    if(L[i]>=alpha | R[i] <= -alpha | (L[i]<=-alpha & R[i] >= alpha)){
      L1[i] = thre_beta_t(theta_grid[i], alph = alpha)-1.96*sd_grid[i]
      R1[i] = thre_beta_t(theta_grid[i], alph = alpha)+1.96*sd_grid[i]
    }
    if(abs(R[i])< alpha & L[i] < - alpha){
      A = (alpha - theta_grid[i])/sd_grid[i]
      B = qnorm(0.05-pnorm(-A),lower.tail = FALSE)
      L1[i] = thre_beta_t(theta_grid[i], alph = alpha)-B*sd_grid[i]
      R1[i] = 0
    }
    if(abs(L[i]) < alpha & R[i]>alpha){
      B = (alpha + theta_grid[i])/sd_grid[i] # 
      A = qnorm(0.05-pnorm(-B),lower.tail = FALSE)
      L1[i] = 0
      R1[i] = thre_beta_t(theta_grid[i], alph = alpha)+A*sd_grid[i]
    }
    if(abs(L[i])<=alpha & abs(R[i]) <= alpha){
      L1[i] = 0
      R1[i] = 0
    }
  }
  return(cbind(L1,R1))
}


STV_cp_all <- function(theta_grid, sd_grid, true_beta, alpha){
  cp_mat = NULL
  zero_mat = NULL
  LR_mat = NULL
  for(j in 1:length(alpha)){
    LR = STV_CI(STV_theta_grid[,j], STV_sd_grid[,j], alpha[j])
    cp_temp = (LR[,1] <= true_beta[,j]) &  (LR[,2] >= true_beta[,j])
    cp_mat = cbind(cp_mat,as.numeric(cp_temp))
    zero_temp = !((LR[,1] <= 0) &  (LR[,2] >= 0))
    zero_mat = cbind(zero_mat, as.numeric(zero_temp))
    LR_mat = cbind(LR_mat,LR)
  }
  return(list(cp_mat=cp_mat,zero_mat=zero_mat))
}




STV_sd <-function(y, z, B, gma, alpha, B_grid, rho=0.001){
  z = as.matrix(z)
  theta = B%*%gma
  BB = t(B)%*%B
  n =  length(y)
  n_grid = nrow(B_grid)
  p = ncol(z)
  zh = rep(0,n)
  U = matrix(0, n, p)
  for(j in 1:p){
    zh = zh + as.vector(z[,j])*h_beta_t(epsi = 0.01, theta[,j], alpha[j])
    U[,j] = as.vector(z[,j])*h1_fun(theta[,j], alpha[j],epsi = 0.01)
  }
  sigma2 = sum((y-zh)^2)/(n-1)
  V = sapply(1:n, function(x) kronecker(U[x,],B[x,]))
  D2 =  (-V%*%t(V) - rho * kronecker(diag(rep(1,p)),BB)) # t(V)%*%V is computationally singular.
  D2_inv = MASS::ginv(D2)
  temp = D2_inv%*%V%*%t(V)%*%t(D2_inv)
  
  sd_all = matrix(0,n_grid,p)
  for(j in 1:p){
    ej = rep(0,p)
    ej[j] = 1
    sd_all[,j] = sqrt(sapply(1:n_grid, function(x) t(kronecker(ej,B_grid[x,]))%*%temp%*%kronecker(ej,B_grid[x,])))
  }
  return(sd_all)
}
regTV_sd <-function(y, z, B, gma, B_grid){
  theta = B%*%gma
  BB = t(B)%*%B
  n =  length(y)
  n_grid = nrow(B_grid)
  p = ncol(z)
  zh = rep(0,n)
  U = matrix(0, n, p)
  for(j in 1:p){
    zh = zh + z[,j]*theta[,j]
    U[,j] = z[,j]
  }
  sigma2 = sum((y-zh)^2)/(n-1)
  V = sapply(1:n, function(x) kronecker(U[x,],B[x,]))
  D2 =  -V%*%t(V)  # t(V)%*%V is computationally singular.
  D2_inv = MASS::ginv(D2)
  sd_all = matrix(0,n_grid,p)
  for(j in 1:p){
    ej = rep(0,p)
    ej[j] = 1
    sd_all[,j] = sqrt(-sapply(1:n_grid, function(x) t(kronecker(ej,B_grid[x,]))%*%D2_inv%*%kronecker(ej,B_grid[x,])))
  }
  return(sd_all)
}

ini_gma_fix_alpha = function(z,y,q_n,alpha){
  p = ncol(z)
  lse = rep(0,p)
  alpha = rep(alpha,p)
  for(i in 1:p){
    lse[i] = sum(z[,i]*y)/(sum((z[,i])^2))
  }
  gma_ini = rep(lse+sign(lse)*alpha,each = q_n)
  return(list(gma_ini = gma_ini, alpha = alpha))
}

ini_gma = function(z,y,q_n){
  p = ncol(z)
  lse = rep(0,p)
  alpha = rep(0,p)
  for(i in 1:p){
    lse[i] = sum(z[,i]*y)/(sum((z[,i])^2))
  }
  alpha = 0.5*abs(lse)
  alpha[alpha<0.5] = 0.5
  gma_ini = rep(lse+sign(lse)*alpha,each = q_n)
  return(list(gma_ini = gma_ini, alpha = alpha))
}


STV_linear_simulation_basis = function(n, p, max_time = 3, err_sd_ratio = 0.2, f_list, f_pos, seed = 12, z_var = 0.1, q_n){
  set.seed(seed)
  #w = runif(n, min = 0, max = max_time)
  z_mat = matrix(z_var, p,p)
  diag(z_mat) = rep(2,p)
  #w = runif(n, min =0, max = max_time)
  w = seq(0,max_time,length.out = n)
  z = mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = z_mat)
  z =  as.matrix(z)
  y = 0
  for(j in 1:length(f_list)){
    y = y + z[, f_pos[j]]*f_list[[j]](w,q_n)
  }
  e = rnorm(n, mean = 0, sd = sqrt(err_sd_ratio*var(y)))
  y = y + e
  y = y[order(w)]
  z = z[order(w),]
  w = w[order(w)]
  return(list(z = z, y = y, w = w))
}

STV_linear_simulation = function(n, p, max_time = 3, err_sd_ratio = 0.2, f_list, f_pos, seed = 12, z_var = 0.1){
  set.seed(seed)
  #w = runif(n, min = 0, max = max_time)
  z_mat = matrix(z_var, p,p)
  diag(z_mat) = rep(2,p)
  w = runif(n, min =0, max = max_time)
  z = mvtnorm::rmvnorm(n, mean = rep(0,p), sigma = z_mat)
  y = 0
  for(j in 1:length(f_list)){
    y = y + z[, f_pos[j]]*f_list[[j]](w)
  }
  e = rnorm(n, mean = 0, sd = sqrt(err_sd_ratio*var(y)))
  y = y + e
  y = y[order(w)]
  z = z[order(w),]
  w = w[order(w)]
  return(list(z = z, y = y, w = w))
}


#generate cubic spline basis
gen_basis_linear<-function(w, q_n, w_bounds = c(0,max(w))){
  N =  length(w)
  knot_set = quantile(seq(w_bounds[1],w_bounds[2],length.out = 100),(1:(q_n-4))/(q_n-3))
  extra = min(diff(knot_set))
  boundary = c(w_bounds[1]-extra,w_bounds[2]+extra)
  B0 = splines::bs(w, knot=knot_set, intercept=TRUE, degree=3,Boundary.knots = boundary)
  B = matrix(B0, nrow=N)
  return(B)
}


reg_linear <- function(p,q_n,z,B){  
  z_ = matrix(0, nrow = n, ncol = p*q_n)
  for(j in 1:p){
    z_[,(1+(j-1)*q_n):(j*q_n)] = z[,j]*B
  }
  gma = solve(t(z_)%*%z_)%*%t(z_)%*%y
  beta_t = NULL
  for(j in 1:p){
    beta_t = cbind(beta_t,B%*%gma[(1+(j-1)*q_n):(j*q_n)])
  }
  return(list(gma = gma, beta_t = beta_t))
} 

plot_all <- function(w,beta_t, title = NULL,y_names = NULL){
  library("reshape2")
  library("ggplot2")
  if(!is.null(y_names) & length(y_names) == ncol(beta_t)){
    colnames(beta_t) = y_names
  }
  test_data = as.data.frame(cbind(beta_t,w))
  test_data_long <- melt(test_data, id="w") 
  ggplot(data=test_data_long,aes(x= w, y=value, colour=variable))+geom_line()+ggtitle(title)
}

plot_all_cp <- function(w,beta_t, title = NULL,y_names = NULL){
  library("reshape2")
  library("ggplot2")
  if(!is.null(y_names) & length(y_names) == ncol(beta_t)){
    colnames(beta_t) = y_names
  }
  test_data = as.data.frame(cbind(beta_t,w))
  test_data_long <- melt(test_data, id="w") 
  ggplot(data=test_data_long,aes(x= w, y=value, colour=variable))+geom_line()+ggtitle(title)+geom_hline(yintercept= 0.95, linetype="dashed", color = "black")+ylim(0,1)
}

plot_true <-function(w,n,p,f_pos,f_list){
  true_beta = matrix(0,nrow = n, ncol = p)
  for(j in 1:length(f_pos)){
    true_beta[,f_pos[j]] = f_list[[j]](w)
  }
  plot_all(w,true_beta,title = "Truth")
}

STV_linear_plot_NR = function(gma,z,y,B,alpha,tau,true_fun){
  beta_t = B%*%gma
  u1 = (beta_t-alpha)/tau
  u2 = (beta_t+alpha)/tau
  h_t = (1+2/pi*atan(u1))*u1*tau/2+(1-2/pi*atan(u2))*u2*tau/2
  plot(z,true_fun(z),type = "l",ylab = "Effect",main = paste("Comparison of true and estimated function \n n = ",length(z),"; q_n =  ",length(gma),"; la =  ",alpha,"; tau =  ",tau,sep = ""),ylim = c(min(h_t,true_fun(z)),max(h_t,true_fun(z))))
  lines(z,h_t,col = 2,lty = 2)
  MSE = sum((y-h_t)^2)
  FP = sum(abs(h_t[true_fun(z) == 0])>10e-3)
  FN = sum(abs(true_fun(z)[abs(h_t)<10e-3])>10e-3)
  legend("topright", legend = c("Truth","Estimate",paste("MSE: ",round(MSE,digits  =  3)),paste("FP: ",round(FP,digits  =  3)),paste("FN: ",round(FP,digits  =  3))), lty = c(1,2,0,0,0),col = 1:2)
  return(c(length(z),length(gma),alpha,tau,MSE,FP,FN))
}

### one variable
STV_linear_deriv<-function(gma,z,y,sum_zh,B,alpha = 1,tau = 0.01){
  p = length(gma)
  beta_t = B%*%as.vector(gma)
  u1 = (beta_t-alpha)/tau
  u2 = (beta_t+alpha)/tau
  h_t = (1+2/pi*atan(u1))*u1*tau/2+(1-2/pi*atan(u2))*u2*tau/2
  f1 = (u1/(1+u1^2)+atan(u1)-u2/(1+u2^2)-atan(u2)+pi)/pi
  f2 = 2*(1/(1+u1^2)^2-1/(1+u2^2)^2)/(tau*pi)
  D1 = colSums(2*as.vector((y-sum_zh)*f1)*B)
  D2_pre = 2*(-z*f1^2+(y-sum_zh)*f2)
  n = length(z)
  D2_list = sapply(1:n, function(i){
    B[i,]%*%t(B[i,])*D2_pre[i]
  }
  )    
  D2 = matrix(rowSums(D2_list),nrow  =  p)
  #lik = sum((y-h_t)^2)
  return(list(D1 = D1, D2 = D2, h_t = h_t))
}

h_beta_t0<-function(epsi, beta_t, alph){
  0.5*(1+2/pi*atan((beta_t-alph)/epsi))*(beta_t-alph)+
    0.5*(1-2/pi*atan((beta_t+alph)/epsi))*(beta_t+alph)
}
h_beta_t=Vectorize(h_beta_t0)

stv_beta_t_thre <-function(gma, q_n, alpha, B){
  nc = length(gma)/q_n
  h_t = NULL
  for(j in 1:nc){
    h_t = cbind(h_t,thre_beta_t(B%*%gma[(1+(j-1)*q_n):(j*q_n)], alpha[j]))
  }
  return(h_t)
}


thre_beta_t0<-function( beta_t, alph){
  thre_beta = rep(0,length(beta_t))
  thre_beta[abs(beta_t)>alph] = (sign(beta_t)*(abs(beta_t)-alph))[abs(beta_t)>alph]
  thre_beta
}
thre_beta_t=Vectorize(thre_beta_t0)

STV_linear_NR = function(p,q_n,gma,z,y,B,alpha = 1,tau = 0.01,max_ite = 50,tol = 10e-5){
  ite = 1
  n = nrow(z)
  llk_pre = -10
  repeat{
    for(j in 1:p){
      if(j == 1){ sum_zh = 0 }
      temp = STV_linear_deriv(gma[(1+(j-1)*q_n):(j*q_n)],z[,j],y,sum_zh,B,alpha = alpha[j],tau = 0.01)
      gma[(1+(j-1)*q_n):(j*q_n)] = gma[(1+(j-1)*q_n):(j*q_n)]-MASS::ginv(temp$D2)%*%temp$D1
      beta_t = B%*%gma[(1+(j-1)*q_n):(j*q_n)]
      u1 = (beta_t-alpha[j])/tau
      u2 = (beta_t+alpha[j])/tau
      h_t = (1+2/pi*atan(u1))*u1*tau/2+(1-2/pi*atan(u2))*u2*tau/2
      sum_zh = sum_zh + z[,j]*h_t
    }
    llk_cur = mean((y-sum_zh)^2)
    if(abs(llk_pre-llk_cur) < tol | ite == max_ite ){
      cat(ite)
      break
    }
    ite = ite+1
    llk_pre = llk_cur
  }
  h_t = NULL
  for(j in 1:p){
    h_t = cbind(h_t,h_beta_t(tau, B%*%gma[(1+(j-1)*q_n):(j*q_n)], alpha[j]))
  }
  return(list(gma = gma,lik = llk_cur, beta_t = h_t))
}

llk_linear <-function(par,z,y,q_n,B,alpha,tau){
  nc = length(par)/q_n
  sum_zb = 0
  for(j in 1:nc){
    beta_t = B%*%par[(1+(j-1)*q_n):(j*q_n)]
    u1 = (beta_t-alpha[j])/tau
    u2 = (beta_t+alpha[j])/tau
    h_t = (1+2/pi*atan(u1))*u1*tau/2+(1-2/pi*atan(u2))*u2*tau/2
    sum_zb = sum_zb + z[,j]*h_t
  }
  llk_cur = mean((y-sum_zb)^2)
  return(llk_cur)
}

llk_linear_stv <-function(par,z,y,q_n,B,alpha){
  nc = length(par)/q_n
  sum_zb = 0
  for(j in 1:nc){
    beta_t = B%*%par[(1+(j-1)*q_n):(j*q_n)]
    h_t = thre_beta_t(beta_t,alpha[j])
    sum_zb = sum_zb + z[,j]*h_t
  }
  llk_cur = mean((y-sum_zb)^2)
  return(llk_cur)
}

llk_linear_stv_penalized <-function(par,z,y,q_n,B,alpha,rho){
  nc = length(par)/q_n
  sum_zb = 0
  penalty = 0
  for(j in 1:nc){
    beta_t = B%*%par[(1+(j-1)*q_n):(j*q_n)]
    h_t = thre_beta_t(beta_t,alpha[j])
    sum_zb = sum_zb + z[,j]*h_t
    penalty = penalty + sum((beta_t)^2)
  }
  llk_cur = mean((y-sum_zb)^2) + rho*penalty/length(y)
  return(llk_cur)
}

llk_linear_stv_penalized2 <-function(par,z,y,q_n,B,alpha,penal){
  nc = length(par)/q_n
  sum_zb = 0
  penalty = 0
  for(j in 1:nc){
    beta_t = B%*%par[(1+(j-1)*q_n):(j*q_n)]
    h_t = thre_beta_t(beta_t,alpha[j])
    sum_zb = sum_zb + z[,j]*h_t
    penalty = penalty + sum(abs(beta_t[abs(beta_t)<alpha[j]]))
  }
  llk_cur = mean((y-sum_zb)^2) + penal*penalty/length(y)
  return(llk_cur)
}

penalty_stv <-function(par,z,y,q_n,B,alpha,penal){
  nc = length(par)/q_n
  penalty = 0
  for(j in 1:nc){
    beta_t = B%*%par[(1+(j-1)*q_n):(j*q_n)]
    penalty = penalty + sum((beta_t)^2)
  }
  final = penal*penalty/length(y)
  return(final)
}

penalty_stv2 <-function(par,z,y,q_n,B,alpha,penal){
  nc = length(par)/q_n
  penalty = 0
  for(j in 1:nc){
    beta_t = B%*%par[(1+(j-1)*q_n):(j*q_n)]
    penalty = penalty + sum(abs(beta_t[abs(beta_t)<alpha[j]]))
  }
  final = penal*penalty/length(y)
  return(final)
}

h1_fun <- function(z,alpha,epsi){
  a = (z-alpha)/epsi
  b = (z+alpha)/epsi
  1/pi*a/(1+a^2)+0.5*(1+2/pi*atan(a))-1/pi*b/(1+b^2)+0.5*(1-2/pi*atan(b))
}

h1_fun_simple <- function(z,alpha){
  h1 = rep(0,length(z))
  h1[abs(z)>alpha] = 1
}


llk_linear_reg <-function(par,z,y,q_n,B,alpha,tau){
  nc = length(par)/q_n
  sum_zb = 0
  for(j in 1:nc){
    beta_t = B%*%par[(1+(j-1)*q_n):(j*q_n)]
    sum_zb = sum_zb + z[,j]*beta_t
  }
  llk_cur = mean((y-sum_zb)^2)
  return(llk_cur)
}
stv_beta_t <-function(gma, q_n, alpha, tau, B){
  nc = length(gma)/q_n
  h_t = NULL
  for(j in 1:nc){
    h_t = cbind(h_t,h_beta_t(tau, B%*%gma[(1+(j-1)*q_n):(j*q_n)], alpha[j]))
  }
  return(h_t)
}



sum_zh <-function(parall,index_remain,z,q_n,B,alpha,tau){
  sum_zh = 0
  for(j in index_remain){
    beta_t = B%*%parall[(1+(j-1)*q_n):(j*q_n)]
    u1 = (beta_t-alpha[j])/tau
    u2 = (beta_t+alpha[j])/tau
    h_t = (1+2/pi*atan(u1))*u1*tau/2+(1-2/pi*atan(u2))*u2*tau/2
    sum_zh = sum_zh + z[,j]*h_t
  }
  return(sum_zh)
}


llk_linear_coordinate <-function(par, index_update,z,y,sum_zh_v,q_n,B,alpha,tau){
  j = index_update
  beta_t = B%*%par
  u1 = (beta_t-alpha[j])/tau
  u2 = (beta_t+alpha[j])/tau
  h_t = (1+2/pi*atan(u1))*u1*tau/2+(1-2/pi*atan(u2))*u2*tau/2
  sum_zh = sum_zh_v + z[,j]*h_t
  llk_cur = mean((y-sum_zh)^2)
  return(llk_cur)
}

norm_constrain <-function(par,llk_optim, p,q_n,B,alpha,tau){
  rho = par[1]
  gma = par[-1]
  for(j in 1:p){
    if(j == 1){ sum_zh = 0 }
    beta_t = B%*%gma[(1+(j-1)*q_n):(j*q_n)]
    u1 = (beta_t-alpha[j])/tau
    u2 = (beta_t+alpha[j])/tau
    h_t = (1+2/pi*atan(u1))*u1*tau/2+(1-2/pi*atan(u2))*u2*tau/2
    sum_zh = sum_zh + z[,j]*h_t
  }
  llk_cur = mean((y-sum_zh)^2)
  norm_con = t(gma)%*%gma + rho*(llk_cur-llk_optim)
  return(norm_con)
}
#updated 11/07
get_true_theta<-function(z,B,alpha,tau,true_fun){
  y = true_fun(z)
  gma = rep(0,ncol(B))
  NR = STV_linear_NR(gma,z,y,B,alpha = alpha,tau = tau,max_ite = 100,tol = 10e-7)
  r1 = STV_linear_plot_NR(gma = NR$gma,z,y,B,alpha = alpha,tau = tau,true_fun = true_fun)
  return(as.vector(NR$gma))
}

## coordinate descent, with all variables updated, set gma = 0s, if signal too weak
stv_coord_all<-function(gma_ini = gma_ini, z=z, y= y, p = p, q_n = q_n, B = B, alpha = alpha, tau = tau){
  bar_value = 10e-2
  gma_all = matrix(gma_ini,nrow = q_n,ncol = p)
  remain_ind = rep(1,p)
  old_obj = -10
  repeat{
    index_remain = (1:p)[remain_ind==1]
    for(jj in index_remain){
      index_update = jj
      index_remain = (1:p)[remain_ind==1]
      sum_zh_v = sum_zh(parall=c(gma_all),index_remain[-jj],z,q_n,B = B,alpha,tau)
      par_up_ini = gma_all[,jj]
      rslt_optim3 = optim(par = par_up_ini, llk_linear_coordinate, index_update = index_update,z=z,y=y,sum_zh_v=sum_zh_v, q_n=q_n,B = B,alpha= alpha, tau = tau,method ="BFGS",control = list(maxit=200) )
      par_up = rslt_optim3$par
      gma_all[,index_update] = par_up
      beta_t_up = stv_beta_t(par_up ,q_n = q_n, alpha = alpha[jj],tau = tau , B = B)
      if(mean(abs(beta_t_up)) < bar_value){
        gma_all[,index_update] = rep(0,q_n)
      }
    }
    if(abs(old_obj-rslt_optim3$value)<10e-5) break
    old_obj = rslt_optim3$value
    cat(old_obj,"\n")
  } 
  return(gma_all)
}
## coordinate descent, update all variables, have inference for all variables.
stv_coord_all2<-function(gma_ini = gma_ini, z=z, y= y, p = p, q_n = q_n, B = B, alpha = alpha, tau = tau){
  bar_value = 10e-2
  gma_all = matrix(gma_ini,nrow = q_n,ncol = p)
  remain_ind = rep(1,p)
  old_obj = -10
  repeat{
    index_remain = (1:p)[remain_ind==1]
    for(jj in index_remain){
      index_update = jj
      index_remain = (1:p)[remain_ind==1]
      sum_xh_v = sum_xh(parall=c(gma_all),index_remain[-jj],z,q_n,B = B,alpha,tau)
      par_up_ini = gma_all[,jj]
      rslt_optim3 = optim(par = par_up_ini, llk_linear_coordinate, index_update = index_update,z=z,y=y,sum_zh_v=sum_zh_v, q_n=q_n,B = B,alpha= alpha, tau = tau,method ="BFGS",control = list(maxit=200) )
      par_up = rslt_optim3$par
      gma_all[,index_update] = par_up
    }
    if(abs(old_obj-rslt_optim3$value)<10e-5) break
    old_obj = rslt_optim3$value
    cat(old_obj,"\n")
  } 
  return(gma_all)
}

###########function list for simulation
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

f1b<-function(x, q_n){
  y1 = func_basis(x,fun = function(x) (-x^2+3), q_n = q_n)
  y=y1*(x<=sqrt(3))
  return(y)
}

f2b<-function(x, q_n){
  y1 = func_basis(x,fun = function(x) 2*log(x+0.01), q_n = q_n)
  y=y1*(x>=1)
  return(y)
}

f3b<-function(x, q_n){
  y1 = func_basis(x,fun = function(x) (-3/(x+1)+1), q_n = q_n)
  y=y1*(x<=2)
  return(y)
}
func_basis <- function(w, fun, q_n,w.min = 0, w.max = 3 ){
  knot_set = quantile(seq(w.min,w.max,length.out = 100),(1:(q_n-4))/(q_n-3))
  extra = min(diff(knot_set))
  boundary = c(w.min-extra,w.max+extra)
  B0 = splines::bs(w, degree = 3, intercept = TRUE,knots = knot_set, Boundary.knots = boundary)
  B = matrix(B0,ncol = q_n)
  y = fun(w)
  coef = solve(t(B)%*%B)%*%t(B)%*%y
  y_hat = B%*%coef
  return(y_hat)
}

#shift == 1: shift up
#shift == -1: shift down
func_basis_coef <- function(w, fun, q_n,w.min = 0, w.max = 3, shift = 1 ){
  knot_set = quantile(seq(w.min,w.max,length.out = 100),(1:(q_n-4))/(q_n-3))
  extra = min(diff(quantile(seq(w.min,w.max,length.out = 100),(1:(q_n-4))/(q_n-3))))
  boundary = c(w.min-extra,w.max+extra)
  B0 = splines::bs(w, degree = 3, intercept = TRUE,knots = knot_set, Boundary.knots = boundary)
  B = matrix(B0,ncol = q_n)
  y = fun(w)
  coef = solve(t(B)%*%B)%*%t(B)%*%y
  y_hat = B%*%coef
  if(shift == 1){
    alpha = ceiling(abs(min(range(y_hat))))
    coef_shift = coef + alpha
  }else{
    alpha = ceiling(abs(max(range(y_hat))))
    coef_shift = coef - alpha
  }
  return(list(coef=coef_shift,alpha= alpha,B = B))
}


f1c<-function(x, q_n){
  fit = func_basis_coef(x,fun = function(x) (-x^2+3), q_n = q_n,shift = 1)
  y= thre_beta_t(fit$B%*%fit$coef,fit$alpha)
  return(y)
}

f2c<-function(x, q_n){
  fit = func_basis_coef(x,fun = function(x) 2*log(x+0.01), q_n = q_n,shift = 1)
  y= thre_beta_t(fit$B%*%fit$coef,fit$alpha)
  return(y)
}

f3c<-function(x, q_n){
  fit = func_basis_coef(x,fun = function(x) (-3/(x+1)+1), q_n = q_n,shift = -1)
  y= thre_beta_t(fit$B%*%fit$coef,fit$alpha)
  return(y)
}