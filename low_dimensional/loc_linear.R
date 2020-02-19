
loc_poly_one <-function(z,w,y,h, w_grid){
  p = ncol(z)
  e_mat = matrix(0,nr = p, nc = 2*p)
  for(j in 1:p){
    e_mat[j,j*2-1] = 1
  }
  u0 = w_grid
  z0 = NULL
  for(j in 1:p){
    z0 = cbind(z0,z[,j],z[,j]*(w-u0))
  }
  wt = dnorm((w-u0)/h)/h
  w0 = diag(wt)
  beta = e_mat%*%(MASS::ginv(t(z0)%*%w0%*%z0))%*%t(z0)%*%w0%*%y
  return(beta)
}



loc_linear_fit_no_cv <- function(x,y,w,w_grid,h){
  N = length(y)
  p = ncol(x)
  L = length(h)
  mse = rep(0,L)
  h_opt = h[1]
  K = length(w_grid)
  # beta_grid = matrix(0,nrow = K,ncol = p)
  beta_grid2 = matrix(0,nrow = K,ncol = p)
  for(k in 1:K){
    #  opt_par = optim(par = rep(0,p*2),loc_lik,u0 = w_grid[k],x=x, y=y,w=w, h=h_opt)$par
    beta = loc_poly_one(z = x,w = w,y=y,h = h_opt, w_grid = w_grid[k])
    beta_grid2[k,] = beta
    #    beta_grid[k,] = opt_par[1:p]
  }
  return(list(beta_grid = beta_grid2, h_opt = h_opt))
}


loc_linear_fit <- function(x,y,w,w_grid, n_folds = 5){
  N = length(y)
  folds_i = sample(rep(1:n_folds, length.out = N))
  h = c(N^{-1/18},N^{-1/9},N^{-2/9},N^{-3/9})
  p = ncol(x)
  L = length(h)
  mse = rep(0,L)
  cv_tmp = matrix(NA,nrow = n_folds, ncol = length(h))
  for(fold in 1:n_folds){
    test_i = which(folds_i == fold)
    x_train = x[-test_i,]
    x_test = x[test_i,]
    y_train = y[-test_i]
    y_test = y[test_i]
    w_train = w[-test_i]
    w_test = w[test_i]
    K = length(w_test)
    for(l in 1:L){
      beta_grid = matrix(0,nrow = K,ncol = p)
      sse = 0
      for(k in 1:K){
        u0 = w_test[k]
        #opt_par = optim(par = rep(0,p*2),loc_lik,u0 = w_test[k],x=x_train, y=y_train,w=w_train, h=h[l])$par
        beta = loc_poly_one(z = x_train,w = w_train,y=y_train,h = h[l], w_grid = w_test[k])
        sse = sse + (y_test[k]-x_test[k,]%*%beta)^2
      }
      cv_tmp[fold,l] = sse/K
    }
  }
  cv <- colMeans(cv_tmp)
  h_opt = h[cv == min(cv)]
  K = length(w_grid)
  beta_grid = matrix(0,nrow = K,ncol = p)
  for(k in 1:K){
    opt_par = optim(par = rep(0,p*2),loc_lik,u0 = w_grid[k],x=x, y=y,w=w, h=h_opt)$par
    beta = loc_poly_one(z = x,w = w,y=y,h = h_opt, w_grid = w_grid[k])
    beta_grid[k,] = opt_par[1:p]
  }
  return(list(beta_grid = beta_grid, h_opt = h_opt))
}


loc_lik <- function(par,u0,x,y,w,h){
  p = ncol(x)
  n = nrow(x)
  coeff = matrix(par,ncol=2, nrow = p)
  U_mat = rbind(rep(1,n),w-u0)
  wt = dnorm((w-u0)/h)/h
  ssek = 0
  for(i in 1:n){
    ssek = ssek + (y[i] - x[i,]%*%coeff%*%U_mat[,i])^2*wt[i]
  }
  return(ssek)
}


loc_linear_optim <- function(x,y,w,w_grid, n_folds = 5){
  N = length(y)
  n_train = N
  folds_i = sample(rep(1:n_folds, length.out = n_train))
  h = c(N^{-1/9},N^{-2/9},N^{-3/9})/2
  p = ncol(x)
  L = length(h)
  mse = rep(0,L)
  cv_tmp = matrix(NA,nrow = n_folds, ncol = length(h))
  for(fold in 1:n_folds){
    test_i = which(folds_i == fold)
    x_train = x[-test_i,]
    x_test = x[test_i,]
    y_train = y[-test_i]
    y_test = y[test_i]
    w_train = w[-test_i]
    w_test = w[test_i]
    K = length(w_test)
    for(l in 1:L){
      beta_grid = matrix(0,nrow = K,ncol = p)
      sse = 0
      for(k in 1:K){
        u0 = w_test[k]
        opt_par = optim(par = rep(0,p*2),loc_lik,u0 = w_test[k],x=x_train, y=y_train,w=w_train, h=h[l])$par
        sse = sse + (y_test[k]-x_test[k,]%*%opt_par[1:p])^2
      }
      cv_tmp[fold,l] = sse/K
    }
  }
  cv <- colMeans(cv_tmp)
  h_opt = h[cv == min(cv)]
  K = length(w_grid)
  beta_grid = matrix(0,nrow = K,ncol = p)
  for(k in 1:K){
    u0 = w_grid[k]
    opt_par = optim(par = rep(0,p*2),loc_lik,u0 = w_grid[k],x=x, y=y,w=w, h=h_opt)$par
    beta_grid[k,] = opt_par[1:p]
  }
  return(list(beta_grid = beta_grid, h_opt = h_opt))
}



sigma_loc_poly <-function(z,w,y,h, w_grid){
  p = ncol(z)
  e_mat = matrix(0,nr = p, nc = 2*p)
  for(j in 1:p){
    e_mat[j,j*2-1] = 1
  }
  sigma2_all = matrix(0,nrow = length(w_grid),nc = p)
  for(i in 1:length(w_grid)){
    u0 = w_grid[i]
    z0 = NULL
    for(j in 1:p){
      z0 = cbind(z0,z[,j],z[,j]*(w-u0))
    }
    wt = dnorm((w-u0)/h)/h
    w0 = diag(wt)
    zwt = t(z0)%*%w0%*%z0
    zwt2 = t(z0)%*%w0^2%*%z0
    A_inv = MASS::ginv(zwt)
    y_hat = z0%*%(A_inv)%*%t(z0)%*%(wt*y)
    sigma2_hat = 1/(sum(diag(w0))-sum(diag(A_inv%*%zwt2)))*sum((y-y_hat)^2*wt)
    sigma2_beta = diag(e_mat%*%(A_inv)%*%zwt2%*%(A_inv)%*%t(e_mat)*sigma2_hat)
    sigma2_all[i,] = sigma2_beta
  }
  return(sqrt(sigma2_all))
}

# 
# h = 0.1
# n = 100
# p = 2
# 
# x = matrix(rnorm(n*p),ncol = p)
# w = runif(n,min = 0, max=5)
# y = x[,1]*sin(w)+x[,2]*cos(w)
# 
# l = 1
# K = 100
# w_grid = seq(0,5,length.out = K)
# L = 10  #number of candidate h.
# h = seq(0.03,by=(0.15-0.03)/(L-1),length=L) # bandwidth h
# true_m = cbind(sin(w_grid),cos(w_grid))
#  
# loc_fit = loc_linear(x,y,w,w_grid, h, true_m)
