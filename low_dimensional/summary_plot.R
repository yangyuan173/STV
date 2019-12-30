## this file is to generate figures for paper use.
n = 500
out_folder = "./figures"
############################
## Results Summary 
############################
source('STV_linear_functions.R')
n_grid =  100
w_grid = seq(0,3,length.out = n_grid)

prefix =  "STV"
STV_beta_all = read.table(paste(prefix,"_beta_grid.txt",sep=""))
STV_cp_all = read.table(paste(prefix,"_cp.txt",sep=""))
STV_sd_all = read.table(paste(prefix,"_sd_grid.txt",sep=""))
STV_cp1 =  apply(matrix(STV_cp_all[,1], nr = n_grid),1, mean)
STV_cp2 =  apply(matrix(STV_cp_all[,2], nr = n_grid),1, mean)
STV_cp3 =  apply(matrix(STV_cp_all[,3], nr = n_grid),1, mean)
STV1 = STV_beta_all[1:n_grid,] # one simulation
STVsd1 = STV_sd_all[1:n_grid,]
STV_m1 = apply(matrix(STV_beta_all[,1], nr = n_grid),1, mean)
STV_m2 = apply(matrix(STV_beta_all[,2], nr = n_grid),1, mean)
STV_m3 = apply(matrix(STV_beta_all[,3], nr = n_grid),1, mean)
STV_m10 = apply(matrix(STV_beta_all[,1], nr = n_grid) == 0,1, mean)
STV_m20 = apply(matrix(STV_beta_all[,2], nr = n_grid) == 0,1, mean)
STV_m30 = apply(matrix(STV_beta_all[,3], nr = n_grid) == 0,1, mean)

prefix =  "regTV"
regTV_beta_all = read.table(paste(prefix,"_beta_grid.txt",sep=""))
regTV_cp_all = read.table(paste(prefix,"_cp.txt",sep=""))
regTV_sd_all = read.table(paste(prefix,"_sd_grid.txt",sep=""))
regTV_cp1 =  apply(matrix(regTV_cp_all[,1], nr = n_grid),1, mean)
regTV_cp2 =  apply(matrix(regTV_cp_all[,2], nr = n_grid),1, mean)
regTV_cp3 =  apply(matrix(regTV_cp_all[,3], nr = n_grid),1, mean)
regTV1 = regTV_beta_all[1:n_grid,] # one simulation
regTVsd1 = regTV_sd_all[1:n_grid,]
regTV_m1 = apply(matrix(regTV_beta_all[,1], nr = n_grid),1, mean)
regTV_m2 = apply(matrix(regTV_beta_all[,2], nr = n_grid),1, mean)
regTV_m3 = apply(matrix(regTV_beta_all[,3], nr = n_grid),1, mean)
regTV_m10 = apply(matrix(regTV_beta_all[,1], nr = n_grid) == 0,1, mean)
regTV_m20 = apply(matrix(regTV_beta_all[,2], nr = n_grid) == 0,1, mean)
regTV_m30 = apply(matrix(regTV_beta_all[,3], nr = n_grid) == 0,1, mean)

prefix =  "loc"
loc_beta_all = read.table(paste(prefix,"_beta_grid.txt",sep=""))
loc_cp_all = read.table(paste(prefix,"_cp.txt",sep=""))
loc_sd_all = read.table(paste(prefix,"_sd_grid.txt",sep=""))
loc_cp1 =  apply(matrix(loc_cp_all[,1], nr = n_grid),1, mean)
loc_cp2 =  apply(matrix(loc_cp_all[,2], nr = n_grid),1, mean)
loc_cp3 =  apply(matrix(loc_cp_all[,3], nr = n_grid),1, mean)
loc1 = loc_beta_all[1:n_grid,] # one simulation
locsd1 = loc_sd_all[1:n_grid,]
loc_m1 = apply(matrix(loc_beta_all[,1], nr = n_grid),1, mean)
loc_m2 = apply(matrix(loc_beta_all[,2], nr = n_grid),1, mean)
loc_m3 = apply(matrix(loc_beta_all[,3], nr = n_grid),1, mean)
loc_m10 = apply(matrix(loc_beta_all[,1], nr = n_grid) == 0,1, mean)
loc_m20 = apply(matrix(loc_beta_all[,2], nr = n_grid) == 0,1, mean)
loc_m30 = apply(matrix(loc_beta_all[,3], nr = n_grid) == 0,1, mean)

w_grid

fa = STV_m1
fb = regTV_m1
fc = loc_m1

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

f0 = f1(w_grid)

#### Figure 4 ####
plot_CP3 <-function(w_grid,truth,STV_cp1,ylab="Coverage Probability",ylab2,main="Method",cex_value=1.4){
  par(mar = c(5,5,2,5))
  plot(w_grid,STV_cp1,ylim = c(0,1),type = "l",lwd = 2, lty = 1,ylab = ylab,xlab = "w",las=1,main = main,cex.main=cex_value,cex.lab=cex_value,cex.axis=cex_value)
  abline(h = 0.95,lty=4)
  par(new = T)
  plot(w_grid, truth, type="l",lwd = 2, col="grey", axes=F, xlab=NA, ylab=NA, cex=1.2)
  axis(side = 4,las=1,cex.axis=cex_value)
  mtext(side = 4, line = 3, ylab2,cex = cex_value)
}
plot_CP3(w_grid,truth = f1(w_grid),STV_cp1,ylab2 = expression(b[1]),main = "STV")

pdf(file = paste(out_folder,"/n",n,"cpB1stv.pdf",sep = ""),width = 6,height = 6)
plot_CP3(w_grid,truth = f1(w_grid),STV_cp1,ylab2 = expression(b[1]),main = "STV")
dev.off()
pdf(file = paste(out_folder,"/n",n,"cpB1regtv.pdf",sep = ""),width = 6,height = 6)
plot_CP3(w_grid,truth = f1(w_grid),regTV_cp1,ylab2 = expression(b[1]),main = "B-spline")
dev.off()
pdf(file = paste(out_folder,"/n",n,"cpB1loc.pdf",sep = ""),width = 6,height = 6)
plot_CP3(w_grid,truth = f1(w_grid),loc_cp1,ylab2 = expression(b[1]),main = "Local polynomial")
dev.off()




#### Figure 3 ####
#### #### #### #### #### #### 
#### plot median of the simulation
#### #### #### #### #### #### 

STV_beta1 = matrix(STV_beta_all[,1], nr = n_grid)
STV_beta2 = matrix(STV_beta_all[,2], nr = n_grid)
STV_beta3 = matrix(STV_beta_all[,3], nr = n_grid)
regTV_beta1 = matrix(regTV_beta_all[,1], nr = n_grid)
regTV_beta2 = matrix(regTV_beta_all[,2], nr = n_grid)
regTV_beta3 = matrix(regTV_beta_all[,3], nr = n_grid)
loc_beta1 = matrix(loc_beta_all[,1], nr = n_grid)
loc_beta2 = matrix(loc_beta_all[,2], nr = n_grid)
loc_beta3 = matrix(loc_beta_all[,3], nr = n_grid)

betas = STV_beta1
f0 = f1(w_grid)
plot_md <-function(w_grid,f0,betas,main = NULL,ylab = "Coefficient function",xlab = "w",cex_value=1.5){
  ylim = range(f0,betas)+c(-0.25,0.25)
  md =  apply(betas, 1, median)
  plot(w_grid,f0,main = main,ylim = ylim,type = "l",ylab = ylab,xlab = xlab,lty = 1,lwd=3,col="grey",las=1,cex.main=cex_value,cex.axis=cex_value,cex.lab=cex_value) 
  
  for (j  in 1:ncol(betas)) {
    lines(w_grid,betas[,j],lty=1,lwd=0.5,col="grey80")
  }
  lines(w_grid,f0,lty = 1,lwd=2,col='red')
  lines(w_grid,md,lty = 1,lwd=2)
}

pdf(file = paste(out_folder,"/n",n,"B1STVmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f1(w_grid),STV_beta1,main = NULL,ylab = expression(b[1]),xlab = "w")
dev.off()

pdf(file = paste(out_folder,"/n",n,"B2STVmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f2(w_grid),STV_beta2,main = NULL,ylab = expression(b[2]),xlab = "w")
dev.off()

pdf(file = paste(out_folder,"/n",n,"B3STVmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f3(w_grid),STV_beta3,main = NULL,ylab = expression(b[3]),xlab = "w")
dev.off()



pdf(file = paste(out_folder,"/n",n,"B1regTVmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f1(w_grid),regTV_beta1,main = NULL,ylab = expression(b[1]),xlab = "w")
dev.off()

pdf(file = paste(out_folder,"/n",n,"B2regTVmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f2(w_grid),regTV_beta2,main = NULL,ylab = expression(b[2]),xlab = "w")
dev.off()

pdf(file = paste(out_folder,"/n",n,"B3regTVmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f3(w_grid),regTV_beta3,main = NULL,ylab = expression(b[3]),xlab = "w")
dev.off()


pdf(file = paste(out_folder,"/n",n,"B1locmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f1(w_grid),loc_beta1,main = NULL,ylab = expression(b[1]),xlab = "w")
dev.off()

pdf(file = paste(out_folder,"/n",n,"B2locmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f2(w_grid),loc_beta2,main = NULL,ylab = expression(b[2]),xlab = "w")
dev.off()

pdf(file = paste(out_folder,"/n",n,"B3locmd.pdf",sep = ""),width = 6,height = 6)
par(mar = c(4.0,4.6,2,2))
plot_md(w_grid,f3(w_grid),loc_beta3,main = NULL,ylab = expression(b[3]),xlab = "w")
dev.off()

