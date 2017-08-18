# Examine the EnKF's true & ensemble B.
#
# M Tsyrulnikov
# 17 Aug 2017

library(mixAK)
library(MASS)
library(stats)
library(plot3D)
library(psych)

source('functions_prior.R')
source('mx_fft.R')
source('fft_Rspe_rearrange.R')
source('mx_fft_Rspe_rearrange.R')
source('create_LocMx.R')
source('Cov2VarCor.R')
source('ndiag.R')


file=paste0("./filter.RData") 
load(file) # contains filter$...
#                            parameters
#                            xi=xi_sparse
#                            enkf_one_world$
#                                           X_f, X_a, B_arr, S_nonloc_arr
#                            S_enkf_nonloc_arr_mean (Sigma)
#                            Cov_mat_enkf
#                            rho,nu,U,sigma

n=    dim(filter$Cov_mat_enkf)[1]
ntime_full=dim(filter$Cov_mat_enkf)[3]  # nu of the asml cycles from the first cycle

seed_for_secondary_fields=filter$parameters$seed_for_secondary_fields
seed_for_filters=filter$parameters$seed_for_filters
seed_for_secondary_fields
seed_for_filters

# skip initial transient in the time series

ini=1 #  >0 !,  # 30
ini1= ini +1
ntime=ntime_full - ini

B=filter$Cov_mat_enkf[,,ini1:ntime_full]  # true spatial fcst-err cvms for all t
#S=enkf_B_arr_mean[,,ini+1:ntime]
S=filter$enkf_S_nonloc_arr_mean[,,ini1:ntime_full]  # ==Sigma  cvms for all t

#!!!!!!!! tmp
#B=S
#!!!!!!!!

S1=filter$enkf_one_world$S_nonloc_arr[,,ini1:ntime_full] # [ix,iy,it]  one-world S cvms for all t


N=filter$parameters$N
M=filter$parameters$M

rem=6.37*10^6 # m
rekm=rem/10^3 # km

mesh <- 2*pi/n  # rad
mesh_km <- mesh * rekm # grid spacing, km
grid_km      <- c(0:(n-1))*mesh_km
step_h=filter$parameters$delta_t / 3600
tgrid_h=c(1:ntime)*step_h
tgrid=c(1:ntime) # time steps
dt_mdl_h=filter$parameters$delta_t /3600  # model time integration step, h
Ta_h=dt_mdl_h * filter$parameters$stride # asml cycle

V <- matrix(0, nrow=n, ncol=ntime)  # FG error variances 

for (i in (1:n)){
  V[i,]=B[i,i,]
}

# replace zero variances w epsilon

eps=0.001
ind_small=which(V < eps, arr.ind=TRUE)
V[ind_small]=eps

std_st = sqrt(V)                    # FG error st.devs.

macroscale_st_km  <- matrix(0, nrow=n, ncol=ntime)
microscale_st_km  <- matrix(0, nrow=n, ncol=ntime)

CRM=B  # initialize

for (t in (1:ntime)){
  CRM[,,t]=diag(1/std_st[,t]) %*% B[,,t] %*% diag(1/std_st[,t])  # Crl Mx
  for (ii in (1:n)){
    macroscale_st_km[ii, t]  = sum(CRM[ii,,t])*mesh_km/2
    iim1=ii-1
    if(ii==1) iim1=n
    iip1=ii+1
    if(ii==n) iip1=1
    microscale_st_km[ii, t] = 1/sqrt( (-CRM[ii,iim1,t] + 2*CRM[ii,ii,t] - CRM[ii,iip1,t])/mesh^2 ) * rekm
  }
}  

Lambda=macroscale_st_km
lambda=microscale_st_km

# S

V_S <-  matrix(0, nrow=n, ncol=ntime)  # FG <S> error variances 
V_S1 <- matrix(0, nrow=n, ncol=ntime)  # FG S in one world: error variances 

macroscale_S_st_km  <- matrix(0, nrow=n, ncol=ntime)
microscale_S_st_km  <- matrix(0, nrow=n, ncol=ntime)

for (i in (1:n)){
  V_S[i,] =S[i,i,]
  V_S1[i,]=S1[i,i,]
}

ind_small=which(V_S < eps, arr.ind=TRUE)
V_S[ind_small]=eps

std_st_S = sqrt(V_S)

CRM_S=S  # initialize

for (t in (1:ntime)){
  CRM_S[,,t]=diag(1/std_st_S[,t]) %*% S[,,t] %*% diag(1/std_st_S[,t])  # Crl Mx
  for (ii in (1:n)){
    macroscale_S_st_km[ii, t]  = sum(CRM_S[ii,,t])*mesh_km/2
    iim1=ii-1
    if(ii==1) iim1=n
    iip1=ii+1
    if(ii==n) iip1=1
    microscale_S_st_km[ii, t] = 1/sqrt( (-CRM_S[ii,iim1,t] + 2*CRM_S[ii,ii,t] - CRM_S[ii,iip1,t])/mesh^2 ) * rekm
  }
} 

Lambda_S=macroscale_S_st_km
lambda_S=microscale_S_st_km

plot(V[1,1:(ntime/2)], main="V (circles) & V_S")
lines(V_S[1,1:(ntime/2)])


#--------------------
# Selected i,t stats

numplot=3
for (iplot in c(1:numplot)) {
  
  nt_loc=15
  tt <- sample(c(1:ntime), nt_loc, replace = FALSE, prob = NULL) # select random time instant
  t=tt[1]
  
  # Diagnostics at the selected t
  
  image2D(B[,,t], main=paste("B. t=",t))
  image2D(CRM[,,t], main=paste("CRM. t=",t))
  
  ii <- sample(c((n/4):(3*n/4)), nt_loc, replace = FALSE, prob = NULL)  # select random spatial point
  i=ii[1]
  
  mx <- max(B[i,(i-14):(i+14),t])
  plot(grid_km[(i-14):(i+14)]-grid_km[i], B[i,(i-14):(i+14),t], 
       type="l", ylim=c(-mx,mx), main=paste("B. t=",t," i=",i))
  
  mx <- max(CRM[i,(i-14):(i+14),t])
  plot(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], 
       type="p", ylim=c(-mx,mx), main=paste("CRM. t=",t," i=",i))
  
  
  pngname=paste0("SpaCrl_t",t,"i",i,".png")
  png(pngname, width=7.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(grid_km[(i-14):(i+14)]-grid_km[i],    CRM[i,(i-14):(i+14),1], 
       xlab="Distance, km", ylab="Correlation", 
       type="l", ylim=c(0,mx), 
       cex.main=1.7, cex.axis=1.3, cex.lab=1.6,
       main=paste("Spatial fcst-err correl") )
  i=ii[2]
  t=tt[2]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[3]
  t=tt[3]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[4]
  t=tt[4]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[5]
  t=tt[5]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[6]
  t=tt[6]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[7]
  t=tt[7]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[8]
  t=tt[8]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[9]
  t=tt[9]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[10]
  t=tt[10]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[11]
  t=tt[11]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[12]
  t=tt[12]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[13]
  t=tt[13]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[14]
  t=tt[14]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  i=ii[15]
  t=tt[15]
  lines(grid_km[(i-14):(i+14)]-grid_km[i], CRM[i,(i-14):(i+14),t], type="l")
  abline(h=0)
  dev.off()

} # End for(iplot in c(1:numplot))  
  
#------------------------------------------------
#------------------------------------------------
# True vs. ensemble covs: 
# Plots: V, Lambda, lambda

plots=FALSE
if(plots){

pngname=paste0("V.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(V, x=grid_km, y=tgrid, xlab="Space, km", ylab="Time, h",
        main=paste('True forecast error vaiances V'),
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("VS.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(V_S, x=grid_km, y=tgrid, xlab="Space, km", ylab="Time, h",
        main=paste('Ensemble forecast error vaiances V_S'),
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("dV_f.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(V_S-V, x=grid_km, y=tgrid, xlab="Space, km", ylab="Time, h",
        main=paste('Errors V_S - V'),
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


pngname=paste0("Lambda.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(Lambda, x=grid_km, y=tgrid, xlab="Space, km", ylab="Time, h",
        main=paste("Local macroscale Lambda"), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


pngname=paste0("Lambda_mic.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(lambda, x=grid_km, y=tgrid, xlab="Space, km", ylab="Time, h",
        main=paste("Local microscale lambda"), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

#------------------------------------------------
# Fit the SDE SPPL order-1 model (switched off)


fit_SPPL=FALSE

                 if(fit_SPPL){
dim_phy_space = 1
P_order=1
nsearch1=41 # exhaust search 

crf_tru = CRM[1,,t]

# symmetrization is needed

for (i in c( 2:((n/2)+1) )) {
  i_symm=n+2-i
  x=(crf_tru[i] + crf_tru[i_symm]) /2
  crf_tru[i_symm]=x
  crf_tru[i]=     x
  
}

PMSDE <- FitHomogenSPPL(dim_phy_space, n, crf_tru, P_order, nsearch1)
crf_mod=PMSDE$crf_mod

nx=20

pngname=paste0("enkf_sppl_approx",t,".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(crf_tru[1:nx], x=grid_km[1:nx], xlab="Space, km", ylab="Correlation",
     main=paste("EnKF crf (circles) vs SPPL model"),
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(crf_mod[1:nx], x=grid_km[1:nx], col="red")
dev.off()

                 } # end if(fit_SPPL)

#------------------------------------------------
# Compare B=B_true and S=S_worldmean

nx=14

numplot=1
for (iplot in c(1:numplot)) {
  
  nt_loc=15
  tt <- sample(c(1:ntime), nt_loc, replace = FALSE, prob = NULL) # select random time instant
  t=tt[1]
  
  # Diagnostics at the selected t
  
  image2D(B[,,t], main=paste("B. t=",t))
  image2D(CRM[,,t], main=paste("CRM. t=",t))
  
  ii <- sample(c((n/4):(3*n/4)), nt_loc, replace = FALSE, prob = NULL)  # select random spatial point
  i=ii[1]
  
  mx <- max( max(B[i,(i-nx):(i+nx),t]) , max(S[i,(i-nx):(i+nx),t]) )
  mn <- min( min(B[i,(i-nx):(i+nx),t]) , min(S[i,(i-nx):(i+nx),t]) )
  
  pngname=paste0("B_and_S",t,".png")
  png(pngname, width=7.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(grid_km[(i-nx):(i+nx)]-grid_km[i], B[i,(i-nx):(i+nx),t],  ylim=c(mn,mx),
       xlab="Space, km", ylab="Covariance", type="p", pch=1, col="blue",
       main=paste("Prior covariances"),
       cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  lines(grid_km[(i-nx):(i+nx)]-grid_km[i], S[i,(i-nx):(i+nx),t], 
        type="l", lwd=1)
  
  leg.txt=c( 'B', as.expression( bquote(Sigma) ) )
  leg.col<-c("blue", "black")
  legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(NA,1), 
         lty=c(NA,1), pt.lwd=2.0, pch=c(1,NA),
         cex=1.4, pt.cex=1, bg="white")
  dev.off()
  
  
  mx <- max( max(CRM[i,(i-nx):(i+nx),t]) , max(CRM_S[i,(i-nx):(i+nx),t]) )
  mn <- min( min(CRM[i,(i-nx):(i+nx),t]) , min(CRM_S[i,(i-nx):(i+nx),t]) )
  
  pngname=paste0("B_and_S_correl",t,".png")
  png(pngname, width=7.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(grid_km[(i-nx):(i+nx)]-grid_km[i], CRM[i,(i-nx):(i+nx),t], ylim=c(mn,mx),
       xlab="Space, km", ylab="Correlation",type="p", pch=1, col="blue",
       main=paste("Prior correlations"),
       cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  lines(grid_km[(i-nx):(i+nx)]-grid_km[i], CRM_S[i,(i-nx):(i+nx),t], 
        type="l", col="black", lwd=1)
  
  leg.txt=c( 'B', as.expression( bquote(Sigma) ) )
  leg.col<-c("blue", "black")
  legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(NA,1), 
         lty=c(NA,1), pt.lwd=2.0, pch=c(1,NA),
         cex=1.4, pt.cex=1, bg="white")
  dev.off()
}  # End for (numplot ... N2

#------------------------------------------------
# Len scales

mx_mac <- max(macroscale_st_km[,t])
mx_std=max(std_st[,t])

plot (macroscale_st_km[,t]/mx_mac, ylim=c(0,1.1), type="l", col="green",
      main=("GREEN: Rel. Len.scale | RED: Rel. std"))
lines(std_st[,t]/mx_std, type="l", col="red")


plot (microscale_st_km[,t]/macroscale_st_km[,t], ylim=c(0.3,1.5), type="l", col="green",
      main=("Rel. curv. radius at the origin"))

#------------------------------------------------

ind_max=which(Lambda == max(Lambda), arr.ind=TRUE)
plot(B[,ind_max[1],ind_max[2]], main="B at point of max Lambda")

ind_min=which(Lambda == min(Lambda), arr.ind=TRUE)
plot(B[,ind_min[1],ind_min[2]], main="B at point of min Lambda")

#------------------------------------------------
# cross-fld crl

crl_Lambda_V=cor(as.vector(Lambda), as.vector(V), method="spearman")
crl_Lambda_V
crl_lambda_V=cor(as.vector(lambda), as.vector(V), method="spearman")
crl_lambda_V

micro_d_macro_scale_mean=mean(microscale_st_km[,10] / Lambda[,10])  # 10: to avoid t=0...9, where B=0
micro_d_macro_scale_mean

crl_macroS_macroB=cor(as.vector(macroscale_S_st_km), as.vector(macroscale_st_km), method="spearman")
crl_macroS_macroB

crl_microS_microB=cor(as.vector(microscale_S_st_km), as.vector(microscale_st_km), method="spearman")
crl_microS_microB

mean_macro=mean(macroscale_st_km)
mean_macro
bias_macro=mean(macroscale_S_st_km - macroscale_st_km)
bias_macro

mean_micro=mean(microscale_st_km)
mean_micro
bias_micro=mean(microscale_S_st_km - microscale_st_km)
bias_micro

#--------------------------------------------------------
# Scatterplots
# 
# lambda bias

est_err <- c()
true_val <- c()
for(step in 1:ntime){
  est_err  <- c(est_err,  microscale_S_st_km[,step]-microscale_st_km[,step])
  true_val <- c(true_val,                           microscale_st_km[,step])
}
ymean_err <- mean(est_err)
xmean_err <- mean(true_val)
ym_e <- round(ymean_err,0)
xm_e <- round(xmean_err,0)

pngname=paste0("Scatt_microscale_bias.png")
png(pngname, width=7.1, height=7.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = expression(lambda[true]), 
     ylab = bquote( paste("<",lambda, "> - ",lambda[true]) ), 
     main = bquote(paste('Bias ',lambda, ', ymean=', .(ym_e), ' xmean=', .(xm_e))),
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
dev.off()

#---------------------------------------------
# Lambda bias

est_err <- c()
true_val <- c()
for(step in 1:ntime){
  est_err  <- c(est_err,  macroscale_S_st_km[,step]-macroscale_st_km[,step])
  true_val <- c(true_val,                           macroscale_st_km[,step])
}
ymean_err <- mean(est_err)
xmean_err <- mean(true_val)
ym_e <- round(ymean_err,0)
xm_e <- round(xmean_err,0)

pngname=paste0("Scatt_macroscale_bias.png")
png(pngname, width=7.1, height=7.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = expression(Lambda[true]), 
     ylab = bquote( paste("<",Lambda, "> - ",Lambda[true]) ), 
     main = bquote(paste('Bias ',Lambda, ', ymean=', .(ym_e), ' xmean=', .(xm_e))), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
#lines(lowess(true_val, est_err, f = 0.3), col = 'lightcoral', lwd=3)
dev.off()


#---------------------------------------------
# Var bias as a fu of B

est_err <- c()
true_val <- c()
for(step in 1:ntime){
  est_err <- c(est_err,   V_S[,step]-V[,step])
  true_val <- c(true_val,            V[,step])
}
ymean_err <- mean(est_err)
xmean_err <- mean(true_val)
ym_e <- round(ymean_err,2)
xm_e <- round(xmean_err,1)

pngname=paste0("Scatt_S_bias.png")
png(pngname, width=7.1, height=7.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, 
     xlab = bquote(paste('B[ii](t)')), 
     ylab = bquote(paste(Sigma[ii](t) - B[ii](t))), 
     main = paste("Bias in sample variances"),
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
dev.off()


#---------------------------------------------
# Crl bias as a fu of CRM_true
# as a fu of shift in space, grid spacings

n_ds_max =10

crf_mean     = c(1:n_ds_max) 
crf_S_mean   = c(1:n_ds_max)

crf_mean[]=0
crf_mean[1]=1

crf_S_mean[]=0
crf_S_mean[1]=1

CRL_ds   = matrix(0, nrow=n, ncol=ntime)  # CRL, shift=ds
CRL_S_ds = matrix(0, nrow=n, ncol=ntime) 

for (ds in 2:n_ds_max) {
  
  for (i in (1:n)){
    ip=i+ds
    if(ip > n) {ip=ip-n}
    CRL_ds[i,]=  CRM[i,ip,]
    CRL_S_ds[i,]=CRM_S[i,ip,]
    
    crf_mean[ds]=  crf_mean[ds]   + mean(CRM[i,ip,])
    crf_S_mean[ds]=crf_S_mean[ds] + mean(CRM_S[i,ip,])
  }
  
  crf_mean[ds]=crf_mean[ds]/n
  crf_S_mean[ds]=crf_S_mean[ds]/n
  
  est_err <- c()
  true_val <- c()
  for(step in 1:ntime){
    est_err <- c(est_err,   CRL_S_ds[,step]-CRL_ds[,step])
    true_val <- c(true_val,                 CRL_ds[,step])
  }
  ymean_err <- mean(est_err)
  xmean_err <- mean(true_val)
  ym_e <- round(ymean_err,2)
  xm_e <- round(xmean_err,1)
  
  pngname=paste0("Scatt_CRLshift",ds,"_bias.png")
  png(pngname, width=7.1, height=7.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(true_val, est_err, 
       xlab = paste("True spatial correlations"), 
       ylab = paste("Correlation biases"), 
       main = paste("Bias in sample correl, distance=", 
                    round(ds*mesh_km,-1),'km'),
       cex.main=1.5, cex.axis=1.2, cex.lab=1.4)
  lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
  dev.off()
  
}

#---------------------------------------------
# Time and space mean correlation

mn <- min( min(crf_S_mean[1:n_ds_max]) , min(crf_S_mean[1:n_ds_max]) )

pngname=paste0("Mean_correl.png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(grid_km[1:n_ds_max]-grid_km[1], crf_mean[1:n_ds_max], 
     ylim=c(mn,1), 
     xlab="Distance, km", 
     ylab="Correlation", type="p", pch=1, col="blue", 
     main=paste("Mean prior correlations"),
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(grid_km[1:n_ds_max]-grid_km[1], crf_S_mean[1:n_ds_max], 
      type="l", lty=1, lwd=1)

leg.txt=c( 'B', as.expression( bquote(Sigma) ) )
leg.col<-c("blue", "black")
legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(NA,1), 
       lty=c(NA,1), pt.lwd=2.0, pch=c(1,NA),
       cex=1.3, pt.cex=1, bg="white")
dev.off()

#--------------------------------------------
# Var bias time ser

nt=min(ntime,50)
ix=41

image2D(V)
image2D(V_S)
image2D(V_S-V)

it1=1
it2=nt

plot(V[ix,it1:it2])
lines(V_S[ix,it1:it2])

plot(V_S[ix,it1:it2]-V[ix,it1:it2])
lines(-V[ix,it1:it2]/20)

plot((V_S[ix,it1:it2]-V[ix,it1:it2])/V[ix,it1:it2], type="p", main="Rel. var bias tim ser")
   
ind_max=which(V == max(V), arr.ind=TRUE)

ix=ind_max[1]  
it_max=ind_max[2] 

it1=max(1, it_max - 10)
it2=min(ntime, it_max + 10)

plot(V_S[ix,it1:it2]-V[ix,it1:it2], main="Abs. var bias tim ser")
  
ix1=max(1, ix - 10)
ix2=min(n, ix + 10)

image2D(V[ix1:ix2, it1:it2], main="V near the point of max V")
image2D(V_S[ix1:ix2, it1:it2], main="V_S near the point of max V")
image2D(V_S[ix1:ix2, it1:it2] - V[ix1:ix2, it1:it2], main="V_S-V near the point of max V")  
plot(apply(V, 1, mean), type="l", main="Time-mean V")
plot(apply(V, 2, mean), type="l", main="Space-mean V")

#--------------------------------------------
# Max CRM_S error

cLoc=filter$parameters$L_loc /1000  # km
cLoc_rad=cLoc / rekm
cLoc_spacings=cLoc_rad *n/(2*pi)

LocMx_ = create_LocMx(n, cLoc_spacings, "GC5")

CRM_loc=CRM # init
CRM_S_loc=CRM

for (t in (1:ntime)){
  CRM_loc[,,t]   = CRM[,,t]   * LocMx_
  CRM_S_loc[,,t] = CRM_S[,,t] * LocMx_
}

ind_max=which(abs(CRM_S_loc-CRM_loc) == max(abs(CRM_S_loc-CRM_loc)), arr.ind=TRUE)
ix1=ind_max[1]
ix2=ind_max[2]
it=ind_max[3]

plot(CRM[ix1,,it], main="CRM, CRM_S (max dCRM_loc)")
lines(CRM_S[ix1,,it])

# look at superdiagonals

p=8  # nu of superdiagonals

dCRM_sdiag=array(0, dim=c(n,p,ntime))

for (t in (1:ntime)){
  dCRM_sdiag[,,t]=
    extract_cyclic_superdiagonals(CRM_S[,,t]-CRM[,,t], p)$superdiagonals
}

max(abs(dCRM_sdiag[,,]))
ind_max=which(abs(dCRM_sdiag) == max(abs(dCRM_sdiag)), arr.ind=TRUE)
ix1=ind_max[1]
ix2=ind_max[2]
it=ind_max[3]
dCRM_sdiag[ix1,ix2,it]

plot(CRM[ix1,,it], main="CRM, CRM_S (max dCRM p sdiags)")
lines(CRM_S[ix1,,it])


ind_max=which(abs(lambda-lambda_S) == max(abs(lambda-lambda_S)), arr.ind=TRUE)
ix=ind_max[1]
it=ind_max[2]

plot(CRM[ix,,it], main="CRM, CRM_S (max dlambda)")
lines(CRM_S[ix,,it])


ind_max=which(abs(Lambda-Lambda_S) == max(abs(Lambda-Lambda_S)), arr.ind=TRUE)
ix=ind_max[1]
it=ind_max[2]

plot(CRM[ix,,it], main="CRM, CRM_S (max dLambda)")
lines(CRM_S[ix,,it])


} # end if(plots) (BIG if)

#------------------------------------------------
#------------------------------------------------
# Compare space-time-mean rms errors in Xf & Xa (one world)

sdXf = sqrt( sum(filter$enkf_one_world$X_f[,ini1:ntime_full]^2) / (n*ntime))
sdXf

sdXamXf = sqrt( sum((filter$enkf_one_world$X_a[,ini1:ntime_full] - 
                     filter$enkf_one_world$X_f[,ini1:ntime_full])^2) / (n*ntime))
sdXamXf

X_true=filter$xi_sparse
sdXfmX = sqrt( sum((filter$enkf_one_world$X_f[,ini1:ntime_full] - X_true[,ini1:ntime_full])^2) / (n*ntime))
sdXamX = sqrt( sum((filter$enkf_one_world$X_a[,ini1:ntime_full] - X_true[,ini1:ntime_full])^2) / (n*ntime))

sdXfmX
sdXamX


#------------------------------------------------
#------------------------------------------------
# Fourier space calcul

# Time mean B

B_clim=matrix(0, nrow=n, ncol=n)

for (i in 1:n) {
  for (j in 1:n){
    B_clim[i,j]=mean(B[i,j,])
  }
}

# Lclz radius

mult_loc=1
if(N == 10){mult_loc=1.0}
if(N == 20){mult_loc=1.2}
 
cLoc=filter$parameters$L_loc /1000 * mult_loc # km
cLoc_rad=cLoc / rekm
cLoc_spacings=cLoc_rad *n/(2*pi)

LocMx = create_LocMx(n, cLoc_spacings, "GC5")

# Mx FFT

FBclim=mx_fft(B_clim, "f")
#tmp=mx_fft(FBclim, "b")
#max(abs(tmp - B_clim)) # test OK.

FBclim_CRM = Cov2VarCor(FBclim)$C
FBclim_var = Cov2VarCor(FBclim)$v

# How diagonal is FBclim?

ndiag_FBclim=ndiag(FBclim)
ndiag_FBclim

plot(abs(FBclim[1,]), main="abs(FBclim[1,])")
plot(abs(FBclim[10,]), main="abs(FBclim[10,])")
plot(abs(FBclim[30,]), main="abs(FBclim[30,])")

# ==> quite diagonal..

# neglect non-diagonal entries of FBclim 
#  (as B_clim should be homogeneous in this model)

FBclim = diag(diag(Re(FBclim)))

FB       = array(0+0i, dim=c(n,n,ntime)) # F(B)
FFB      = array(0+0i, dim=c(n,n,ntime)) # F(B)

FS1loc   = array(0+0i, dim=c(n,n,ntime)) # FFT(S1_loc)
FBclim_err = array(0+0i, dim=c(n,n,ntime)) # FFT(B_clim-B)
FS1loc_err = array(0+0i, dim=c(n,n,ntime)) # FFT(S1_loc-B)

maFB      =matrix(0, nrow=n, ncol=n)
mae_FBclim=matrix(0, nrow=n, ncol=n)
mae_FS1loc=matrix(0, nrow=n, ncol=n)
dS1loc=c(1:ntime)

for (t in 1:ntime){
  S1_loc=S1[,,t] *LocMx
  FB      [,,t] = mx_fft(B[,,t], "f")  # \hat B =0
  FBclim_err[,,t] = FBclim - FB[,,t] # mx_fft(B_clim - B[,,t], "f")
  FS1loc[,,t] = mx_fft(S1_loc, "f")
  FS1loc_err[,,t] = FS1loc[,,t] - FB[,,t] # mx_fft(S1_loc  - B[,,t], "f")

  maFB      =maFB       + abs(FB      [,,t])
  mae_FBclim=mae_FBclim + abs(FBclim_err[,,t])
  mae_FS1loc=mae_FS1loc + abs(FS1loc_err[,,t])
  
  dS1loc[t] = norm(abs(S1_loc  - B[,,t]), "F") / 
              norm(abs(          B[,,t]), "F")
}

rmse_S1loc_rel  = mean(dS1loc)
rmse_S1loc_rel

t=ntime/2.5
tmp=Re(mx_fft_Rspe_rearrange(FB[,,t]))
image2D(tmp, main="FB")

tmp=Re(mx_fft_Rspe_rearrange(FBclim))
image2D(tmp, main="FBclim")

tmp=Re(mx_fft_Rspe_rearrange(FS1loc[,,t]))
image2D(tmp, main="FS1loc")

tmp=Re(mx_fft_Rspe_rearrange(FBclim_err[,,t]))
image2D(tmp, main="FBclim_err")

tmp=Re(mx_fft_Rspe_rearrange(FS1loc_err[,,t]))
image2D(tmp, main="FS1loc_err")

maFB      =maFB       /ntime
mae_FBclim=mae_FBclim /ntime
mae_FS1loc=mae_FS1loc /ntime

dFa=mae_FS1loc - mae_FBclim
dFB=maFB - mae_FBclim

mean(dFB)
max(dFB)
min(dFB)

# ==> dFB>0 ==> Bclim is always better than 0.


wvns=c((-n/2):(n/2))
wvns_posit=c(0:(n/2))

          FBclimPlots=FALSE
          if(FBclimPlots){

pngname=paste0("dFB",mult_loc,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(Re(mx_fft_Rspe_rearrange(dFB)),
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="<0: 0 bett,  >0: Bclim bett",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("dFa_loc",mult_loc,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(Re(mx_fft_Rspe_rearrange(dFa)),
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="<0: S bett,  >0: Bclim bett",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("maeFBclim",mult_loc,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(Re(mx_fft_Rspe_rearrange(mae_FBclim)),
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="Mean abs err FBclim",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("maeFS1loc",mult_loc,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(Re(mx_fft_Rspe_rearrange(mae_FS1loc)),
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="Mean abs err FS1loc",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

dFa_sign=dFa
dFa_sign[,]=0
mx=max(abs(dFa))
d=0 # mx/20

ind_Sworse=which(dFa > d, arr.ind = TRUE)
dFa_sign[ind_Sworse]=1

ind_Sbett=which(dFa< -d, arr.ind = TRUE)
dFa_sign[ind_Sbett]=-1

pngname=paste0("dFa_sign.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(dFa_sign, xlab="Wavenumber", ylab="Wavenumber",
        x=wvns, y=wvns,
        main=paste("1:CLIM better,-1:S better"),
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

sum(dFa[ind_Sbett])
sum(dFa[ind_Sworse])

                } # End if(FBclimPlots)

#-----------------------------------------------
#-----------------------------------------------
# Compute Spectral Variances.
# NB after mx_fft_Rspe_rearrange, the wvns go as follows:
# f(n/2),f(-n/2 +1), ..., f(-1), f(0), f(1),..., f(n/2)
# so that m=0 is at the (n/2+1)st place in the (n+1)-vector
# (for the even n).
#
# NB: (1) Spectral Variances are purely real (for both real & cplx 
#         phy-space processes).
#     (2) v(-m)=v(m) (for a real phy-space process).
#--
# NB:
# im=n/2 + m +1
# m=im - n/2 -1

np1=n+1
nd2=n/2
nd2p1=nd2 +1

# NB: In what follows, the v-matrices/vectors correspond to only non-negative wvns

v_true =matrix(0, nrow=nd2p1, ncol=ntime)  # v(m,t)
v_S1loc=matrix(0, nrow=nd2p1, ncol=ntime)

tmp=diag( Re(FBclim) ) 
v_clim=tmp[1:nd2p1]

for (t in 1:ntime){
  tmp=diag( Re(FB[,,t]) ) 
  v_true[,t]=tmp[1:nd2p1]
  
  tmp=diag( Re(FS1loc[,,t]) ) 
  v_S1loc[,t]=tmp[1:nd2p1]
}

v_S1loc_full =matrix(0, nrow=n, ncol=ntime)
v_true_full  =matrix(0, nrow=n, ncol=ntime)

for (t in 1:ntime){
  v_S1loc_full[,t] =c(v_S1loc [,t], rev(v_S1loc[2:nd2,t]))
  v_true_full[,t]  =c(v_true[,t],   rev(v_true [2:nd2,t]))
}

#-----------------------------------------------
#-----------------------------------------------
# Examine Spectral Correlations

FCRM_e  =FS1loc # init
FCRM_t  =FB     # init
FFCRM_t =FB     # init

for (t in 1:ntime) {
  std=sqrt(diag(Re(FB[,,t])))      # NB: spe VARs are real
  FCRM_t[,,t]= diag(1/std) %*% FB[,,t] %*% diag(1/std)
  
  std=sqrt(diag(Re(FS1loc[,,t])))  # NB: spe VARs are real
  FCRM_e[,,t]= diag(1/std) %*% FS1loc[,,t] %*% diag(1/std)
  
  FFCRM_t[,,t]=mx_fft(FCRM_t[,,t], "b")
}

t=ntime/3.5

image2D(Re(FCRM_t[,,t]) - diag(n), main="Re(FCRM_t[,,t]) - I")
image2D(Im(FCRM_t[,,t]), main="Im(FCRM_t[,,t])")

image2D(Re(FFCRM_t[,,t]), main="Re(FFCRM_t[,,t])")  # ~diag!
#image2D(Im(FFCRM_t[,,t]), main="Im(FFCRM_t[,,t])")  # =0!
#max(abs(Im(FFCRM_t[,,t])))
# ==> FFCRM_t is real =>

FFCRM_t=Re(FFCRM_t)

i=n/3.4
plot(FFCRM_t[i,,t])
plot(B[i,,t])

norm( Re(FCRM_t[,,t])- diag(n), "F" )
norm( Im(FCRM_t[,,t]), "F" )

D = diag( diag(FFCRM_t[,,t]) ) # leave only the diagonal in FFCRM_t
norm(                D, "F" )
norm( FFCRM_t[,,t],     "F" )
norm( FFCRM_t[,,t] - D, "F" )

norm( Re(FCRM_e[,,t] - FCRM_t[,,t]), "F" )
norm( Im(FCRM_e[,,t] - FCRM_t[,,t]), "F" )

image2D(Re(FCRM_e[,,t]) - diag(n), main="Re(FCRM_e[,,t]) - I")
image2D(Im(FCRM_e[,,t]), main="Im(FCRM_e[,,t])")

#-----------------------------------------------
#-----------------------------------------------
# Fit and evaluate SSCM 

 B_est      = array(0, dim=c(n,n,ntime))
FB_est      = array(0, dim=c(n,n,ntime))
C_tilde_mod = array(0, dim=c(n,n,ntime))

uu = matrix(0,nrow=n, ncol=ntime)

tru_norm=c(1:ntime)
dmod_norm=c(1:ntime)
rel_dmod=c(1:ntime)

for (t in 1:ntime) {
  D = mx_fft(FCRM_t[,,t], "b")
  
  # Fit the SSCM
  
  M=mx_fft(diag(v_true_full[,t]), "b") /n
  b=diag(M %*% D %*% M)
  N=M*t(M)
  u=solve(N, b)  # real
  uu[,t] = Re(u)
  
  U = diag(Re(u))
    
  C_tilde_mod[,,t] = mx_fft(U, "f")
  
  #image2D(Re(FCRM_t[,,t]))
  #image2D(Re(C_tilde_mod[,,t]))
  #image2D(Im(FCRM_t[,,t]))
  #image2D(Im(C_tilde_mod[,,t]))
  #plot(Re(u))
  #lines(Re(diag(FFCRM_t[,,t])))
  
  
  #FB_est[,,t] = diag(sqrt(v_opt_full[,t])) %*% C_tilde[,,t] %*% 
  #              diag(sqrt(v_opt_full[,t]))
  
  # v_true
  FB_est[,,t] = diag(sqrt(v_true_full[,t])) %*% C_tilde_mod[,,t] %*% 
                diag(sqrt(v_true_full[,t]))
  
  
  # transplant v_true & the diagonal C+: poor ==> diagonal C is not good..
  #D = diag( diag(mx_fft(FCRM_t[,,t], "b")) )
  #C_tilde[,,t] = mx_fft(D, "f")
  #FB_est[,,t] = diag(sqrt(v_true_full[,t])) %*% C_tilde[,,t] %*% 
  #              diag(sqrt(v_true_full[,t]))
  
  # take FCRM_t and v_S1loc:  error ~ 30%
  #FB_est[,,t] = diag(sqrt(v_S1loc_full[,t])) %*% FCRM_t[,,t] %*% 
  #              diag(sqrt(v_S1loc_full[,t]))
  
  # take FCRM_t and v_opt: small error ~15%
  #FB_est[,,t] = diag(sqrt(v_opt_full[,t])) %*% FCRM_t[,,t] %*% 
  #              diag(sqrt(v_opt_full[,t]))
  
  # take FCRM_t and v_true: OK, exact
  #FB_est[,,t] = diag(sqrt(v_true_full[,t])) %*% FCRM_t[,,t] %*% 
  #              diag(sqrt(v_true_full[,t]))
  
  B_est[,,t] = mx_fft(FB_est[,,t], "b")
  
  #image2D(B[,,t])
  #image2D(Re(B_est[,,t]))
  
  #tru_norm[t]  = ( norm(abs(               FB[,,t]), "F") )
  #dmod_norm[t] = ( norm(abs(FB_est[,,t]  - FB[,,t]), "F") )
  tru_norm[t]  = norm(abs(              B[,,t]), "F")
  dmod_norm[t] = norm(abs(B_est[,,t]  - B[,,t]), "F")
  rel_dmod[t] = dmod_norm[t] / tru_norm[t]
}

rms_tru  = mean(tru_norm)
rms_dmod = mean(dmod_norm)
rms_tru
rms_dmod
rel_misfit=rms_dmod / rms_tru
rel_misfit
mean_rel_misfit=mean(rel_dmod)
mean_rel_misfit

S1loc_rel_misfit=rmse_S1loc_rel

# Select the typical (as measured by the misfit) t, and plot

rrel=dmod_norm / tru_norm
t=which(abs(rrel[]-mean_rel_misfit) == min(abs(rrel[]-mean_rel_misfit)), arr.ind=TRUE)

pngname=paste0("B_tru_se",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(B[,,t], x=grid_km, y=grid_km, xlab="Space, km", ylab="Space, km",
        main=paste('Typical true B'),
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

B_mod=Re( B_est[,,t] )

pngname=paste0("B_mod_se",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(B_mod, x=grid_km, y=grid_km, xlab="Space, km", ylab="Space, km",
        main=paste('Typical model B'),
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

S_loc=S1[,,t] * LocMx

pngname=paste0("S_loc_se",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(S_loc, x=grid_km, y=grid_km, xlab="Space, km", ylab="Space, km",
        main=paste('Typical localized S'),
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


RR=Re(mx_fft_Rspe_rearrange(FB[,,t]))
pngname=paste0("Re_SpeCVM_tru_se",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(RR[2:np1,2:np1],   # avoid wvn=-n/2
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="Re(Spectral CVM)",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


RR=Re(mx_fft_Rspe_rearrange(FCRM_t[,,t]))
pngname=paste0("Re_SpeCRM_tru_se",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(RR[2:np1,2:np1],   # avoid wvn=-n/2
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="Re(Spectral CRM)",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


pngname=paste0("Re_SpeCRM-I_tru_se",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(RR[2:np1,2:np1] - diag(n),  # avoid wvn=-n/2
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="Re(Spectral CRM) - I",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

II=Im(mx_fft_Rspe_rearrange(FCRM_t[,,t]))
pngname=paste0("Im_SpeCRM_tru_se",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(II[2:np1,2:np1],
        x=wvns, y=wvns,
        xlab="Wavenumber", ylab="Wavenumber",
        main="Im(Spectral CRM)",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


pngname=paste0("uu",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(Re(uu),
        x=grid_km, y=tgrid_h,
        xlab="Time, h", ylab="Space, km",
        main="u-vectors",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


# Select the "typical" pair (t,i)  with the closest relative Quadratic ROW norms to 
# both B_est and S1_loc

rel_dmod_row=matrix(0, nrow=n, ncol=ntime)
rel_dS1_row =matrix(0, nrow=n, ncol=ntime)

for (t in 1:ntime){
  S1loc=S1[,,t]*LocMx
  for(i in 1:n){
    rel_dmod_row[i,t] = sqrt( sum( (Re(B_est[i,,t])  - B[i,,t])^2 ) / sum( B[i,,t]^2 ) )
    rel_dS1_row[i,t]  = sqrt( sum( (S1loc[i,]        - B[i,,t])^2 ) / sum( B[i,,t]^2 ) )
  }
}

it=which( abs(rel_dmod_row - mean_rel_misfit) + abs(rel_dS1_row - S1loc_rel_misfit) == 
      min(abs(rel_dmod_row - mean_rel_misfit) + abs(rel_dS1_row - S1loc_rel_misfit) ), 
      arr.ind=TRUE)
i=it[1]
t=it[2]
i
t
rel_dmod_row[i,t]
rel_dS1_row[i,t]

B_mod=Re( B_est[,,t] )
S_loc=S1[,,t] * LocMx

nx=10
ix1=i-nx
ix2=i+nx

# extend the rows up and down using periodicity

circ_km=2*pi*rekm
ggrid=c(grid_km[] - circ_km-grid_km[i], grid_km[]-grid_km[i], grid_km[] + circ_km-grid_km[i])
bb=c(B[i,,t], B[i,,t], B[i,,t])
bbm=c(B_mod[i,], B_mod[i,], B_mod[i,])
bbe=c(S_loc[i,], S_loc[i,], S_loc[i,])

ixx1=ix1+n
ixx2=ix2+n

mn <- min( min(bb[ixx1:ixx2]), min(bbm[ixx1:ixx2]), min(bbe[ixx1:ixx2]) )
mx <- max( max(bb[ixx1:ixx2]), max(bbm[ixx1:ixx2]), max(bbe[ixx1:ixx2]) )

pngname=paste0("Typic3Covs_se",seed_for_secondary_fields,".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(ggrid[ixx1:ixx2], bb[ixx1:ixx2],  ylim=c(mn,mx),
     xlab="Space, km", ylab="Covariance", type="p", pch=1, col="black",
     main=paste("Prior covariances"),
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(ggrid[ixx1:ixx2], bbm[ixx1:ixx2], 
      type="l", lwd=1, lty=1)
lines(ggrid[ixx1:ixx2], bbe[ixx1:ixx2], 
      type="l", lwd=1, lty=2)

leg.txt=c( 'True B',  'Model B', 'Ensemble B')
leg.col<-c("black", "black", "black")
legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(NA,1,1), 
       lty=c(NA,1,2), pt.lwd=2.0, pch=c(1,NA,NA),
       cex=1.4, pt.cex=1, bg="white")
dev.off()


#--------------------
# uu: 

u_tmean=apply(uu, 1, mean)
u_tsd  =apply(uu, 1, sd)

mx=max(u_tmean, u_tsd)
pngname=paste0("u_tmean_tsd",seed_for_secondary_fields,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(u_tmean, ylim=c(0,mx), type="l", lty=1,
        main="Time mean, sd u-vectors",
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(u_tsd, lty=2)

leg.txt=c( 'Time mean u',  'Time std u')
leg.col<-c("black", "black")
legend("bottom", inset=0, leg.txt, col=leg.col, lwd=c(1,1), 
       lty=c(1,2), pt.lwd=2.0, pch=c(NA,NA),
       cex=1.4, pt.cex=1, bg="white")
dev.off()


# Selected i,t stats

numplot=3
for (iplot in c(1:numplot)) {
  
  nt_loc=15
  tt <- sample(c(1:ntime), nt_loc, replace = FALSE, prob = NULL) # select random time instant
  t=tt[1]
  
  ii <- sample(c((n/4):(3*n/4)), nt_loc, replace = FALSE, prob = NULL)  # select random spatial point
  i=ii[1]
  
  mx=max(uu)
  
  pngname=paste0("u_",t,"i",i,".png")
  png(pngname, width=7.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(grid_km, uu[,t],
       xlab="Grid, km", ylab="u", 
       type="l", ylim=c(0,mx), 
       cex.main=1.7, cex.axis=1.3, cex.lab=1.6,
       main=paste("u-vectors at various t") )
  i=ii[2]
  t=tt[2]
  lines(grid_km, uu[,t], type="l", lty=2)
  i=ii[3]
  t=tt[3]
  lines(grid_km, uu[,t], type="l", lty=3)

  abline(h=0)
  dev.off()
  
} # End for(iplot in c(1:numplot))  

#------------------------------------------------
# Write key findings to the file

outfile=paste0("out_inpSeed",seed_for_secondary_fields,".txt")
unlink(outfile)
sink(outfile)

cat("max(V) / min(V)=")
print(max(V)/min(V))

cat("max(Lambda) / min(Lambda)=")
print(max(Lambda)/min(Lambda))


cat("max(lambda) / min(lambda)=")
print(max(lambda)/min(lambda))

cat("V_decile_ratio=")
V_decile_ratio=quantile(V[, 10:ntime], probs=0.9) / quantile(V[, 10:ntime], probs=0.1) # 10: to avoid ini transient
print(V_decile_ratio)

cat("Lambda_decile_ratio=")
Lambda_decile_ratio=quantile(Lambda[, 10:ntime], probs=0.9) / quantile(Lambda[, 10:ntime], probs=0.1)
print(Lambda_decile_ratio)

cat("lambda_decile_ratio=")
lambda_decile_ratio=quantile(lambda[, 10:ntime], probs=0.9) / quantile(lambda[, 10:ntime], probs=0.1)
print(lambda_decile_ratio)


cat("\n")

cat("SSCM_mean_rel_misfit=")
print(mean_rel_misfit)

cat("S1loc_rel_misfit=")
print(S1loc_rel_misfit)



sink()
#------------------------------------------------
