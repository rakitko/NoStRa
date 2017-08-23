# Examine the EnKF's true & ensemble B.
#
# M Tsyrulnikov
# 23 Aug 2017

library(mixAK)
library(MASS)
library(stats)
library(plot3D)
library(psych)

source('functions_prior.R')
source('est_crf_timser.R')
source('mx_fft.R')
source('fft_Rspe_rearrange.R')
source('mx_fft_Rspe_rearrange.R')
source('create_LocMx.R')
source('Diag_flt.R')
source('Cov2VarCor.R')
source('ndiag.R')


SelfClim = FALSE # TRUE FALSE  # if TRUE, WRITE the clim stats to the
# B_tmean_$seed_for_secondary_fields.RData file, 
# where  seed_for_secondary_fields is taken from the input filter.RData file.
# otherwise, READ several B_tmean matrices from the B_tmean_$seed_clim.RData files.
ExtClim = !SelfClim  


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

ini=30 #  >0 !,  # 30
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

#!!!!!
#B=CRM  # TMP
#!!!!!

#---------------------------------------------
# Lclz radius

mult_loc=1
#if(N == 10){mult_loc=1.0}
#if(N == 20){mult_loc=1.2}

cLoc=filter$parameters$L_loc /1000 * mult_loc # km
cLoc_rad=cLoc / rekm
cLoc_spacings=cLoc_rad *n/(2*pi)

LocMx = create_LocMx(n, cLoc_spacings, "GC5")


CRM_loc=CRM # init
CRM_S_loc=CRM

for (t in (1:ntime)){
  CRM_loc[,,t]   = CRM[,,t]   * LocMx
  CRM_S_loc[,,t] = CRM_S[,,t] * LocMx
}

#--------------------
# Selected i,t stats

SpaCrl_plots=FALSE
if(SpaCrl_plots){
  
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

}  # End if(SpaCrl_plots)

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
# Time mean B for current seed --- B_clim

B_tmean=matrix(0, nrow=n, ncol=n)

for (i in 1:n) {
  for (j in 1:n){
    B_tmean[i,j]=mean(B[i,j,])
  }
}

B_clim=B_tmean
B_tmean_loc=B_tmean

if(SelfClim){  # Write  B_tmean  to the file
  save(B_tmean, file=paste0("B_tmean_", seed_for_secondary_fields,".RData"))
  # use  B_clim=B_tmean  assigned above
  nclim=1
  
}else{
  # Read from several B_tmean files & average over them
  
  B_clim[,]=0
  nclim=0
  
  for(seed in 
      c(10500,
        10510,
        10530,
        10590,
        101520,
        101570,
        101610
              ) ){
    
    file=paste0("./B_tmean_", seed, ".RData") 
    load(file)
    B_clim=B_clim + B_tmean
    nclim=nclim +1
  }
  nclim
  B_clim = B_clim / nclim
}


plot(diag(B_clim), type="l", lwd=2, main=paste0("diag(B_clim, B_tmean_loc)", nclim))
lines(diag(B_tmean_loc), type="l")

image2D(B_tmean_loc, main="B_tmean_loc")  
image2D(B_clim, main="B_clim")  

# Spe space clim

FB_clim=mx_fft(B_clim, "f")
image2D(Re(FB_clim), main="Re(FB_clim)")
ndiag(FB_clim)
FB_clim_diag=diag(diag(FB_clim))  # the statio component

#------------------------------------------------
# Perform time filtering of S1 (a single-world realization) (non-lclz):
# T1(1)=S1(1),
# T1(k) = mu*T1(k-1) + (1-mu)*S1(k)
# mu=1-eps  =>  tau=1/eps (time scale in time steps)

mult_loc=1.0
#if(N == 10){mult_loc=1.0}
#if(N == 20){mult_loc=1.2}

cLoc=filter$parameters$L_loc /1000 * mult_loc # km
cLoc_rad=cLoc / rekm
cLoc_spacings=cLoc_rad *n/(2*pi)
LocMx = create_LocMx(n, cLoc_spacings, "GC5")


mult_loc_T=1.7
cLoc=filter$parameters$L_loc /1000 * mult_loc_T # km
cLoc_rad=cLoc / rekm
cLoc_spacings=cLoc_rad *n/(2*pi)

LocMx_T = create_LocMx(n, cLoc_spacings, "GC5")

if(Ta_h == 6)  w_EVS=0.5
if(Ta_h == 12) w_EVS=0.6

w_EVT=0.7
# 0.8: 74.35
# 0.75: 73.8
# 0.7: 73.7

#-------------------------------------
# FILTERING

#----------------
simple_FLT=TRUE  # TRUE  FALSE

if(simple_FLT == TRUE){
  
  if(Ta_h == 6)  eps=0.08
  if(Ta_h == 12) eps=0.35
  # Ta=6: eps=0.08 
  # Ta=12: eps=0.3-0.35
  # Ta=24: eps=0.45
  mu=1-eps
  T1=S1 # init
  
  for (t in 2:ntime) {
    T1[,,t]=mu*T1[,,t-1] + eps* S1[,,t]
  }
  
}

#----------------
# Use a spatial flt in the cov fcst.

spe_FLT=!simple_FLT
if(spe_FLT == TRUE){

  
if(Ta_h == 6)  eps=0.08
if(Ta_h == 12) eps=0.3
  # Ta=6: eps=0.08 
  # Ta=12: eps=0.3
  # Ta=24: eps=0.45 
mu=1-eps  
  
  # Define the spe-space filter
  # 
  # tranfu(m) = 1/(1 + (m/m0)^2),
  # 
  # whr m0 is the spe scale such that lambda=1/m0
  # 
  # R spe vector :  f0, f1, ..., f(n/2),f(-n/2 +1),..., f(-1)
   
tranfu=c(1:n) # the filter transfer fu
m_max=n/2

coef_lambda=0.4

flt_length=mesh * coef_lambda
m0=1/flt_length

tranfu[]=0

m=0
im=m+1
tranfu[im]=1 # m=0

m=m_max
im=m+1
tranfu[im]=1/(1 + (m/m0)^2)  # m=m_max

for (m in 1:(m_max-1)){
  im=m+1
  tranfu[im] = 1/(1 + (m/m0)^2)
  tranfu[n-m+1] = tranfu[im]
}
plot(tranfu, main="tranfu")
#tranfu

respfu=fft(tranfu, inverse=TRUE)
#max(Im(respfu))
plot(Re(respfu[1:10]), main="respfu[1:10]")

SFLT=diag(tranfu)  # diag mx

# To maintain the expectation of T, which is reduced at the wvn m
# by the factor of tranfu(n)^2,
# we add the respective B_clim wvn=m component multiplied by (1 - tranfu(m)^2),
# ie pre & post multiply FB_clim by diag(sqrt(1 - tranfu(:)^2))

SFLT_clim=diag(sqrt(1-tranfu^2))

T1=S1 # init
for (t in 2:ntime) {
  T=T1[,,t-1]
  fT=mx_fft(T, "f")
  fT_flt=SFLT %*% fT %*% SFLT
  T_flt=Re( mx_fft(fT_flt, "b") )
  T1[,,t]=mu*T_flt + eps* S1[,,t] #*LocMx
  
  FB_clim_flt = SFLT_clim %*% FB_clim %*% SFLT_clim
  B_clim_flt=Re( mx_fft(FB_clim_flt, "b") )
  T1[,,t]=mu*(T_flt + B_clim_flt) + eps* S1[,,t] #*LocMx
}

} # End if(spe_FLT
   
#-------------------------------------End Filtering
# Evaluate the results

F_B=c(1:ntime)
F_Bclim_mB=c(1:ntime)
F_S1loc_mB=c(1:ntime)
F_Tloc_mB=c(1:ntime)
F_EV_mB=c(1:ntime)
F_EVT_mB=c(1:ntime)

Tr_B=c(1:ntime)
Tr_climBmB=c(1:ntime)
Tr_S1mB=c(1:ntime)
Tr_EVmB=c(1:ntime)
Tr_T1mB=c(1:ntime)

for (t in 1:ntime) {
  # traditional hybrid EnVar
  S1loc=S1[,,t] * LocMx
  T1loc=T1[,,t] * LocMx_T
  B_EV =w_EVS * S1loc + (1-w_EVS)*B_clim
  B_EVT=w_EVT * T1loc + (1-w_EVT)*B_clim
  
  F_B[t]         =norm(         B[,,t], "F")
  F_Bclim_mB[t]  =norm(B_clim - B[,,t], "F")
  F_S1loc_mB[t]  =norm(S1loc  - B[,,t], "F")
  F_Tloc_mB[t]   =norm(T1loc  - B[,,t], "F")
  F_EV_mB[t]     =norm(B_EV   - B[,,t], "F")
  F_EVT_mB[t]    =norm(B_EVT  - B[,,t], "F")

  #Tr_B[t]      =sqrt( sum( (diag( (         B[,,t]) ))^2 ) )
  #Tr_climBmB[t]=sqrt( sum( (diag( (B_clim - B[,,t]) ))^2 ) )
  #Tr_S1mB[t]   =sqrt( sum( (diag( (S1loc  - B[,,t]) ))^2 ) )
  #Tr_EVmB[t]   =sqrt( sum( (diag( (B_EV   - B[,,t]) ))^2 ) )
  #Tr_T1mB[t]   =sqrt( sum( (diag( (T1loc  - B[,,t]) ))^2 ) )
  
  Tr_B[t]      =( sum( (diag( (         B[,,t]) )) ) ) # biases
  Tr_climBmB[t]=( sum( (diag( (B_clim - B[,,t]) )) ) )
  Tr_S1mB[t]   =( sum( (diag( (S1loc  - B[,,t]) )) ) )
  Tr_EVmB[t]   =( sum( (diag( (B_EV   - B[,,t]) )) ) )
  Tr_T1mB[t]   =( sum( (diag( (T1loc  - B[,,t]) )) ) )
}


t=141
S1loc=S1[,,t] * LocMx
B_EV=w_EVS * S1loc + (1-w_EVS)*B_clim
T1loc=T1[,,t]*LocMx_T
B_EVT=w_EVT * T1loc + (1-w_EVT)*B_clim
  
plot(diag(B[,,t]), main=paste0("B at t=", t))
lines(diag(B_clim))
lines(diag(S1[,,t] * LocMx), col="green")
lines(diag(B_EV), col="blue")
lines(diag(T1[,,t]*LocMx_T), col="red")

image2D(B_clim, main="B")
image2D(B[,,t], main="B")
image2D(B_EV, main="B_EV")
image2D(T1[,,t]*LocMx_T, main="T")


F_B_mean=mean(F_B)
F_Bclim_mB_mean=mean(F_Bclim_mB)
F_S1loc_mB_mean=mean(F_S1loc_mB)
F_Tloc_mB_mean=mean(F_Tloc_mB)
F_EV_mB_mean=mean(F_EV_mB)
F_EVT_mB_mean=mean(F_EVT_mB)

Tr_B_mean=mean(Tr_B)
Tr_Bclim_mB_mean=mean(Tr_climBmB)
Tr_S1loc_mB_mean=mean(Tr_S1mB)
Tr_EV_mB_mean=mean(Tr_EVmB)
Tr_T1_Tloc_mB_mean=mean(Tr_T1mB)

F_B_mean
F_Bclim_mB_mean
F_S1loc_mB_mean
F_Tloc_mB_mean
F_EVT_mB_mean
F_EV_mB_mean

Tr_B_mean
Tr_Bclim_mB_mean
Tr_S1loc_mB_mean
Tr_EV_mB_mean
Tr_T1_Tloc_mB_mean

#------------
# Plots

ntime_plot=ntime/2
mn=min(
  min(F_B         [1:ntime_plot]),
  min(F_S1loc_mB  [1:ntime_plot]),
  min(F_Tloc_mB   [1:ntime_plot]), 
  min(F_EV_mB     [1:ntime_plot])
) *1.0

mx=max(
  max(F_B         [1:ntime_plot]),
  max(F_S1loc_mB  [1:ntime_plot]),
  max(F_Tloc_mB[1:ntime_plot]), 
  max(F_EV_mB     [1:ntime_plot])
) *1.1

pngname=paste0("BEvST_Fnorm.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))

plot(F_B[1:ntime_plot], 
     ylim=c(mn,mx), 
     xlab="Time step", ylab="Frobenius Norm", type="p", pch=1, col="black",
     main=paste("Errors in covariances"),
     cex.main=1.5, cex.axis=1.2, cex.lab=1.4)

lines(F_EV_mB[1:ntime_plot], type="l", lwd=1.5, lty=3, col="black")
lines(F_S1loc_mB[1:ntime_plot], type="l", lwd=1, col="black")
lines(F_Tloc_mB[1:ntime_plot], type="l", lwd=2.5, col="black")

leg.txt=c(as.expression( bquote( "|B|" ) ), 
          as.expression( bquote( "|" ~ "B"[hybr] ~ "-B|" ) ) ,
          as.expression( bquote( "|" ~ "S"[loc] ~ "-B|" ) ), 
          as.expression( bquote( "|" ~ tilde("S") ~ "-B|" ) ) )
leg.col<-c("black", "black", "black", "black")
legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(NA,1.5,1,2.5), 
       lty=c(NA,3,1,1), pt.lwd=1.0, pch=c(1,NA,NA,NA),
       cex=1.0, pt.cex=1, bg="white")
dev.off()


ntime_plot=ntime/2
mn=min(
  min(F_S1loc_mB  [1:ntime_plot]),
  min(F_EVT_mB    [1:ntime_plot]), 
  min(F_EV_mB     [1:ntime_plot])
) *1.0

mx=max(
  max(F_S1loc_mB  [1:ntime_plot]),
  max(F_EVT_mB    [1:ntime_plot]), 
  max(F_EV_mB     [1:ntime_plot])
) *1

pngname=paste0("SEET_Fnorm.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))

plot(F_S1loc_mB[1:ntime_plot], 
     ylim=c(mn,mx), 
     xlab="Time step", ylab="Frobenius Norm",  type="l", lwd=1.5, lty=3, col="black",
     main=paste("Errors in covariances"),
     cex.main=1.5, cex.axis=1.2, cex.lab=1.4)

lines(F_EV_mB [1:ntime_plot], type="l", lwd=1, col="black")
lines(F_EVT_mB[1:ntime_plot], type="l", lwd=2.5, col="black")

leg.txt=c(
  as.expression( bquote( "|" ~ "S"[loc] ~ "-B|" ) ),
  as.expression( bquote( "|" ~ "B"[hybr] ~ "-B|" ) ) ,
  as.expression( bquote( "|" ~ "B"[hybr-T] ~ "-B|" ) ) 
)

leg.col<-c("black", "black", "black")
legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(1.5,1,2.5), 
       lty=c(3,1,1), pt.lwd=1.0, pch=c(NA,NA,NA),
       cex=1.0, pt.cex=1, bg="white")
dev.off()



ntime_plot=ntime/2
mn=min(
  min(F_Bclim_mB  [1:ntime_plot]), 
  min(F_S1loc_mB  [1:ntime_plot]),
  min(F_Tloc_mB   [1:ntime_plot])
) *1.0

mx=max(
  max(F_Bclim_mB  [1:ntime_plot]), 
  max(F_S1loc_mB  [1:ntime_plot]),
  max(F_Tloc_mB   [1:ntime_plot])
  ) *1

pngname=paste0("CST_Fnorm.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))

plot(F_Bclim_mB[1:ntime_plot], 
     ylim=c(mn,mx), 
     xlab="Time step", ylab="Frobenius Norm",  type="p", col="black",
     main=paste("Errors in covariances"),
     cex.main=1.5, cex.axis=1.2, cex.lab=1.4)

lines(F_S1loc_mB [1:ntime_plot], type="l", lwd=1.5, lty=3, col="black")
lines(F_Tloc_mB[1:ntime_plot], type="l", lwd=2, col="black")

leg.txt=c(
          as.expression( bquote( "|" ~ "B"[clim] ~ "-B|" ) ),
          as.expression( bquote( "|" ~ "S"[loc]  ~ "-B|" ) ),
          as.expression( bquote( "|" ~ "T"[loc]  ~ "-B|" ) )
          )
          
leg.col<-c("black", "black", "black")
legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(NA,1.5,2), 
       lty=c(1,3,1), pt.lwd=1.0, pch=c(1,NA,NA),
       cex=1.0, pt.cex=1, bg="white")
dev.off()


errvar=FALSE
if(errvar == TRUE){
  
pngname=paste0("BClimS1T1_trace.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(Tr_B[1:ntime_plot], ylim=c(0,max(Tr_B[1:ntime_plot])),
     xlab="Time step", ylab="Rms, diagonal entries", type="p", pch=1,
     main=paste("Errors in variances.Ta=", Ta_h,"h"),
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(Tr_S1mB[1:ntime_plot], type="l", lwd=1, col="blue")
lines(Tr_T1mB[1:ntime_plot], type="l", lwd=2, col="black")
lines(Tr_climBmB[1:ntime_plot], type="l", lwd=3, col="green")
leg.txt<-c('|B|', '|S-B|', '|T-B|')
leg.col<-c("black", "blue", "black")
legend("topright", inset=0, leg.txt, col=leg.col, lwd=c(NA,1,2), 
       lty=c(NA,1,1), pt.lwd=2.0, pch=c(1,NA,NA),
       cex=1.4, pt.cex=1, bg="white")
dev.off()

}  # End if(errvar == TRUE)

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

cat("F_B_mean=")
print(F_B_mean)

cat("F_Bclim_mB_mean=")
print(F_Bclim_mB_mean)

cat("F_S1loc_mB_mean=")
print(F_S1loc_mB_mean)

cat("F_Tloc_mB_mean=")
print(F_Tloc_mB_mean)

cat("F_EV_mB_mean=")
print(F_EV_mB_mean)

cat("F_EVT_mB_mean=")
print(F_EVT_mB_mean)




sink()
#------------------------------------------------
