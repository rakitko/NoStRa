# Examine the output of the DSM on S1:
# 1) Plot: 
#    xi,
#    the secondary field that are perturbed, 
#    V, 
#    Lambda
#    lambda.
# 2) Calc: 
#    Variability (max/min and the ratio 90% / 10% quantile) in
#      V, Lambda, microscale
#    crl(theta, V)
#    crl(V, Lambda)
#
# M Tsyrulnikov
# 13 Aug 2017

library(mixAK)
library(MASS)
library(stats)
library(plot3D)
library(psych)
library(nortest)
library(moments)
source('est_crf_timser.R')


universe= "All" # Rho  Nu  U  Sigma All

datafile=paste0("./model.RData")


load(datafile) # contains:
#           parameters:
#             dim=n, delta_t=delta_t_model, 
#             stride defines the asml time step: Ta=delta_t * stride
#             ...
#           all scnd fields: 
#             Rho, Nu, U, Sigma
#           Cov_mat=CVM_truth[1:n,1:n, 1:ntime] 


n=filter$parameters$dim
ntime_model=filter$parameters$time
ntime=ntime_model
M=filter$parameters$M
seed_for_secondary_fields=filter$parameters$seed_for_secondary_fields
seed_for_secondary_fields

xi=filter$xi_full
rho=filter$Rho
nu=filter$Nu
sigma=filter$Sigma
U=filter$U


if(universe == "Rho"){
  theta_name="rho" 
  theta=rho
  
} else if (universe == "Nu"){
  theta_name="nu"
  theta=nu
  
} else if (universe == "U"){
  theta_name="U"
  theta=U
  
} else if (universe == "Sigma"){
  theta_name="sigma"
  theta=sigma
  
} else if (universe == "All"){
  theta_name="all"
  theta=rho # !!!!!!! just to fill in smth
} 


B=filter$Cov_mat

rem=6.37*10^6 # m
rekm=rem/10^3 # km

mesh <- 2*pi/n  # rad
mesh_km <- mesh * rekm # grid spacing, km
grid_km      <- c(0:(n-1))*mesh_km
step_h=filter$parameters$delta_t / 3600
tgrid_h=c(1:ntime)*step_h

V <- matrix(0, nrow=n, ncol=ntime)

for (i in (1:n)){
  V[i,]=B[i,i,]
}
std_st = sqrt(V)

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

#------------------------------------------------
# Plots: theta, xi, V, Lambda, lambda

                  if(theta_name == "all"){

pngname=paste0("rho.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(rho, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste("Rho"), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("nu.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(nu, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste("Nu"), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6) 
dev.off()

pngname=paste0("sigma.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(sigma, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste("Sigma"), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6) 
dev.off()

pngname=paste0("U.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(U, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste("U"), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6) 
dev.off()

                  }else{

pngname=paste0(theta_name, ".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(theta, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste(theta_name), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6) 
dev.off()

                  }


pngname=paste0("xi.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(xi, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=as.expression( bquote(xi) ), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


pngname=paste0("V.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(log(V), x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste('log(Var xi)'), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


pngname=paste0("Lambda.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(Lambda, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=as.expression( bquote(Lambda) ), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()


#pngname=paste0("Lambda_mic.png")
#png(pngname, width=5.1, height=5.1, units = "in", res=300)
#par(mai=c(1.2,1.2,0.7,0.7))
#image2D(lambda, x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
#        main=paste("lambda"), 
#        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
#dev.off()


#--------------------------------------------------
# Selected i,t stats

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


pngname=paste0("SpaCrl_Pert_", i, ".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(grid_km[(i-14):(i+14)]-grid_km[i],    CRM[i,(i-14):(i+14),1], 
     xlab="Distance, km", ylab="Correlation", 
     type="l", ylim=c(0,mx), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6,
     main=paste("Spatial field correlations") )
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


mx_mac <- max(macroscale_st_km[,t])
mx_std=max(std_st[,t])

plot (macroscale_st_km[,t]/mx_mac, ylim=c(0,1.1), type="l", col="green",
      main=("GREEN: Rel. Len.scale | RED: Rel. std"))
lines(std_st[,t]/mx_std, type="l", col="red")


plot (microscale_st_km[,t]/macroscale_st_km[,t], type="l", col="green",
      main=("Rel. curv. radius at the origin"))

#------------------------------------------------
# Magnitudes of the spa-tim variability

outfile=paste0("out_",theta_name,".txt")
unlink(outfile)
sink(outfile)


cat("ntime=")
print(ntime)

cat("M=")
print(M)

cat("max(theta) / min(theta)=")
print(max(theta) / min(theta))

cat("max(V) / min(V)=")
print(max(V)/min(V))

cat("max(Lambda) / min(Lambda)=")
print(max(Lambda)/min(Lambda))


cat("max(lambda) / min(lambda)=")
print(max(lambda)/min(lambda))

V_decile_ratio=quantile(V[, 1:ntime], probs=0.9) / quantile(V[, 1:ntime], probs=0.1)
print(V_decile_ratio)

Lambda_decile_ratio=quantile(Lambda[, 1:ntime], probs=0.9) / quantile(Lambda[, 1:ntime], probs=0.1)
print(Lambda_decile_ratio)


sink()

#------------------------------------------------

ind_max=which(Lambda == max(Lambda), arr.ind=TRUE)
plot(B[,ind_max[1],ind_max[2]])

ind_min=which(Lambda == min(Lambda), arr.ind=TRUE)
plot(B[,ind_min[1],ind_min[2]])

#------------------------------------------------
# cross-fld crl

# time shift
# V should be delayed wrt theta

                                if(theta_name != "all"){
nsh=20
ccor_theta_V = c(1:nsh)
ccor_theta_V[]=0

for (sh in (1:20)) {
  ntime_minus_sh=ntime - sh
  ccor_theta_V[sh] = cor(as.vector(V[,(sh+1):ntime]), as.vector(theta[,1:ntime_minus_sh]),
                     method="spearman") # spearman better than pearson,
               
}
ccor_theta_V
crl_theta_V_abs_max=max(abs(ccor_theta_V))
shift_theta_V_opt=which(abs(ccor_theta_V) == abs(crl_theta_V_abs_max))
crl_theta_V_max=ccor_theta_V[shift_theta_V_opt]
shift_theta_V_opt
crl_theta_V_max

crl_Lambda_V=cor(as.vector(Lambda), as.vector(V), method="spearman")
crl_lambda_V=cor(as.vector(lambda), as.vector(V), method="spearman")

micro_d_macro_scale_mean=mean(microscale_st_km / Lambda)

                                }
#------------------------------------------------
# In the U-universe, check if V & Lambda_xi are correlated with 
# dU/dx (x is  an integer in the centered finite difference)
# dU/dt (t is  an integer in the centered finite difference)

dUdx = matrix(0, nrow=n, ncol=ntime)
dUdt = matrix(0, nrow=n, ncol=ntime)
             
if(universe == "U") {

  for (t in c(2:(ntime-1))) {
    for (i in c(1:n)) {
      im1=i-1
      ip1=i+1
      if(i == 1) im1=n
      if(i == n) ip1=1
      dUdx[i,t]=(theta[ip1,t] - theta[im1,t])/2
      dUdt[i,t]=(theta[i,t+1] - theta[i,t-1])/2
    }
  }
  
# smooth dUdt
  
  nsmoo=64 # >0 (16 is good, 64 is even better - for Lambda and V)
  dUdx_=dUdx # initialize

  for (smoo in c(1:nsmoo)){
    dUdx=dUdx_
  
    for (t in c(2:(ntime-1))) {
      for (i in c(1:n)) {
        im1=i-1
        ip1=i+1
        if(i == 1) im1=n
        if(i == n) ip1=1
        dUdx_[i,t]=( dUdx[i,t]+ 
                  0.5*(dUdx[ip1,t]+dUdx[im1,t]+dUdx[i,t-1]+dUdx[i,t+1]) )/3
      }
    }
  }
 
  image2D(dUdx_)
  crl_Lambda_dUdx_smoothed =
    cor( as.vector(Lambda          [,2:(ntime-1)]) , 
         as.vector(dUdx_[,2:(ntime-1)]), method="pearson" ) # spearman  pearson
  crl_lambda_dUdx_smoothed =
    cor( as.vector(microscale_st_km[,2:(ntime-1)]) , 
         as.vector(dUdx_[,2:(ntime-1)]), method="pearson" )
  crl_V_dUdx_smoothed =
    cor( as.vector(V               [,2:(ntime-1)]) , 
         as.vector(dUdx_[,2:(ntime-1)]), method="pearson" )
  
}

# ==> Both V & Lambda positively correlate with a smoothed dU/dx,
#     with nsmoo=64, the crls are both =0.6.
#------------------------------------------------
#------------------------------------------------
# Append the output file

sink(outfile, append=TRUE)
cat("\n")

                  if(theta_name != "all"){
cat("shift_theta_V_opt=")
print(shift_theta_V_opt)

cat("crl_theta_V_max=")
print(crl_theta_V_max)

cat("crl_Lambda_V=")
print(crl_Lambda_V)

cat("crl_lambda_V=")
print(crl_lambda_V)

cat("micro_d_macro_scale_mean=")
print(micro_d_macro_scale_mean)

                  }
cat("\n")


if(universe == "U") {
  cat("nsmoo=")
  print(nsmoo)
  
  cat("crl_Lambda_dUdx_smoothed=")
  print(crl_Lambda_dUdx_smoothed)
  
  cat("crl_lambda_dUdx_smoothed=")
  print(crl_lambda_dUdx_smoothed)
  
  cat("crl_V_dUdx_smoothed=")
  print(crl_V_dUdx_smoothed)
}


sink()
print("Normal finish")
#------------------------------------------------
