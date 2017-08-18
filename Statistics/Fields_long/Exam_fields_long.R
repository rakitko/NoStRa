# Examine the output from one world of the DSM on S1:
# 1) Calc time crf of xi
# 3) Gau q-q plot for xi (unconditionally on the scnd flds)
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

datafile=paste0("./fields_long.RData")


load(datafile) # contains:
#           parameters:
#             dim=n, delta_t=delta_t_model, 
#             stride defines the asml time step: Ta=delta_t * stride
#             ...
#           all scnd fields: 
#             Rho, Nu, U, Sigma
#           xi=xi_full

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

rem=6.37*10^6 # m
rekm=rem/10^3 # km

mesh <- 2*pi/n  # rad
mesh_km <- mesh * rekm # grid spacing, km
grid_km      <- c(0:(n-1))*mesh_km
step_h=filter$parameters$delta_t / 3600
ntime4plots_steps=600 # floor(4000/step_h)
ntime4plots_steps=min(ntime,ntime4plots_steps)
tgrid_h=c(1:ntime)*step_h

#------------------------------------------------
# Plots: xi and the scnd flds

nt1=5
if(ntime < nt1*ntime4plots_steps) nt1=1
dt1=ntime/nt1

for (it1 in 1:nt1){
  t1=(it1 -1)*dt1 +1 
  t2=t1 + ntime4plots_steps -1
  
  pngname=paste0("xi_",t1, ".png")
  png(pngname, width=5.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(xi[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
          xlab="Space, km", ylab="Time, h",
          main=as.expression( bquote(xi) ), 
          cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  dev.off()
  
  pngname=paste0("rho_",t1, ".png")
  png(pngname, width=5.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(rho[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
          xlab="Space, km", ylab="Time, h",
          main=as.expression( bquote(rho) ), 
          cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  dev.off()
  
  
  pngname=paste0("nu_",t1, ".png")
  png(pngname, width=5.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(nu[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
          xlab="Space, km", ylab="Time, h",
          main=as.expression( bquote(nu) ), 
          cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  dev.off()
  
  pngname=paste0("sigma_",t1, ".png")
  png(pngname, width=5.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(sigma[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
          xlab="Space, km", ylab="Time, h",
          main=as.expression( bquote(sigma) ), 
          cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  dev.off()
  
  pngname=paste0("U_",t1, ".png")
  png(pngname, width=5.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  image2D(U[,t1:t2], x=grid_km, y=tgrid_h[t1:t2], 
          xlab="Space, km", ylab="Time, h",
          main=paste('U'), 
          cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
  dev.off()
  
}

#------------------------------------------------
# NGaussianity of  xi

xx=as.vector(xi)
nn=length(xx)
rarefy=16 # >1 to save CPU time while uniformly sampling the whole spacetimea 64 is fast & OK.
ind=seq(from=1, to=nn, by=rarefy)
yy=xx[ind]

pngname=paste0("QQplot_xi_perturb_",theta_name,".png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
qqnorm(yy, main = paste("Normal q-q plot: xi"),
       xlab="Theoretical quantile", ylab="Sample quantile", 
       cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
qqline(yy)
dev.off()

kurtosis_xi=kurtosis(yy)
skewness_xi=skewness(xx)

#------------------------------------------------
# Time crf estm

dtmax=240
ttcrf=matrix(0, nrow=n, ncol=dtmax+1)

nt1=2
dt1=ntime/nt1

for (it1 in 1:nt1){
  t1=(it1 -1)*dt1 +1 
  t2=t1 + ntime/nt1 -1
  for (i in 1:n){
    ttcrf[i,]=est_crf_timser(xi[i,t1:t2], dtmax)$crf
  }
  tcrf=colMeans(ttcrf)
  
  pngname=paste0("TimCrl_mean_",t1, ".png")
  png(pngname, width=7.1, height=5.1, units = "in", res=300)
  par(mai=c(1.2,1.2,0.7,0.7))
  plot(tgrid_h[1:(dtmax+1)] - grid_km[1], tcrf,  ylim=c(0,1),
       xlab="Time, h", ylab="Correlation",
       main="Mean temporal crf")
  dev.off()
  
}

t1=2
t2=ntime

for (i in 1:n){
  ttcrf[i,]=est_crf_timser(xi[i,t1:t2], dtmax)$crf
}
tcrf=colMeans(ttcrf)

pngname=paste0("TimCrl_mean_", ".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(tgrid_h[1:(dtmax+1)] - grid_km[1], tcrf,  ylim=c(0,1),
     xlab="Time, h", ylab="Correlation",
     main="Mean temporal correlation")
dev.off()

it_1h=1/step_h +1
tcrf[it_1h]

it_6h=6/step_h +1
tcrf[it_6h]

it_12h=12/step_h +1
tcrf[it_12h]

it_24h=24/step_h +1
tcrf[it_24h]

it_48h=48/step_h +1
tcrf[it_48h]

it_96h=96/step_h +1
tcrf[it_96h]

it_120h=120/step_h +1
tcrf[it_120h]

it_144h=144/step_h +1
tcrf[it_144h]

#------------------------------------------------
# Spatial crf estm

dsmax=20
nt=ntime /10
sscrf=matrix(0, nrow=nt, ncol=dsmax+1)

for (t in 1:nt){
  sscrf[t,]=est_crf_timser(xi[,t], dsmax)$crf
}
scrf=colMeans(sscrf)

pngname=paste0("SpaCrl_mean_",theta_name, ".png")
png(pngname, width=7.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(grid_km[1:(dsmax+1)] - grid_km[1], scrf, ylim=c(-0.1,1),
     xlab="Distance, km", ylab="Correlation",
     main="Mean spatial correlation")
dev.off()

#------------------------------------------------
#------------------------------------------------
# Write output file

outfile=paste0("out_xi_",theta_name,".txt")
unlink(outfile)
sink(outfile, append=TRUE)

cat("\n")

cat("kurtosis_xi=")
print(kurtosis_xi)

cat("skewness_xi=")
print(skewness_xi)

cat("\n")

cat("time crf at 6, 12, 24, 48, 96, 120, 144 h")

cat("\n")

it_6h=6/step_h +1
tcrf[it_6h]

it_12h=12/step_h +1
tcrf[it_12h]

it_24h=24/step_h +1
tcrf[it_24h]

it_48h=48/step_h +1
tcrf[it_48h]

it_96h=96/step_h +1
tcrf[it_96h]

it_120h=120/step_h +1
tcrf[it_120h]

it_144h=144/step_h +1
tcrf[it_144h]

cat("tcrf[1:25]=")
print(tcrf[1:25])


sink()
#------------------------------------------------
