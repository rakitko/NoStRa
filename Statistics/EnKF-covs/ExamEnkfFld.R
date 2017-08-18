# Examine 1 realization of the t-x rand.fld & the EnKF analysis field
#
# M Tsyrulnikov
# 16 May 2017

library(mixAK)
library(MASS)
library(stats)
library(plot3D)
library(psych)
library(binhf) # shift

load('./truth_5000.RData') #  X_true[x,t]
dims <- dim(X_true)
dims
n=dims[1]
ntime=dims[2]

load('./enkf_5000.RData') # enkf_res, contains:
#   X_a, X_f, ...

X_enkf_f=enkf_res$X_f
X_enkf_a=enkf_res$X_a

rem=6.37*10^6 # m
rekm=rem/10^3 # km

mesh <- 2*pi/n  # rad
mesh_km <- mesh * rekm # grid spacing, m
grid_km      <- c(0:(n-1))*mesh_km
grid_km_plus <- c( grid_km, (grid_km[n] + mesh_km) ) # periodicity
k_max=n/2
step_h=1 # h

#-------------------------------------------------------
#-------------------------------------------------------
# Plot a fragment of the field

ntime_short=990
t_start=11
if(t_start + ntime_short -1 > ntime) {t_start = ntime - ntime_short}
t_fin=t_start + ntime_short -1
tgrid_h=c(t_start:t_fin)

pngname=paste0("X_true.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(X_true[,t_start:t_fin], x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste('X_true'), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("X_enkf_f.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(X_enkf_f[,t_start:t_fin], x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste('X_enkf_f'), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()

pngname=paste0("X_enkf_a.png")
png(pngname, width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
image2D(X_enkf_a[,t_start:t_fin], x=grid_km, y=tgrid_h, xlab="Space, km", ylab="Time, h",
        main=paste('X_enkf_a'), 
        cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
dev.off()
