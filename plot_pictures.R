#setwd('/mnt/dev/rakitko/hbef/')
name <- 'thin_6'
path <- paste0('./Images/',name,'/')

variable <- paste0('_',name)


def_plt <- par('plt')

load(paste0(path,'/parameters.Rdata'))
dim           <- parameters$dim
time          <- parameters$time / parameters$thin_time_for_filer
delta_t       <- parameters$delta_t
Re            <- 6370000  
#m             <- 4 #OBS GRID MESH SIZE
TOTAL         <- parameters$M # number of total repeats (=L, in old versions)


load(paste0(path,'diag_mean.RData'))
load(paste0(path,'spectr_mean.RData'))
load(paste0(path,'Cov_mat_enkf.RData'))
load(paste0(path,'macroscale_st_km_mean.RData'))
load(paste0(path,'microscale_st_km_mean.RData'))
load(paste0(path,'macroscale_arr.RData'))
load(paste0(path,'microscale_arr.RData'))
load(paste0(path,'diag_arr.RData'))


rem <- 6.37*10^6 # m
rekm <- rem/10^3 # km
mesh <- 2*pi/dim  # rad
mesh_km <- mesh * rekm # grid spacing, km

CRM <- array(NA, dim = c(dim, dim ,time))
macroscale_st_km  <- matrix(0, nrow=dim, ncol=time)
microscale_st_km  <- matrix(0, nrow=dim, ncol=time)
std_st            <- matrix(0, nrow=dim, ncol=time)

for (t in (1:time)){
  for (ii in (1:dim)){
    std_st[ii, t] = sqrt(Cov_mat_enkf[ii,ii,t])      # B is the Cov Mx
  }
}

for (t in (1:time)){
  CRM[,,t]=diag(1/std_st[,t]) %*% Cov_mat_enkf[,,t] %*% diag(1/std_st[,t])      # Crl Mx
  for (ii in (1:dim)){
    macroscale_st_km[ii, t]  = sum(CRM[ii,,t])*mesh_km/2
    iim1=ii-1
    if(ii==1) iim1=dim
    iip1=ii+1
    if(ii==dim) iip1=1
    microscale_st_km[ii, t] = 1/sqrt( (-CRM[ii,iim1,t] + 2*CRM[ii,ii,t] - CRM[ii,iip1,t])/mesh^2 ) * rekm
  }
}  
microscale_st_km_true <- microscale_st_km
macroscale_st_km_true <- macroscale_st_km

est_err <- c()
true_val <- c()
for(step in 51:time){
  est_err <- c(est_err, microscale_st_km_mean[,step]-microscale_st_km_true[,step])
  true_val <- c(true_val, microscale_st_km_true[,step])
}
mean_err <- mean(est_err)

m_e <- round(mean_err,0)
png(paste0(path,'microscale',variable,'.png'), width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = expression(lambda[true]), ylab = bquote(paste('<',lambda,'> - ',lambda[true])), main = bquote(paste('(c)  Bias ',lambda, ', mean = ', .(m_e))), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
dev.off()




#---------------------------------------------

est_err <- c()
true_val <- c()
for(step in 51:time){
  est_err <- c(est_err, macroscale_st_km_mean[,step]-macroscale_st_km_true[,step])
  true_val <- c(true_val, macroscale_st_km_true[,step])
}
mean_err <- mean(est_err)
m_e <- round(mean_err,0)

png(paste0(path,'/macroscale',variable,'.png'), width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = expression(Lambda[true]), ylab = bquote(paste('<',Lambda,'> - ',Lambda[true])), main = bquote(paste('(b)  Bias ',Lambda, ', mean = ', .(m_e))), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
#lines(lowess(true_val, est_err, f = 0.3), col = 'lightcoral', lwd=3)
dev.off()

#---------------------------------------------


sd_micro <- apply(microscale_arr, c(1,2), sd)
est_err <- c()
true_val <- c()
for(step in 51:time){
  est_err <- c(est_err, sd_micro[,step])
  true_val <- c(true_val, microscale_st_km_true[,step])
}
mean_err <- mean(est_err)
m_e <- round(mean_err,0)
png(paste0(path,'/sd_microscale',variable,'.png'), width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = expression(lambda[true]), ylab = bquote(paste('SD(',lambda,')')), main = bquote(paste('Std.dev ',lambda, ', mean = ', .(m_e))), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
#lines(lowess(true_val, est_err, f = 0.3), col = 'lightcoral', lwd=3)
dev.off()

#---------------------------------------------

sd_macro <- apply(macroscale_arr, c(1,2), sd)
est_err <- c()
true_val <- c()
for(step in 51:time){
  est_err <- c(est_err, sd_macro[,step])
  true_val <- c(true_val, macroscale_st_km_true[,step])
}
mean_err <- mean(est_err)
m_e <- round(mean_err,0)
png(paste0(path,'/sd_macroscale',variable,'.png'), width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = expression(Lambda[true]), ylab = bquote(paste('SD(',Lambda,')')), 
     main = bquote(paste('Std.dev ',Lambda, ', mean = ', .(m_e))), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
dev.off()

#---------------------------------------------

est_err <- c()
true_val <- c()
for(step in 50:time){
  est_err <- c(est_err, diag_mean[,step]-diag(Cov_mat_enkf[,,step]))
  true_val <- c(true_val, diag(Cov_mat_enkf[,,step]))
}
mean_err <- mean(est_err)
m_e <- round(mean_err,2)
png(paste0(path,'/B',variable,'.png'), width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = paste0('B'), ylab = '<S> - B', main = bquote(paste('(a)  Bias S;  mean = ', .(m_e))), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
dev.off()

#---------------------------------------------

sd_B <- apply(diag_arr, c(1,2), sd)
est_err <- c()
true_val <- c()
for(step in 50:time){
  est_err <- c(est_err, sd_B[,step])
  true_val <- c(true_val, diag(Cov_mat_enkf[,,step]))
}
mean_err <- mean(est_err)
m_e <- round(mean_err,2)
png(paste0(path,'/sd_B',variable,'.png'), width=5.1, height=5.1, units = "in", res=300)
par(mai=c(1.2,1.2,0.7,0.7))
plot(true_val, est_err, xlab = paste0('B'), ylab = 'SD(S)', main = bquote(paste('St. dev S;  mean = ', .(m_e))), 
     cex.main=1.7, cex.axis=1.3, cex.lab=1.6)
lines(lowess(true_val, est_err, f = 0.3), col = 'aliceblue', lwd=3)
dev.off()

