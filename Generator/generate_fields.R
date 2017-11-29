#setwd('/mnt/dev/rakitko/hbef/')



## sources
source("./functions_0.6.R")
library(Rcpp)
library(RcppArmadillo)
library(plot3D)
Rcpp::sourceCpp('generate_fields_functions.cpp')

path                            <- './Out'

# -----------

config <- read.table('./config.txt', sep = ';')


# NB: time_filter is the number of assimilation steps (analyses), whereas
#     time is the number of model time steps
#
#     delta_t is the model time step, so that
#     the time interval between the cosecutive analyses is delta_t*stride 
#     (stride is the synonym for thin_time_for_filter)

dim                             <- config[config$V1 == "dim", 2]
thin_time_for_filter             <- config[config$V1 == "stride", 2]   # model time steps in one asml cycle
time                            <- config[config$V1 == "time_filter", 2] * thin_time_for_filter  # nu of model time steps
delta_t                         <- config[config$V1 == "delta_t", 2]
Re                              <- config[config$V1 == "Re", 2]
u_mean                          <- config[config$V1 == "u_mean", 2]
u_char                          <- config[config$V1 == "u_char", 2]
sd_X                            <- config[config$V1 == "sd_X", 2]
sqrt_R                          <- config[config$V1 == "sqrt_R", 2]
RHO_kappa                       <- config[config$V1 == "RHO_kappa", 2]
RHO_pi                          <- config[config$V1 == "RHO_pi", 2]
Nu_kappa                        <- config[config$V1 == "Nu_kappa", 2]
Nu_pi                           <- config[config$V1 == "Nu_pi", 2]
Sigma_kappa                     <- config[config$V1 == "Sigma_kappa", 2]
sd_U                            <- config[config$V1 == "sd_U", 2]
m                               <- config[config$V1 == "m", 2]                 #OBS GRID MESH SIZE
TOTAL                           <- config[config$V1 == "M", 2]             # number of total repeats (=L, in old versions)
N                               <- config[config$V1 == "N", 2]
localization_mult               <- config[config$V1 == "localization_mult", 2]
L_mult                          <- config[config$V1 == "L_mult", 2]
L_perturb_mult                  <- config[config$V1 == "L_perturb_mult", 2]
inflation_enkf                  <- config[config$V1 == "inflation_enkf", 2]
inflation_var                   <- config[config$V1 == "inflation_var", 2]
seed_for_secondary_fields       <- config[config$V1 == "seed_for_secondary_fields", 2]
seed_for_filters                <- config[config$V1 == "seed_for_filters", 2]
perform_enkf                    <- config[config$V1 == "perform_enkf", 2]
perform_VAR                     <- config[config$V1 == "perform_VAR", 2]
compute_field_true_covs         <- config[config$V1 == "compute_field_true_covs", 2]

stride=thin_time_for_filter
ntime_cycles=time/stride
ntime_model_time_steps=time
# -----------


## calculated parameters
a                               <- 1 / sqrt(2*pi*Re)
delta_x                         <- 2*pi*Re/dim
L                               <- L_mult*delta_x
L_perturb                       <- L*L_perturb_mult
T_sc                            <- L / u_char
RHO_mean                        <- 1/T_sc * sum(1/(((seq(-dim/2,dim/2,1))*L/Re)^2+1)^2) / sum(1/(((seq(-dim/2,dim/2,1))*L/Re)^2+1))
Nu_mean                         <- L^2*RHO_mean 
R                               <- rep(sqrt_R^2, dim%/%m+min(1,dim%%m))

L_loc                           <- localization_mult*L
C_enkf                          <- create_loc_matrix(dim, rep(1,dim), L_loc)
inflation_enkf_bounds           <- c(1.07, 2.72)
inflation_enkf_coef             <- c(sqrt(1),sqrt(1),sqrt(1))*sqrt(inflation_enkf)
inflation_beta                  <- 0.4*0
# -----------

# initialize arrays and dirs
enkf_S_nonloc_arr_mean          <- array(0, dim = c(dim, dim, time/thin_time_for_filter))
filter                          <- list()
dir.create(path)
dir.create(paste0(path,'/DATA'))


# generate secondary fields
params_for_field                <- get_params_for_field(L_perturb, u_char, RHO_kappa, Nu_kappa, sd_U, Sigma_kappa)

sd_RHO                          <- log(RHO_kappa)
sd_Nu                           <- log(Nu_kappa)
sd_Sigma                        <- log(Sigma_kappa)

# create list for parameters (for output file)
parameters                      <- list()
parameters$dim                  <- dim
parameters$time                 <- time
parameters$delta_t              <- delta_t  # model time integration step, s
parameters$u_mean               <- u_mean
parameters$u_char               <- u_char
parameters$sd_X                 <- sd_X
parameters$Rho_kappa            <- RHO_kappa
parameters$Rho_pi               <- RHO_pi
parameters$Nu_kappa             <- Nu_kappa
parameters$Nu_pi                <- Nu_pi
parameters$sd_U                 <- sd_U
parameters$sqrt_R               <- sqrt_R
parameters$L_loc                <- L_loc
parameters$sigma_kappa          <- Sigma_kappa
parameters$M                    <- TOTAL
parameters$mesh_obs             <- m
parameters$L_mean               <- L
parameters$L_perturb            <- L_perturb
parameters$N                    <- N
parameters$inflation_enkf_coef  <- inflation_enkf 
parameters$stride               <- thin_time_for_filter 
parameters$seed_for_secondary_fields <- seed_for_secondary_fields
parameters$seed_for_filters     <- seed_for_filters


U_sigma_mean                    <- (sd_U)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
Sigma_mean                      <- (sd_X)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
Sigma_sigma_mean                <- (sd_Sigma)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)

RHO_epilon                      <- exp(log(RHO_kappa)*qnorm(RHO_pi)) / (1 - exp(log(RHO_kappa)*qnorm(RHO_pi)))
RHO_sigma_mean                  <- (sd_RHO)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)

Nu_epilon                       <- exp(log(Nu_kappa)*qnorm(Nu_pi)) / (1 - exp(log(Nu_kappa)*qnorm(Nu_pi)))
Nu_sigma_mean                   <- (sd_Nu)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)



set.seed(seed_for_secondary_fields)

#-----------------------------------------------
# scnd fields
#    generating RHO
if(sd_RHO != 0){
  RHO_start <- start_value(a, params_for_field$RHO_sigma_mean, params_for_field$RHO_mean, params_for_field$Nu_mean, Re, dim)
  RHO       <- generate_field_cpp(time, dim, RHO_start, delta_t, matrix(u_mean, nrow = dim, ncol = time), 
                                  matrix(params_for_field$RHO_mean, nrow = dim, ncol = time), matrix(params_for_field$RHO_sigma_mean, nrow = dim, ncol = time), 
                                  matrix(params_for_field$Nu_mean, nrow = dim, ncol = time), 1, create_cov_matrix)
  RHO       <- RHO_mean * ModifLogist(RHO, 1, params_for_field$RHO_epsilon)
}else{
  RHO       <- matrix(RHO_mean, nrow = dim, ncol = time)
}

#    generating u
if(sd_U != 0){
  U_start <- start_value(a, U_sigma_mean, params_for_field$RHO_mean, params_for_field$Nu_mean, Re, dim)
  U       <- generate_field_cpp(time, dim, U_start, delta_t, matrix(u_mean, nrow = dim, ncol = time), 
                                matrix(params_for_field$RHO_mean, nrow = dim, ncol = time), matrix(params_for_field$U_sigma_mean, nrow = dim, ncol = time), 
                                matrix(params_for_field$Nu_mean, nrow = dim, ncol = time), 1, create_cov_matrix)
  U       <- u_mean + U
}else{
  U       <- matrix(u_mean, nrow = dim, ncol = time)
}

#    generating nu
if(sd_Nu != 0){
  Nu_start <- start_value(a, params_for_field$Nu_sigma_mean, params_for_field$RHO_mean, params_for_field$Nu_mean, Re, dim)
  Nu       <- generate_field_cpp(time, dim, Nu_start, delta_t, matrix(u_mean, nrow = dim, ncol = time), 
                                 matrix(params_for_field$RHO_mean, nrow = dim, ncol = time), matrix(params_for_field$Nu_sigma_mean, nrow = dim, ncol = time), 
                                 matrix(params_for_field$Nu_mean, nrow = dim, ncol = time), 1, create_cov_matrix)
  Nu       <- Nu_mean * ModifLogist(Nu, 1, params_for_field$Nu_epsilon)
}else{
  Nu       <- matrix(Nu_mean, nrow = dim, ncol = time)
}

#    generating Sigma
if(sd_Sigma != 0){
  Sigma_start <- start_value(a, params_for_field$Sigma_sigma_mean, params_for_field$RHO_mean, params_for_field$Nu_mean, Re, dim)
  Sigma       <- generate_field_cpp(time, dim, Sigma_start, delta_t, matrix(u_mean, nrow = dim, ncol = time), 
                                    matrix(params_for_field$RHO_mean, nrow = dim, ncol = time), matrix(params_for_field$Sigma_sigma_mean, nrow = dim, ncol = time), 
                                    matrix(params_for_field$Nu_mean, nrow = dim, ncol = time), 1, create_cov_matrix)
  Sigma       <- Sigma_mean * ModifLogist(Sigma, 1, 0)
}else{
  Sigma       <- matrix(Sigma_mean, nrow = dim, ncol = time)
}

# Estimation of B_start
cat('Estimating of B_start')
ind_obs     <- seq(1,dim,m)
ind_time    <- seq(1,time,thin_time_for_filter)

start_X     <- start_value(a, Sigma_mean, RHO_mean, Nu_mean, Re, dim)
start_A     <- diag(rep(mean(R), dim))
X_temp      <- generate_field_cpp(time, dim, start_X, delta_t, U, RHO, Sigma, Nu, 1, create_cov_matrix)
OBS_NOISE   <- matrix(rnorm((dim%/%m+min(1,dim%%m))*time,0,1),nrow=dim%/%m+min(1,dim%%m),ncol=time)
OBS         <- X_temp[ind_obs,] +sqrt(R)*OBS_NOISE
X_true      <- X_temp[,ind_time]
kalman_res  <- kalman(time, dim, delta_t, U, RHO, Sigma, Nu, L_C, R, m, OBS, start_X, start_A, thin_time_for_filter)

# # временный вывод 
# plot(kalman_res$B_arr[,1,1])
# #View(kalman_res$X_a)
# #View(kalman_res$X_f)
# tmom <- 1
# rmse(kalman_res$X_a[tmom,ind_time], X_temp[tmom,ind_time])
# rmse(kalman_res$X_f[tmom,ind_time], X_temp[tmom,ind_time])
# rmse(kalman_res$X_a[tmom,ind_time+1], X_temp[tmom,ind_time+1])
# rmse(kalman_res$X_f[tmom,ind_time+1], X_temp[tmom,ind_time+1])

B_kalman_mean <- apply(kalman_res$B_arr,c(2,3),function(x) mean(x, na.rm = TRUE))
#------------------------------------------------------
# Worlds

X             <- list()
set.seed(seed_for_filters)
print('Generating worlds')
for(iter in 1:TOTAL){  
  cat("\r",paste0(round(iter/TOTAL*100,0),'%'))
  ind_obs     <- seq(1,dim,m)
  ind_time    <- seq(1,time,thin_time_for_filter)
  
  start_X     <- start_value(a, Sigma_mean, RHO_mean, Nu_mean, Re, dim)
  X[[iter]]   <- generate_field_cpp(time, dim, start_X, delta_t, U, RHO, Sigma, Nu, 1, create_cov_matrix)
  X_true                   <- X[[iter]][,ind_time]
  save(X_true, file = paste0(path,'/DATA/truth_',iter, '.Rdata'))
  
  if(iter == TOTAL){
    filter$xi_full   <- X[[iter]]
    filter$xi_sparse <- X[[iter]][,ind_time]
  }
  
  OBS_NOISE   <- matrix(rnorm((dim%/%m+min(1,dim%%m))*time,0,1),nrow=dim%/%m+min(1,dim%%m),ncol=time)
  OBS         <- X[[iter]][ind_obs,] +sqrt(R)*OBS_NOISE
  
  start_X     <- X[[iter]][,1]
  
  if(perform_enkf == 1){
    enkf_res    <- ensembl_kalman_cpp(time, dim, delta_t, U, RHO, Sigma, Nu,
                                 1, R, ind_obs-1, OBS, start_X,
                                 create_cov_matrix, N, inflation_enkf_bounds, inflation_enkf_coef, C_enkf, thin_time_for_filter)

    enkf_S_nonloc_arr_mean   <- enkf_S_nonloc_arr_mean + enkf_res$S_nonloc_arr
    if(iter == TOTAL){
      enkf_one_world <- enkf_res
    }
    enkf_res$S_nonloc_arr    <- NULL
    enkf_res$B_arr    <- NULL
    save(enkf_res, file = paste0(path,'/DATA/enkf_',iter, '.Rdata'))
  }
  
  if(perform_VAR == 1){
    VAR_res    <- VAR(time, dim, delta_t, U, RHO, Sigma, Nu, R, m, OBS, start_X, B_kalman_mean, inflation_var)
    save(VAR_res, file = paste0(path,'/DATA/VAR_',iter, '.Rdata'))
  }
  
}

#-----------------------------------------------

print('Covariance matrix of xi')

if(compute_field_true_covs == 1){
  X_arr <- array(NA, dim = c(TOTAL, dim, time))
  
  for(iter in 1:TOTAL){
    cat("\r",paste0(round(iter/TOTAL*100,0),'%'))
    X_arr[iter,,] <- X[[iter]]
  }
  
  Cov_mat<- array(NA, dim = c(dim, dim, time))

  for(t in 1:time){
    Cov_mat[,,t] <- cov(X_arr[,,t])
  }
  filter$Cov_mat  <- Cov_mat
}

#-----------------------------------------------
print('EnKF')

if(perform_enkf == 1){
  enkf_S_nonloc_arr_mean        <- enkf_S_nonloc_arr_mean / TOTAL
  
  Cov_mat_enkf                  <- array(0, dim = c(dim, dim, time / thin_time_for_filter))
  
  for(iter in 1:TOTAL){
    load(paste0(path,'/DATA/truth_',iter,'.Rdata'))  
    load(paste0(path,'/DATA/enkf_',iter,'.Rdata'))  
    for(step in 1:(time/thin_time_for_filter)){
      Cov_mat_enkf[,,step] <- Cov_mat_enkf[,,step] + ((enkf_res$X_f[,step] - X_true[,step])) %*% t(enkf_res$X_f[,step] - X_true[,step])
    }
    cat("\r",paste0(round(iter/TOTAL*100,0),'%'))
  }
  
  Cov_mat_enkf                  <- Cov_mat_enkf / TOTAL
  filter$enkf_one_world         <- enkf_one_world
  filter$enkf_S_nonloc_arr_mean <- enkf_S_nonloc_arr_mean
  filter$Cov_mat_enkf           <- Cov_mat_enkf
}


if(perform_enkf == 1){
  if(TOTAL == 1){  # evaluate time mean fcst rmse
    print('One-world fcst-err stats')
    
    nn=dim*ntime_cycles
    fc_RMSE=sqrt( sum( (enkf_one_world$X_f - X_true)^2 ) /nn )
    print("fc_RMSE=")
    print(fc_RMSE)
    fe_max=max(enkf_res$X_f - X_true)
    print("fe_max=")
    print(fe_max)
    #xtmax=which(enkf_res$X_f - X_true == fe_max, arr.ind=TRUE)
    #xmax=xtmax[1]
    #tmax=xtmax[2]
    #plot(X_true[xmax, (tmax-100):(tmax+100)], main="X_true near fe_max")
    #lines(enkf_res$X_f[xmax, (tmax-100):(tmax+100)])
  }
} 
#-----------------------------------------------

filter$parameters             <- parameters
filter$Rho                    <- RHO
filter$Nu                     <- Nu
filter$U                      <- U
filter$Sigma                  <- Sigma[,1:time]

#-----------------------------------------------
# write output (normally, either model or filter output is computed)

if(compute_field_true_covs == 1) {                     # model output
  save(filter, file = paste0(path,'/model.RData'))
}

if(perform_enkf == 1) {                                 # filter output
  save(filter, file = paste0(path,'/filter.RData'))
}

if(compute_field_true_covs != 1 & perform_enkf != 1) {  # fields_long only
  save(filter, file = paste0(path,'/fields_long.RData'))
}
#--------------------------------------------
