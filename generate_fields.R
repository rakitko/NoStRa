#setwd('/mnt/dev/rakitko/hbef/')
path                            <- './test'
config <- read.table('./config.txt', sep = ';')



## sources
source("./functions_0.6.R")
library(Rcpp)
library(RcppArmadillo)
library(plot3D)
Rcpp::sourceCpp('generate_fields_functions.cpp')
# -----------
dim                             <- config[config$V1 == "dim", 2]
thin_time_for_filer             <- config[config$V1 == "stride", 2]
time                            <- config[config$V1 == "time", 2] * thin_time_for_filer
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
TOTAL                           <- config[config$V1 == "TOTAL", 2]             # number of total repeats (=L, in old versions)
N                               <- config[config$V1 == "N", 2]
localization_mult               <- config[config$V1 == "localization_mult", 2]
L_mult                          <- config[config$V1 == "L_mult", 2]
L_perturb_mult                  <- config[config$V1 == "L_perturb_mult", 2]
inflation_enkf                  <- config[config$V1 == "inflation_enkf", 2]
seed_for_secondary_fields       <- config[config$V1 == "seed_for_secondary_fields", 2]
seed_for_filters                <- config[config$V1 == "seed_for_filters", 2]
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
enkf_S_nonloc_arr_mean          <- array(0, dim = c(dim, dim, time/thin_time_for_filer))
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
parameters$delta_t              <- delta_t
parameters$u_mean               <- u_mean
parameters$u_char               <- u_char
parameters$sd_X                 <- sd_X
parameters$Rho_kappa            <- RHO_kappa
parameters$Rho_pi               <- RHO_pi
parameters$Nu_kappa             <- Nu_kappa
parameters$Nu_pi                <- Nu_pi
parameters$sd_U                 <- sd_U
parameters$sdR                  <- sqrt_R
parameters$L_loc                <- L_loc
parameters$sigma_kappa          <- Sigma_kappa
parameters$sd_Sigma             <- sd_Sigma
parameters$M                    <- TOTAL
parameters$mesh_obs             <- m
parameters$L_mean               <- L
parameters$L_perturb            <- L_perturb
parameters$N                    <- N
parameters$C_enkf               <- C_enkf
parameters$inflation_enkf_coef  <- inflation_enkf_coef 
parameters$stride               <- thin_time_for_filer



U_sigma_mean                    <- (sd_U)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
Sigma_mean                      <- (sd_X)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
Sigma_sigma_mean                <- (sd_Sigma)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)

RHO_epilon                      <- exp(log(RHO_kappa)*qnorm(RHO_pi)) / (1 - exp(log(RHO_kappa)*qnorm(RHO_pi)))
RHO_sigma_mean                  <- (sd_RHO)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)

Nu_epilon                       <- exp(log(Nu_kappa)*qnorm(Nu_pi)) / (1 - exp(log(Nu_kappa)*qnorm(Nu_pi)))
Nu_sigma_mean                   <- (sd_Nu)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)



set.seed(seed_for_secondary_fields)

###################################
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

X             <- list()
set.seed(seed_for_filters)
print('Generating worlds')
for(iter in 1:TOTAL){  
  cat("\r",paste0(round(iter/TOTAL*100,0),'%'))
  start_X     <- start_value(a, Sigma_mean, RHO_mean, Nu_mean, Re, dim)
  X[[iter]]   <- generate_field_cpp(time, dim, start_X, delta_t, U, RHO, Sigma, Nu, 1, create_cov_matrix)

  ind_obs     <- seq(1,dim,m)
  ind_time    <- seq(1,time,thin_time_for_filer)
  
  OBS_NOISE   <- matrix(rnorm((dim%/%m+min(1,dim%%m))*time,0,1),nrow=dim%/%m+min(1,dim%%m),ncol=time)
  OBS         <- X[[iter]][ind_obs,] +sqrt(R)*OBS_NOISE

  start_X     <- X[[iter]][,1]
  enkf_res    <- ensembl_kalman_cpp(time, dim, delta_t, U, RHO, Sigma, Nu,
                                 1, R, ind_obs-1, OBS, start_X,
                                 create_cov_matrix, N, inflation_enkf_bounds, inflation_enkf_coef, C_enkf, thin_time_for_filer)

  enkf_S_nonloc_arr_mean   <- enkf_S_nonloc_arr_mean + enkf_res$S_nonloc_arr
  enkf_res$S_nonloc_arr    <- NULL
  X_true                   <- X[[iter]][,ind_time]
  save(enkf_res, file = paste0(path,'/DATA/enkf_',iter, '.Rdata'))
  save(X_true, file = paste0(path,'/DATA/truth_',iter, '.Rdata'))
}

enkf_S_nonloc_arr_mean        <- enkf_S_nonloc_arr_mean / TOTAL

Cov_mat_enkf                  <- array(0, dim = c(dim, dim, time / thin_time_for_filer))

for(iter in 1:TOTAL){
  load(paste0(path,'/DATA/truth_',iter,'.Rdata'))  
  load(paste0(path,'/DATA/enkf_',iter,'.Rdata'))  
  for(step in 1:(time/thin_time_for_filer)){
    Cov_mat_enkf[,,step] <- Cov_mat_enkf[,,step] + ((enkf_res$X_f[,step] - X_true[,step])) %*% t(enkf_res$X_f[,step] - X_true[,step])
  }
  cat("\r",paste0(round(iter/TOTAL*100,0),'%'))
}

Cov_mat_enkf                  <- Cov_mat_enkf / TOTAL
filter                        <- list()
filter$parameters             <- parameters
filter$enkf_S_nonloc_arr_mean <- enkf_S_nonloc_arr_mean
filter$Cov_mat_enkf           <- Cov_mat_enkf
save(filter, file = paste0(path,'/filter.RData'))
