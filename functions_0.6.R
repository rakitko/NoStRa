get_params_for_field <- function(L, u_char, RHO_kappa, Nu_kappa, sd_U, Sigma_kappa){
  T_sc          <- L / u_char
  RHO_mean      <- 1/T_sc * sum(1/(((seq(-dim/2,dim/2,1))*L/Re)^2+1)^2) / sum(1/(((seq(-dim/2,dim/2,1))*L/Re)^2+1))
  Nu_mean       <- L^2*RHO_mean 
  
  sd_RHO   <- log(RHO_kappa)
  sd_Nu   <- log(Nu_kappa)
  sd_Sigma   <- log(Sigma_kappa)
  
  Sigma_sigma_mean <- (sd_Sigma)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
  
  RHO_epsilon <- exp(log(RHO_kappa)*qnorm(RHO_pi)) / (1 - exp(log(RHO_kappa)*qnorm(RHO_pi)))
  RHO_sigma_mean <- (sd_RHO)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
  
  Nu_epsilon <- exp(log(Nu_kappa)*qnorm(Nu_pi)) / (1 - exp(log(Nu_kappa)*qnorm(Nu_pi)))
  Nu_sigma_mean <- (sd_Nu)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
  
  
  U_sigma_mean <- (sd_U)*sqrt(2)/a*(sum(1/((seq(-dim/2,dim/2,1))^2 * Nu_mean/(Re^2) + RHO_mean)))^(-1/2)
  
  res <- list()
  res$RHO_mean <- RHO_mean
  res$Nu_mean <- Nu_mean
  res$sd_RHO <- sd_RHO
  res$sd_Nu <- sd_Nu
  res$RHO_epsilon <- RHO_epsilon
  res$RHO_sigma_mean <- RHO_sigma_mean
  res$Nu_epsilon <- Nu_epsilon
  res$Nu_sigma_mean <- Nu_sigma_mean
  res$U_sigma_mean <- U_sigma_mean
  res$Sigma_sigma_mean <- Sigma_sigma_mean
  return(res)
}

solve_ctrimat <- function(a,b,c,d){
  solution <- rep(NA, length(a)+1)
  c_pr <- rep(NA, length(c))
  d_pr <- rep(NA, length(d))
  for(i in 1:length(a)){
    if(i == 1){
      c_pr[i] <- c[i]/b[i]
      d_pr[i] <- d[i]/b[i]
    }else{
      c_pr[i] <- c[i]/(b[i] - a[i-1]*c_pr[i-1])
      d_pr[i] <- (d[i] - a[i-1]*d_pr[i-1])/(b[i]-a[i-1]*c_pr[i-1])
    }
  }
  n <- length(a)+1
  d_pr[n] <- (d[n] - a[n-1]*d_pr[n-1])/(b[n]-a[n-1]*c_pr[n-1])
  solution[n] <- d_pr[n]
  for(i in (n-1):1){
    solution[i] <- d_pr[i] - c_pr[i]*solution[i+1]
  }
  return(solution)
}

ModifLogist <- function(x,b,eps) (1+exp(b)) * (1+eps) * exp(x-b) /(1 + exp(x-b)) - eps

solve_system <- function(a_0,a,b,c,c_0,d){
  dim <- length(d)
   X_1 <- solve_ctrimat_cpp(a[1:(dim-2)],b[1:(dim-1)],c[1:(dim-2)],d[1:(dim-1)])
   X_2 <- solve_ctrimat_cpp(a[1:(dim-2)],b[1:(dim-1)],c[1:(dim-2)],c(-c_0,rep(0,(dim-3)),-c[dim-1]))
  #X_1 <- Solve.tridiag(a[1:(dim-2)],b[1:(dim-1)],c[1:(dim-2)],d[1:(dim-1)])
  #X_2 <- Solve.tridiag(a[1:(dim-2)],b[1:(dim-1)],c[1:(dim-2)],c(-c_0,rep(0,(dim-3)),-c[dim-1]))
  
  
  x_last <- (d[dim]-a_0*X_1[1]-a[dim-1]*X_1[dim-1])/(a_0*X_2[1]+a[dim-1]*X_2[dim-1]+b[dim])
  result <- X_1 + x_last*X_2
  result <- c(result, x_last)
  return(result)
}

generate_parameter <- function(time, dim, delta_t, u_mean, u_par, L_C_par, L){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  rho <- u_par/Re
  L_mu <- L
  nu <- u_par*L_mu*(1-L_mu/Re)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C_par)
  c             <- C[1,]
  cm <- fft(c, inverse = FALSE)/dim
  am <- 2*pi*(Re(cm))^2*Re
  am <- c(am[1:(dim/2)],am[dim/2+1]/2,am[dim/2+1]/2,am[(dim/2+2):dim])
  vec <- c(0:(dim/2),(-dim/2):(-1))
  sigma <- sqrt(2)*(sum(am/(rho+nu/Re^2*vec^2)))^(-1/2)
  bm <- sigma^2/2*am/(rho + nu/Re^2*vec^2)
  bm <- c(bm[1:(dim/2)],bm[(dim/2+1)]+bm[(dim/2+2)],bm[(dim/2+3):(dim+1)])
  sum(bm)
  xi_re <- rep(0,dim)
  xi_im <- rep(0,dim)
  xi_re[1] <- rnorm(1,0,1)
  xi_re[2:(dim/2+1)] <- rnorm(dim/2,0,1/2)
  xi_re[(dim/2+1):dim] <- xi_re[(dim/2+1):dim] + rev(xi_re[2:(dim/2+1)])
  xi_im[2:(dim/2+1)] <- rnorm(dim/2,0,1/sqrt(2))
  xi_im[(dim/2+1):dim] <- xi_im[(dim/2+1):dim] - rev(xi_im[2:(dim/2+1)])
  start <- Re(fft(sqrt(bm)*(xi_re+1i*xi_im), inverse = TRUE))
  X <- matrix(NA, nrow = dim, ncol = time)
  X[,1] <- start
  # Diagnostics
  B_spa <- Re(fft(bm, inverse = TRUE))
  
  for(i in 2:time){
    a_0   <- u_mean/(2*delta_x)-nu/(delta_x^2)
    a     <- rep(-u_mean/(2*delta_x)-nu/(delta_x^2),dim-1)
    b     <-rep(rho+1/delta_t+2*nu/(delta_x^2), dim)
    c     <- rep(u_mean/(2*delta_x)-nu/(delta_x^2),dim-1)
    c_0   <- -u_mean/(2*delta_x)-nu/(delta_x^2)
    d     <- X[,i-1]/delta_t + sigma*C%*%rnorm(dim,0,sqrt(delta_x)/(sqrt(delta_t)))
    X[,i] <- solve_system(a_0,a,b,c,c_0,d)
  }
  return(list(X=X, B_spa=B_spa))
}

generate_field <- function(time, dim, start_X, delta_t, U, RHO, Sigma, Nu, L_C){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  X <- matrix(NA, nrow = dim, ncol = time)
  X[,1] <- start_X
  for(i in 2:time){
    print((diag(Sigma[,i])%*%C%*%rnorm(dim,0,sqrt(delta_x)/sqrt(delta_t)))[1])
    X[,i] <- predict_cpp(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], X[,i-1], diag(Sigma[,i])%*%C%*%rnorm(dim,0,sqrt(delta_x)/sqrt(delta_t)))  
    #predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], X[,i-1], rep(0,dim))    
  }
  return(X)
}

lorenz <- function(time, dim, start_X, delta_t, F_X){
  SPCmod_1 <- function(t,X,params)  {
    for_1 <- (1:dim)%%dim+1
    rec_1 <- (1:dim-2)%%dim+1
    rec_2 <- (1:dim-3)%%dim+1
    res <- -X[rec_2]*X[rec_1] + X[rec_1]*X[for_1]- X+F_X
    list(res)
  }
  
  
  xstart <- rnorm(dim,0, sd_X)
  times <- seq(0, time*4, length=time)
  times <- seq(0, delta_t*time/(3600*24*5), by = delta_t/(3600*24*5))
  
  params <- c(1)
  
  ## classical Runge-Kutta 4th order
  out <- rk(start_X, times, SPCmod_1,params,hini = 1/240, 
            method = "rk4")
  
  return(t(as.matrix(out)))
}

kalman <- function(time, dim, delta_t, U, RHO, Sigma, Nu, L_C, R, m, OBS_NOISE, start_X, start_A){
  X_a           <- matrix(NA, nrow = dim, ncol = time)
  X_f           <- matrix(NA, nrow = dim, ncol = time)
  B_arr         <- array(NA,c(time,dim,dim))
  A_arr         <- array(NA,c(time,dim,dim))
  A_arr_mod     <- array(NA,c(time,dim,dim))
  X_f[,1]       <- start_X
  X_a[,1]       <- start_X
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  
  A_arr[1,,]    <- start_A
  PHI           <- apply(start_A,2, function(x) predict(dim, delta_x, delta_t, U[,2], Nu[,2], RHO[,2], x, rep(0,dim))) 
  B             <- apply(t(PHI),2, function(x) predict(dim, delta_x, delta_t, U[,2], Nu[,2], RHO[,2], x, rep(0,dim))) 
  B             <- B + diag(Sigma[,2])%*%C%*%C%*%diag(Sigma[,2])*delta_x*delta_t
  
  ind_obs       <- seq(1,dim,m)
  
  for(i in 2:time){
    B_arr[i,,]    <- B
    
    y             <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]              
    BH            <- B[,ind_obs]
    HBH           <- B[ind_obs, ind_obs]
    X_f[,i]       <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], X_a[,i-1], rep(0,dim))  
    z             <- solve(HBH + diag(R), y-X_f[ind_obs,i])
    X_a[,i]       <- X_f[,i] + BH %*% z
    #A             <- B - BH %*% solve(HBH+diag(R)) %*% t(BH)
    K  <- BH %*% solve(HBH + diag(R))
    #A             <- (diag(rep(1,length(ind_obs))) - K[,ind_obs])%*% B %*% t(diag(rep(1,length(ind_obs))) - K[,ind_obs]) + K%*%diag(R)%*%t(K)
    #A             <- (diag(rep(1,dim)) - K,ind_obs])%*% B %*% t(diag(rep(1,dim)) - K[,ind_obs]) + K%*%diag(R)%*%t(K)
    A     <- B - K%*%B[ind_obs,]
    #A             <- B - K %*% t(BH) - BH%*% t(K) +K %*% HBH %*% t(K) + K%*%diag(R)%*%t(K)
    A_arr[i,,]    <- A
    if(i < time){
      PHI           <- apply(A,2, function(x) predict(dim, delta_x, delta_t, U[,i+1], Nu[,i+1], RHO[,i+1], x, rep(0,dim))) 
      B             <- apply(t(PHI),2, function(x) predict(dim, delta_x, delta_t, U[,i+1], Nu[,i+1], RHO[,i+1], x, rep(0,dim))) 
      B             <- B + diag(Sigma[,i+1])%*%C%*%C%*%diag(Sigma[,i+1])*delta_x*delta_t
      #!!!
      #B             <- A + diag(Sigma[,i+1])%*%C%*%C%*%diag(Sigma[,i+1])
      #B             <- t(B)
    }
    
    #eig_decomp <- eigen(B)
    #print(sum(diag(B)))
    #print(c(i,eig_decomp$values[30]))
  }
  return(list(X_a=X_a, X_f=X_f, B_arr = B_arr, A_arr = A_arr, A_arr_mod = A_arr_mod))
}


VAR <- function(time, dim, delta_t, U, RHO, Sigma, Nu, R, m, OBS_NOISE, start_X, start_B, inflation_var){
  X_a           <- matrix(NA, nrow = dim, ncol = time)
  X_f           <- matrix(NA, nrow = dim, ncol = time)
  X_a[,1]       <- start_X
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  B             <- start_B * inflation_var
  ind_obs       <- seq(1,dim,m)
  BH            <- B[,ind_obs]
  HBH           <- B[ind_obs, ind_obs]  
  for(i in 2:time){
    y             <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]              
    X_f[,i]       <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], X_a[,i-1], rep(0,dim))  
    z             <- solve(HBH + diag(R), y-X_f[ind_obs,i])
    X_a[,i]       <- X_f[,i] + BH %*% z
  }
  return(list(X_a=X_a))
}


create_cov_matrix <- function(dim, sigma, L){
    return(diag(rep(1,dim)))
  #   Re            <- 6370000  
#   delta_x       <- 2*pi*Re/dim
#   result <- matrix(NA, ncol = dim, nrow = dim)
#   for(i in 1:dim){
#     for(j in 1:dim){
#       if(abs(i-j)<dim/2){
#         z <- abs(i-j)
#       }else{
#         z <- dim - abs(i-j)
#       }
#       z <- 2*sin(z*delta_x/Re/2)*Re
#       result[i,j] <- (1 + z/L)*exp(-z/L)*sigma[i]*sigma[j]
#       #result[i,j] <- exp(-z/L)
#       #result[i,j] <- max(0, 1-z/L)
#     }
#   }
#   return(result)
}

predict <- function(dim, delta_x, delta_t, u, nu, rho, analysis, noise){
  for_1 <- (1:dim)%%dim+1
  rec_1 <- (1:dim-2)%%dim+1
  a_0   <- u[dim]/(2*delta_x)-nu[dim]/(delta_x^2)
  a     <- -u[2:dim]/(2*delta_x) - nu[2:dim]/(delta_x^2)
  b     <- rep(1/delta_t, dim) + (u[for_1]-u[rec_1])/(2*delta_x)*0 + rho + 2*nu/(delta_x^2)
  c     <- u[1:(dim-1)]/(2*delta_x) - nu[1:(dim-1)]/(delta_x^2)
  c_0   <- -u[1]/(2*delta_x)-nu[1]/(delta_x^2)
  d     <- analysis/delta_t + noise
  
  result <- solve_system_cpp(a_0,a,b,c,c_0,d)
  return(result)
  #return(analysis*0.994)
}


enkf <- function(time, dim, delta_t, N, R, m, inflation_enkf, c_enkf, L_C, U, Nu, RHO, Sigma, OBS_NOISE, start_X, start_A, X, type){
  if(type == 'posterior'){
    Re            <- 6370000  
    delta_x       <- 2*pi*Re/dim
    ind_obs       <- seq(1,dim,m)
    C             <- create_cov_matrix(dim, rep(1,dim), L_C)
    C_loc <- create_loc_matrix(dim, delta_x, c_enkf)
    x_a <- matrix(NA, ncol = time, nrow = dim)
    x_f <- matrix(NA, ncol = time, nrow = dim)
    S_arr <- array(NA, dim=c(time,dim,dim))
    B_arr <- array(NA, dim=c(time,dim,dim))
    B_mean <- matrix(0, ncol=dim, nrow = dim)
    temp <- c(1:time)
    mean_diag_B_f <- c(1:time)
    x_a[,1]       <- start_X
    Q <- C%*%C
    #x_en_ape    <- t(mvrnorm(N,X[,1], start_A))
    x_en_a    <- t(mvrnorm(N,X[,1], start_A))
  }else{
    
    
    
    
    Re            <- 6370000  
    delta_x       <- 2*pi*Re/dim
    ind_obs       <- seq(1,dim,m)
    C             <- create_cov_matrix(dim, rep(1,dim), L_C)
    C_loc <- create_loc_matrix(dim, delta_x, c_enkf)
    x_a <- matrix(NA, ncol = time, nrow = dim)
    x_f <- matrix(NA, ncol = time, nrow = dim)
    x_a[,1]       <- start_X
    x_f[,1]       <- start_X  
    x_en_a    <- t(mvrnorm(N,X[,1], start_A))
    S_arr <- array(NA, dim=c(time,dim,dim))
    B_arr <- array(NA, dim=c(time,dim,dim))
    B_mean <- matrix(0, ncol=dim, nrow = dim)
    temp <- c(1:time)
    mean_diag_B_f <- c(1:time)
    Q <- C%*%C
    
  }
  for(i in (2:time)){
    
    if(type=='posterior'){
      x_en_f  <- apply(x_en_a,2,function(x) {
        predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
      }) 
      
      x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
      mean_en <- rep(rowSums(x_en_f)/N, N)
      x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation_enkf
      x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
      x_me <- (mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
      ###!!!!
      M <- x_en_f + t(x_me) - (matrix(rep(x_f[,i],N), nrow=dim))
      
      #B_f  <- cov(t(x_en_f)) + cov(x_me)
      
      M <- x_en_f  - (matrix(rep(x_f[,i],N), nrow=dim))
      
      S_me     <- 1/N*t(x_me)%*%(x_me)
      S_me     <- S_me * C_loc
      S_pe     <- 1/N*M%*%t(M)
      S_pe      <- S_pe * C_loc
      Q_a <- S_me
      P_a <- S_pe
      B_f     <- P_a + Q_a
      
      
      ### !!!
      #B_f <- 1/N*M%*%t(M)
      
      #S_arr[i,,] <- B_f
      #B_f <- B_f * C_loc
      #B_arr[i,,] <- B_f
      K   <- B_f[,ind_obs] %*% solve(B_f[ind_obs,ind_obs] + diag(R))
      
      #!!!!!!
      temp[i] <- t(x_obs-mean_en[ind_obs])%*%(x_obs-mean_en[ind_obs])
      mean_diag_B_f[i] <-mean(diag(B_f))
      B_mean <- B_mean + B_f
      #?????????????????????????????
      x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i]) 
      A <- B_f - K%*%B_f[ind_obs,]
      #x_en_a  <- t(mvrnorm(N,x_a[,i],A))
      
      A_decomp <- eigen(A)
      D <- diag(sqrt(sapply(A_decomp$values, function(x) max(x,1e-16))))
      L <- A_decomp$vectors
      A_sqrt <- L%*%D%*%t(L)
      Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
      #Samp_matr <- scale(Samp_matr, TRUE, FALSE)
      #Samp_matr <- Samp_matr %*% svd(Samp_matr, nu = 0)$v
      #Samp_matr <- scale(Samp_matr, FALSE, TRUE)
      x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N) + A_sqrt%*%t(Samp_matr)
      
      # when N < dim we can't rescale sample matrix
      # so we take the first N columns of L%*%D
      #M <- (L %*% D ) [,1:N] * sqrt(N)
      #image2D(cov(t(M)),main='cov')
      #image2D(A, main = 'A')
      #image2D(M%*%t(M), main='A_eigencut')
      #x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N)+ M
      
      
    }else{
      x_en_f  <- apply(x_en_a,2,function(x) {
        #predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, diag(Sigma[,i])%*%C%*%rnorm(dim,0,sqrt(delta_x)/sqrt(delta_t)))
        predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
      }) 
      
      x_me <- diag(Sigma[,i])%*%C%*%matrix(rnorm(N*dim,0,sqrt(delta_x)*sqrt(delta_t)), nrow=dim)
      x_en_f <- x_en_f + x_me
      x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
      mean_en <- rep(rowSums(x_en_f)/N, N)
      x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation_enkf
      x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
      M <- x_en_f - (matrix(rep(x_f[,i],N), nrow=dim))
      B_f <- 1/N*M%*%t(M)
      
      #B_f  <- cov(t(x_en_f))
      S_arr[i,,] <- B_f
      B_f <- B_f * C_loc
      B_arr[i,,] <- B_f
      K   <- B_f[,ind_obs] %*% solve(B_f[ind_obs,ind_obs] + diag(R))
      
      #!!!!!!
      temp[i] <- t(x_obs-mean_en[ind_obs])%*%(x_obs-mean_en[ind_obs])
      mean_diag_B_f[i] <-mean(diag(B_f))
      B_mean <- B_mean + B_f
      
      x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i]) 
      x_en_a <- apply(x_en_f,2, function(x){x + K%*%(x_obs + sqrt(R)*rnorm(dim%/%m+min(1,dim%%m),0,1)-x[ind_obs])} )
      
      #x_a[,i] <- rowSums(x_en_a)/N
      
    }
  }
  return(list(X_a=x_a, B_mean=B_mean/time, omf = temp/dim, mean_diag_B_f = mean_diag_B_f, S_arr =S_arr, B_arr=B_arr, X_f=x_f))
}

henkf <- function(chi, phi, beta, theta_p, theta_q,time, dim, delta_t, N, R, m, inflation_enkf, c_enkf, L_C, U, Nu, RHO, Sigma, OBS_NOISE, start_X, start_A, start_B, X, include_obs){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_enkf)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  x_a[,1]       <- start_X
  x_f[,1]       <- start_X  
  x_en_a    <- t(mvrnorm(N,X[,1], start_A))
  S_arr <- array(NA, dim=c(time,dim,dim))
  B_arr <- array(NA, dim=c(time,dim,dim))
  B_mean <- matrix(0, ncol=dim, nrow = dim)
  temp <- c(1:time)
  mean_diag_B_f <- c(1:time)
  Q <- C%*%C
  Q_a <- mean(Sigma^2)*Q*delta_t*delta_x
  
  start_P <- start_B - mean(Sigma^2)*Q*delta_t*delta_x
  P_a <- start_P
  P_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr[1,,] <- Q_a
  P_a_arr[1,,] <- P_a
  N_mat <- (P_a+Q_a)[ind_obs,ind_obs] + diag(R)
  
  for(i in (2:time)){
    x_pe  <- apply(x_en_a,2,function(x) {
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    
    x_me <- diag(Sigma[,i])%*%C%*%matrix(rnorm(N*dim,0,sqrt(delta_x)*sqrt(delta_t)), nrow=dim)
    mean_en <- rep(rowSums(x_pe)/N, N)
    x_pe  <- matrix(mean_en, nrow = dim) + (x_pe - matrix(mean_en, nrow = dim))*inflation_enkf
    
    #x_me  <- t(mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    
    
    
    
    Q_f <- Q_a
    P_f <- P_a
    
    
    S_me     <- x_me %*% t(x_me)/N
    S_me     <- S_me * C_loc
    
    S_pe     <- (x_pe - matrix(rep(x_f[,i],N), nrow = dim))%*%t(x_pe - matrix(rep(x_f[,i],N), nrow = dim))/N
    S_pe      <- S_pe * C_loc
    
    Q_a <- (chi*Q_f + N*S_me)/(chi+N)
    P_a <- (phi * P_f + N * S_pe)/(phi + N)
    
    
    #image2D(N*(S_me+S_pe))
    Q_a_arr[i,,] <- Q_a
    P_a_arr[i,,] <- P_a
    x_en_f <- x_pe + x_me
    
    
    if(all(include_obs == TRUE)){
      v <- (x_obs - x_f[ind_obs,i])
      N_tilde <- (P_a+Q_a)[ind_obs,ind_obs] + diag(R)
      N_mat <- beta * N_mat + (1-beta) * v%*%t(v)
      #N_decomp <- eigen(N_mat)
      #D <- diag((sapply(N_decomp$values, function(x) max(x,1e-16))))
      #L <- N_decomp$vectors
      #N_mat <- L%*%D%*%t(L)
      
      
      #delta_P <- P_a[,ind_obs] %*% solve(N_mat) %*% (v%*%t(v)-N_mat) %*% solve(N_mat) %*% P_a[ind_obs,]
      #delta_Q <- Q_a[,ind_obs] %*% solve(N_mat) %*% (v%*%t(v)-N_mat) %*% solve(N_mat) %*% Q_a[ind_obs,]
      
      delta_P <- P_a[,ind_obs] %*% solve(N_tilde) %*% (N_mat-N_tilde) %*% solve(N_tilde) %*% P_a[ind_obs,]
      delta_Q <- Q_a[,ind_obs] %*% solve(N_tilde) %*% (N_mat-N_tilde) %*% solve(N_tilde) %*% Q_a[ind_obs,]
      Q_a <- Q_a + 1/theta_q*delta_Q
      P_a <- P_a + 1/theta_p*delta_P
    }
    
    
    
    
    B_f <- Q_a + P_a
    #B_decomp <- eigen(B_f)
    #D <- diag((sapply(B_decomp$values, function(x) max(x,1e-16))))
    #L <- B_decomp$vectors
    #B_f <- L%*%D%*%t(L)
    
    
    S_arr[i,,] <- B_f
    B_f <- B_f * C_loc
    B_arr[i,,] <- B_f
    K   <- B_f[,ind_obs] %*% solve(B_f[ind_obs,ind_obs] + diag(R))
    
    #!!!!!!
    temp[i] <- t(x_obs-mean_en[ind_obs])%*%(x_obs-mean_en[ind_obs])
    mean_diag_B_f[i] <-mean(diag(B_f))
    B_mean <- B_mean + B_f
    
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i]) 
    x_en_a <- apply(x_en_f,2, function(x){x + K%*%(x_obs + sqrt(R)*rnorm(dim%/%m+min(1,dim%%m),0,1)-x[ind_obs])} )
    
    #x_a[,i] <- rowSums(x_en_a)/N
    #print(max(abs(x_a[,i])))
    
    
  }
  return(list(X_a=x_a, B_mean=B_mean/time, omf = temp/dim, mean_diag_B_f = mean_diag_B_f, S_arr =S_arr, B_arr=B_arr, X_f=x_f))
}


thbef <- function(chi, phi,mu_hbef, time, dim, delta_t, N, R, m, inflation_thbef, c_thbef, L_C, U, Nu, RHO, Sigma, OBS_NOISE, start_X, start_A, start_B, X, include_obs){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_thbef)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  x_a[,1]       <- start_X
  x_f[,1]       <- start_X  
  x_en_a    <- t(mvrnorm(N,X[,1], start_A))
  S_arr <- array(NA, dim=c(time,dim,dim))
  B_arr <- array(NA, dim=c(time,dim,dim))
  B_mean <- matrix(0, ncol=dim, nrow = dim)
  temp <- c(1:time)
  mean_diag_B_f <- c(1:time)
  Q <- C%*%C
  Q_a <- mean(Sigma^2)*Q*delta_t*delta_x
  
  start_P <- start_B - mean(Sigma^2)*Q*delta_t*delta_x
  P_a <- start_P
  P_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr[1,,] <- Q_a
  P_a_arr[1,,] <- P_a
  
  
  for(i in (2:time)){
    #print(i)
    x_pe  <- apply(x_en_a,2,function(x) {
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    
    x_me <- diag(Sigma[,i])%*%C%*%matrix(rnorm(N*dim,0,sqrt(delta_x)*sqrt(delta_t)), nrow=dim)
    mean_en <- rep(rowSums(x_pe)/N, N)
    x_pe  <- matrix(mean_en, nrow = dim) + (x_pe - matrix(mean_en, nrow = dim))*inflation_thbef
    
    #x_me  <- t(mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    
    
    
    
    Q_f <- Q_a
    P_f <- P_a
    
    
    S_me     <- x_me %*% t(x_me)/N
    S_me     <- S_me * C_loc
    
    S_pe     <- (x_pe - matrix(rep(x_f[,i],N), nrow = dim))%*%t(x_pe - matrix(rep(x_f[,i],N), nrow = dim))/N
    S_pe      <- S_pe * C_loc
    
    Q_a <- (chi*Q_f + N*S_me)/(chi+N)
    P_a <- (phi * P_f + N * S_pe)/(phi + N)
    
    
    #image2D(N*(S_me+S_pe))
    Q_a_arr[i,,] <- Q_a
    P_a_arr[i,,] <- P_a
    x_en_f <- x_pe + x_me
    
    if(all(include_obs == TRUE, chi!=0, phi!=0)){
      v <- (x_obs - x_f[ind_obs,i])
      N_mat <- (P_a+Q_a)[ind_obs,ind_obs] + diag(R)
      delta_P <- P_a[,ind_obs] %*% solve(N_mat) %*% (v%*%t(v)-N_mat) %*% solve(N_mat) %*% P_a[ind_obs,]
      delta_Q <- Q_a[,ind_obs] %*% solve(N_mat) %*% (v%*%t(v)-N_mat) %*% solve(N_mat) %*% Q_a[ind_obs,]
      Q_a <- Q_a + 1/chi*delta_Q
      P_a <- P_a + 1/phi*delta_P
    }
    
    
    
    
    B_f <- Q_a + P_a
    S_arr[i,,] <- B_f
    B_f <- B_f * C_loc
    B_arr[i,,] <- B_f
    K   <- B_f[,ind_obs] %*% solve(B_f[ind_obs,ind_obs] + diag(R))
    
    
    #!!!!!!
    temp[i] <- t(x_obs-mean_en[ind_obs])%*%(x_obs-mean_en[ind_obs])
    mean_diag_B_f[i] <-mean(diag(B_f))
    B_mean <- B_mean + B_f
    
    
    
    
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i]) 
    
    #!!!!!!!!
    #mu_hbef <- 0.0001
    A     <- B_f - K%*%B_f[ind_obs,]
    E <- ((x_en_f - (matrix(rep(x_f[,i],N), nrow=dim)))/sqrt(N))
    
    #E_1 <- (t(E)%*%(E))
    E_1 <- (t(E)%*%C_loc%*%(E))
    
    #E_1 <- solve(t(E)%*%(E))
    E_1_decomp <- eigen(E_1)
    L <- E_1_decomp$vectors
    #M <- diag((E_1_decomp$values^2+mu_hbef*rep(1,N))^(-1))%*%t(L)%*%t(E)%*%(C_loc*A)%*%E%*%L
    # <- diag((E_1_decomp$values^2+mu_hbef*rep(1,N))^(-1))%*%t(L)%*%t(E)%*%(A)%*%E%*%L
    J <- t(L) %*% t(E) %*% (A*C_loc) %*% E %*% L
    #J <- t(L) %*% t(E) %*% (A) %*% E %*% L
    M <- matrix(NA,N,N)
    for(row in 1:N){
      for(col in 1:N){
        M[row,col] <- J[row,col] / (E_1_decomp$values[row]*E_1_decomp$values[col] + mu_hbef)
      }
    }
    
    W <- L%*%M%*%t(L)
    W_decomp <- eigen(W)
    L <- W_decomp$vectors
    D <- diag(sqrt(sapply(W_decomp$values, function(x) max(x,1e-16))))
    W_sqrt <- L%*%D%*%t(L)
    
    x_en_a <- x_a[,i] + E%*%W_sqrt*sqrt(N)
    #x_en_a <- apply(x_en_f,2, function(x){x + K%*%(x_obs + sqrt(R)*rnorm(dim%/%m+min(1,dim%%m),0,1)-x[ind_obs])} )
    #image2D((E%*%W_sqrt)%*%t(E%*%W_sqrt))
    #image2D(cov(t(x_en_a)))
    #image2D(A)
    #image2D(abs(A-cov(t(x_en_a))))
    
    #x_a[,i] <- rowSums(x_en_a)/N
    #print(i)
    
  }
  return(list(X_a=x_a, B_mean=B_mean/time, omf = temp/dim, mean_diag_B_f = mean_diag_B_f, S_arr =S_arr, B_arr=B_arr, X_f=x_f))
}


tenkf <- function(time, dim, delta_t, N, R, m, inflation_tenkf, c_tenkf, L_C, U, Nu, RHO, Sigma, OBS_NOISE, start_X, start_A, X, type){
  if(N != dim) {print('Alarm! N!=dim'); return(-1)}
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_tenkf)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  S_arr <- array(NA, dim=c(time,dim,dim))
  B_arr <- array(NA, dim=c(time,dim,dim))
  B_mean <- matrix(0, ncol=dim, nrow = dim)
  temp <- c(1:time)
  mean_diag_B_f <- c(1:time)
  x_a[,1]       <- start_X
  Q <- C%*%C
  x_en_a    <- t(mvrnorm(N,X[,1], start_A))
  for(i in (2:time)){
    x_en_f  <- apply(x_en_a,2,function(x) {
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    mean_en <- rep(rowSums(x_en_f)/N, N)
    x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation_tenkf
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    x_me <- (mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    M <- x_en_f + t(x_me) - (matrix(rep(x_f[,i],N), nrow=dim))
    B_f <- 1/N*M%*%t(M)
    
    S_arr[i,,] <- B_f
    B_f <- B_f * C_loc
    B_arr[i,,] <- B_f
    K   <- B_f[,ind_obs] %*% solve(B_f[ind_obs,ind_obs] + diag(R))
    
    #!!!!!!
    temp[i] <- t(x_obs-mean_en[ind_obs])%*%(x_obs-mean_en[ind_obs])
    mean_diag_B_f[i] <-mean(diag(B_f))
    B_mean <- B_mean + B_f
    #?????????????????????????????
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i]) 
    A <- B_f - K%*%B_f[ind_obs,]
    #x_en_a  <- t(mvrnorm(N,x_a[,i],A, empirical = TRUE))
    
    A_decomp <- eigen(A)
    D <- diag(sqrt(sapply(A_decomp$values, function(x) max(x,1e-16))))
    L <- A_decomp$vectors
    A_sqrt <- L%*%D%*%t(L)
    #Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
    #Samp_matr <- scale(Samp_matr, TRUE, FALSE)
    #Samp_matr <- Samp_matr %*% svd(Samp_matr, nu = 0)$v
    #Samp_matr <- scale(Samp_matr, FALSE, TRUE)
    #x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N) + A_sqrt%*%t(Samp_matr)
    
    Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
    decomp <- eigen(N*solve(Samp_matr)%*%A%*%t(solve(Samp_matr)))
    T <- decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
    Samp_matr <- Samp_matr %*% T
    1/N*Samp_matr%*%t(Samp_matr)-A
    # when N < dim we can't rescale sample matrix
    # so we take the first N columns of L%*%D
    #M <- (L %*% D ) [,1:N] * sqrt(N)
    #image2D(cov(t(M)),main='cov')
    #image2D(A, main = 'A')
    #image2D(M%*%t(M), main='A_eigencut')
    
    x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N)+ Samp_matr
  }
  return(list(X_a=x_a, B_mean=B_mean/time, omf = temp/dim, mean_diag_B_f = mean_diag_B_f, S_arr =S_arr, B_arr=B_arr, X_f=x_f))
}


thbef_ex <- function(time, dim, delta_t, N, R, m, inflation_thbef, c_thbef, L_C, U, Nu, RHO, Sigma, OBS_NOISE, start_X, start_A, X, type){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_thbef)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  S_arr <- array(NA, dim=c(time,dim,dim))
  B_arr <- array(NA, dim=c(time,dim,dim))
  B_mean <- matrix(0, ncol=dim, nrow = dim)
  temp <- c(1:time)
  mean_diag_B_f <- c(1:time)
  x_a[,1]       <- start_X
  Q <- C%*%C
  x_en_a    <- t(mvrnorm(N,X[,1], start_A))
  for(i in (2:time)){
    x_en_f  <- apply(x_en_a,2,function(x) {
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    mean_en <- rep(rowSums(x_en_f)/N, N)
    x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation_thbef
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    x_me <- (mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    M <- x_en_f + t(x_me) - (matrix(rep(x_f[,i],N), nrow=dim))
    B_f <- 1/N*M%*%t(M)
    
    S_arr[i,,] <- B_f
    B_f <- B_f * C_loc
    B_arr[i,,] <- B_f
    K   <- B_f[,ind_obs] %*% solve(B_f[ind_obs,ind_obs] + diag(R))
    
    #!!!!!!
    temp[i] <- t(x_obs-mean_en[ind_obs])%*%(x_obs-mean_en[ind_obs])
    mean_diag_B_f[i] <-mean(diag(B_f))
    B_mean <- B_mean + B_f
    #?????????????????????????????
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i]) 
    A <- B_f - K%*%B_f[ind_obs,]
    #x_en_a  <- t(mvrnorm(N,x_a[,i],A, empirical = TRUE))
    
    
    
    
    M_1 <- solve(t(M)%*%M)
    W <- M_1%*%t(M)%*%A%*%M%*%M_1
    W_decomp <- eigen(W)
    D <- diag(sqrt(sapply(W_decomp$values, function(x) max(x,1e-16))))
    L <- W_decomp$vectors
    W_sqrt <- L%*%D%*%t(L)
    E <- E%*%W_sqrt
    #A_decomp <- eigen(A)
    #D <- diag(sqrt(sapply(A_decomp$values, function(x) max(x,1e-16))))
    #L <- A_decomp$vectors
    #A_sqrt <- L%*%D%*%t(L)
    #Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
    #Samp_matr <- scale(Samp_matr, TRUE, FALSE)
    #Samp_matr <- Samp_matr %*% svd(Samp_matr, nu = 0)$v
    #Samp_matr <- scale(Samp_matr, FALSE, TRUE)
    #x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N) + A_sqrt%*%t(Samp_matr)
    
    Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
    decomp <- eigen(N*solve(Samp_matr)%*%A%*%t(solve(Samp_matr)))
    T <- decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
    Samp_matr <- Samp_matr %*% T
    1/N*Samp_matr%*%t(Samp_matr)-A
    # when N < dim we can't rescale sample matrix
    # so we take the first N columns of L%*%D
    #M <- (L %*% D ) [,1:N] * sqrt(N)
    #image2D(cov(t(M)),main='cov')
    #image2D(A, main = 'A')
    #image2D(M%*%t(M), main='A_eigencut')
    
    x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N)+ Samp_matr
  }
  return(list(X_a=x_a, B_mean=B_mean/time, omf = temp/dim, mean_diag_B_f = mean_diag_B_f, S_arr =S_arr, B_arr=B_arr, X_f=x_f))
}



hbef <- function(time, dim, delta_t, N, R, m, inflation_hbef, c_hbef, L_C, U, Nu, RHO, Sigma, OBS_NOISE, chi, phi, start_X, start_B, start_A, X){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_hbef)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  x_a[,1]       <- start_X
  Q <- C%*%C
  #Q_a <- cov(mvrnorm(N,rep(0,dim),mean(Sigma^2)*Q*delta_t*delta_x))
  Q_a <- mean(Sigma^2)*Q*delta_t*delta_x
  
  start_P <- start_B - mean(Sigma^2)*Q*delta_t*delta_x
  cond <- 1000
  vec <- eigen(start_P)$values
  decomp_mat <- eigen(start_P)$vectors
  for(i in 1:dim){
    if(vec[i] < max(vec)/cond) vec[i] <- max(vec)/cond
  }
  start_P <- (decomp_mat)%*%diag(vec)%*%t(decomp_mat)
  #P_a <- cov(mvrnorm(50*N,rep(0,dim),start_P))
  P_a <- start_P
  
  #x_en_ape    <- t(mvrnorm(N,X[,1], start_P))
  x_en_ape    <- t(mvrnorm(N,X[,1], start_A))
  
  
  S_arr <- array(NA,dim=c(time,dim,dim))
  P_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr[1,,] <- Q_a
  P_a_arr[1,,] <- P_a
  temp <- c()
  temp_me <- c()
  temp_f <- c()
  temp_B_a <- matrix(0,ncol=dim, nrow=dim)
  
  for(i in (2:time)){
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    x_pe  <- apply(x_en_ape,2,function(x) {
      #predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, diag(Sigma[,i])%*%C%*%rnorm(dim,0,1/sqrt(delta_t)))
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    x_me  <- t(mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    
    Q_f <- Q_a
    P_f <- P_a
    # P_f <- P_a *inflation_hbef
    
    B_f <- chi/(chi+N)*Q_f + phi/(phi+N)*P_f
    #x_en_f <- N/(chi+N)*x_me + N/(phi+N)*x_pe + t(mvrnorm(N,rep(0,dim),B_f))
    mean_en <- rep(rowSums(x_pe)/N, N)
    x_pe  <- matrix(mean_en, nrow = dim) + (x_pe - matrix(mean_en, nrow = dim))*inflation_hbef
    
    M <- x_pe  - (matrix(rep(x_f[,i],N), nrow=dim))
    
    S_me     <- 1/N*x_me%*%t(x_me)
    S_me    <- cov(t(x_me))
    S_me     <- S_me * C_loc
    S_pe     <- 1/N*M%*%t(M)
    S_pe      <- S_pe * C_loc
    S_arr[i,,] <- S_pe
    Q_a <- (chi*Q_f + N*S_me)/(chi+N)
    P_a <- (phi * P_f + N * S_pe)/(phi + N)
    B_a     <- P_a + Q_a
    #B_a <- P_a+ diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x
    
    Q_a_arr[i,,] <- Q_a
    P_a_arr[i,,] <- P_a
    #!!!!!!!!!!!!!!!!!!!
    #x_en_f <- x_me+x_pe
    #mean_en <- rep(rowSums(x_en_f)/N, N)
    #x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation
    #B_a <- cov(t(x_en_f))
    #B_a <- B_a * C_loc
    
    K   <- B_a[,ind_obs] %*% solve(B_a[ind_obs,ind_obs] + diag(R))
    A <- B_a - K%*%B_a[ind_obs,]
    #A             <- B_a - K %*% t(B_a[,ind_obs]) - B_a[,ind_obs]%*% t(K) +K %*% B_a[ind_obs,ind_obs] %*% t(K) + K%*%diag(R)%*%t(K)
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i])
    
    A_decomp <- eigen(A)
    D <- diag(sqrt(sapply(A_decomp$values, function(x) max(x,1e-16))))
    L <- A_decomp$vectors
    A_sqrt <- L%*%D%*%t(L)
    Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
    Samp_matr <- A_sqrt %*% t(Samp_matr)
    #decomp <- eigen(N*solve(Samp_matr)%*%A%*%t(solve(Samp_matr)))
    #T <- decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
    #Samp_matr <- Samp_matr %*% T
    #1/N*Samp_matr%*%t(Samp_matr)-A
    
    x_en_ape <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N)+ Samp_matr
    
    
    
    
    M <- diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x
    temp <- c(temp,base::norm((M-Q_a), type = 'F'))
    temp_me <- c(temp_me,base::norm(M-S_me, type = 'F'))
    temp_f <- c(temp_f,base::norm(M-Q_f, type = 'F'))
    temp_B_a <- temp_B_a+B_a
  }
  
  return(list(Q_a=Q_a, P_a=P_a, S_me=S_me,S_pe=S_pe, S_arr=S_arr, Q_a_arr = Q_a_arr,P_a_arr=P_a_arr, Q_f=Q_f,P_f=P_f,norm_Q_a=temp,norm_Q_f = temp_f,norm_S_me = temp_me, X_a=x_a, X_f=x_f))
  
}

shbef <- function(time, dim, delta_t, N, R, m, inflation_shbef, c_shbef, L_C, U, Nu, RHO, Sigma, OBS_NOISE, chi, phi, start_X, start_B, start_A, X){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_shbef)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  x_a[,1]       <- start_X
  Q <- C%*%C
  Q_a <- mean(Sigma^2)*Q*delta_t*delta_x
  P_a <- start_P
  
  #x_en_ape    <- t(mvrnorm(N,X[,1], start_P))
  x_en_ape    <- t(mvrnorm(N,X[,1], start_A))
  
  
  S_arr <- array(NA,dim=c(time,dim,dim))
  P_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr[1,,] <- Q_a
  P_a_arr[1,,] <- P_a
  temp <- c()
  temp_me <- c()
  temp_f <- c()
  temp_B_a <- matrix(0,ncol=dim, nrow=dim)
  
  for(i in (2:time)){
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    x_pe  <- apply(x_en_ape,2,function(x) {
      #predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, diag(Sigma[,i])%*%C%*%rnorm(dim,0,1/sqrt(delta_t)))
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    x_me  <- t(mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    
    Q_f <- Q_a
    #P_f <- P_a
    #!!!!!!!!
    P_f <- start_A
    
    
    # P_f <- P_a *inflation_hbef
    
    #B_f <- chi/(chi+N)*Q_f + phi/(phi+N)*P_f
    #x_en_f <- N/(chi+N)*x_me + N/(phi+N)*x_pe + t(mvrnorm(N,rep(0,dim),B_f))
    mean_en <- rep(rowSums(x_pe)/N, N)
    x_pe  <- matrix(mean_en, nrow = dim) + (x_pe - matrix(mean_en, nrow = dim))*inflation_shbef
    
    M <- x_pe  - (matrix(rep(x_f[,i],N), nrow=dim))
    
    S_me     <- 1/N*x_me%*%t(x_me)
    #S_me    <- cov(t(x_me))
    S_me     <- S_me * C_loc
    S_pe     <- 1/N*M%*%t(M)
    S_pe      <- S_pe * C_loc
    S_arr[i,,] <- S_pe
    Q_a <- (chi*Q_f + N*S_me)/(chi+N)
    P_a <- (phi * P_f + N * S_pe)/(phi + N)
    B_a     <- P_a + Q_a
    
    
    Q_a_arr[i,,] <- Q_a
    P_a_arr[i,,] <- P_a
    #!!!!!!!!!!!!!!!!!!!
    ind <- sample(1:N, floor(N^2/(phi+N)))
    if(length(ind) == N){      
      x_pe_1 <- x_pe[,ind]
      x_en_f <- x_me+x_pe_1 + matrix(rep(x_f[,i],N), nrow = dim)
    }
    if(length(ind) == 0){
      x_pe_2 <- t(mvrnorm(N - length(ind), rep(0,dim), P_f))
      x_en_f <- x_me+x_pe_2 + matrix(rep(x_f[,i],N), nrow = dim)
    }
    if(length(ind)>0 && length(ind)<N){
      x_pe_1 <- x_pe[,ind]
      x_pe_2 <- t(mvrnorm(N - length(ind), rep(0,dim), P_f))
      x_en_f <- x_me+cbind(x_pe_1, x_pe_2) + matrix(rep(x_f[,i],N), nrow = dim)
      #x_en_f <- x_me+cbind(x_pe_1, matrix(x_pe_2, nrow=dim)) + matrix(rep(x_f[,i],N), nrow = dim)
    }
    
    
    
    
    #mean_en <- rep(rowSums(x_en_f)/N, N)
    #x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation
    #B_a <- cov(t(x_en_f))
    #B_a <- B_a * C_loc
    
    K   <- B_a[,ind_obs] %*% solve(B_a[ind_obs,ind_obs] + diag(R))
    A <- B_a - K%*%B_a[ind_obs,]
    #A             <- B_a - K %*% t(B_a[,ind_obs]) - B_a[,ind_obs]%*% t(K) +K %*% B_a[ind_obs,ind_obs] %*% t(K) + K%*%diag(R)%*%t(K)
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i])
    
    
    x_en_ape <- apply(x_en_f,2, function(x){x + K%*%(x_obs + sqrt(R)*rnorm(dim%/%m+min(1,dim%%m),0,1)-x[ind_obs])} )
    
    
    
    M <- diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x
    temp <- c(temp,base::norm((M-Q_a), type = 'F'))
    temp_me <- c(temp_me,base::norm(M-S_me, type = 'F'))
    temp_f <- c(temp_f,base::norm(M-Q_f, type = 'F'))
    temp_B_a <- temp_B_a+B_a
  }
  
  return(list(Q_a=Q_a, P_a=P_a, S_me=S_me,S_pe=S_pe, S_arr=S_arr, Q_a_arr = Q_a_arr,P_a_arr=P_a_arr, Q_f=Q_f,P_f=P_f,norm_Q_a=temp,norm_Q_f = temp_f,norm_S_me = temp_me, X_a=x_a, X_f=x_f))
  
}








thbef_ex <- function(time, dim, delta_t, N, R, m, inflation_thbef, c_thbef, L_C, U, Nu, RHO, Sigma, OBS_NOISE, chi, phi, start_X, start_B, start_A, X){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_thbef)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  x_a[,1]       <- start_X
  Q <- C%*%C
  Q_a <- mean(Sigma^2)*Q*delta_t*delta_x
  
  #############
  #start_P <- start_B - mean(Sigma^2)*Q*delta_t*delta_x
  #cond <- 1000
  #vec <- eigen(start_P)$values
  #decomp_mat <- eigen(start_P)$vectors
  #for(i in 1:dim){
  #  if(vec[i] < max(vec)/cond) vec[i] <- max(vec)/cond
  #}
  #start_P <- (decomp_mat)%*%diag(vec)%*%t(decomp_mat)
  ##P_a <- cov(mvrnorm(50*N,rep(0,dim),start_P))
  #P_a <- start_P
  ##x_en_ape    <- t(mvrnorm(N,X[,1], start_P))
  #x_en_ape    <- t(mvrnorm(N,X[,1], start_A))
  ###############
  start_P <- start_B - mean(Sigma^2)*Q*delta_t*delta_x
  P_a <- start_P
  
  x_en_ape <- t(mvrnorm(N,X[,1], start_A))
  
  
  P_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr <- array(NA,dim=c(time,dim,dim))
  Q_a_arr[1,,] <- Q_a
  P_a_arr[1,,] <- P_a
  temp <- c()
  temp_me <- c()
  temp_f <- c()
  temp_B_a <- matrix(0,ncol=dim, nrow=dim)
  norm_f_en <- c()
  norm_w_en <- c()
  norm_u_en <- c()
  norm_v_en <- c()
  for(i in (2:time)){
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    x_pe  <- apply(x_en_ape,2,function(x) {
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    mean_en <- rep(rowSums(x_pe)/N, N)
    x_pe  <- matrix(mean_en, nrow = dim) + (x_pe - matrix(mean_en, nrow = dim))*inflation_thbef
    
    x_me  <- t(mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    
    Q_f <- Q_a
    #!!!!!!!!!!!!!!!!!
    P_f <- start_P
    #P_f <- P_a
    
    B_f <- chi/(chi+N)*Q_f + phi/(phi+N)*P_f
    #!!!!!!!!!!!!!!!!!!!!!!
    #B_f <- chi/(chi+N)*mean(Sigma^2)*Q*delta_t*delta_x + phi/(phi+N)*start_P
    
    
    S_me     <- x_me %*% t(x_me)/N
    S_me     <- S_me * C_loc
    
    S_pe     <- (x_pe - matrix(rep(x_f[,i],N), nrow = dim))%*%t(x_pe - matrix(rep(x_f[,i],N), nrow = dim))/N
    #S_pe    <- cov(t(x_pe))
    S_pe      <- S_pe * C_loc
    
    Q_a <- (chi*Q_f + N*S_me)/(chi+N)
    P_a <- (phi * P_f + N * S_pe)/(phi + N)
    
    #image2D(N*(S_me+S_pe))
    Q_a_arr[i,,] <- Q_a
    P_a_arr[i,,] <- P_a
    
    B_a <- P_a + Q_a  
    
    
    
    B_a <- B_a * C_loc
    
    #######
    #B_a <- P_a+Q_a
    K   <- B_a[,ind_obs] %*% solve(B_a[ind_obs,ind_obs] + diag(R))
    
    T <- diag(rep(1,N)) + t(E[ind_obs,]) %*% diag(R^(-1)) %*% E[ind_obs,]
    T_decomp <- eigen(T)
    T <- T_decomp$vectors %*% diag((T_decomp$values)^(-1/2)) %*% t(T_decomp$vectors)
    
    
    A <- B_a - K%*%B_a[ind_obs,]
    #A             <- B_a - K %*% t(B_a[,ind_obs]) - B_a[,ind_obs]%*% t(K) +K %*% B_a[ind_obs,ind_obs] %*% t(K) + K%*%diag(R)%*%t(K)
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i])
    x_en_ape <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N)+ sqrt(N)*E %*% T
    #print(B_a[1,1])
    
    
    
    M <- diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x
    temp <- c(temp,base::norm((M-Q_a), type = 'F'))
    temp_me <- c(temp_me,base::norm(M-S_me, type = 'F'))
    temp_f <- c(temp_f,base::norm(M-Q_f, type = 'F'))
    temp_B_a <- temp_B_a+B_a
  }
  
  return(list(Q_a=Q_a, P_a=P_a, S_me=S_me,S_pe=S_pe, Q_a_arr = Q_a_arr,P_a_arr=P_a_arr, Q_f=Q_f,P_f=P_f,
              norm_Q_a=temp,norm_Q_f = temp_f,norm_S_me = temp_me, X_a=x_a, X_f=x_f,
              norm_f_en=norm_f_en))
  
}



etkf <- function(time, dim, delta_t, N, R, m, inflation_enkf, c_enkf, L_C, U, Nu, RHO, Sigma, OBS_NOISE, start_X, start_A, X, type){
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  ind_obs       <- seq(1,dim,m)
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc <- create_loc_matrix(dim, delta_x, c_enkf)
  x_a <- matrix(NA, ncol = time, nrow = dim)
  x_f <- matrix(NA, ncol = time, nrow = dim)
  S_arr <- array(NA, dim=c(time,dim,dim))
  B_arr <- array(NA, dim=c(time,dim,dim))
  B_mean <- matrix(0, ncol=dim, nrow = dim)
  temp <- c(1:time)
  mean_diag_B_f <- c(1:time)
  x_a[,1]       <- start_X
  Q <- C%*%C
  #x_en_ape    <- t(mvrnorm(N,X[,1], start_A))
  x_en_a    <- t(mvrnorm(N,X[,1], start_A))
  for(i in (2:time)){
    x_en_f  <- apply(x_en_a,2,function(x) {
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    }) 
    
    x_f[,i] <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x_a[,i-1], rep(0,dim))
    mean_en <- rep(rowSums(x_en_f)/N, N)
    x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation_enkf
    x_obs <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]
    x_me <- (mvrnorm(N,rep(0,dim),diag(Sigma[,i])%*%Q%*%diag(Sigma[,i])*delta_t*delta_x))
    
    if(type == 'ensemble'){
      mean_en <- rep(rowSums(x_en_f + t(x_me))/N, N)
      E <- (x_en_f + t(x_me) - (matrix(mean_en, nrow=dim)))/sqrt(N-1)
    }
    
    if(type == 'deterministic'){
      mean_en <- rep(x_f[,i], N)
      E <- (x_en_f + t(x_me) - (matrix(mean_en, nrow=dim)))/sqrt(N)
    }
    
    B_f  <- E %*% t(E)
    S_arr[i,,] <- B_f
    B_f <- B_f * C_loc
    B_arr[i,,] <- B_f
    K   <- B_f[,ind_obs] %*% solve(B_f[ind_obs,ind_obs] + diag(R))
    A <- B_f - K%*%B_f[ind_obs,]
    
    T <- diag(rep(1,N)) + t(E[ind_obs,]) %*% diag(R^(-1)) %*% E[ind_obs,]
    T_decomp <- eigen(T)
    T <- T_decomp$vectors %*% diag((T_decomp$values)^(-1/2)) %*% t(T_decomp$vectors)
    
    #!!!!!!
    temp[i] <- t(x_obs-mean_en[ind_obs])%*%(x_obs-mean_en[ind_obs])
    mean_diag_B_f[i] <-mean(diag(B_f))
    B_mean <- B_mean + B_f
    #?????????????????????????????
    
    x_a[,i] <- x_f[,i] + K%*%(x_obs - x_f[ind_obs,i]) 
    A <- B_f - K%*%B_f[ind_obs,]
    
    if(type == 'deterministic'){
      x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N)+ sqrt(N)*E %*% T
    }
    if(type == 'ensemble'){
      x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N)+ sqrt(N-1)*E %*% T
    }
  }
  return(list(X_a=x_a, B_mean=B_mean/time, omf = temp/dim, mean_diag_B_f = mean_diag_B_f, S_arr =S_arr, B_arr=B_arr, X_f=x_f))
}

kalman_test <- function(time, dim, delta_t, U, RHO, Sigma, Nu, L_C, R, m, OBS_NOISE, start_X, start_A, N, phi, c_kalman, inflation_kalman){
  X_a           <- matrix(NA, nrow = dim, ncol = time)
  X_a_test           <- matrix(NA, nrow = dim, ncol = time)
  X_f           <- matrix(NA, nrow = dim, ncol = time)
  B_arr         <- array(NA,c(time,dim,dim))
  A_arr         <- array(NA,c(time,dim,dim))
  A_arr_mod     <- array(NA,c(time,dim,dim))
  X_a[,1]       <- start_X
  X_a_test[,1]       <- start_X
  Re            <- 6370000  
  delta_x       <- 2*pi*Re/dim
  C             <- create_cov_matrix(dim, rep(1,dim), L_C)
  C_loc         <- create_loc_matrix(dim, delta_x, c_hbef)
  
  A_arr[1,,]    <- start_A
  PHI           <- apply(start_A,2, function(x) predict(dim, delta_x, delta_t, U[,2], Nu[,2], RHO[,2], x, rep(0,dim))) 
  B             <- apply(t(PHI),2, function(x) predict(dim, delta_x, delta_t, U[,2], Nu[,2], RHO[,2], x, rep(0,dim))) 
  B             <- B + diag(Sigma[,2])%*%C%*%C%*%diag(Sigma[,2])*delta_x*delta_t
  x_en_a        <- t(mvrnorm(N,X[,1], start_A))
  ind_obs       <- seq(1,dim,m)
  temp_pf <- c()
  temp_pa <- c()
  temp_s <- c()
  for(i in 2:time){
    x_en_f  <- apply(x_en_a,2,function(x) {
      predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], x, rep(0,dim))
    })
    mean_en <- rep(rowSums(x_en_f)/N, N)
    x_en_f  <- matrix(mean_en, nrow = dim) + (x_en_f - matrix(mean_en, nrow = dim))*inflation_kalman
    
    B_arr[i,,]    <- B
    
    y             <- X[ind_obs,i] + sqrt(R)*OBS_NOISE[,i]              
    BH            <- B[,ind_obs]
    HBH           <- B[ind_obs, ind_obs]
    X_f[,i]       <- predict(dim, delta_x, delta_t, U[,i], Nu[,i], RHO[,i], X_a[,i-1], rep(0,dim))  
    z             <- solve(HBH + diag(R), y-X_f[ind_obs,i])
    X_a[,i]       <- X_f[,i] + BH %*% z
    
    if(i==2){
      P_f           <- B_arr[i,,] - diag(Sigma[,i])%*%C%*%C%*%diag(Sigma[,i])*delta_x*delta_t
    }else{
      P_f           <- B_arr[i-1,,] - diag(Sigma[,i-1])%*%C%*%C%*%diag(Sigma[,i-1])*delta_x*delta_t
    }
    
    S             <- cov(t(x_en_f))
    B_test        <- (phi*P_f + N*S*C_loc)/(phi+N) + diag(Sigma[,i])%*%C%*%C%*%diag(Sigma[,i])*delta_x*delta_t
    BH_test            <- B_test[,ind_obs]
    HBH_test           <- B_test[ind_obs, ind_obs]
    X_a_test[,i]  <- X_f[,i] + BH_test %*% z
    
    temp_pf <- c(temp_pf,base::norm(P_f - B + diag(Sigma[,i])%*%C%*%C%*%diag(Sigma[,i])*delta_x*delta_t, type = 'F'))
    temp_s <- c(temp_pf,base::norm(S*C_loc - B + diag(Sigma[,i])%*%C%*%C%*%diag(Sigma[,i])*delta_x*delta_t, type = 'F'))
    temp_pa <- c(temp_pf,base::norm(B_test - B, type = 'F'))
    
    temp_pf <- sqrt(mean((diag(P_f - B + diag(Sigma[,i])%*%C%*%C%*%diag(Sigma[,i])*delta_x*delta_t))^2))
    temp_s <- sqrt(mean((diag(S*C_loc - B + diag(Sigma[,i])%*%C%*%C%*%diag(Sigma[,i])*delta_x*delta_t))^2))
    temp_pa <- sqrt(mean((diag(B_test - B))^2))
    
    #A             <- B - BH %*% solve(HBH+diag(R)) %*% t(BH)
    K  <- BH %*% solve(HBH + diag(R))
    #A             <- (diag(rep(1,length(ind_obs))) - K[,ind_obs])%*% B %*% t(diag(rep(1,length(ind_obs))) - K[,ind_obs]) + K%*%diag(R)%*%t(K)
    #A             <- (diag(rep(1,dim)) - K,ind_obs])%*% B %*% t(diag(rep(1,dim)) - K[,ind_obs]) + K%*%diag(R)%*%t(K)
    A     <- B - K%*%B[ind_obs,]
    #A             <- B - K %*% t(BH) - BH%*% t(K) +K %*% HBH %*% t(K) + K%*%diag(R)%*%t(K)
    A_arr[i,,]    <- A
    A_decomp <- eigen(A)
    D <- diag(sqrt(sapply(A_decomp$values, function(x) max(x,1e-16))))
    L <- A_decomp$vectors
    A_sqrt <- L%*%D%*%t(L)
    #Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
    #Samp_matr <- scale(Samp_matr, TRUE, FALSE)
    #Samp_matr <- Samp_matr %*% svd(Samp_matr, nu = 0)$v
    #Samp_matr <- scale(Samp_matr, FALSE, TRUE)
    #x_en_a <- matrix(rep(x_a[,i],N),nrow=dim,ncol=N) + A_sqrt%*%t(Samp_matr)
    
    Samp_matr <- t(matrix(rnorm(N*dim), nrow=dim, ncol = N))
    decomp <- eigen(N*solve(Samp_matr)%*%A%*%t(solve(Samp_matr)))
    T <- decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
    Samp_matr <- Samp_matr %*% T
    1/N*Samp_matr%*%t(Samp_matr)-A
    # when N < dim we can't rescale sample matrix
    # so we take the first N columns of L%*%D
    #M <- (L %*% D ) [,1:N] * sqrt(N)
    #image2D(cov(t(M)),main='cov')
    #image2D(A, main = 'A')
    #image2D(M%*%t(M), main='A_eigencut')
    
    x_en_a <- matrix(rep(X_a[,i],N),nrow=dim,ncol=N)+ Samp_matr
    #x_en_a  <- t(mvrnorm(N,X_a[,i],A))
    
    
    if(i < time){
      PHI           <- apply(A,2, function(x) predict(dim, delta_x, delta_t, U[,i+1], Nu[,i+1], RHO[,i+1], x, rep(0,dim))) 
      B             <- apply(t(PHI),2, function(x) predict(dim, delta_x, delta_t, U[,i+1], Nu[,i+1], RHO[,i+1], x, rep(0,dim))) 
      B             <- B + diag(Sigma[,i+1])%*%C%*%C%*%diag(Sigma[,i+1])*delta_x*delta_t
      #!!!
      #B             <- A + diag(Sigma[,i+1])%*%C%*%C%*%diag(Sigma[,i+1])
      #B             <- t(B)
    }
    
    #eig_decomp <- eigen(B)
    #print(sum(diag(B)))
    #print(c(i,eig_decomp$values[30]))
  }
  return(list(X_a=X_a, X_f=X_f, B_arr = B_arr, A_arr = A_arr, A_arr_mod = A_arr_mod, X_a_test=X_a_test,
              temp_pf=temp_pf,temp_pa=temp_pa,temp_s=temp_s))
}




create_loc_matrix <- function(dim, delta_x,c){
  M <- matrix(NA, ncol = dim, nrow = dim)
  Re            <- 6370000
  delta_x       <- 2*pi*Re/dim
  for(i in 1:dim){
    for(j in 1:dim){
      if(abs(i-j)<dim/2){
        z <- abs(i-j)
      }else{
        z <- dim - abs(i-j)
      }
      z <- 2*sin(z*delta_x/(Re*2))*Re
      #M[i,j] <- z
      #z <- z*delta_x
      #      M[i,j] <- exp(-(z/c)^2)
      if(z>=0 && z<=c) M[i,j] <- -1/4*(z/c)^5+1/2*(z/c)^4+5/8*(z/c)^3-5/3*(z/c)^2+1
      if(z>=c && z<=2*c) M[i,j] <- 1/12*(z/c)^5-1/2*(z/c)^4+5/8*(z/c)^3+5/3*(z/c)^2-5*(z/c)+4-2/3*c/z
      if(z>=2*c) M[i,j] <- 0
    }
  }
  return(M)
}

start_value <- function(a, sigma_mean, RHO_mean, Nu_mean, Re, dim){
  bm <- a^2*sigma_mean^2/2*1/(RHO_mean + Nu_mean/(Re^2) * (seq(-dim/2,dim/2,1))^2)
  z_re <- rep(0,dim+1)
  z_im <- rep(0,dim+1)
  z_re[dim/2+1] <- rnorm(1,0,sqrt(bm[dim/2+1]))
  z_re[1] <- rnorm(1,0,sqrt(bm[1]))
  for(i in 2:(dim/2)){
    z_re[i] <- rnorm(1,0, sqrt(bm[i]/2))
    z_im[i] <- rnorm(1,0, sqrt(bm[i]/2))
  }
  z_re[(dim/2+2):(dim)] <- rev(z_re[2:(dim/2)])
  z_im[(dim/2+2):(dim)] <- -rev(z_im[2:(dim/2)])
  
  start_X_function <- function(x){
    sum((z_re+1i*z_im) %*% exp(1i*(seq(-dim/2,dim/2,1)) * x / Re ))
  }
  return(Re(sapply((1:dim)*delta_x, start_X_function)))
}


logit_gen_inv <- function(x,b){
  return((1+exp(b))*exp(x-b)/(1+exp(x-b)))
}

logit_gen <- function(x,b){
  return(b-log((1+exp(b))/x-1))
}

rmse <- function(x,y){
  temp <- as.vector(x-y)
  temp <- temp^2
  res <- sqrt(mean(temp))
  return(res)
}