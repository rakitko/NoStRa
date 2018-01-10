lorenz05_2 <- function(ntime, n, X1, dt, F_Lorenz, J, sd_noise){
  
  # Lorenz-2005 Model-II system
  # Integrate using RK2
  #----------
  # Lorenz-05_2:
  # 
  # Set K<<n. Let K be Odd. Then, with J=(K-1)/2, the model is:
  # 
  # dx/dt=f(x), where x and f are the n-vectors
  #***********************************************************************
  # f(x)_i = - W(i-2*K) * W(i-K)  +  1/K * sum_{j=-J}^J  W(i-K+j)*X(i+K+j),  
  # where
  # W(i) = 1/K * sum_{j=-J}^J  X(i-j)
  #*********************************************************************** 
  # (i is the space index, all spatial indices are modulo n)
  # 
  # Time scale: dt=0.05 corresponds to 6h in the atmosphere (Lorenz&Emanuel 1998)
  #---------
  # RK2 (2nd-order Runge-Kutta, the half-step version):
  # 
  # (1) x[k+1/2]=x[k] + dt/2 * f(x[k])
  # (2) x[k+1] = x[k] + dt * f(x[k+1/2])
  # (k is the time index)
  #
  # Arguments:
  #
  # ntime is the number of time steps (inclu the starting step 1)
  # n is the dim-ty of the state vector X
  # X1=X[,1] is the initial condition
  # dt is the time step (unitless, dt=0.05 corresponds to 6h in the atmosphere)
  # F_Lorenz is the constant
  # J (or, equiv-ly, K) - integer, defines the spatial length scale (in the range 0...3)
  # sd_noise - St.dev of the white noise added to the model field each time step (=0 normally)
  # 
  # Return the whole history X[:,:]
  #  
  # M Tsyrulnikov
  # Dec 2017
  
  X=matrix(0, nrow=n, ncol=ntime)
  X[,1] = X1  # ini cond
  
  K=2*J+1
  
  noise=rnorm(n, 0, sd_noise)
  
  ## classical Runge-Kutta 2nd order, the half-step version
 
  for(k in 2:ntime){
    km1=k-1
    x=X[,km1] # starting x at current time step
    x_half = x + (dt/2)*Oprt_Lorenz05(x,      n, F_Lorenz, J,K)
    x_full = x +  dt   *Oprt_Lorenz05(x_half, n, F_Lorenz, J,K)
    X[,k] = x_full  + noise
  }
  
  return(X)
}


Oprt_Lorenz05 <- function(x, n, F_Lorenz, J,K)  {
  
  # The rhs operator of the Lorenz-05 Model-II system
  # x is the n-vector of the current state of the system
  # F_Lorenz is the constant.
  # J,K are two integer constants
  
  # First, compute the spatially smoothed X_i, i.e. W_i:
  # W(i) = 1/K * sum_{j=-J}^J  X(i-j)
  
  W=c(1:n)
  Oprt_Lorenz05=c(1:n)

  for(i in 1:n){
    ind=((i-J-1):(i+J-1)) %%n +1
    W[i]=mean(x[ind])
  }
  
  # Compute the rhs:
  # Oprt_Lorenz05[i]= - W(i-2*K) * W(i-K)  +  1/K * sum_{j=-J}^J  W(i-K+j)*X(i+K+j)

  for(i in 1:n){
    ind1=(i-(2*K)-1) %%n +1  # i-2K modulo n
    ind2=(i-   K -1) %%n +1  # i-K  modulo n
    
    ind3=((i-K-J-1):(i-K+J-1)) %%n +1  # i-K+j,  j \in [-J,J]
    ind4=((i+K-J-1):(i+K+J-1)) %%n +1  # i+K+j,  j \in [-J,J]
    
    Oprt_Lorenz05[i]= - W[ind1]*W[ind2] +  mean( W[ind3] * x[ind4]) - x[i] + F_Lorenz
  }
  
  return(Oprt_Lorenz05)
}

