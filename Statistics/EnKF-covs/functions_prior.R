create_cvm <- function(n,L){
  Re=6.37*10^6 # m
  cvm  <- matrix(ncol = n, nrow = n)
  ds=2*pi*Re/n
  
  for(i in (1:n)) {
    for(j in (1:n)) {
      d=abs(i-j)
      if ((d> floor(n/2))) {d <- n-d}
      rho <- ds * d
      cvm[i,j] <- (1 + rho/L)*exp(-rho/L)
      #cvm[i,j] <- exp(-rho/L)
    }
  }
  return(cvm)
}



negLaplMx <- function(n){
  Re=6.37*10^6 # m
  negLaplMx  <- matrix(0, ncol = n, nrow = n)
  ds=2*pi/n # *Re
  a <- -1 #/ds^2
  a2 <- 2 #/ds^2
  
  for(i in (1:n)) {
    im1 <- i-1
    ip1 <- i+1
    if ((i==1)) {im1 <- n}
    if ((i==n)) {ip1 <- 1}
    negLaplMx[i,i] <- a2
    negLaplMx[i,im1] <- a
    negLaplMx[i,ip1] <- a
  }
  return(negLaplMx)
}




SkewRegNegLaplMx <- function(n, ttheta, epsilon){
  
    # ttheta[1:n] - pointwise degree of skewness (left-right contrast)
  #         (does not violate diagonal dominance)
  # epsilon - small positive rregularizing number to be added to the  main diagonal
  #           (ensures diagonla dominance)
  
  SkewRegNegLaplMx  <- matrix(0, ncol = n, nrow = n)
  
  ds   <- 2*pi/n 
  mult <- 1/ds^2
  diag_entry <- (2 + epsilon) *mult
  
  for(i in (1:n)) {
    im1 <- i-1
    ip1 <- i+1
    if ((i==1)) {im1 <- n}
    if ((i==n)) {ip1 <- 1}
    SkewRegNegLaplMx[i,i]   <- diag_entry
    SkewRegNegLaplMx[i,im1] <- (-1 - ttheta[i])  *mult
    SkewRegNegLaplMx[i,ip1] <- (-1 + ttheta[i])  *mult
  }
  return(SkewRegNegLaplMx)
}



create_cvm_loc <- function(n,L,Lloc){
  
  # Create CVM & Localiz. mx.
  
  cvm    <- matrix(ncol = n, nrow = n)
  LocMx  <- matrix(ncol = n, nrow = n)
  
  Re=6.37*10^6 # m
  Lrad=L / Re
  Llocrad=Lloc / Re
  ds=2*pi/n
  
  for(i in (1:n)) {
    for(j in (1:n)) {
      d=abs(i-j)
      if ((d> floor(n/2))) {d <- n-d}
      rho <- ds * d
      r <- 2*sin(rho/2)  # chordal distance
      cvm[i,j]   <- (1 + r/Lrad)*exp(-r/Lrad)
      LocMx[i,j] <- exp(-0.5*(r/Llocrad)^2)
      #cvm[i,j] <- exp(-rho/L)
    }
  }
  return(list("cvm" = cvm, "LocMx" = LocMx))
}





extract_cyclic_superdiagonals <- function(A,p){
  
  # extract the main diag & p super-diagonals from mx A.
  # Do this CYCLICALLY, i.e. each super-diagonal has the same length n
  # and continues as on the circle
  # p is the nu of super-diags
  # returns digonals (an n*(p+1) mx) including the main diagonal.
  
  n <- dim(A)[1]  # dim-ty
  sdiagonals      <- matrix(nrow=n, ncol=p+1)
  superdiagonals  <- matrix(nrow=n, ncol=p)
  
  # Extract the main diagonal
  
  #diagonals[,1] <- diag(A)
  maindiag <- diag(A)
  
  # Extract the super-diagonals
  
  for (i in 1:p){
    superdiagonals[1:(n-i),   i] <- 
      diag(matrix(as.vector(A[-(n-i+1):-n, -1:-i]), nrow=n-i, ncol=n-i))  # diagonals
    superdiagonals[(n-i+1):n, i] <- 
      diag(matrix(as.vector(A[-1:-(n-i), -(i+1):-n]), nrow=i, ncol=i))  # corners
  }
  
  sdiagonals <- cbind(maindiag, superdiagonals)
  
  return(list("sdiagonals"=sdiagonals, "superdiagonals"=superdiagonals))
}





construct_symm_mx_from_cyclic_superdiagonals <- function(sdiagonals){
  
  # From the main diag & p super-diagonals (all put in 'sdiagonals')
  #  of a symmetric mx A, build the whole A.
  # NB: each super-diagonal has the same length n and continues as on the circle.
  # returns A.
  
  n <- dim(sdiagonals)[1]       # dim-ty
  p <- dim(sdiagonals)[2] -1    # nu of superdiags
  
  A <- matrix(0, nrow=n, ncol=n) # initialize by 0's

  # Fill in the super-diagonals
  
  for (i in 1:p){
    a <- matrix(0, nrow=n-i, ncol=n-i) # submx with the MAIN part of the
                                       # i-th A-superdiag on its main diag
    diag(a) <- sdiagonals[1:(n-i),   i+1]  
    A[-(n-i+1):-n, -1:-i] <- A[-(n-i+1):-n, -1:-i] + a # add the MAIN part of the diagonal
    
    b <- matrix(0, nrow=i, ncol=i)     # submx with the CORNER part of the
                                       # i-th A-superdiag on its main diag
    diag(b) <- sdiagonals[(n-i+1):n, i+1]
    A[-1:-(n-i), -(i+1):-n] <- A[-1:-(n-i), -(i+1):-n] + b       # corners
  }
  
  # Fill in the subdiagonals by symmetry
  
  A <- A + t(A)
  
  # Finally, fill in the main diagonal
  
  diag(A) <- sdiagonals[,1]
  
  return(A)
}





Diag_flt <- function(S){
  
  n <- dim(S)[1]  # dim-ty
  
  S_flt <- matrix(ncol = n, nrow = n)
  
  Inm1 <- diag(n-1) # the Id mx of size (n-1)*(n-1)
  
  Ushift   <- matrix(ncol = n, nrow = n)
  Dshift   <- matrix(ncol = n, nrow = n)
  
  colnm1 <- cbind(rep(0, n-1))
  rown   <- rbind(rep(0, n))
  
  # Mx with one sub-diagonal: shift DOWN 
  
  Dshift   <- cbind(Inm1, colnm1)
  Dshift   <- rbind(rown, Dshift)
  Dshift[1,n] <- 1
  
  # Mx with one super-diagonal: shift UP 
  
  Ushift   <- cbind(colnm1, Inm1)
  Ushift   <- rbind(Ushift, rown)
  Ushift[n,1] <- 1
  
  # S_flt = sum [ w(l) S(l) ],
  # where the weights w(l) sum up to 1,
  # l=0,1,.., lmax
  #
  # S(l) = Ushift^l * S *  Ushift^l +
  #        Dshift^l * S *  Dshift^l.

  #w <- c(1/2, 1/4)
  w <- c(1/3, 1/3)
  
  S_flt <- w[1] * S +
           w[2] * Ushift %*% S %*% t(Ushift) +
           w[2] * Dshift %*% S %*% t(Dshift)
 
  return(S_flt)
}




w2fields_S1 <- function(n, m, eps){
  
# weights w(x) on the grid with n points on S1 - for the mdl:
#
# xi(x) := w(x) xi_1(x)  + sqrt(1-w^2(x)) xi_2(x),
# where xi_1(x) and xi_2(x) are homogeneous and
#       
#  w(x) := (1 + eps + sin(mx)) / (2 + 2*eps),
#
# m - integer, eps>0. 
# Small eps => w_min close to 0, w_max close to 1.
  
  w <- c(1:n)
  
  for(i in (1:n)) {
    x <- 2*pi/n * (i-1)
    w[i] <- (1 + eps + sin(m*x)) / (2 + 2*eps)
  }
  return(w)
}


