mx_fft_Rspe_rearrange = function(A) {
  
  # For plotting purposes, rearrange the the ordering of rows & cols
  # of a cplx n by n mx A
  # computed with the R's fft, so that the wvns go as follows:
  # f(n/2),f(-n/2 +1), ..., f(-1), f(0), f(1),..., f(n/2)
  #
  # NB: n needs to be even.
  #
  # M Tsy 3 Jul 2017
  
  source("fft_Rspe_rearrange.R")
  
  n = dim(A)[1]
  np1=n+1
  A_symm=matrix(0+0i, nrow=np1, ncol=np1)
  
  for (m in 1:n){ # loop over rows of A
    A_symm[m,]=fft_Rspe_rearrange(A[m,])
  }
  
  for (m in 1:np1){ # loop over cols of A_symm
    A_symm[,m]=fft_Rspe_rearrange(A_symm[(1:n),m])
  }
  
  return(A_symm)
}
  
  
  
  