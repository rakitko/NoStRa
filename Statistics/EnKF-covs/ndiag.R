ndiag = function(A) {
  
  # portion of the non-diag elements in the Frobenius norm:
  # fA=sum_i (A_ij^2)
  # fdiag=sum_i A_ii^2
  # fnondiag=fA - fdiag
  # ndiag=sqrt(fnondiag / fA)
  
  fA=sum(abs(A)^2)
  fdiag=sum(diag(abs(A)^2))
  fnondiag=fA - fdiag
  ndiag=sqrt(fnondiag / fA)
  
  return(ndiag)
}
