Gaspari <- function(r,c){

  # 5-th order positive definite picewise rational function, 
  # Gaspari & Cohn, QJ 1999, p.748, Eq.(4.10)
  
  Gaspari=0
  x=r / c
  if(r < c)             Gaspari = -1/4*x^5 + 1/2*x^4 + 5/8*x^3 - 5/3*x^2 + 1
  if(r >= c && r < 2*c) Gaspari = 1/12*x^5 - 1/2*x^4 + 5/8*x^3 + 5/3*x^2 - 5*x + 4 - 2/(3*x)

  return(Gaspari)
}