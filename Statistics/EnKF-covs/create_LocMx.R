create_LocMx <- function(n, L_spacings, type){
  
  # Create the Localiz. mx on S1.
  # L_spacings is the len scale measured  in grid spacings.
  # type = "Gau" or "GC5" (Gaspari-Cohn)
  # M Tsy Jul 2017
  
  source('Gaspari.R')
  
  LocMx  <- matrix(ncol = n, nrow = n)
  
  sp2rad=2 *pi /n
  Lrad=L_spacings * sp2rad # 2pi radians ~ n grid spacings   =>  1 spacing = 2pi/n rad.
 
  if(type != "Gau" & type != "GC5"){
    message("Wrong type=")
    print(type)  
    stop("create_LocMx. Set type = Gau or GC5")
  }
  
  for(i in (1:n)) {
    for(j in (1:n)) {
      d=abs(i-j)    # spacings
      if ((d> floor(n/2))) {d <- n-d}
      rho <-  d * sp2rad
      r <- 2*sin(rho/2)  # chordal distance
      if(type == "Gau") {
        LocMx[i,j] <- exp(-0.5*(r/Lrad)^2)
      } else if(type == "GC5") {
        LocMx[i,j] <- Gaspari(r, Lrad)     
      }
    }
  }
  return(LocMx)
} 
