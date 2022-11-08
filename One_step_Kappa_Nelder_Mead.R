nll.kappa.4p.nm <- function(kappa.pars,data,lower.kappa=c(-Inf,1e-6,-Inf,-Inf),upper.kappa=c(Inf,Inf,Inf,Inf),na.var = NA){
 
  ## Import extremes package 
  #require(extRemes)

  len.kap.par <- length(kappa.pars)
  if(len.kap.par != 4){stop('Need exactly 4 parameters in kappa.pars')}
  
  # Extrat data
  x <- data$maxima
  n.count <- data$n.count
  
  x <- x[n.count!=0]
  n.count <- n.count[n.count!=0]
  
  # Check bounds
  between <- (kappa.pars >= lower.kappa & kappa.pars <= upper.kappa)
  if(sum(between)!=4){ return(na.var) }
  
  # Assign params, Coles terminology
  xi <- kappa.pars[1]
  
  alpha <- kappa.pars[2]
  
  k <- kappa.pars[3]
  
  h <- kappa.pars[4]


  # Check bounds as per Hosking 1994
  if(k != 0){
	if(k < 0 & (xi-(alpha/k)) <= max(x)){ return(na.var) }
  
	if(h > 0 & ((xi-(alpha/k)*(1-(h^k))) >= min(x))){ return(na.var) }
  
	if(h <= 0 & k > 0 & ((xi - alpha/k) >= min(x))){ return(na.var) }
  } else if(k == 0){
	
	if((h > 0)&((xi +alpha*log(h)) >= min(x))){return(na.var)}
	
  }

  
  # Avoid edge cases by setting h and k to zero if very small in magnitude
  if(abs(k) < 0.001){ k = 0 }
  
  if(abs(h) < 0.001){ h = 0 }
  
  
  # Calculate likelihood
  if(k==0  && h ==0){
    # gumbel  
    
    x.adj <- (x - xi)/alpha
    nll <- sum(log(alpha) + x.adj+ exp(-x.adj))
    
  }else if(h==0){
    # gev  
    
    x.adj <- (1/k)*log(1+(k/alpha)*(x-xi))
    
    nll <- sum(log(alpha) + (1+k)*x.adj + exp(-x.adj))
    
  }else if(k==0){
    # kappa k=0  
    
    
    y <- (x - xi)/alpha
    nll <- sum(log(alpha)+y+(1-(1/h))*log(1-h*exp(-y)))
    
    
  } else {
    # kappa
    
    y <- (x - xi)/alpha
    y <- 1 + k * y
    Fx <- (1-h*((y)^(-1/k)))^(1/h)
    nll <- sum(log(alpha)+(1+(1/k))*log(y)+(1-(1/h))*log((1-h*((y)^(-1/k)))))
    
    
  }
  
  return(nll)
}

one.step.kappa.mle.fit.ext.nm <- function(data.in,lower.kappa,upper.kappa,num.it=10){
  # data.in must be a dataframe with $maxima and $n.count
  
  
  require(lmom)
  
  data.in.npos <- data.in
  
  # Test if Lmoments are defined for the Kappa
  # If not defined set to GEV
  err <- try(pelkap(samlmu(data.in.npos$maxima)),TRUE)
  
  if(typeof(err)=="double"){
    
    init.pars <- pelkap(samlmu(data.in.npos$maxima))
    init.pars[3] <- -init.pars[3]
    
  } else {
    
    init.pars <- c(pelgev(samlmu(data.in.npos$maxima)),0.001)
    init.pars[3] <- -init.pars[3]
    
  }
  
  sign.p3 <- sign(init.pars[3])
  test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
  
  if(!is.finite(test.lmom)){
 
    between <- (init.pars >= lower.kappa & init.pars <= upper.kappa)
    init.pars[!between] <- 0.5*(lower.kappa[!between] + upper.kappa[!between])
    test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
    
  }
  
  
  
  if(!is.finite(test.lmom)){
    
    init.pars[1] <- 0.99*min(data.in.npos$maxima)
    test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
    
  }
  
 
  
  if(!is.finite(test.lmom)){
    
    init.pars[3] <- sign.p3*0.1
    test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
    
  }

  
  
  if(!is.finite(test.lmom)){
    
    init.pars[3] <- -sign.p3*0.1
    
  }

  
  
  init.pars.4kap <- list("location" = init.pars[1], "scale" = init.pars[2], "shape" = init.pars[3], "shape2" = init.pars[4])
  
  for(ii in 1:num.it){
        
    
    if(ii == 1){
      
      res.kap.temp <- optim(par = init.pars.4kap, fn = nll.kappa.4p.nm,data = data.in.npos, lower.kappa=lower.kappa, upper.kappa=upper.kappa, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      
      res.kap <- res.kap.temp
      init.pars.keep <- init.pars.4kap
      
    } else {
      
      
      test.guess <- NA
      init.pars.v1 <- init.pars
      while (!is.finite(test.guess)) {
        
        init.pars <- c(init.pars.v1[1:2], runif(1,min=lower.kappa[3],upper.kappa[3]), runif(1,min=lower.kappa[4],upper.kappa[4]))
        #init.pars <- c(runif(1,min=lower.kappa[1],upper.kappa[1]), runif(1,min=lower.kappa[2],upper.kappa[2]), runif(1,min=lower.kappa[3],upper.kappa[3]), runif(1,min=lower.kappa[4],upper.kappa[4]))
        test.guess <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
         
      }
      
      init.pars.4kap <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3], "shape2" = init.pars[4])
      
      
      res.kap.temp <- optim(par = init.pars.4kap, fn = nll.kappa.4p.nm,data = data.in.npos, lower.kappa=lower.kappa, upper.kappa=upper.kappa, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      
      if(res.kap.temp$value < res.kap$value){
        
        res.kap <- res.kap.temp
        init.pars.keep <- init.pars.4kap
        
      }
      
    }
    
    
    
    
  }
  
  return(list("Kappa.params.coles" = res.kap$par, "Kappa.out" = res.kap, "Kappa.init.pars"= init.pars.keep))
  
}
