nll.kappa.4p.nm <- function(kappa.pars,data,lower.kappa=c(-Inf,1e-6,-Inf,-Inf),upper.kappa=c(Inf,Inf,Inf,Inf),na.var = NA){
 
  # Negative log likelihood function for the four parameter Kappa distribution. 
  #
  # Kappa CDF takes the form (note this is the same sign on the shape parameter using in Coles (2001) for the GEV, this is the negative sign of what is presented by Hosking and Wallis (1997)):
  # F(x) = [1 - h{1+k(x-xi)/alpha}^(-1/k)]^(1/h)
  #
  # Inputs
  # kappa.pars <- 4 x 1 vector of Kappa parameters c(xi, alpha, k, h)
  #
  # data <- Data frame of n years, must have at least two columns. Must have one column labelled "maxima" with the annual maxima used for fitting. Another column must be called "n.count" and indicates the number of threshold exceedances in that year. 
  # other columns allowed but not used
  #
  # lower.kappa & upper.kappa <- respectively 4 x 1 vectors of the lower and upper bounds of parameter values of the Kappa distribution allowed 
  #
  # na.var <- return value of the function in the case that the nll is not defined or the parameter bounds are exceeded
  #
  # Outputs
  #
  # nll <- the value of the negative log likleihood function for the specified parameter. Or na.var in the event that nll is not defined or the parameter bounds are exceeded

  
  # Check input length
  len.kap.par <- length(kappa.pars)
  if(len.kap.par != 4){stop('Need exactly 4 parameters in kappa.pars')}
  
  # Extrat data
  x <- data$maxima
  n.count <- data$n.count
  
  # Remove annual maxima if there are no threshold exceedances that year
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


  # Check bounds as per Hosking 1994 (adjusted for Coles terminology)
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
    # GUMBEL
    
    x.adj <- (x - xi)/alpha
    nll <- sum(log(alpha) + x.adj+ exp(-x.adj))
    
  }else if(h==0){
    # GEV
    
    x.adj <- (1/k)*log(1+(k/alpha)*(x-xi))
    
    nll <- sum(log(alpha) + (1+k)*x.adj + exp(-x.adj))
    
  }else if(k==0){
    # Kappa k=0  
    
    
    y <- (x - xi)/alpha
    nll <- sum(log(alpha)+y+(1-(1/h))*log(1-h*exp(-y)))
    
    
  } else {
    # Kappa
    
    y <- (x - xi)/alpha
    y <- 1 + k * y
    Fx <- (1-h*((y)^(-1/k)))^(1/h)
    nll <- sum(log(alpha)+(1+(1/k))*log(y)+(1-(1/h))*log((1-h*((y)^(-1/k)))))
    
    
  }
  
  return(nll)
}

one.step.kappa.mle.fit.ext.nm <- function(data.in,lower.kappa,upper.kappa,num.it=10){
  # Wrapper code to implement the calibration of the one-step Kappa
  #
  # Input
  # data.in <- must be a dataframe, each row represents a separate year. One column must be called "maxima" and contain the annual maxima of that year. Once column must be called "n.count"
  # lower.kappa & upper.kappa <- respectively 4 x 1 vectors of the lower and upper bounds of parameter values of the Kappa distribution allowed 
  # num.it <- the number of Nelder - Mead iterations to run
  #
  # Output
  # List object with three levels:
  # 1) Kappa.params.coles <- 4 x 1 vector of Kappa parameters c(xi, alpha, k, h)
  # 2) Kappa.out <- output from running of the optim function (see ?optim for more details)
  # 3) Initial parameters <- list of the initial parameters used to give the Nelder Mead optimisation (these are the parameters that given the best optimisation from num.it runs) 
   
  # Load L-moment package
  require(lmom)
  
  data.in.npos <- data.in
  
  # Calculate initial parameters
  # Test if L moments are defined for the Kappa
  # If not defined then calculate GEV parameters (negative k parameter values is to switch from Hosking and Wallis notation to Coles)
  err <- try(pelkap(samlmu(data.in.npos$maxima)),TRUE)
  
  if(typeof(err)=="double"){
    
    init.pars <- pelkap(samlmu(data.in.npos$maxima))
    init.pars[3] <- -init.pars[3]
    
  } else {
    
    init.pars <- c(pelgev(samlmu(data.in.npos$maxima)),0.001)
    init.pars[3] <- -init.pars[3]
    
  }
  
  # take sign of shape parameter k
  sign.p3 <- sign(init.pars[3])
  
  # Test if nll is defined for l moment estimate of initial parameters
  test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
  
  # If nll is defined for L moment estimate of initial parameters, set any parameter which fall outside the parameter bounds to the middle of the bound, again test if nll is defined
  if(!is.finite(test.lmom)){
 
    between <- (init.pars >= lower.kappa & init.pars <= upper.kappa)
    init.pars[!between] <- 0.5*(lower.kappa[!between] + upper.kappa[!between])
    test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
    
  }
  
  
  # If nll is still not defined then set the location parameter to be less than the smallest maxima
  if(!is.finite(test.lmom)){
    
    init.pars[1] <- 0.99*min(data.in.npos$maxima)
    test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
    
  }
  
 
  # If nll is still not defined, try setting the shape parameter k to +/- 0.1
  if(!is.finite(test.lmom)){
    
    init.pars[3] <- sign.p3*0.1
    test.lmom <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
    
  }

  if(!is.finite(test.lmom)){
    
    init.pars[3] <- -sign.p3*0.1
    
  }

  
  # Set up initial parameters list
  init.pars.4kap <- list("location" = init.pars[1], "scale" = init.pars[2], "shape" = init.pars[3], "shape2" = init.pars[4])
  
  # Four the number of set iterations, optimise the kappa parameters. In the first instance use the initial parameters estimated above, in subsequent instances use random values from the parameter space. Keep track of the lowest nll value, this becomes our optimised parameter values
  for(ii in 1:num.it){
        
    # use above initial parameters in first pass
    if(ii == 1){
      
      res.kap.temp <- optim(par = init.pars.4kap, fn = nll.kappa.4p.nm,data = data.in.npos, lower.kappa=lower.kappa, upper.kappa=upper.kappa, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      
      res.kap <- res.kap.temp
      init.pars.keep <- init.pars.4kap
      
    } else {
      
      # For subsequent runs use random parameter starting values
      test.guess <- NA
      init.pars.v1 <- init.pars
      while (!is.finite(test.guess)) {
        
        init.pars <- c(init.pars.v1[1:2], runif(1,min=lower.kappa[3],upper.kappa[3]), runif(1,min=lower.kappa[4],upper.kappa[4]))
        test.guess <- nll.kappa.4p.nm(init.pars,data.in.npos,lower.kappa = lower.kappa,upper.kappa=upper.kappa,na.var = NA)
         
      }
      
      init.pars.4kap <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3], "shape2" = init.pars[4])
      
      
      res.kap.temp <- optim(par = init.pars.4kap, fn = nll.kappa.4p.nm,data = data.in.npos, lower.kappa=lower.kappa, upper.kappa=upper.kappa, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      
	  # If the new nll is less than the previous best, then set this runs results as the new best
      if(res.kap.temp$value < res.kap$value){
        
        res.kap <- res.kap.temp
        init.pars.keep <- init.pars.4kap
        
      }
      
    }
    
    
    
    
  }
  
  return(list("Kappa.params.coles" = res.kap$par, "Kappa.out" = res.kap, "Kappa.init.pars"= init.pars.keep))
  
}
