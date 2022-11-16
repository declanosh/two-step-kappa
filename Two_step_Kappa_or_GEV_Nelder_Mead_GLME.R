source("./One_step_Kappa_Nelder_Mead.R")

kappa.pars.5.to.4.Coles <- function(gpa.pars,np.pars){
  # Convert Binomial frequency and GPA intensity to Kappa
  #
  # Input
  # gpa.pars <- GPA parameters c(location, scale, shape)
  # np.pars <- Binomial parameters, not that we use expectation np rather than p c(np, n)
  #
  # Outputs
  # out <- vector of Kappa parameters (location = xi, scale = alpha, shape = k, 2nd shape = h)
  
  
  # Extract GPA parameters
  xi.gp <- gpa.pars[1]
  alpha.gp <- gpa.pars[2]
  k.gp <- gpa.pars[3]
  
  # Extract Binomial parameters
  np <- np.pars[1]
  n <- np.pars[2]
  
  # Set Kappa shape parameter
  k <- k.gp
  
  # Avoid edge case for GPA shape
  if(k.gp == 0){k.gp = 0.001}
  
  # Calculate Kappa parameters
  xi <- xi.gp - (alpha.gp/k.gp) + (alpha.gp/k.gp)/(np^(-k.gp))
  
  alpha <- (alpha.gp)/(np^(-k.gp))
  
  h <- 1/n
  
  # Output kappa parameters
  out <- c("xi" = xi,"alpha" = alpha,"k"= k,"h"= h)
  names(out) <- c("xi", "alpha", "k", "h")
  return(out)
}

gp.shape.prior <- function(k,p=9,q=6){
  # Calculate the prior function for Kappa shape k from Martins and Stedinger (2000)
  
  pi <- ((0.5 + k)^(p-1))*((0.5-k)^(q-1))/beta(p,q)
  
  return(pi)
}

kappa.h.prior <- function(h,p=1,q=3){

   # Calculate the prior function for Kappa 2nd shape h, adapted from Martins and Stedinger (2000)
  
  pi <- ((h)^(p-1))*((1-h)^(q-1))/beta(p,q)
 
  
  return(pi)
}

nll.gpd.nm <- function(pars,data,lower.gp=c(-Inf,1e-6,-0.5),upper.gp=c(Inf,Inf,0.5),na.var=NA,prior.k=F){
  
  # Negative log likelihood function for the Genearlised Pareto distribution (GPD). 
  #
  # GPD quantile takes the form (note this is the same sign on the shape parameter using in Coles (2001) for the GEV, this is the negative sign of what is presented by Hosking and Wallis (1997)):
  # x(F) = xi + (alpha/k){1-((-log(F))^-k)}
  #
  # Inputs
  # pars <- 3 x 1 vector of GPD parameters c(xi, alpha, k)
  #
  # data <- n x 1 vector of POT maxima
  #
  # lower.gp & upper.gp <- respectively 3 x 1 vectors of the lower and upper bounds of parameter values of the GPD allowed 
  #
  # na.var <- return value of the function in the case that the nll is not defined or the parameter bounds are exceeded
  #
  # prior.k <- logical, if true applies the GLME prior function to the likelihood estimate
  #
  # Outputs
  # nll <- the value of the negative log likleihood function for the specified parameter. Or na.var in the event that nll is not defined or the parameter bounds are exceeded

  # Load extRemes package
  require(extRemes)
  
  # Assign data
  x <- data
  
  # check that parameter input vector is the right length
  len.par <- length(pars)
  if(len.par != 3){stop('Need exactly 3 parameters for Pareto')}
  
  
  # Check bounds on parameter vector
  between <- (pars >= lower.gp & pars <= upper.gp)
  if(sum(between)!=length(pars)){ return(na.var) }
  

  
  # Assign parameters
  xi <- pars[1]
  alpha <- pars[2]
  k <- pars[3]

  
  # Abort function if scale is negative
  if(alpha <= 0){return(na.var)}
  
  # Abort function if any maxima are too close to the bounds 
  tol.ub1 <- (xi-(alpha/k)) - max(x)
  
  if(k < 0 & tol.ub1 < 0.001 ){ return(na.var) }
  
  tol.lb1 <- min(x) - xi
  if(tol.lb1 < 0.001){  return(na.var) }
  
  # Avoid edge cases by setting  k to zero if very small in magnitude
  if(abs(k) < 0.001){ k = 0 }
  
  # Calculate likelihood
  
  if(k==0){
    # Exponential case (k=0)
    x.adj <- (x-xi)/alpha
    
    gpd.F <- 1 - exp(-x.adj)
    
  
    nll <- sum(-devd(x,xi,alpha,k,log = T,type="GP"))
    
	
  } else {
    # GPA case
	
    x.adj <- (1/k)*log(1+(k/alpha)*(x-xi))
    
    gpd.F <- 1 - exp(-x.adj)
    

    nll <- sum(-devd(x,xi,alpha,k,log = T,type="GP"))
  }
  
  # Add prior function if using GLME appraoch 
  if(prior.k ==T){
    
    pi <- gp.shape.prior(k)
    
    if(pi>0){
      
      nll <- nll-log(pi)
      
    } else {
      
      return(na.var)
      
    }
    
  }
  
  # Return value of nll
  return(nll)
}

nll.kappa.np.h.nm <- function(np.pars, gpa.pars, data, lower.np = c(0,0),upper.np = 
                           c(Inf,1),na.var=NA, fix.h.zero = F, prior.h = F){
  
  # Negative log likelihood function for Kappa distribution. 
  # This function calculates the NLL over AM data for the 4p Kappa distribution. Input are estimates of the Binomial parameters and the GPA parameters (pre calibrated from a concurrent POT series). 
  # These two sets of parameters are used to calculate the parameters of the four parameter Kappa and then subsequently calculate the nll.
  # Kappa CDF takes the form (note this is the same sign on the shape parameter using in Coles (2001) for the GEV, this is the negative sign of what is presented by Hosking and Wallis (1997):
  # F(x) = [1 - h{1+k(x-xi)/alpha}^(-1/k)]^(1/h)
  #
  # Inputs
  # np.pars <- 2 x 1 vector of the binomial quantities to be estimated in the form c(np, 1/n = h). Rather than the individual parameters of the Binomial distribution n and p, for ease of application in 
  # estimating the Kappa parameters we focus on np and 1/n = h instead. In the event that fix.h.zero == T, this should be a 1x1 scalar representing the expectation of the binomial parameter from infinite trials (alternatively the Poisson rate parameter lambda)
  #
  # gpa.pars <- 3 x 1 vector of Generalised Pareto parameters c(location, scale, shape)
  #
  # data <- Data frame of n years, must have at least two columns. Must have one column labelled "maxima" with the annual maxima used for fitting. Another column must be called "n.count" and indicates the number of threshold exceedances in that year. 
  # other columns allowed but not used
  #
  # lower.np & upper.np <- respectively 3 x 1 vectors of the lower and upper bounds of parameter values of the Binomial quantities c(np, 1/n) allowed 
  #
  # na.var <- return value of the function in the case that the nll is not defined or the parameter bounds are exceeded
  #
  # fix.h.zero <- logical, if TRUE set h = 0 (i.e. GEV) and overwrite 2nd input in np.pars
  #
  # prior.h <- logical, if TRUE applies the GLME prior function to the likelihood estimate
  #
  # Outputs
  # nll <- the value of the negative log likleihood function for the specified parameter. Or na.var in the event that nll is not defined or the parameter bounds are exceeded
  
  # Load extRemes package
  require(extRemes)
  
  # Check length of input vectors
  len.par <- length(np.pars)
  if((len.par != 2)&(fix.h.zero==F)){ stop('Need exactly 2 parameters: np and h') }
  if((len.par != 1)&(fix.h.zero==T)){ stop('Need exactly 1 parameters: np / lambda') } 

  
  # Set 1/n = h = 0 if fix.h.zero=TRUE
  if(fix.h.zero==T){
  
	  np.pars[2] <- 0
    lower.np[2] <- 0
	  
  }
  
  # Check parameter bounds 
  between <- (np.pars >= lower.np & np.pars <= upper.np)
  if(sum(between)!=length(np.pars)){ return(na.var) }
  
  # Assign GPD parameters
  xi.gp <- gpa.pars[1]
  alpha.gp <- gpa.pars[2]
  k.gp <- gpa.pars[3]
  
  # Assign binomial parameters
  np <- np.pars[1]
  h <- np.pars[2]
  
  
  #Assign Kappa parameters  Coles terminology
  xi <- xi.gp - (alpha.gp/k.gp) + (alpha.gp/k.gp)/(np^(-k.gp))
  alpha <- (alpha.gp)/(np^(-k.gp))
  k <- k.gp
  kappa.pars <- c(xi, alpha, k, h)

  # Calculate nll
  nll <- nll.kappa.4p.nm(kappa.pars,data=data,na.var=na.var)
  
  # Apply GLME prior function if prior.h = TRUE
  if(prior.h ==T){
    
    pi <- kappa.h.prior(h)
    
    if(pi>0){
      
      nll <- nll-log(pi)
      
    } else {
      
      return(na.var)
      
    }
    
  }
  
  # Return  values of nll
  return(nll)
  
  
}

two.step.kappa.gmle.k.fit.POT.nm <- function(data.in,lower.gp, upper.gp, lower.np, upper.np,fix.h.zero = F,prior.h = F, prior.k = F, num.it=2){
  # Wrapper code to implement the calibration of the two-step Kappa. The code executes two mle optimisations. The first estimates Generalised Pareto parameters from POT data, the second estimates the Binomial parameters from AMS
  #
  # Input
  # data.in <- list with two attributes. 1) AM a dataframe, each row represents a separate year. One column must be called "maxima" and contain the annual maxima of that year. Once column must be called "n.count" and contain the 
  # number of threshold exceedances in each year. 
  #
  # lower.gp & upper.gp <- respectively 3 x 1 vectors of the lower and upper bounds of parameter values of the GPD allowed 
  #
  # lower.np & upper.np <- respectively 2 x 1 vectors of the lower and upper bounds of parameter values of the binomial quantities to be estimated in the form c(np, 1/n = h). Rather than the individual parameters of the Binomial distribution n and p, for ease of application in 
  # estimating the Kappa parameters we focus on np and 1/n = h instead. In the event that fix.h.zero == T, this should be a 1x1 scalar representing the expectation of the binomial parameter from infinite trials (alternatively the Poisson rate parameter lambda)
  #
  # fix.h.zero <- logical, if TRUE set h = 0 (i.e. GEV) and overwrite 2nd input in np.pars
  #
  # prior.h <- logical, if TRUE applies the GLME prior function to the likelihood estimate
  #
  # num.it <- the number of Nelder - Mead iterations to run
  #
  # Output
  # List object with five levels:
  # 1) Kappa.params.coles <- 4 x 1 vector of Kappa parameters c(xi, alpha, k, h)
  # 2) GPA.out <- output from running of the optim function for the POT/GPA likelihood function (see ?optim for more details)
  # 3) Binom.out <- output from running of the optim function for the AM/Kappa likelihood function (see ?optim for more details)
  # 4) gpa.init.pars <- list of the initial gpa parameters used in the the Nelder Mead optimisation (these are the parameters that given the best optimisation from num.it runs) 
  # 5) binom.init.pars <- list of the initial Binomial parameters used in the the Nelder Mead optimisation (these are the parameters that given the best optimisation from num.it runs) 

  
    
  # Require the L-moments package
  require(lmom)
  
  # Call error if prior.h and fix.h.zero are both TRUE (prior on h with fixed h serves no purpose)
  if((prior.h==T)&(fix.h.zero==T)){ stop("h is fixed at zero, can't be used with h prior")}
  
  # Set optimiser based on how many parameters are optimised in the 2nd likelihood step
  if(fix.h.zero==T){
  
	optimiser = "Brent"
  
  } else {
  
	optimiser="Nelder-Mead"
  
  }
  
  # Assign variables
  data.in.AM <- data.in$AM
  data.in.AM.npos <- data.in.AM[data.in.AM$n.count > 0,] # Subset AM to only include years where there threshold exceedances
  data.in.POT <- data.in$POT
  
  ## First Likelihood step - Generalised Pareto / POT 
  
  # Estimate initial parameters for the Generalised Pareto using l-moments
  init.pars <- pelgpa(samlmu(data.in.POT))
  init.pars[3] <- -init.pars[3]
  sign.p3 <- sign(init.pars[3])
  
  # Test l-moment estimate, if not defined then set threshold to 99% of lowest POT maxima
  test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
  
  if(!is.finite(test.lmom)){
    
    init.pars[1] <- 0.99*min(data.in.POT)
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  # If this estimate is undefined then try the centre of the bounds for any previous initial parameter value estimates which exceeds the parameter bounds
  if(!is.finite(test.lmom)){
    
    between <- (init.pars >= lower.gp & init.pars <= upper.gp)
    init.pars[!between] <- 0.5*(lower.kappa[!between] + upper.kappa[!between])
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  # If this estimate is undefined then try in turn setting the shape parameter to +/- 0.1
  if(!is.finite(test.lmom)){

    init.pars[3] <- sign.p3*0.1
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  if(!is.finite(test.lmom)){
    
    init.pars[3] <- -sign.p3*0.1
    
  }
  
  # Define first initial parameter vector
  init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
  
  # For each iteration of the Nelder-Mead optimisation
  for(ii in 1:num.it){
   
    # For the first iteration try the "best estimate" initial parameters calculated above
    if(ii == 1){
      
      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = Inf, control=list(trace=0,maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      res.gpa <- res.gpa.temp
      init.pars.gpa.keep <- init.pars
      
	# For subsequent estimates try random values of the shape parameter k (using values of location and scale from previous estimate). Once a set of initial conditions which is gives a defined nll is obtained then rerun the Nelder-Mead optimisation
    } else {
      
      test.guess <- NA
      init.pars.v1 <- init.pars
      while (!is.finite(test.guess)) {
        
        init.pars <- c(init.pars.v1[1:2], runif(1,min=lower.gp[3],upper.gp[3]))
        test.guess <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k = prior.k, na.var = NA)
        
      }
      
      init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      
      
	  # If the new optimisation gives a better nll value than the previous best case then set this one as the new best case
      if(res.gpa.temp$value < res.gpa$value){
        
        res.gpa <- res.gpa.temp
        init.pars.gpa.keep <- init.pars
      }
      
    } 
    
      
  
    
  }
  
  ## 2nd likelihood step - Kappa / AM
  
  # Assign Generalised Pareto parameters obtained from 1st likelihood step
  gpa.pars <- res.gpa$par
  
  
  # If h is fixed at zero, then estimate np as the mean of the number of threshold exceedances each year in the sample
  if(fix.h.zero==T){
	
	init.pars <- c(mean(data.in.AM$n.count))
	init.pars.np <- list("np" = mean(data.in.AM$n.count))
  
   # Else estimate np as the mean of the number of threshold exceedances each year in the sample and n as the maximum number of exceedances in one year of the sample
  } else {
  
    init.pars <- c(mean(data.in.AM$n.count), 1/max(data.in.AM$n.count))
	init.pars.np <- list("np" = mean(data.in.AM$n.count), "h" = 1/max(data.in.AM$n.count))
  
  }

  
  # For each iteration of the Nelder-Mead optimisation
  for(ii in 1:num.it){
  
	 # For the first iteration try the "best estimate" initial parameters calculated above
    if(ii == 1){
      
		  if(fix.h.zero==T){
			  
		      res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf,method=optimiser, lower = lower.np[1], upper = upper.np[1], control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		  
			 } else {
			  
		    res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf,method=optimiser, control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		  
		  }
        
		
		
		res.binom <- res.binom.temp
    init.pars.np.keep <- init.pars
      
	# For subsequent estimates... 
    } else {
      
      test.guess <- NA
      init.pars.v1 <- init.pars.np
      
      while (!is.finite(test.guess)) {
          
			# If h is fixed randomly guess np between its lower bound and the minimum of its upper bound and twice the mean from the data sample
		    if(fix.h.zero==T){
	
			    init.pars.temp <- c(runif(1,min=lower.np[1],min(c(upper.np[1],2*mean(data.in.AM$n.count)))))
  
			# If h is not fixed randomly np and h between their lower and upper bounds
		    } else {
  
			    init.pars.temp <- c(runif(1,min=lower.np[1],upper.np[1]), runif(1,min=lower.np[2],upper.np[2]))
  
		    }
          
		    test.guess <- nll.kappa.np.h.nm(init.pars.temp,gpa.pars=gpa.pars, data.in.AM.npos,lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = NA)
          
      }
        
	# Assign initial parameter estimates as lists
    if(fix.h.zero==T){

			init.pars.np <- list("np" = init.pars.temp[1])
		
		}else{
		
			init.pars.np <- list("np" = init.pars.temp[1], "h"= init.pars.temp[2])
		
		}
        
		
		# Run Nelder-Mead optimisation 
		if(fix.h.zero==T){
			res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf, method = optimiser, lower=lower.np[1],upper=upper.np[1], control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		} else {
			res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf, method = optimiser, control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		}
        
          
          
      }
	  
      # If the new optimisation gives a better nll value than the previous best case then set this one as the new best case
      if(res.binom.temp$value < res.binom$value){
      
        res.binom <- res.binom.temp
        init.pars.np.keep <- init.pars
        
      }
    
	
      
    }
    
	# Ensure h= 0 in output for fixed h case
	if(fix.h.zero == T){ res.binom$par[2] <- 0 }
  temp <- res.binom$par
  temp[2] <- 1/temp[2]
  
  # Convert optimised GPD and Binomial parameters to 4 parameter Kappa
  kappa.pars <- kappa.pars.5.to.4.Coles(res.gpa$par,temp)
 
  #Return outputs
  return(list("Kappa.params.coles" = kappa.pars,"GPA.out" = res.gpa, "Binom.out" = res.binom, "gpa.init.pars"= init.pars.gpa.keep, "binom.init.pars" = init.pars.np.keep))
  
}

