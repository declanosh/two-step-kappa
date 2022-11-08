source("./One_step_Kappa_Nelder_Mead.R")

kappa.pars.5.to.4.Coles <- function(gpa.pars,np.pars){
  
  xi.gp <- gpa.pars[1]
  alpha.gp <- gpa.pars[2]
  k.gp <- gpa.pars[3]
  
  np <- np.pars[1]
  n <- np.pars[2]
  
  
  k <- k.gp
  
  if(k.gp == 0){k.gp = 0.001}
  xi <- xi.gp - (alpha.gp/k.gp) + (alpha.gp/k.gp)/(np^(-k.gp))
  
  alpha <- (alpha.gp)/(np^(-k.gp))
  
  h <- 1/n
  
  out <- c("xi" = xi,"alpha" = alpha,"k"= k,"h"= h)
  names(out) <- c("xi", "alpha", "k", "h")
  return(out)
}

gp.shape.prior <- function(k,p=9,q=6){
  
  pi <- ((0.5 + k)^(p-1))*((0.5-k)^(q-1))/beta(p,q)
  
  return(pi)
}

kappa.h.prior <- function(h,p=1,q=3){
  
  #p=1.5,q=6
  
  pi <- ((h)^(p-1))*((1-h)^(q-1))/beta(p,q)
 
  
  return(pi)
}

nll.gpd.nm <- function(pars,data,lower.gp=c(-Inf,1e-6,-0.5),upper.gp=c(Inf,Inf,0.5),na.var=NA,prior.k=F){
  
  
  len.par <- length(pars)
  if(len.par != 3){stop('Need exactly 3 parameters for Pareto')}
  
  
  # Check bounds
  between <- (pars >= lower.gp & pars <= upper.gp)
  if(sum(between)!=length(pars)){ return(na.var) }
  
  require(extRemes)
  x <- data
  
  xi <- pars[1]
  alpha <- pars[2]
  k <- pars[3]

  # need to change this k <- pars[3]
  
  if(alpha <= 0){return(na.var)}
  
  tol.ub1 <- (xi-(alpha/k)) - max(x)
  
  if(k < 0 & tol.ub1 < 0.001 ){ return(na.var) }
  
  tol.lb1 <- min(x) - xi
  if(tol.lb1 < 0.001){  return(na.var) }
  
  # Avoid edge cases by setting  k to zero if very small in magnitude
  if(abs(k) < 0.001){ k = 0 }
  
  if(k==0){
    
    x.adj <- (x-xi)/alpha
    
    gpd.F <- 1 - exp(-x.adj)
    
    #nll <- sum(log(alpha) + x.adj  - log(gpd.F))
    #nll <- sum(log(alpha) + x.adj)
    nll <- sum(-devd(x,xi,alpha,k,log = T,type="GP"))
    
  } else {
    
    x.adj <- (1/k)*log(1+(k/alpha)*(x-xi))
    
    gpd.F <- 1 - exp(-x.adj)
    

    nll <- sum(-devd(x,xi,alpha,k,log = T,type="GP"))
  }
  
  if(prior.k ==T){
    
    pi <- gp.shape.prior(k)
    
    if(pi>0){
      
      nll <- nll-log(pi)
      
    } else {
      
      return(na.var)
      
    }
    
  }
  
  
  return(nll)
}

nll.kappa.np.h.nm <- function(np.pars, gpa.pars, data, lower.np = c(0,0),upper.np = 
                           c(Inf,Inf),na.var=NA, fix.h.zero = F, prior.h = F){
  

  
  # Check params are in bounds
  # Check bounds
  
  
  len.par <- length(np.pars)
  if((len.par != 2)&(fix.h.zero==F)){ stop('Need exactly 2 parameters: np and h') }
  if((len.par != 1)&(fix.h.zero==T)){ stop('Need exactly 1 parameters: np / lambda') } 

  if(fix.h.zero==T){
  
	  np.pars[2] <- 0
    lower.np[2] <- 0
	  
  }
  
  between <- (np.pars >= lower.np & np.pars <= upper.np)
  if(sum(between)!=length(np.pars)){ return(na.var) }
  
  require(extRemes)
  
  xi.gp <- gpa.pars[1]
  alpha.gp <- gpa.pars[2]
  k.gp <- gpa.pars[3]
  
  np <- np.pars[1]
  
  h <- np.pars[2]
  
  
  #Coles terminology
  xi <- xi.gp - (alpha.gp/k.gp) + (alpha.gp/k.gp)/(np^(-k.gp))
  
  alpha <- (alpha.gp)/(np^(-k.gp))
  
  k <- k.gp
  
  #if((n < np) | (n <= 0) | (np <= 0)){ return(na.var)}
  #if((h < 0) | (h > 1) | (np <= 0)){ return(na.var) }
  
  kappa.pars <- c(xi, alpha, k, h)

  nll <- nll.kappa.4p.nm(kappa.pars,data=data,na.var=na.var)
  
  if(prior.h ==T){
    
    pi <- kappa.h.prior(h)
    
    if(pi>0){
      
      nll <- nll-log(pi)
      
    } else {
      
      return(na.var)
      
    }
    
  }
  
  return(nll)
  
  
}

two.step.kappa.gmle.k.fit.POT.nm <- function(data.in,lower.gp, upper.gp, lower.np, upper.np,fix.h.zero = F,prior.h = F, prior.k = F, num.it=2){
  # data.in must be a dataframe with $maxima and $n.count
  
  
  require(lmom)
  
  if((prior.h==T)&(fix.h.zero==T)){ stop("h is fixed at zero, can't be used with h prior")}
  
  if(fix.h.zero==T){
  
	optimiser = "Brent"
  
  } else {
  
	optimiser="Nelder-Mead"
  
  }
  
  data.in.AM <- data.in$AM
  data.in.AM.npos <- data.in.AM[data.in.AM$n.count > 0,]
  data.in.POT <- data.in$POT
  
  init.pars <- pelgpa(samlmu(data.in.POT))
  init.pars[3] <- -init.pars[3]
  sign.p3 <- sign(init.pars[3])
  
  test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
  
  if(!is.finite(test.lmom)){
    
    init.pars[1] <- 0.99*min(data.in.POT)
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  if(!is.finite(test.lmom)){
    
    between <- (init.pars >= lower.gp & init.pars <= upper.gp)
    init.pars[!between] <- 0.5*(lower.kappa[!between] + upper.kappa[!between])
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  if(!is.finite(test.lmom)){

    init.pars[3] <- sign.p3*0.1
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  if(!is.finite(test.lmom)){
    
    init.pars[3] <- -sign.p3*0.1
    
  }
  
  init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
  
  for(ii in 1:num.it){
   
    
    if(ii == 1){
      
      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = Inf, control=list(trace=0,maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      res.gpa <- res.gpa.temp
      init.pars.gpa.keep <- init.pars
      
    } else {
      
      test.guess <- NA
      init.pars.v1 <- init.pars
      while (!is.finite(test.guess)) {
        
        init.pars <- c(init.pars.v1[1:2], runif(1,min=lower.gp[3],upper.gp[3]))
        test.guess <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k = prior.k, na.var = NA)
        
      }
      
      init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      
      
      if(res.gpa.temp$value < res.gpa$value){
        
        res.gpa <- res.gpa.temp
        init.pars.gpa.keep <- init.pars
      }
      
    } 
    
      
  
    
  }
  
  gpa.pars <- res.gpa$par
  
  if(fix.h.zero==T){
	
	init.pars <- c(mean(data.in.AM$n.count))
	init.pars.np <- list("np" = mean(data.in.AM$n.count))
  
  } else {
  
    init.pars <- c(mean(data.in.AM$n.count), max(data.in.AM$n.count))
	init.pars.np <- list("np" = mean(data.in.AM$n.count), "h" = 1/max(data.in.AM$n.count))
  
  }

  
  
  for(ii in 1:num.it){
  
    if(ii == 1){
      
		  if(fix.h.zero==T){
			  
		      res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf,method=optimiser, lower = lower.np[1], upper = upper.np[1], control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		  
			 } else {
			  
		    res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf,method=optimiser, control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		  
		  }
        
		
		
		res.binom <- res.binom.temp
    init.pars.np.keep <- init.pars
      
    } else {
      
      test.guess <- NA
      init.pars.v1 <- init.pars.np
      
      while (!is.finite(test.guess)) {
          
		    if(fix.h.zero==T){
	
			    init.pars.temp <- c(runif(1,min=lower.np[1],min(c(upper.np[1],2*mean(data.in.AM$n.count)))))
  
		    } else {
  
			    init.pars.temp <- c(runif(1,min=lower.np[1],upper.np[1]), runif(1,min=lower.np[2],upper.np[2]))
  
		    }
          
		    test.guess <- nll.kappa.np.h.nm(init.pars.temp,gpa.pars=gpa.pars, data.in.AM.npos,lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = NA)
          
      }
        
		
    if(fix.h.zero==T){

			init.pars.np <- list("np" = init.pars.temp[1])
		
		}else{
		
			init.pars.np <- list("np" = init.pars.temp[1], "h"= init.pars.temp[2])
		
		}
        
		
		
		if(fix.h.zero==T){
			res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf, method = optimiser, lower=lower.np[1],upper=upper.np[1], control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		} else {
			res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf, method = optimiser, control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		}
        
          
          
      }
      
      if(res.binom.temp$value < res.binom$value){
      
        res.binom <- res.binom.temp
        init.pars.np.keep <- init.pars
        
      }
    
	
      
    }
    
	if(fix.h.zero == T){ res.binom$par[2] <- 0 }
  temp <- res.binom$par
  temp[2] <- 1/temp[2]
  
  kappa.pars <- kappa.pars.5.to.4.Coles(res.gpa$par,temp)
 
  
  return(list("Kappa.params.coles" = kappa.pars,"GPA.out" = res.gpa, "Binom.out" = res.binom, "gpa.init.pars"= init.pars.gpa.keep, "binom.init.pars" = init.pars.np.keep))
  
}

two.step.kappa.gmle.k.fit.POT.ObsData.nm <- function(data.in,lower.gp, upper.gp, lower.np, upper.np,fix.h.zero = F,prior.h = F, prior.k = F, num.it=2){
  # data.in must be a dataframe with $maxima and $n.count
  
  
  require(lmom)
  
  if((prior.h==T)&(fix.h.zero==T)){ stop("h is fixed at zero, can't be used with h prior")}
  
  if(fix.h.zero==T){
  
	optimiser = "Brent"
  
  } else {
  
	optimiser="Nelder-Mead"
  
  }
  
  data.in.AM <- data.in$AM
  data.in.AM.npos <- data.in.AM[data.in.AM$n.count > 0,]
  data.in.POT <- data.in$POT
  
  init.pars <- pelgpa(samlmu(data.in.POT))
  init.pars[3] <- -init.pars[3]
  sign.p3 <- sign(init.pars[3])
  
  test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
  
  if(!is.finite(test.lmom)){
    
    init.pars[1] <- 0.99*min(data.in.POT)
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  if(!is.finite(test.lmom)){
    
    between <- (init.pars >= lower.gp & init.pars <= upper.gp)
    init.pars[!between] <- 0.5*(lower.kappa[!between] + upper.kappa[!between])
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  if(!is.finite(test.lmom)){

    init.pars[3] <- sign.p3*0.1
    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = NA)
    
  }
  
  if(!is.finite(test.lmom)){
    
    init.pars[3] <- -sign.p3*0.1
    
  }
  
  init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
  
  for(ii in 1:num.it){
   
    
    if(ii == 1){
      
      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = Inf, control=list(trace=0,maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      res.gpa <- res.gpa.temp
      init.pars.gpa.keep <- init.pars
      
    } else {
      
      test.guess <- NA
      init.pars.v1 <- init.pars
      while (!is.finite(test.guess)) {
        
        init.pars <- c(init.pars.v1[1:2], runif(1,min=lower.gp[3],upper.gp[3]))
        test.guess <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, prior.k = prior.k, na.var = NA)
        
      }
      
      init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, prior.k =prior.k, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
      
      
      if(res.gpa.temp$value < res.gpa$value){
        
        res.gpa <- res.gpa.temp
        init.pars.gpa.keep <- init.pars
      }
      
    } 
    
      
  
    
  }
  
  gpa.pars <- res.gpa$par
  
  if(fix.h.zero==T){
	
	init.pars <- c(mean(data.in.AM$n.count))
	init.pars.np <- list("np" = mean(data.in.AM$n.count))
  
  } else {
  
    init.pars <- c(mean(data.in.AM$n.count), max(data.in.AM$n.count))
	init.pars.np <- list("np" = mean(data.in.AM$n.count), "h" = 1/max(data.in.AM$n.count))
  
  }

  
  
  for(ii in 1:num.it){
  
    if(ii == 1){
      
		  if(fix.h.zero==T){
			  
		      res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf,method=optimiser, lower = lower.np[1], upper = upper.np[1], control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		  
			 } else {
			  
		    res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf,method=optimiser, control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		  
		  }
        
		
		
		res.binom <- res.binom.temp
    init.pars.np.keep <- init.pars
      
    } else {
      
      test.guess <- NA
      init.pars.v1 <- init.pars.np
      
      while (!is.finite(test.guess)) {
          
		    if(fix.h.zero==T){
	
			    init.pars.temp <- c(runif(1,min=lower.np[1],min(c(upper.np[1],2*mean(data.in.AM$n.count)))))
  
		    } else {
  
			    init.pars.temp <- c(runif(1,min=lower.np[1],upper.np[1]), runif(1,min=lower.np[2],upper.np[2]))
  
		    }
          
		    test.guess <- nll.kappa.np.h.nm(init.pars.temp,gpa.pars=gpa.pars, data.in.AM.npos,lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = NA)
          
      }
        
		
    if(fix.h.zero==T){

			init.pars.np <- list("np" = init.pars.temp[1])
		
		}else{
		
			init.pars.np <- list("np" = init.pars.temp[1], "h"= init.pars.temp[2])
		
		}
        
		
		
		if(fix.h.zero==T){
			res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf, method = optimiser, lower=lower.np[1],upper=upper.np[1], control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		} else {
			res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.h.nm,gpa.par=gpa.pars ,data = data.in.AM.npos, lower.np=lower.np, upper.np=upper.np,fix.h.zero = fix.h.zero, prior.h = prior.h, na.var = Inf, method = optimiser, control=list(trace = 0, maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
		}
        
          
          
      }
      
      if(res.binom.temp$value < res.binom$value){
      
        res.binom <- res.binom.temp
        init.pars.np.keep <- init.pars
        
      }
    
	
      
    }
    
	if(fix.h.zero == T){ res.binom$par[2] <- 0 }
  temp <- res.binom$par
  temp[2] <- 1/temp[2]
  
  kappa.pars <- kappa.pars.5.to.4.Coles(res.gpa$par,temp)
 
  
  return(list("Kappa.params.coles" = kappa.pars,"GPA.out" = res.gpa, "Binom.out" = res.binom, "gpa.init.pars"= init.pars.gpa.keep, "binom.init.pars" = init.pars.np.keep))
  
}

#two.step.kappa.mle.fit.POT.ObsData.nm <- function(data.in,lower.gp, upper.gp, lower.np, upper.np, num.it=5){
#  # data.in must be a dataframe with $maxima and $n.count
#  
#  
#  require(lmom)
#  #require(hydromad)
#  
#  data.in.AM <- data.in$AM
#  #data.in.AM.npos <- data.in.AM[data.in.AM$n.count > 0,]
#  data.in.POT <- data.in$POT
#  
#  init.pars <- pelgpa(samlmu(data.in.POT))
#  init.pars[3] <- -init.pars[3]
#  sign.p3 <- sign(init.pars[3])
#  
#  test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, na.var = NA)
#  
#  if(!is.finite(test.lmom)){
#    
#    init.pars[1] <- 0.99*min(data.in.POT)
#    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, na.var = NA)
#    
#  }
#  
#  if(!is.finite(test.lmom)){
#    
#    between <- (init.pars >= lower.gp & init.pars <= upper.gp)
#    init.pars[!between] <- 0.5*(lower.kappa[!between] + upper.kappa[!between])
#    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, na.var = NA)
#    
#  }
#  
#  if(!is.finite(test.lmom)){
#    
#    init.pars[3] <- sign.p3*0.1
#    test.lmom <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, na.var = NA)
#    
#  }
#  
#  if(!is.finite(test.lmom)){
#    
#    init.pars[3] <- -sign.p3*0.1
#    
#  }
#  
#  init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
#  
#  for(ii in 1:num.it){
#    
#    
#    if(ii == 1){
#      
#      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
#      res.gpa <- res.gpa.temp
#      init.pars.gpa.keep <- init.pars
#      
#    } else {
#      
#      test.guess <- NA
#      init.pars.v1 <- init.pars
#      while (!is.finite(test.guess)) {
#        
#        init.pars <- c(init.pars.v1[1:2], runif(1,min=lower.gp[3],upper.gp[3]))
#        test.guess <- nll.gpd.nm(init.pars,data.in.POT,lower.gp=lower.gp, upper.gp=upper.gp, na.var = NA)
#        
#      }
#      
#      init.pars.gp <- list("location" = init.pars[1], "scale"= init.pars[2], "shape" = init.pars[3])
#      res.gpa.temp <-  optim(par = init.pars.gp, fn = nll.gpd.nm,data = data.in.POT, lower.gp=lower.gp, upper.gp=upper.gp, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
#      
#      
#      if(res.gpa.temp$value < res.gpa$value){
#        
#        res.gpa <- res.gpa.temp
#        init.pars.gpa.keep <- init.pars
#      }
#      
#    } 
#    
#    
#    
#    
#  }
#  
#  gpa.pars <- res.gpa$par
#  init.pars <- c(mean(data.in.AM$n.count), max(data.in.AM$n.count))
#  init.pars.np <- list("np" = mean(data.in.AM$n.count), "n" = max(data.in.AM$n.count))
#  
#  
#  for(ii in 1:num.it){
#    
#    if(ii == 1){
#      
#      res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.nm,gpa.par=gpa.pars ,data = data.in.AM, lower.np=lower.np, upper.np=upper.np, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
#      res.binom <- res.binom.temp
#      init.pars.np.keep <- init.pars
#      
#    } else {
#      
#      test.guess <- NA
#      init.pars.v1 <- init.pars.np
#      while (!is.finite(test.guess)) {
#        
#        init.pars.temp <- c(runif(1,min=lower.np[1],upper.np[1]), runif(1,min=lower.np[2],upper.np[2]))
#        test.guess <- nll.kappa.np.nm(init.pars.temp,gpa.pars=gpa.pars, data.in.AM,lower.np=lower.np, upper.np=upper.np, na.var = NA)
#        
#      }
#      
#      init.pars.np <- list("np" = init.pars.temp[1], "n"= init.pars.temp[2])
#      res.binom.temp <- optim(par = init.pars.np, fn = nll.kappa.np.nm,gpa.par=gpa.pars ,data = data.in.AM, lower.np=lower.np, upper.np=upper.np, na.var = Inf, control=list(maxit = 1000000, reltol = sqrt(.Machine$double.eps)))
#      
#      
#      
#    }
#    
#    if(res.binom.temp$value < res.binom$value){
#      
#      res.binom <- res.binom.temp
#      init.pars.np.keep <- init.pars
#      
#    }
#    
#    
#  }
#  
#  temp <- res.binom$par
#  temp[2] <- 1/temp[2]
#
#  
#  kappa.pars <- kappa.pars.5.to.4.Coles(res.gpa$par,temp)
#  
#  
#  return(list("Kappa.params.coles" = kappa.pars,"GPA.out" = res.gpa, "Binom.out" = res.binom, "gpa.init.pars"= init.pars.gpa.keep, "binom.init.pars" = init.pars.np.keep))
#  
#}
