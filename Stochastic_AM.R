require(lmom)
Gen.AM.Kappa.ns.shape2 <- function(nyears=100, start.year = 0, kappa.para = c(loc = 0, scale = 1, shape = 0, shape2.int = 0, shape2.slope = 0 )){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  annual.max <- rep(0,nyears)
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  kappa.para[3] <- -kappa.para[3]
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    para.temp <- c(kappa.para[1:3],kappa.para[4]+kappa.para[5]*(start.year+ii-1))
    annual.max[ii] <- quakap(runif(1),para.temp)
    
  }
  #annual.max <- annual.max[annual.max!=0]
  
  
  
  return(annual.max)
  
}

Gen.AM.POT.nsBinom <- function(nyears=100, start.year = 0, GPA.para = c(loc = 0, scale = 1, shape = 0), BN.para = c(p = 1, n.int = 10, n.grad = 0)){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  annual.max <- rep(0,nyears)
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  GPA.para[3] <- -GPA.para[3]
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    n.temp <- round(BN.para[2] + BN.para[3]*(start.year+ii-1))
    n.events <-  qbinom(runif(1),n.temp, BN.para[1])
    annual.max[ii] <- max(quagpa(runif(n.events),GPA.para))
       
       
  }
  annual.max[annual.max == -Inf] = NA
  #annual.max <- annual.max[annual.max!=0]
  return(annual.max)
  
}

Gen.AM.POT.nsBinom.n.return <- function(nyears=100, start.year = 0, GPA.para = c(loc = 0, scale = 1, shape = 0), BN.para = c(p = 1, n.int = 10, n.grad = 0)){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  annual.max <- rep(0,nyears)
  n.events.store <- rep(0,nyears)
  POT <- vector(mode = "numeric",length = 0)
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  GPA.para[3] <- -GPA.para[3]
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    n.temp <- round(BN.para[2] + BN.para[3]*(start.year+ii-1))
    n.events <-  qbinom(runif(1),n.temp, BN.para[1])
    n.events.store[ii]<- n.events
    
    
    if(n.events != 0){ 
      
      max.temp <- quagpa(runif(n.events),GPA.para) 
      annual.max[ii] <- max(max.temp,GPA.para[1])
    
    } else {
      
      annual.max[ii] <- NA 
      
      }
    
    
    
    if(n.events != 0){ POT[(length(POT)+1):(length(POT)+n.events)] <- max.temp }
    
    
  }
  
  #annual.max <- annual.max[annual.max!=0]
  annual.max[annual.max == -Inf] = NA
  out <- list("AM" = data.frame(year = seq(0,nyears-1),n.count = n.events.store, maxima = annual.max),"POT"= POT )
  return(out)
  
}

Gen.AM.nsPOT.nsBinom.n.return <- function(nyears=100, start.year = 0, GPA.para = c(loc.int = 0, loc.grad = 0, scale.int = 1, scale.grad = 0, shape.int = 0, shape.grad = 0), BN.para = c(np.int = 1,np.grad = 0, n.int = 10, n.grad = 0),np.constant = T, np.scale.factor=0.5){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  
  if(BN.para[2]!=0 & BN.para[4] != 0){ stop("One of np and n must have zero gradient")}
  if( np.scale.factor < 0 || np.scale.factor > 1){ stop("np scale factor must be between 0 and 1 (it is effectively p)")}
  
  
  annual.max <- rep(0,nyears)
  n.events.store <- rep(0,nyears)
  POT <- vector(mode = "numeric",length = 0)
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  
  years <- seq(start.year,start.year+nyears-1,by=1)
  gp.loc <- GPA.para[1]+GPA.para[2]*years
  gp.scale <- GPA.para[3]+GPA.para[4]*years
  gp.shape <- GPA.para[5]+GPA.para[6]*years
  
  gp.shape <- -gp.shape
  
  gp.pars <- cbind(gp.loc,gp.scale,gp.shape)
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    n.temp <- round(BN.para[3] + BN.para[4]*(start.year+ii-1))
    
    if(np.constant==T){
      
      np.temp <- BN.para[1]
      
    } else if (BN.para[2]!=0){
      
      np.temp <- BN.para[1] + BN.para[2]*(start.year+ii-1)
      
    } else if (BN.para[4]!=0){
      
      np.temp <- n.temp * np.scale.factor
      
    }
    
    
    
    
    
    p.temp <- np.temp/n.temp
    n.events <-  qbinom(runif(1),n.temp, p.temp)
    n.events.store[ii]<- n.events
    
    
    if(n.events != 0){ 
      
      max.temp <- quagpa(runif(n.events),gp.pars[ii,]) 
      annual.max[ii] <- max(max.temp,gp.pars[ii,1])
      
    } else {
      
      annual.max[ii] <- NA 
      
    }
    
    
    
    if(n.events != 0){ POT[(length(POT)+1):(length(POT)+n.events)] <- max.temp }
    
    
  }
  
  #annual.max <- annual.max[annual.max!=0]
  annual.max[annual.max == -Inf] = NA
  out <- list("AM" = data.frame(year = seq(0,nyears-1),n.count = n.events.store, maxima = annual.max),"POT"= POT )
  return(out)
  
}


Gen.AM.POT.ns.Pois.n.return <- function(nyears=100, start.year = 0, GPA.para = c(loc = 0, scale = 1, shape = 0), Pois.para = c(np.int = 10, np.grad = 0)){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  annual.max <- rep(0,nyears)
  n.events.store <- rep(0,nyears)
  POT <- vector(mode = "numeric",length = 0)
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  GPA.para[3] <- -GPA.para[3]
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    np.temp <- round(Pois.para[1] + Pois.para[2]*(start.year+ii-1))
    n.events <-  qpois(runif(1),np.temp)
    n.events.store[ii]<- n.events
    
    
    if(n.events != 0){ 
      
      max.temp <- quagpa(runif(n.events),GPA.para) 
      annual.max[ii] <- max(max.temp,GPA.para[1])
      
    } else {
      
      annual.max[ii] <- NA 
      
    }
    
    
    
    if(n.events != 0){ POT[(length(POT)+1):(length(POT)+n.events)] <- max.temp }
    
    
  }
  
  #annual.max <- annual.max[annual.max!=0]
  annual.max[annual.max == -Inf] = NA
  out <- list("AM" = data.frame(year = seq(0,nyears-1),n.count = n.events.store, maxima = annual.max),"POT"= POT )
  return(out)
  
}



Gen.AM.POT.ns.np.allp.Binom.n.return <- function(nyears=100, start.year = 0, GPA.para = c(loc = 0, scale = 1, shape = 0), BN.para = c(np.int = 10, np.grad = 0, n.int = 10)){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  np.int <- BN.para[1]
  np.grad <- BN.para[2]
  n.int <- BN.para[3]
  
  annual.max <- rep(0,nyears)
  n.events.store <- rep(0,nyears)
  POT <- vector(mode = "numeric",length = 0)
  
  
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  GPA.para[3] <- -GPA.para[3]
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    np.temp <- np.int + np.grad*(start.year+ii-1)
    n.events <-  qbinom(runif(1), n.int, np.temp/n.int)
    n.events.store[ii]<- n.events
    
    
    if(n.events != 0){ 
      
      max.temp <- quagpa(runif(n.events),GPA.para) 
      annual.max[ii] <- max(max.temp,GPA.para[1])
      
    } else {
      
      annual.max[ii] <- NA 
      
    }
    
    
    
    if(n.events != 0){ POT[(length(POT)+1):(length(POT)+n.events)] <- max.temp }
    
    
  }
  
  #annual.max <- annual.max[annual.max!=0]
  out <- list("AM" = data.frame(year = seq(0,nyears-1),n.count = n.events.store, maxima = annual.max),"POT"= POT )
  return(out)
  
}



Gen.AM.POT.ns.np.alln.Binom.n.return <- function(nyears=100, start.year = 0, GPA.para = c(loc = 0, scale = 1, shape = 0), BN.para = c(np.int = 10, np.grad = 0, n.int = 10)){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  
  np.int <- BN.para[1]
  np.grad <- BN.para[2]
  n.int <- BN.para[3]
  
  annual.max <- rep(0,nyears)
  n.events.store <- rep(0,nyears)
  POT <- vector(mode = "numeric",length = 0)
  p.int <- np.int/n.int
  
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  GPA.para[3] <- -GPA.para[3]
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    np.temp <- np.int + np.grad*(start.year+ii-1)
    n.events <-  qbinom(runif(1), round(np.temp/p.int), p.int)
    n.events.store[ii]<- n.events
    
    
    if(n.events != 0){ 
      
      max.temp <- quagpa(runif(n.events),GPA.para) 
      annual.max[ii] <- max(max.temp,GPA.para[1])
      
    } else {
      
      annual.max[ii] <- NA 
      
    }
    
    

    
    if(n.events != 0){ POT[(length(POT)+1):(length(POT)+n.events)] <- max.temp }
    
    
  }
  
  #annual.max <- annual.max[annual.max!=0]
  out <- list("AM" = data.frame(year = seq(0,nyears-1),n.count = n.events.store, maxima = annual.max),"POT"= POT )
  return(out)
  
}

Gen.AM.Kappa.ns.loc <- function(nyears=100, start.year = 0, kappa.para = c(loc.int = 0, loc.slope = 0, scale = 1, shape = 0, shape2 = 0 )){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters c(location int, location slope, scale, shape, 2nd shape)
  #
  annual.max <- rep(0,nyears)
  
  # Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
  kappa.para[4] <- -kappa.para[4]
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    para.temp <- c(kappa.para[1]+kappa.para[2]*(start.year+ii-1),kappa.para[3:5])
    annual.max[ii] <- quakap(runif(1),para.temp)
    
  }
  #annual.max <- annual.max[annual.max!=0]
  return(annual.max)
  
}


Gen.AM <- function(nyears=100, af.para = c(5,10),dist.type="uniform",kappa.para=c(0,1,0,1)){
  # Generate AM series
  # nyears - number of years
  # af.para - parameters of the distribution to use  for arrival frequency (order of params as below):
  #       - fixed: c(value)
  #       - uniform: c(lower, upper)
  #       - binomial: c(n,p)
  #       - poisson: c(rate)
  #       - geometric: c(k,p)
  #
  #
  # dist.type - distribution type:
  #       - fixed
  #       - uniform (N.B. discrete uniform is used)
  #       - binomial
  #       - poisson
  #       - geometric
  #       - linear.normal (i.e. reproduce linear regression with a normal dist)
  #
  # kappa.para - Kappa parent parameters c(location, scale, shape, 2nd shape)
  #
  annual.max <- rep(0,nyears)
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    # Generate number of events by poisson or fixed to rate parameter depending on the switch apply.pois
    if(dist.type=="fixed"){
      
      number.events <- af.para[1]
      
    }else if(dist.type=="uniform"){
      
      number.events <- af.para[1] + floor((af.para[2] - af.para[1] + 1)*runif(1))
      
    }else if(dist.type=="binomial"){
      
      number.events <- qbinom(runif(1),af.para[1],af.para[2])
      
    }else if(dist.type=="poisson"){
      
      number.events <- qpois(runif(1),af.para[1])
      
    }else if(dist.type=="geometric"){
      
      number.events <- qgeom(runif(1),af.para[1],af.para[2])
      
    }else if (dist.type=="linear.normal"){
      
      number.events <-  max(round(qnorm(runif(1),af.para[1],af.para[2])),0)
      
    }
    
    # If randomly generate no events then set annual max to zero
    if(number.events==0){
      
      annual.max[ii] <- NA
      
    }else{
      
      random.f <- runif(number.events)
      peaks <- quakap(f = random.f,para = kappa.para)
      annual.max[ii] <- max(peaks)
      
    }
    
  }
  #annual.max <- annual.max[annual.max!=0]
  return(annual.max)
  
}


Repeat.AM.Simulations <- function(num.sim=100, nyears=100, af.para = c(5,10),dist.type="uniform",kappa.para=c(0,1,0,1)){
  # Genearte multiple AM series and calculate lmoments
  
  # num.sim - number of simulations
  # nyears - number of years
  # af.para - parameters of the distribution to use  for arrival frequency (order of params as below):
  #       - fixed: c(value)
  #       - uniform: c(lower, upper)
  #       - binomial: c(n,p)
  #       - poisson: c(rate)
  #       - geometric: c(k,p)
  #
  #
  # dist.type - distribution type:
  #       - fixed
  #       - uniform (N.B. discrete uniform is used)
  #       - binomial
  #       - poisson
  #       - geometric
  #
  # kappa.para - Kappa parent parameters c(location, scale, shape, 2nd shape)
  
  annual.max.lmom <- data.frame(L1 = rep(0,num.sim), L2 = rep(0,num.sim), Lskew = rep(0,num.sim), Lkurt = rep(0,num.sim))
  
  for(jj in 1:num.sim){
    
    annual.max.rep <- Gen.AM(nyears=nyears, af.para = af.para,dist.type=dist.type,kappa.para = kappa.para)
    annual.max.lmom[jj ,] <- samlmu(annual.max.rep,nmom = 4)
    
  }
  
  return(annual.max.lmom)
  
}   


Plot.Convergence <- function(num.sim=100, nyears=100, af.para = c(5,10), dist.type="uniform", kappa.para=c(0,1,0,1),Suppress.Plot = F){
  # Plot convergence to GEV from GPA
  # num.sim - number of simulations
  # nyears - number of years
  # af.para - parameters of the distribution to use  for arrival frequency (order of params as below):
  #       - fixed: c(value)
  #       - uniform: c(lower, upper)
  #       - binomial: c(n,p)
  #       - poisson: c(rate)
  #       - geometric: c(k,p)
  #
  #
  # dist.type - distribution type:
  #       - fixed
  #       - uniform (N.B. discrete uniform is used)
  #       - binomial
  #       - poisson
  #       - geometric
  #
  # kappa.para - Kappa parent parameters c(location, scale, shape, 2nd shape)
  # Suppress.Plot - If T do not plot the output
  
  AM.sim.lmom <- Repeat.AM.Simulations(num.sim=num.sim, nyears=nyears, af.para=af.para, dist.type=dist.type, kappa.para=kappa.para)
  kap.par <- pelkap(c(1,1,mean(AM.sim.lmom[,3]),mean(AM.sim.lmom[,4])))
  
  
  # Plot Lmoment diagram with GLO, GEV, GPA and Kappa (h=0.1)
  
  for(ii in 1:length(af.para)){
    
    if(ii==1){
      para.text <-af.para[ii]
    }else{
      para.text <- paste(para.text,af.para[ii],sep=", ")  
    }
  }
  
  if(Suppress.Plot==T){
    
    title0 <- paste0("GPA with ",dist.type," arrivals (",para.text,")"," converges to Kappa h = ",round(kap.par[4],3),"\n Num sim: ",num.sim," Num year per sim: ",nyears)
    info <- lmrd(distributions = "GLO GEV GPA",legend.lmrd = F,ylim = c(-0.1,0.4),main=title0)
    
    # Add Lmoments for 100 replicates
    kappa.lmom <- lmrkap(kappa.para, nmom=4)
    lmrdpoints(kappa.lmom[3], kappa.lmom[4], col="blue")
    
    lmrdpoints(AM.sim.lmom[,3],AM.sim.lmom[,4],col="black")
    lmrdpoints(mean(AM.sim.lmom[,3]),mean(AM.sim.lmom[,4]),col="red",pch=16)
    
    
    
  }
  
  return(c(mean(AM.sim.lmom[,3]),mean(AM.sim.lmom[,4])))
  
}

Gen.AM.POT.ns.np.alln.Binom.n.return.v2 <- function(nyears=100, start.year = 0, GPA.para = c(loc = 0, scale = 1, shape = 0), BN.para = c(np.int = 10, np.grad = 0, n.int = 10),Coles=T,do.POT.year = F){
  # Generate AM series
  # nyears - number of years
  # start.year -first year for regression purposes
  #
  # kappa.para - Kappa parent parameters in Coles notation c(location, scale, shape, 2nd shape intercept, 2nd shape slope)
  #
  
  np.int <- BN.para[1]
  np.grad <- BN.para[2]
  n.int <- BN.para[3]
  
  annual.max <- rep(0,nyears)
  n.events.store <- rep(0,nyears)
  POT <- vector(mode = "numeric",length = 0)
  
  if(do.POT.year==T){
  
	POT.year <- vector(mode = "numeric",length = 0)
  
  }
  
  
  
  p.int <- np.int/n.int
  
  if(Coles==T){
  
	# Invert the shape parameter to match Hoskign and Wallis notation (inputs are in Coles notation)
	GPA.para[3] <- -GPA.para[3]
  
  }
  
  
  # For each year generate n, the number of events (using Poisson distribution). THen generate n peaks using the GPA dist
  for(ii in 1:nyears){
    
    
    np.temp <- np.int + np.grad*(start.year+ii-1)
    n.events <-  qbinom(runif(1), round(np.temp/p.int), p.int)
    n.events.store[ii]<- n.events
    
    
    if(n.events != 0){ 
      
      max.temp <- quagpa(runif(n.events),GPA.para) 
      annual.max[ii] <- max(max.temp,GPA.para[1])
      
    } else {
      
      annual.max[ii] <- NA 
      
    }
    
    

    
    if(n.events != 0){ 
		POT[(length(POT)+1):(length(POT)+n.events)] <- max.temp 
		
		if(do.POT.year==T){
			POT.year[(length(POT.year)+1):(length(POT.year)+n.events)] <- rep(start.year+ii-1,n.events)
		}
		
		}
    
    
  }
  
  #annual.max <- annual.max[annual.max!=0]
  out <- list("AM" = data.frame(year = seq(0,nyears-1),n.count = n.events.store, maxima = annual.max),"POT"= data.frame(year=POT.year,POT=POT)) 
  return(out)
  
}