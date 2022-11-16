Gen.AM.POT.nsBinom.n.return <- function(nyears=100, start.year = 0, GPA.para = c(loc = 0, scale = 1, shape = 0), BN.para = c(p = 1, n.int = 10, n.grad = 0)){
  # Generate AM series from a coupled Generalised Pareto and Binomial parent
  #
  # Inputs
  # nyears - number of years
  # start.year -first year for regression purposes
  # GPA.para - GPA parent parameters in Coles notation c(location, scale, shape)
  # BN.para - Binomial parent parameters, n can be a linear function c(p, n.intercept, n.gradient)
  #
  # Outputs
  # annual.max - an nyears x 1 vector of stochastically generated maxima
  
  # Initialise output
  annual.max <- rep(0,nyears)
  n.events.store <- rep(0,nyears)
  POT <- vector(mode = "numeric",length = 0)
  
  # Invert the shape parameter to match Hosking and Wallis notation (inputs are in Coles notation)
  GPA.para[3] <- -GPA.para[3]
  
  # For each year generate n.events, the number of events in that year (using Binomial distribution). Then generate n.events peaks using the GPA distribution. Store the maximum of that year as in the annual maxima vector.
  # Also store all the generated peaks in the POT vector
 for(ii in 1:nyears){
    n.temp <- round(BN.para[2] + BN.para[3]*(start.year+ii-1))
    n.events <-  qbinom(runif(1),n.temp, BN.para[1])
    n.events.store[ii]<- n.events
    
    # If no events in a year, then set annual max to NA
    if(n.events != 0){ 
      
      max.temp <- quagpa(runif(n.events),GPA.para) 
      annual.max[ii] <- max(max.temp,GPA.para[1])
    
    } else {
      
      annual.max[ii] <- NA 
      
      }

    if(n.events != 0){ POT[(length(POT)+1):(length(POT)+n.events)] <- max.temp }
    
    
  }
  
  # Set negative infinity annual max to NA
  annual.max[annual.max == -Inf] = NA
  
  # Make an output list. Make a dataframe for the year, annual maxima and n.count. Also output the POT peaks in the same list
  out <- list("AM" = data.frame(year = seq(0,nyears-1),n.count = n.events.store, maxima = annual.max),"POT"= POT )
  return(out)
  
}
