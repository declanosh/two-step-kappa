
# 08/11/2022 Declan O'Shea
#We provide this script as a very rough (but working) example of the two-step Kappa prior to our paper resubmission date. Due to time constraints we could not provide a more polished version of the code in time. However, we plan to rectify this in the coming days.

# Users should set the working directory to the file containing the modelling scripts
setwd("C:/Users/declano/OneDrive - The University of Melbourne/Documents/Analysis/Kappa_TwoStep_Calibration/Modelling_Functions/GIT")
source("./Stochastic_AM.R")
source("./One_step_Kappa_Nelder_Mead.R")
source("./Two_step_Kappa_or_GEV_Nelder_Mead_GLME.R")

# Set bounds on Generalised Pareto (GPD) parameters (location, scale, shape), using Coles (2001) parametrisation 
lower.gp <- c(1e-3,1e-3,-10)
upper.gp <- c(1e3,1e3,10)

# Set bounds on the following quantities of the Binomial distribution (np and h=1/n). Note we use h=11/n so that we can more easily implment an GMLE appraoch on h (not used in the paper)
lower.np <- c(1e-3,1e-3)
upper.np <- c(1e3,1)


# Randomly generate AMS and POT from a coupled GPD - Binomial parameter
data <- Gen.AM.POT.nsBinom.n.return(nyears=100, start.year = 0, GPA.para = c(loc = 20, scale = 14, shape = 0.1), BN.para = c(p = 0.5, n.int =  10, n.grad = 0))
  

# Fit two-step GEV
two.step.gev.out <- two.step.kappa.gmle.k.fit.POT.nm(data,lower.gp=lower.gp, upper.gp=upper.gp, lower.np=lower.np, upper.np=upper.np,fix.h.zero = T,num.it=10)
      
# Fit two-step Kappa
 two.step.kappa.out <- two.step.kappa.gmle.k.fit.POT.nm(data,lower.gp=lower.gp, upper.gp=upper.gp, lower.np=lower.np, upper.np=upper.np,num.it=10)
     