
# 16/11/2022 Declan O'Shea
# We provide this script, which is able to reproduce a simplified version of Figure 4d from our paper (but only showing the two-step Kappa and two-step GEV)

# Users should set the working directory to the file containing the modelling scripts
setwd("C:/Users/declano/OneDrive - The University of Melbourne/Documents/Analysis/Kappa_TwoStep_Calibration/Modelling_Functions/GIT")
source("./Stochastic_AM.R")
source("./One_step_Kappa_Nelder_Mead.R")
source("./Two_step_Kappa_or_GEV_Nelder_Mead_GLME.R")

require(lmom)
require(matrixStats)

# Set bounds on Generalised Pareto (GPD) parameters (location, scale, shape), using Coles (2001) parametrisation 
lower.gp <- c(1e-3,1e-3,-10)
upper.gp <- c(1e3,1e3,10)

# Set bounds on the following quantities of the Binomial distribution (np and h=1/n). Note we use h=11/n so that we can more easily implment an GMLE appraoch on h (not used in the paper)
lower.np <- c(1e-3,1e-3)
upper.np <- c(1e3,1)

# Set meta variables
n.it <- 100
nyears <- 100
AEPs <- c(2,5,10,20,50,100,200,500,1000,2000,5000,10000)

Fs <- 1-1/AEPs
tsk.qua <- matrix(NA,n.it,length(AEPs))
colnames(tsk.qua) <- AEPs

tsg.qua <- matrix(NA,n.it,length(AEPs))
colnames(tsg.qua) <- AEPs

for(ii in 1:n.it){
  
  cat("Iteration: ", ii,"\n")
  
  # Randomly generate AMS and POT from a coupled GPD - Binomial parameter
  data <- Gen.AM.POT.nsBinom.n.return(nyears=nyears, start.year = 0, GPA.para = c(loc = 20, scale = 14, shape = 0.1), BN.para = c(p = 0.5, n.int =  10, n.grad = 0))
  
  # Fit two-step GEV and calculate quantiles
  two.step.gev.out <- two.step.kappa.gmle.k.fit.POT.nm(data,lower.gp=lower.gp, upper.gp=upper.gp, lower.np=lower.np, upper.np=upper.np,fix.h.zero = T,num.it=10)
  pars.temp <- two.step.gev.out$Kappa.params.coles
  pars.temp[3] <- -pars.temp[3] # convert from Coles to HW notations
  tsg.qua[ii,] <- quagev(Fs,pars.temp[1:3]) 
  
  # Fit two-step Kappa and calculate quantiles
  two.step.kappa.out <- two.step.kappa.gmle.k.fit.POT.nm(data,lower.gp=lower.gp, upper.gp=upper.gp, lower.np=lower.np, upper.np=upper.np,num.it=10)
  pars.temp <- two.step.kappa.out$Kappa.params.coles
  pars.temp[3] <- -pars.temp[3] # convert from Coles to HW notations
  tsk.qua[ii,] <- quakap(Fs,pars.temp) 

  
  
}

# Calculate Kappa parent parameters
kappa.parent <- kappa.pars.5.to.4.Coles(c(20,14,0.1),c(5,10))
kappa.parent.HW <- kappa.parent
kappa.parent.HW[3] <- -kappa.parent.HW[3]

# Calculate confidence intervals
tsk.qua.ci <- colQuantiles(tsk.qua,probs = c(0.05,0.5,0.95))
tsg.qua.ci <- colQuantiles(tsg.qua,probs = c(0.05,0.5,0.95))

# Set some plotting parameters
ylim.max <- max(tsk.qua.ci)
red <- rgb(213/255,94/255,0/255,1)
l.blue <- rgb(86/255,180/255,233/255,1)

# Generate plot
plot(0,0,xlim = c(2.1,9.5),ylim = c(0,ylim.max),ylab = "Rainfall (mm)",xlab="Annual exceedance probability (1 in X)", xaxt='n')
gum <- -log(-log(Fs))
axis(1,gum,labels=AEPs)
title(paste("Frequency plot:",n.it,"simulations of",nyears,"years"))


# Plot median estimates
lines(gum,quakap(Fs,kappa.parent.HW),lwd = 3,col=black)
lines(gum,tsg.qua.ci[,2],type="o",pch=17,cex=1.2,lwd=3,col=l.blue)
lines(gum,tsk.qua.ci[,2],type="o",pch=18,cex=1.2,lwd=3,col=red)

# Plot confidence intervals
lines(gum,tsg.qua.ci[,1],lwd=1.5,lty=5,col=l.blue)
lines(gum,tsk.qua.ci[,1],lwd=1.5,lty=5,col=red)
lines(gum,tsg.qua.ci[,3],lwd=1.5,lty=5,col=l.blue)
lines(gum,tsk.qua.ci[,3],lwd=1.5,lty=5,col=red)

# Add legend
legend("topleft",inset=0.02,c("Parent","Two-step Kappa","Two-step GEV"),lty=c(1,1,1),pch=c(NA,18,17),lwd=c(2,2,2),col=c(black,red,l.blue))

