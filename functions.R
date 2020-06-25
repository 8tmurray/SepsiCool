### Functions for SepsiCool Trial Design
### Author: Thomas A Murray
### e-mail: murra484@umn.edu


### Hazard estimates for the historical control arm on days t=1,...,60
hazC0 = c(0.0313,0.0900,0.0409,0.0164,0.0227,0.0196,0.0131,0.0196,0.0204,0.0177,
          0.0157,0.0107,0.0066,0.0089,0.0088,0.0100,0.0067,0.0053,0.0086,0.0093,
          0.0053,0.0081,0.0065,0.0065,0.0071,0.0116,0.0070,0.0061,0.0047,0.0058,
          0.0062,0.0075,0.0067,0.0028,0.0015,0.0014,0.0023,0.0022,0.0021,0.0022,
          0.0023,0.0012,0.0010,0.0015,0.0021,0.0031,0.0024,0.0058,0.0038,0.0050,
          0.0030,0.0025,0.0018,0.0011,0.0011,0.0017,0.0020,0.0017,0.0039,0.0017)



### Function to generate a time to death given the hazard for death on days 1,2,... 
# Input: haz - hazards, i.e., Prob(Y = t | Y >= t), for t=1,2,...
# Output: y* where y* = length(haz)+1 implies y* > length(haz)
get.response <- function(haz){
  y.star = which(rbinom(length(haz)+1,1,c(haz,1))==1)[1]
  return(y.star)
}




###Function to Get Data on Concurrent Arms
# Input: nobs - number of observations
#        hazE, hazC - hazards in each concurrent arm
# Output: dat - data frame with variables y*, arm (1=E, 2=C), enroll (enrollment time in calendar day)
library(dplyr)
get.current.data <- function(nobs,hazE,hazC){
  arm = as.vector(sapply(1:(nobs/4),function(p) sample(c("1","1","2","2"))))
  y.star = sapply(arm,function(a) ifelse(a=="1",get.response(hazE),get.response(hazC)))
  enroll = sort(sample(1:855,size=nobs,replace=TRUE))
  dat = tibble(y.star=y.star,arm=arm,enroll=enroll)
  return(dat)
}





### Function to Calculate Log-Rank Test Statistic based on Current Data
library(survival)
get.lrt.stat <- function(dat){
  fit.lrt <- coxph(Surv(y,delta)~arm,data=dat)
  lrt.stat = as.vector(fit.lrt$coef/sqrt(fit.lrt$var))
  return(lrt.stat)
}





# Function to Calculate Prob(RMS(E) > RMS(C) | Curr Data)
library(BayesLogit)
get.curr.only.stats <- function(dat,maxt=60,nsamps=4000,warmup=200){
  ### Data pre-preprocessing
  y = dat$y; delta = dat$delta; arm = dat$arm
  
  # Design matrix (n x p)
  X = rbind(cbind(1,diag(maxt),-0.5,-0.5*diag(maxt)),
            cbind(1,diag(maxt),0.5,0.5*diag(maxt)))
  
  # Number at Risk (n x 1)
  r = rbind(cbind(sapply(1:maxt,function(t) sum(y>=t & arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(y>=t & arm=="2"))))
  
  # Number of Events (n x 1)
  d = rbind(cbind(sapply(1:maxt,function(t) sum(y==t & delta==1 & arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(y==t & delta==1 & arm=="2"))))
  
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(d-0.5*r)
  
  # AR(1) Precision Matrix
  A = diag(rep(0,maxt)); A[1,2] = A[maxt,maxt-1] = 1
  for(t in 2:(maxt-1)) A[t,t-1] = A[t,t+1] = 1
  Q = diag(rowSums(A))-A
  
  # Hyperprior specification, Ghosh defaults
  nu.alpha = nu.a = nu.beta = nu.b = 5 
  s.alpha = s.beta = 5
  s.a = s.b = 2.5
  
  # Storage objects
  rmsE.samps = rmsC.samps = rep(NA,nsamps)
  
  # Initialize Hyperparameter values
  tau.alpha = tau.a = tau.beta = tau.b = 1
  
  # Initialize regression coefficients
  reg.coef = cbind(c(log(mean(d/r))-log((1-mean(d/r))),rep(0,2*maxt+1)))
  
  # Gibbs Sampler    
  for(samp in 1:(nsamps+warmup)){
    # Sample Regression Coefficients
    P0.alpha = cbind(0,rbind(0,tau.a*Q)); P0.alpha[1,1] = tau.alpha
    P0.beta = cbind(0,rbind(0,tau.b*Q)); P0.beta[1,1] = tau.beta
    P0.fill = matrix(0,nrow=nrow(P0.alpha),ncol=ncol(P0.beta))
    P0 = rbind(cbind(P0.alpha,P0.fill),cbind(t(P0.fill),P0.beta)) # Prior Precision
    omega = rpg.devroye(nrow(X),r,X%*%reg.coef) # update auxillary parameters
    V = chol2inv(chol(t(X)%*%(X*omega) + P0)) # Variance matrix
    m = V%*%(t(X)%*%kappa)                    # Mean vector
    reg.coef = m + t(chol(V))%*%rnorm(ncol(X))
    
    # Sample hyperparameters
    alpha = reg.coef[1,]; a = cbind(reg.coef[1:maxt+1,])
    beta = reg.coef[maxt+2,]; b = cbind(reg.coef[1:maxt+maxt+2,])
    tau.alpha = rgamma(1,shape=(nu.alpha+1)/2,rate=(nu.alpha*s.alpha^2+alpha^2)/2)
    tau.a = rgamma(1,shape=(nu.a+maxt-1)/2,rate=(nu.a*s.a^2+t(a)%*%Q%*%a)/2)
    tau.beta = rgamma(1,shape=(nu.beta+1)/2,rate=(nu.beta*s.beta^2+beta^2)/2)
    tau.b = rgamma(1,shape=(nu.b+maxt-1)/2,rate=(nu.b*s.b^2+t(b)%*%Q%*%b)/2)
    
    # transform to difference in rms and save
    haz = 1/(1+exp(-X%*%reg.coef)); hazE = haz[1:maxt]; hazC = haz[1:maxt+maxt]
    if(samp>warmup){
      rmsE.samps[samp-warmup] = -sum(diff(c(1,cumprod(1-hazE),0))*(0:maxt))
      rmsC.samps[samp-warmup] = -sum(diff(c(1,cumprod(1-hazC),0))*(0:maxt))
    } 
  }
  
  # Calculate Prob(RMS(E) > RMS(C) | Data)
  stats = c(mean(rmsE.samps > rmsC.samps),mean(rmsE.samps),mean(rmsC.samps),var(rmsE.samps-rmsC.samps))
  return(stats)
}







# Function to Calculate Prob(RMS(E) > RMS(C) | Curr and Hist Data)
library(BayesLogit)
get.curr.hist.stats <- function(dat,maxt=60,nsamps=4000,warmup=200,spike.gamma=4000,p.gamma=0.5){
  ### Data pre-preprocessing
  y = dat$y; delta = dat$delta; arm = dat$arm
  
  # Design matrix (n x p)
  X = rbind(cbind(1,diag(maxt),1,1*diag(maxt),0),
            cbind(1,diag(maxt),0,0*diag(maxt),0),
            cbind(1,diag(maxt),0,0*diag(maxt),1))
  
  # Number at Risk (n x 1)
  r = rbind(cbind(sapply(1:maxt,function(t) sum(y>=t & arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(y>=t & arm=="2"))),
            cbind(sapply(1:maxt,function(t) sum(y>=t & arm=="3"))))
  
  # Number of Events (n x 1)
  d = rbind(cbind(sapply(1:maxt,function(t) sum(y==t & delta==1 & arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(y==t & delta==1 & arm=="2"))),
            cbind(sapply(1:maxt,function(t) sum(y==t & delta==1 & arm=="3"))))
  
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(d-0.5*r)
  
  # AR(1) Precision Matrix
  A = diag(rep(0,maxt)); A[1,2] = A[maxt,maxt-1] = 1
  for(t in 2:(maxt-1)) A[t,t-1] = A[t,t+1] = 1
  Q = diag(rowSums(A))-A
  
  # Hyperprior specification
  nu.alpha = nu.a = nu.beta = nu.b = nu.gamma = 5 
  s.alpha = 5; s.a = s.b = s.beta = s.gamma = 2.5
  spike.gamma = spike.gamma; p.gamma = p.gamma
  
  # Storage objects
  rmsE.samps = rmsC.samps = iota.gamma.samps = rep(NA,nsamps)
  
  # Initialize Hyperparameter values
  tau.alpha = tau.a = tau.beta = tau.b = tau.gamma = 1
  iota.gamma = rbinom(1,1,p.gamma)
  
  # Initialize regression coefficients
  reg.coef = cbind(c(log(mean(d/r))-log((1-mean(d/r))),rep(0,2*maxt+1),0))
  
  # Gibbs Sampler    
  for(samp in 1:(nsamps+warmup)){
    # Sample Regression Coefficients
    P0.alpha = cbind(0,rbind(0,tau.a*Q)); P0.alpha[1,1] = tau.alpha
    P0.beta = cbind(0,rbind(0,tau.b*Q)); P0.beta[1,1] = tau.beta
    P0.gamma = cbind(iota.gamma*spike.gamma+(1-iota.gamma)*tau.gamma)
    P0.ab.fill = matrix(0,nrow=nrow(P0.alpha),ncol=ncol(P0.beta))
    P0.ag.fill = matrix(0,nrow=nrow(P0.alpha),ncol=ncol(P0.gamma))
    P0.bg.fill = matrix(0,nrow=nrow(P0.beta),ncol=ncol(P0.gamma))
    P0 = rbind(cbind(P0.alpha,P0.ab.fill,P0.ag.fill),
               cbind(t(P0.ab.fill),P0.beta,P0.bg.fill),
               cbind(t(P0.ag.fill),t(P0.bg.fill),P0.gamma))
    omega = rpg.devroye(nrow(X),r,X%*%reg.coef) # update auxillary parameters
    V = chol2inv(chol(t(X)%*%(X*omega) + P0)) # Variance matrix
    m = V%*%(t(X)%*%kappa)                    # Mean vector
    reg.coef = m + t(chol(V))%*%rnorm(ncol(X))
    
    # Sample hyperparameters
    alpha = reg.coef[1,]; a = cbind(reg.coef[1:maxt+1,])
    tau.alpha = rgamma(1,shape=(nu.alpha+1)/2,rate=(nu.alpha*s.alpha^2+alpha^2)/2)
    tau.a = rgamma(1,shape=(nu.a+maxt-1)/2,rate=(nu.a*s.a^2+t(a)%*%Q%*%a)/2)
    
    beta = reg.coef[maxt+2,]; b = cbind(reg.coef[1:maxt+maxt+2,])
    tau.beta = rgamma(1,shape=(nu.beta+1)/2,rate=(nu.beta*s.beta^2+beta^2)/2)
    tau.b = rgamma(1,shape=(nu.b+maxt-1)/2,rate=(nu.b*s.b^2+t(b)%*%Q%*%b)/2)
    
    gamma = reg.coef[3+2*maxt,]
    tau.gamma = rgamma(1,shape=(nu.gamma+iota.gamma-1)/2,rate=(nu.gamma*s.gamma^2+(1-iota.gamma)*gamma^2)/2)
    iota.gamma = rbinom(1,1,p.gamma*dnorm(gamma,0,1/sqrt(spike.gamma))/(p.gamma*dnorm(gamma,0,1/sqrt(spike.gamma)) + (1-p.gamma)*dnorm(gamma,0,1/sqrt(tau.gamma))))
    
    # transform to difference in rms and save
    haz = 1/(1+exp(-X%*%reg.coef)); hazE = haz[1:maxt]; hazC = haz[1:maxt+maxt]
    rmsE = -sum(diff(c(1,cumprod(1-hazE),0))*(0:maxt)); rmsC = -sum(diff(c(1,cumprod(1-hazC),0))*(0:maxt))
    if(samp>warmup){
      rmsE.samps[samp-warmup] = rmsE
      rmsC.samps[samp-warmup] = rmsC
      iota.gamma.samps[samp-warmup] = iota.gamma
    } 
  }
  
  # Calculate Prob(RMS(E) > RMS(C) | Data)
  stats = c(mean(rmsE.samps > rmsC.samps),mean(rmsE.samps),mean(rmsC.samps),var(rmsE.samps-rmsC.samps),mean(iota.gamma.samps))
  return(stats)
}
