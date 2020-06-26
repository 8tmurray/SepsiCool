### Functions for SepsiCool Trial Design
### Written by Thomas A Murray on 12/21/2017
### Last Modified on 01/13/2020
library(BayesLogit); library(dplyr); library(survival); library(Matrix)

### Hazard estimates for the historical control arm from t=1,...,60
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
get.current.data <- function(nobs,hazE,hazC){
  arm = as.vector(sapply(1:(nobs/4),function(p) sample(c("1","1","2","2"))))
  y.star = sapply(arm,function(a) ifelse(a=="1",get.response(hazE),get.response(hazC)))
  enroll = sort(sample(1:855,size=nobs,replace=TRUE))
  dat = tibble(y.star=y.star,arm=arm,enroll=enroll)
  return(dat)
}





### Function to Calculate Log-Rank Test Statistic based on Current Data
get.lrt.stat <- function(dat){
  fit.lrt <- coxph(Surv(y,delta)~arm,data=dat)
  lrt.stat = as.vector(fit.lrt$coef/sqrt(fit.lrt$var))
  return(lrt.stat)
}




### Non-Proportional Hazards/Odds Commensurate Prior Model
npo.comm.sampler <- function(dat,maxt=60,nsamps=2000,warmup=200,p.zeta=0.5,R.zeta=4000,p.z=0.5,R.z=4000){

  # Design matrix (n x p)
  X = rbind(cbind(1,diag(maxt),1,1*diag(maxt),0,0*diag(maxt)), #Current Experimental Arm (arm = "1")
            cbind(1,diag(maxt),0,0*diag(maxt),0,0*diag(maxt)), #Current Control Arm (arm = "2")
            cbind(1,diag(maxt),0,0*diag(maxt),1,1*diag(maxt))) #Historical Control Arm (arm = "3")
  
  # Number at Risk (n x 1)
  r = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="2"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="3"))))
  
  # Number of Events (n x 1)
  d = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="2"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="3"))))
  
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(d-0.5*r)
  
  # AR(1) Precision Matrix for GMRF Priors (Same for all three GMRF priors)
  Q = matrix(0,maxt,maxt); diag(Q) = c(1,rep(2,maxt-2),1)
  Q[cbind(2:maxt,2:maxt-1)] = Q[cbind(2:maxt-1,2:maxt)] = -1
  
  # Constraint Matrix (each GMRF vector must sum to zero)
  A = as.matrix(t(bdiag(c(0,rep(1,maxt)),c(0,rep(1,maxt)),c(0,rep(1,maxt)))))
  
  # Storage object
  samps = matrix(NA,nrow=nsamps,ncol=6)

  # Initial values
  beta = rep(0,ncol(X)); 
  tau.gamma = tau.lambda = tau.zeta = 1
  tau.g = tau.l = tau.z = R.z
  nu.zeta = nu.z = 1
  tau.zeta.star = R.zeta; tau.z.star = R.z

  # Gibbs Sampler    
  for(samp in 1:(nsamps+warmup)){
    # Sample parameters
    omega = rpg.devroye(nrow(X),r,X%*%beta) # Generate Polya-gamma auxillary parameters
    m0 = cbind(c(-5,rep(0,ncol(X)-1))) #Prior Mean
    P0 = as.matrix(bdiag(tau.gamma,tau.g*Q,tau.lambda,tau.l*Q,tau.zeta.star,tau.z.star*Q)) #Prior Precision
    P = crossprod(X,X*omega) + P0 #Full Conditional Posterior Precision
    V = chol2inv(chol(P)) #Full Conditional Posterior Covariance
    m = V%*%(crossprod(X,kappa)+P0%*%m0) # Full Conditional Posterior Mean
    beta.temp = m + t(chol(V))%*%rnorm(ncol(X)) #sample beta without constraints
    beta = beta.temp - V%*%t(A)%*%chol2inv(chol(A%*%V%*%t(A)))%*%(A%*%beta.temp) #apply constraints (pg 37, eqn 2.30 in Rue and Held)
    gamma = beta[1,]; g = cbind(beta[1:maxt+1,]); lambda = beta[1+(maxt+1),]; l = cbind(beta[1:maxt+1+(maxt+1),]); zeta = beta[1+2*(maxt+1),]; z = cbind(beta[1:maxt+1+2*(maxt+1),])

    # Sample hyperparameters
    tau.gamma = rgamma(1,shape=(7+1)/2,rate=(7*5^2+(gamma+5)^2)/2)
    tau.g = rgamma(1,shape=(1+maxt-1)/2,rate=(1*2.5^2/(maxt-1)+crossprod(g,Q%*%g))/2)
    tau.lambda = rgamma(1,shape=(7+1)/2,rate=(7*2.5^2+lambda^2)/2)
    tau.l = rgamma(1,shape=(1+maxt-1)/2,rate=(1*2.5^2/(maxt-1)+crossprod(l,Q%*%l))/2)
    tau.zeta = rgamma(1,shape=(7+(1-nu.zeta))/2,rate=(7*2.5^2+(1-nu.zeta)*zeta^2)/2)
    nu.zeta = rbinom(1,1,1/(1+(1-p.zeta)/p.zeta*exp(log(tau.zeta/R.zeta)/2-(tau.zeta-R.zeta)*zeta^2/2)))
    tau.zeta.star = nu.zeta*R.zeta + (1-nu.zeta)*tau.zeta
    tau.z = rgamma(1,shape=(1+(1-nu.z)*(maxt-1))/2,rate=(1/R.z+(1-nu.z)*crossprod(z,Q%*%z))/2)
    nu.z = rbinom(1,1,1/(1+(1-p.z)/p.z*exp(log(tau.z/R.z)*(maxt-1)/2-(tau.z-R.z)*crossprod(z,Q%*%z)/2)))
    tau.z.star = nu.z*R.z + (1-nu.z)*tau.z

    # Calculate and store posterior summaries
    if(samp>warmup){
      haz = 1/(1+exp(-X%*%beta)); hazE = haz[1:maxt]; hazC = haz[maxt+1:maxt]; hazC0 = haz[2*maxt+1:maxt] #Hazards
      resE = c(-sum(diff(c(1,cumprod(1-hazE)))*(1:maxt-1)),prod(1-hazE))
      resC = c(-sum(diff(c(1,cumprod(1-hazC)))*(1:maxt-1)),prod(1-hazC))
      samps[samp-warmup,] = c(resE,resC,nu.zeta,1*(tau.z.star>=R.z))
    } 
  }
  
  #Return relevant samples of posterior summaries
  return(samps)
}


### Proportional Hazards/Odds Commensurate Prior Model
po.comm.sampler <- function(dat,maxt=60,nsamps=2000,warmup=200,p.zeta=0.5,R.zeta=4000){
  
  # Design matrix (n x p)
  X = rbind(cbind(1,diag(maxt),1,1*diag(maxt),0), #Current Experimental Arm (arm = "1")
            cbind(1,diag(maxt),0,0*diag(maxt),0), #Current Control Arm (arm = "2")
            cbind(1,diag(maxt),0,0*diag(maxt),1)) #Historical Control Arm (arm = "3")
  # Number at Risk (n x 1)
  r = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="2"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="3"))))
  
  # Number of Events (n x 1)
  d = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="2"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="3"))))
  
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(d-0.5*r)
  
  # AR(1) Precision Matrix for GMRF Priors (Same for all three GMRF priors)
  Q = matrix(0,maxt,maxt); diag(Q) = c(1,rep(2,maxt-2),1)
  Q[cbind(2:maxt,2:maxt-1)] = Q[cbind(2:maxt-1,2:maxt)] = -1
  
  # Constraint Matrix (each GMRF vector must sum to zero)
  A = as.matrix(t(bdiag(c(0,rep(1,maxt)),c(0,rep(1,maxt)),0)))[1:2,]
  
  # Storage object
  samps = matrix(NA,nrow=nsamps,ncol=5)
  
  # Initial values
  beta = rep(0,ncol(X)); 
  tau.gamma = tau.lambda = tau.zeta = 1
  tau.g = tau.l = R.zeta
  nu.zeta = 1
  tau.zeta.star = R.zeta
  
  # Gibbs Sampler    
  for(samp in 1:(nsamps+warmup)){
    # Sample parameters
    omega = rpg.devroye(nrow(X),r,X%*%beta) # Generate Polya-gamma auxillary parameters
    m0 = cbind(c(-5,rep(0,ncol(X)-1))) #Prior Mean
    P0 = as.matrix(bdiag(tau.gamma,tau.g*Q,tau.lambda,tau.l*Q,tau.zeta.star)) #Prior Precision
    P = crossprod(X,X*omega) + P0 #Full Conditional Posterior Precision
    V = chol2inv(chol(P)) #Full Conditional Posterior Covariance
    m = V%*%(crossprod(X,kappa)+P0%*%m0) # Full Conditional Posterior Mean
    beta.temp = m + t(chol(V))%*%rnorm(ncol(X)) #sample beta without constraints
    beta = beta.temp - V%*%t(A)%*%chol2inv(chol(A%*%V%*%t(A)))%*%(A%*%beta.temp) #apply constraints (pg 37, eqn 2.30 in Rue and Held)
    gamma = beta[1,]; g = cbind(beta[1:maxt+1,]); lambda = beta[1+(maxt+1),]; l = cbind(beta[1:maxt+1+(maxt+1),]); zeta = beta[1+2*(maxt+1),]

    # Sample hyperparameters
    tau.gamma = rgamma(1,shape=(7+1)/2,rate=(7*5^2+(gamma+5)^2)/2)
    tau.g = rgamma(1,shape=(1+maxt-1)/2,rate=(1*2.5^2/(maxt-1)+crossprod(g,Q%*%g))/2)
    tau.lambda = rgamma(1,shape=(7+1)/2,rate=(7*2.5^2+lambda^2)/2)
    tau.l = rgamma(1,shape=(1+maxt-1)/2,rate=(1*2.5^2/(maxt-1)+crossprod(l,Q%*%l))/2)
    tau.zeta = rgamma(1,shape=(7+(1-nu.zeta))/2,rate=(7*2.5^2+(1-nu.zeta)*zeta^2)/2)
    nu.zeta = rbinom(1,1,1/(1+(1-p.zeta)/p.zeta*exp(log(tau.zeta/R.zeta)/2-(tau.zeta-R.zeta)*zeta^2/2)))
    tau.zeta.star = nu.zeta*R.zeta + (1-nu.zeta)*tau.zeta
    
    # Calculate and store posterior summaries
    if(samp>warmup){
      haz = 1/(1+exp(-X%*%beta)); hazE = haz[1:maxt]; hazC = haz[maxt+1:maxt]; hazC0 = haz[2*maxt+1:maxt] #Hazards
      resE = c(-sum(diff(c(1,cumprod(1-hazE)))*(1:maxt-1)),prod(1-hazE))
      resC = c(-sum(diff(c(1,cumprod(1-hazC)))*(1:maxt-1)),prod(1-hazC))
      samps[samp-warmup,] = c(resE,resC,nu.zeta)
    } 
  }
  
  #Return relevant samples of posterior summaries
  return(samps)
}






### Non-proportional Hazards/Odds Model
npo.trad.sampler <- function(dat,maxt=60,nsamps=2000,warmup=200){
  
  # Design matrix (n x p)
  X = rbind(cbind(1,diag(maxt),1,1*diag(maxt)), #Current Experimental Arm (arm = "1")
            cbind(1,diag(maxt),0,0*diag(maxt))) #Current Control Arm (arm = "2")
  
  # Number at Risk (n x 1)
  r = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="2"))))
  
  # Number of Events (n x 1)
  d = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="2"))))
  
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(d-0.5*r)
  
  # AR(1) Precision Matrix for GMRF Priors (Same for all three GMRF priors)
  Q = matrix(0,maxt,maxt); diag(Q) = c(1,rep(2,maxt-2),1)
  Q[cbind(2:maxt,2:maxt-1)] = Q[cbind(2:maxt-1,2:maxt)] = -1
  
  # Constraint Matrix (each GMRF vector must sum to zero)
  constr = as.matrix(t(bdiag(c(0,rep(1,maxt)),c(0,rep(1,maxt)))))
  
  # Storage object
  samps = matrix(NA,nrow=nsamps,ncol=4)
  
  # Initial values
  beta = rep(0,ncol(X)); 
  tau.gamma = tau.lambda = 1
  tau.g = tau.l = 1
  
  # Gibbs Sampler    
  for(samp in 1:(nsamps+warmup)){
    # Sample parameters
    omega = rpg.devroye(nrow(X),r,X%*%beta) # Generate Polya-gamma auxillary parameters
    m0 = cbind(c(-5,rep(0,ncol(X)-1))) #Prior Mean
    P0 = as.matrix(bdiag(tau.gamma,tau.g*Q,tau.lambda,tau.l*Q)) #Prior Precision
    P = crossprod(X,X*omega) + P0 #Full Conditional Posterior Precision
    V = chol2inv(chol(P)) #Full Conditional Posterior Covariance
    m = V%*%(crossprod(X,kappa)+P0%*%m0) # Full Conditional Posterior Mean
    beta.temp = m + t(chol(V))%*%rnorm(ncol(X)) #sample beta without constraints
    beta = beta.temp - V%*%t(constr)%*%chol2inv(chol(constr%*%V%*%t(constr)))%*%(constr%*%beta.temp) #apply constraints (pg 37, eqn 2.30 in Rue and Held)
    gamma = beta[1,]; g = cbind(beta[1:maxt+1,]); lambda = beta[1+(maxt+1),]; l = cbind(beta[1:maxt+1+(maxt+1),])
    
    # Sample hyperparameters
    tau.gamma = rgamma(1,shape=(7+1)/2,rate=(7*5^2+(gamma+5)^2)/2)
    tau.g = rgamma(1,shape=(1+maxt-1)/2,rate=(1*2.5^2/(maxt-1)+crossprod(g,Q%*%g))/2)
    tau.lambda = rgamma(1,shape=(7+1)/2,rate=(7*2.5^2+lambda^2)/2)
    tau.l = rgamma(1,shape=(1+maxt-1)/2,rate=(1*2.5^2/(maxt-1)+crossprod(l,Q%*%l))/2)
    
    # Calculate and store posterior summaries
    if(samp>warmup){
      haz = 1/(1+exp(-X%*%beta)); hazE = haz[1:maxt]; hazC = haz[maxt+1:maxt]; hazC0 = haz[2*maxt+1:maxt] #Hazards
      resE = c(-sum(diff(c(1,cumprod(1-hazE)))*(1:maxt-1)),prod(1-hazE))
      resC = c(-sum(diff(c(1,cumprod(1-hazC)))*(1:maxt-1)),prod(1-hazC))
      samps[samp-warmup,] = c(resE,resC)
    } 
  }
  
  #Return relevant samples of posterior summaries
  return(samps)
}




#Proportional Hazards/Odds Model
po.trad.sampler <- function(dat,maxt=60,nsamps=2000,warmup=200){
  
  # Design matrix (n x p)
  X = rbind(cbind(1,diag(maxt),1), #Current Experimental Arm (arm = "1")
            cbind(1,diag(maxt),0)) #Current Control Arm (arm = "2")
  
  # Number at Risk (n x 1)
  r = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y>=t & dat$arm=="2"))))
  
  # Number of Events (n x 1)
  d = rbind(cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="1"))),
            cbind(sapply(1:maxt,function(t) sum(dat$y==t & dat$delta==1 & dat$arm=="2"))))
  
  # Pseudo-outcome vector (n x 1)
  kappa = cbind(d-0.5*r)
  
  # AR(1) Precision Matrix for GMRF Priors (Same for all three GMRF priors)
  Q = matrix(0,maxt,maxt); diag(Q) = c(1,rep(2,maxt-2),1)
  Q[cbind(2:maxt,2:maxt-1)] = Q[cbind(2:maxt-1,2:maxt)] = -1
  
  # Constraint Matrix (each GMRF vector must sum to zero)
  constr = rbind(c(0,rep(1,maxt),0))
  
  # Storage object
  samps = matrix(NA,nrow=nsamps,ncol=4)
  
  # Initial values
  beta = rep(0,ncol(X)); 
  tau.gamma = tau.lambda = 1
  tau.g = 1
  
  # Gibbs Sampler    
  for(samp in 1:(nsamps+warmup)){
    # Sample parameters
    omega = rpg.devroye(nrow(X),r,X%*%beta) # Generate Polya-gamma auxillary parameters
    m0 = cbind(c(-5,rep(0,ncol(X)-1))) #Prior Mean
    P0 = as.matrix(bdiag(tau.gamma,tau.g*Q,tau.lambda)) #Prior Precision
    P = crossprod(X,X*omega) + P0 #Full Conditional Posterior Precision
    V = chol2inv(chol(P)) #Full Conditional Posterior Covariance
    m = V%*%(crossprod(X,kappa)+P0%*%m0) # Full Conditional Posterior Mean
    beta.temp = m + t(chol(V))%*%rnorm(ncol(X)) #sample beta without constraints
    beta = beta.temp - V%*%t(constr)%*%chol2inv(chol(constr%*%V%*%t(constr)))%*%(constr%*%beta.temp) #apply constraints (pg 37, eqn 2.30 in Rue and Held)
    gamma = beta[1,]; g = cbind(beta[1:maxt+1,]); lambda = beta[1+(maxt+1),]

    # Sample hyperparameters
    tau.gamma = rgamma(1,shape=(7+1)/2,rate=(7*5^2+(gamma+5)^2)/2)
    tau.g = rgamma(1,shape=(1+maxt-1)/2,rate=(1*2.5^2/(maxt-1)+crossprod(g,Q%*%g))/2)
    tau.lambda = rgamma(1,shape=(7+1)/2,rate=(7*2.5^2+lambda^2)/2)
    
    # Calculate and store posterior summaries
    if(samp>warmup){
      haz = 1/(1+exp(-X%*%beta)); hazE = haz[1:maxt]; hazC = haz[maxt+1:maxt]; hazC0 = haz[2*maxt+1:maxt] #Hazards
      resE = c(-sum(diff(c(1,cumprod(1-hazE)))*(1:maxt-1)),prod(1-hazE))
      resC = c(-sum(diff(c(1,cumprod(1-hazC)))*(1:maxt-1)),prod(1-hazC))
      samps[samp-warmup,] = c(resE,resC)
    }
  }
  
  #Return relevant samples of posterior summaries
  return(samps)
}