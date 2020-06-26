### Code for Simulation to Assess OCs in Treatment Effect and Historical Bias
### By Thomas A Murray
#setwd("Path to Working Directory")#Not run
source('functions.R')

### Load Historical Data
dat.hist = read.csv('Historical-Control-Data.csv',header=TRUE) %>% as_tibble()

### Function to simulate one trial under each of the possible methods
simulate.trial = function(trial.effect=0,treatment.effect=1,nobs=956,analysis.days=c(305,610,915)){

  # Determine hazards in E and C based on trial effect and treatment effect
  hazC = exp(trial.effect)*hazC0/(1-hazC0)/(1+exp(trial.effect)*hazC0/(1-hazC0))
  #trial.effects = 5*60:1/60 - mean(5*60:1/60) #Checking whether GMRF spike and slab can distinguish a strong NPO Historical Bias
  #hazC = exp(trial.effects)*hazC0/(1-hazC0)/(1+exp(trial.effects)*hazC0/(1-hazC0))
  rrr = treatment.effect*c(rep(0.4,15),0.4-0.2*(1:15)/15,rep(0.2,30))
  hazE = hazC0; hazE[1] = (1-rrr[1])*hazC[1]
  for(t in 2:60) hazE[t] = 1-(1-(1-rrr[t])*(1-prod(1-hazC[1:t])))/prod(1-hazE[1:(t-1)])
  rmsE = round(-sum(diff(c(1,cumprod(1-hazE),0))*(0:60)),4)
  rmsC = round(-sum(diff(c(1,cumprod(1-hazC),0))*(0:60)),4)

  # Generate Data from Concurrent Arms
  dat.all.curr = get.current.data(nobs=nobs,hazE=hazE,hazC=hazC)

  #Storage Objects
  lrt.stats = matrix(NA,nrow=1,ncol=length(analysis.days))
  po.trad.stats = npo.trad.stats = matrix(NA,nrow=4,ncol=length(analysis.days))
  po.comm.stats = matrix(NA,nrow=5,ncol=length(analysis.days))
  npo.comm.stats = matrix(NA,nrow=7,ncol=length(analysis.days))

  # Posterior summary at each analysis
  for(analysis in 1:length(analysis.days)){
    dat.curr = filter(dat.all.curr,enroll<=analysis.days[analysis]) %>% mutate(y = ifelse(y.star==61,60,pmin(y.star,analysis.days[analysis] - enroll + 1)),delta = 1*(y.star<=analysis.days[analysis]-enroll+1 & y.star<=60)) %>% select(y,delta,arm)
    dat = rbind(dat.curr,dat.hist)
    
    # Calculate log-rank test statistic (z-score)
    lrt.stats[,analysis] = get.lrt.stat(dat=dat.curr)
    
    # Calculate Prob(rmsE > rmsC | Data) under each Bayesian model
    po.trad.samps <- apply(po.trad.sampler(dat=dat),1,function(x) c(x[1]+60*x[2],x[3]+60*x[4]))
    po.trad.stats[,analysis] = c(mean(po.trad.samps[1,]>po.trad.samps[2,]),var(po.trad.samps[1,]-po.trad.samps[2,]),rowMeans(po.trad.samps))

    npo.trad.samps <- apply(npo.trad.sampler(dat=dat),1,function(x) c(x[1]+60*x[2],x[3]+60*x[4]))
    npo.trad.stats[,analysis] = c(mean(npo.trad.samps[1,]>npo.trad.samps[2,]),var(npo.trad.samps[1,]-npo.trad.samps[2,]),rowMeans(npo.trad.samps))

    po.comm.samps <- apply(po.comm.sampler(dat=dat),1,function(x) c(x[1]+60*x[2],x[3]+60*x[4],x[5]))
    po.comm.stats[,analysis] = c(mean(po.comm.samps[1,]>po.comm.samps[2,]),var(po.comm.samps[1,]-po.comm.samps[2,]),rowMeans(po.comm.samps))

    npo.comm.samps <- apply(npo.comm.sampler(dat=dat),1,function(x) c(x[1]+60*x[2],x[3]+60*x[4],x[5],x[6],x[5]*x[6]))
    npo.comm.stats[,analysis] = c(mean(npo.comm.samps[1,]>npo.comm.samps[2,]),var(npo.comm.samps[1,]-npo.comm.samps[2,]),rowMeans(npo.comm.samps))
  }

  #Compile and Return Results
  res = list(); res$rms = c(rmsE,rmsC); res$lrt.stats = lrt.stats
  res$po.trad.stats = po.trad.stats; res$npo.trad.stats = npo.trad.stats
  res$po.comm.stats = po.comm.stats; res$npo.comm.stats = npo.comm.stats
  return(res)
}

### General Simulation Set-up
niters = 2000 # Number of Simulation Iterations
trial.effects = rep(0,niters) #runif(niters)-0.5
treatment.effects = rep(0,niters) #rep(1,niters) #2*runif(niters)
#Simulate 2000 trials under key null and alternative scenarios with no historical bias
#Simulate 10000 trials under null and alternative scenarios with varying historical bias
#Should take 48 hours to simulate 2000 trials using 24 cores

# Run Simulation in Parallel
library(foreach); library(doParallel); library(jsonlite)
RNGkind("L'Ecuyer"); registerDoParallel(cores=2)
foo <- foreach(iter = 1:niters, .packages=c('dplyr','jsonlite','BayesLogit','survival','INLA')) %dopar% {

  # Simulate trial
  res = simulate.trial(trial.effect=trial.effects[iter],treatment.effect=treatment.effects[iter])
  res$iter = iter
  
  # Write Out Statistics
  write(toJSON(res),file="Results/Stats.json",append=TRUE)
}
stopImplicitCluster()


### Null Simulation With No Historical Control Bias
niters = 2000; trial.effects = rep(0,niters); treatment.effects = rep(0,niters)
library(foreach); library(doParallel); library(jsonlite); RNGkind("L'Ecuyer"); registerDoParallel(cores=2)
foo <- foreach(iter = 1:niters, .packages=c('dplyr','jsonlite','BayesLogit','survival','INLA')) %dopar% {
  
  # Simulate trial
  res = simulate.trial(trial.effect=trial.effects[iter],treatment.effect=treatment.effects[iter])
  res$iter = iter
  
  # Write Out Statistics
  write(toJSON(res),file="Results/Null-Stats.json",append=TRUE)
}
stopImplicitCluster()


### Alternative Simulation With No Historical Control Bias
niters = 2000; trial.effects = rep(0,niters); treatment.effects = rep(1,niters)
library(foreach); library(doParallel); library(jsonlite); RNGkind("L'Ecuyer"); registerDoParallel(cores=2)
foo <- foreach(iter = 1:niters, .packages=c('dplyr','jsonlite','BayesLogit','survival','INLA')) %dopar% {
  
  # Simulate trial
  res = simulate.trial(trial.effect=trial.effects[iter],treatment.effect=treatment.effects[iter])
  res$iter = iter
  
  # Write Out Statistics
  write(toJSON(res),file="Results/Alt-Stats.json",append=TRUE)
}
stopImplicitCluster()


### Null Simulation With Historical Control Bias
niters = 10000; trial.effects = runif(niters)-0.5; treatment.effects = rep(0,niters)
library(foreach); library(doParallel); library(jsonlite); RNGkind("L'Ecuyer"); registerDoParallel(cores=2)
foo <- foreach(iter = 1:niters, .packages=c('dplyr','jsonlite','BayesLogit','survival','INLA')) %dopar% {
  
  # Simulate trial
  res = simulate.trial(trial.effect=trial.effects[iter],treatment.effect=treatment.effects[iter])
  res$iter = iter
  
  # Write Out Statistics
  write(toJSON(res),file="Results/Null-Bias-Stats.json",append=TRUE)
}
stopImplicitCluster()


### Alternative Simulation With Historical Control Bias
niters = 10000; trial.effects = runif(niters)-0.5; treatment.effects = rep(1,niters)
library(foreach); library(doParallel); library(jsonlite); RNGkind("L'Ecuyer"); registerDoParallel(cores=2)
foo <- foreach(iter = 1:niters, .packages=c('dplyr','jsonlite','BayesLogit','survival','INLA')) %dopar% {
  
  # Simulate trial
  res = simulate.trial(trial.effect=trial.effects[iter],treatment.effect=treatment.effects[iter])
  res$iter = iter
  
  # Write Out Statistics
  write(toJSON(res),file="Results/Alt-Bias-Stats.json",append=TRUE)
}
stopImplicitCluster()




##################################
### Analyze Simulation Results ###
##################################
library(gsDesign); library(jsonlite); library(mgcv)
#setwd("Path to Your Working Directory")#Not Run
setwd("C:/Users/Thomas/Dropbox/Research/Active/SepsiCool/R-Software/Public-Files")

# Determine Monitoring Boundaries
des = gsDesign(k=3,test.type=6,alpha=0.025,beta=0.1,sfu=sfHSD,sfupar=-4,sfl=sfHSD,sflpar=-2)
a = des$upper$bound; b = des$lower$bound; p = pnorm(a); q = pnorm(b)

###############
### Table 1 ###
###############
### Null Scenario
stats <- lapply(readLines("Results/Null-Stats.json"),function(x) fromJSON(x))

# Calculate Proportion of Trials Recommending E at Each Analysis
lrt.stat <- sapply(stats,function(x) x$lrt[1,])
lrt.recE <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]>=a[1],z[1]<a[1] & z[2]>=a[2],z[1]<a[1] & z[2]<a[2] & z[3]>=a[3])))
lrt.ess <- sum(956*c(305,610,855)/855*c(lrt.recE[1:2],1-sum(lrt.recE[1:2])))

po.trad.stat <- sapply(stats,function(x) x$po.trad.stats[1,])
po.trad.recE <- rowMeans(apply(po.trad.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
po.trad.ess <- sum(956*c(305,610,855)/855*c(po.trad.recE[1:2],1-sum(po.trad.recE[1:2])))

npo.trad.stat <- sapply(stats,function(x) x$npo.trad.stats[1,])
npo.trad.recE <- rowMeans(apply(npo.trad.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
npo.trad.ess <- sum(956*c(305,610,855)/855*c(npo.trad.recE[1:2],1-sum(npo.trad.recE[1:2])))

po.comm.stat <- sapply(stats,function(x) x$po.comm.stats[1,])
po.comm.recE <- rowMeans(apply(po.comm.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
po.comm.ess <- sum(956*c(305,610,855)/855*c(po.comm.recE[1:2],1-sum(po.comm.recE[1:2])))

npo.comm.stat <- sapply(stats,function(x) x$npo.comm.stats[1,])
npo.comm.recE <- rowMeans(apply(npo.comm.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
npo.comm.ess <- sum(956*c(305,610,855)/855*c(npo.comm.recE[1:2],1-sum(npo.comm.recE[1:2])))

# Table 1 (Null Scenario)
meat <- rbind(lrt.recE, po.trad.recE, npo.trad.recE, po.comm.recE, npo.comm.recE)
tab <- cbind(round(meat,3),round(rowSums(meat),3),round(c(lrt.ess, po.trad.ess, npo.trad.ess, po.comm.ess, npo.comm.ess),1))
colnames(tab) = c("1","2","3","Tot","mSS"); rownames(tab) = c("LRT","Trad PO","Trad NPO","Comm PO","Comm NPO")
write.table(tab,sep=" & ",quote=FALSE)


### Alternative Scenario
stats <- lapply(readLines("Results/Alt-Stats.json"),function(x) fromJSON(x))

# Calculate Proportion of Trials Recommending E at Each Analysis
lrt.stat <- sapply(stats,function(x) x$lrt[1,])
lrt.recE <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]>=a[1],z[1]<a[1] & z[2]>=a[2],z[1]<a[1] & z[2]<a[2] & z[3]>=a[3])))
lrt.ess <- sum(956*c(305,610,855)/855*c(lrt.recE[1:2],1-sum(lrt.recE[1:2])))

po.trad.stat <- sapply(stats,function(x) x$po.trad.stats[1,])
po.trad.recE <- rowMeans(apply(po.trad.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
po.trad.ess <- sum(956*c(305,610,855)/855*c(po.trad.recE[1:2],1-sum(po.trad.recE[1:2])))

npo.trad.stat <- sapply(stats,function(x) x$npo.trad.stats[1,])
npo.trad.recE <- rowMeans(apply(npo.trad.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
npo.trad.ess <- sum(956*c(305,610,855)/855*c(npo.trad.recE[1:2],1-sum(npo.trad.recE[1:2])))

po.comm.stat <- sapply(stats,function(x) x$po.comm.stats[1,])
po.comm.recE <- rowMeans(apply(po.comm.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
po.comm.ess <- sum(956*c(305,610,855)/855*c(po.comm.recE[1:2],1-sum(po.comm.recE[1:2])))

npo.comm.stat <- sapply(stats,function(x) x$npo.comm.stats[1,])
npo.comm.recE <- rowMeans(apply(npo.comm.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
npo.comm.ess <- sum(956*c(305,610,855)/855*c(npo.comm.recE[1:2],1-sum(npo.comm.recE[1:2])))

# Table 1 (Alternative Scenario)
meat <- rbind(lrt.recE, po.trad.recE, npo.trad.recE, po.comm.recE, npo.comm.recE)
tab <- cbind(round(meat,3),round(rowSums(meat),3),round(c(lrt.ess, po.trad.ess, npo.trad.ess, po.comm.ess, npo.comm.ess),1))
colnames(tab) = c("1","2","3","Tot","mSS"); rownames(tab) = c("LRT","Trad PO","Trad NPO","Comm PO","Comm NPO")
write.table(tab,sep=" & ",quote=FALSE)



###########################
### Figure 2 Null Panel ###
###########################
### Null Scenario
stats <- c(lapply(readLines("Results/Null-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Null-Bias-Stats.json"),function(x) fromJSON(x)))

# True RMS
rmsE = sapply(stats,function(x) x$rms[1])
rmsC = sapply(stats,function(x) x$rms[2])
rmsC0 = 39.8068

# Prediction Points for \mu_C
pred.dat <- data.frame(rmsC=seq(min(rmsC),max(rmsC),length.out=100))

# Determine which Trials Recommend E and Estimate then Plot How Power Changes with \mu_C
lrt.recE <- apply(sapply(stats,function(x) x$lrt[1,]),2,function(z) z[1]>=a[1] | z[2]>=a[2] | z[3]>=a[3])
lrt.pred <- predict(gam(lrt.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

po.trad.recE <- apply(sapply(stats,function(x) x$po.trad.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
po.trad.pred <- predict(gam(po.trad.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

npo.trad.recE <- apply(sapply(stats,function(x) x$npo.trad.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
npo.trad.pred <- predict(gam(npo.trad.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

po.comm.recE <- apply(sapply(stats,function(x) x$po.comm.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
po.comm.pred <- predict(gam(po.comm.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

npo.comm.recE <- apply(sapply(stats,function(x) x$npo.comm.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
npo.comm.pred <- predict(gam(npo.comm.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

# Null Panel
pdf(file="Figure2-Null.pdf", width=3.5, height=3.5, family="serif")
par(mfrow=c(1,1),mar=c(5,4,2,1)+0.1)
plot(pred.dat$rmsC,lrt.pred$fit,type="l",lwd=3,las=1,ylim=c(0,0.2),ylab="Probability of Recommending E",xlab=expression(mu[C]),main="Null")
lines(pred.dat$rmsC,po.trad.pred$fit,lwd=3,lty=2,col=2)
lines(pred.dat$rmsC,npo.trad.pred$fit,lwd=3,lty=3,col=3)
lines(pred.dat$rmsC,po.comm.pred$fit,lwd=3,lty=4,col=4)
lines(pred.dat$rmsC,npo.comm.pred$fit,lwd=3,lty=5,col=5)
abline(h=c(0.025,0.90),lty=3); abline(v=rmsC0,lty=3)
legend("topleft",c("LRT","Trad PO","Trad NPO","Comm PO","Comm NPO"),lty=1:5,col=1:5,seg.len = 2,lwd=3,bty="n")
dev.off()



################
### Figure 4 ###
################
# Calculate RIG and estimate then plot how it changes with \mu_C
rig <- sapply(stats,function(x) x$npo.trad.stats[2,3]/x$npo.comm.stats[2,3]-1)
rig.pred <- predict(gam(rig~s(rmsC)),newdata=pred.dat,type="response",se.fit=TRUE)

# Get Prob(nu = 1 | Data) and estimate then plot how it changes with \mu_C
ppe <- sapply(stats,function(x) x$npo.comm.stats[5,3])
ppe.pred <- predict(gam(ppe~s(rmsC)),newdata=pred.dat,type="response",se.fit=TRUE)

rig.ppe.pred <- predict(gam(rig~s(ppe)),newdata=data.frame(ppe=seq(0,1,by=0.01)),type="response")

# Create Figure 4
pdf(file="Figure4.pdf", width=6.5, height=2.5, family="serif")
par(mfrow=c(1,3),mar=c(5,4,1,1)+0.1)

plot(rmsC,rig,pch=20,las=1,col='grey60',cex=0.2,ylab="RIG",xlab=expression(mu[C]))
lines(pred.dat$rmsC,rig.pred$fit,lwd=3)
abline(h=0,lty=3);abline(v=rmsC0,lty=3)

plot(rmsC,ppe,pch=20,las=1,col='grey60',cex=0.2,ylab=expression(paste("PPE for ",zeta,sep="")),xlab=expression(mu[C]))
lines(pred.dat$rmsC,ppe.pred$fit,lwd=3)
abline(v=rmsC0,lty=3)

plot(ppe,rig,pch=20,las=1,col='grey60',cex=0.2,ylab="RIG",xlab=expression(paste("PPE for ",zeta,sep="")))
lines(seq(0,1,by=0.01),rig.ppe.pred,lwd=3)
abline(h=0,lty=3);abline(v=rmsC0,lty=3)

dev.off()


##################################
### Figure 2 Alternative Panel ###
##################################
### Alternative Scenario
stats <- c(lapply(readLines("Results/Alt-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Alt-Bias-Stats.json"),function(x) fromJSON(x)))

# True RMS
rmsE = sapply(stats,function(x) x$rms[1])
rmsC = sapply(stats,function(x) x$rms[2])
rmsC0 = 39.8068

# Prediction Points for \mu_C
pred.dat <- data.frame(rmsC=seq(min(rmsC),max(rmsC),length.out=100))

# Determine which Trials Recommend E and Estimate then Plot How Power Changes with \mu_C
lrt.recE <- apply(sapply(stats,function(x) x$lrt[1,]),2,function(z) z[1]>=a[1] | z[2]>=a[2] | z[3]>=a[3])
lrt.pred <- predict(gam(lrt.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

po.trad.recE <- apply(sapply(stats,function(x) x$po.trad.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
po.trad.pred <- predict(gam(po.trad.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

npo.trad.recE <- apply(sapply(stats,function(x) x$npo.trad.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
npo.trad.pred <- predict(gam(npo.trad.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

po.comm.recE <- apply(sapply(stats,function(x) x$po.comm.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
po.comm.pred <- predict(gam(po.comm.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

npo.comm.recE <- apply(sapply(stats,function(x) x$npo.comm.stats[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
npo.comm.pred <- predict(gam(npo.comm.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

# Alternative Panel
pdf(file="Figure2-Alt.pdf", width=3.5, height=3.5, family="serif")
par(mfrow=c(1,1),mar=c(5,4,2,1)+0.1)
plot(pred.dat$rmsC,lrt.pred$fit,type="l",lwd=3,las=1,ylim=c(0,1),ylab="Probability of Recommending E",xlab=expression(mu[C]),main="Alternative")
lines(pred.dat$rmsC,po.trad.pred$fit,lwd=3,lty=2,col=2)
lines(pred.dat$rmsC,npo.trad.pred$fit,lwd=3,lty=3,col=3)
lines(pred.dat$rmsC,po.comm.pred$fit,lwd=3,lty=4,col=4)
lines(pred.dat$rmsC,npo.comm.pred$fit,lwd=3,lty=5,col=5)
abline(h=c(0.025,0.90),lty=3); abline(v=rmsC0,lty=3)
legend("bottomleft",c("LRT","Trad PO","Trad NPO","Comm PO","Comm NPO"),lty=1:5,col=1:5,seg.len = 2,lwd=3,bty="n")
dev.off()

################
### Figure 3 ###
################
# Get E[\mu_C | Data] and estimate then plot how (E[\mu_C | D] - \mu_C) changes with \mu_C
po.trad.rmsC <- sapply(stats,function(x) x$po.trad.stats[4,3])
po.trad.pred <- predict(gam((po.trad.rmsC-rmsC)~s(rmsC)),pred.dat,type="response",se.fit=TRUE)

npo.trad.rmsC <- sapply(stats,function(x) x$npo.trad.stats[4,3])
npo.trad.pred <- predict(gam((npo.trad.rmsC-rmsC)~s(rmsC)),pred.dat,type="response",se.fit=TRUE)

po.comm.rmsC <- sapply(stats,function(x) x$po.comm.stats[4,3])
po.comm.pred <- predict(gam((po.comm.rmsC-rmsC)~s(rmsC)),pred.dat,type="response",se.fit=TRUE)

npo.comm.rmsC <- sapply(stats,function(x) x$npo.comm.stats[4,3])
npo.comm.pred <- predict(gam((npo.comm.rmsC-rmsC)~s(rmsC)),pred.dat,type="response",se.fit=TRUE)

# Create Figure
pdf(file="Figure3.pdf", width=6.5, height=2.5, family="serif")
par(mfrow=c(1,3),mar=c(5,4.3,2,1)+0.1)

plot(pred.dat$rmsC,po.trad.pred$fit,lwd=3,las=1,ylab=expression(paste(E,"[",hat(mu)[C]-mu[C]," | ",mu[c],"]",sep="")),xlab=expression(mu[C]),main=" ",cex=0.2, ylim=c(-1.5,1.5),lty=2,col=2)
lines(pred.dat$rmsC,npo.trad.pred$fit,lwd=3,lty=3,col=3)
lines(pred.dat$rmsC,po.comm.pred$fit,lwd=3,lty=4,col=4)
lines(pred.dat$rmsC,npo.comm.pred$fit,lwd=3,lty=5,col=5)
abline(h=0,lty=3); abline(v=rmsC0,lty=3)
legend("bottomleft",c("Trad PO","Trad NPO","Comm PO","Comm NPO"),lwd=3,lty=c(1,3:5),col=2:5,seg.len = 2,bty="n")

plot(po.comm.rmsC,npo.trad.rmsC,pch=20,col="grey60",las=1,ylab=expression(paste(hat(mu)[C]," from Comm PO Model",sep="")),xlab=expression(paste(hat(mu)[C]," from Trad NPO Model",sep="")),main=" ",cex=0.2)
abline(h=rmsC0,lty=3); abline(v=rmsC0,lty=3); lines(c(0,55),c(0,55),lty=3)

plot(npo.comm.rmsC,npo.trad.rmsC,pch=20,col="grey60",las=1,ylab=expression(paste(hat(mu)[C]," from Comm NPO Model",sep="")),xlab=expression(paste(hat(mu)[C]," from Trad NPO Model",sep="")),main=" ",cex=0.2)
abline(h=rmsC0,lty=3); abline(v=rmsC0,lty=3); lines(c(0,55),c(0,55),lty=3)

dev.off()

