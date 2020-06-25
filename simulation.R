### Computer simulation to assess operating characteristics of the design
### Reproduces Table 1 and Figures 2-4
### Author: Thomas A Murray
### e-mail: murra484@umn.edu
setwd("your-working-directory")
source('functions.R')

### Load Historical Data
library(dplyr)
dat.hist = read.csv('Historical-Control-Data.csv',header=TRUE) %>% as_tibble() %>% mutate(arm=3)

###########################
###########################
### Computer Simulation ###
###########################
###########################
### To reproduce Table 1, this code must be run once under the null and once under the altertive
### To reproduce Figures 2-4, this code must be run again once under the null and once under the alternative
### See annotations in the code for details on how to modify the code to reproduce each element
nobs = 956 # Sample Size for the Current Trial
analysis.days = c(305,610,915) # Anaysis Days
niters = 2000 # Number of Simulation Iterations (2000 for Table 1, 10000 for Figures 2-4)

### Run Simulation (We used 'doParallel' and 'foreach' to run the simulations in parallel)
library(jsonlite)
for(iter in 1:niters){
  #Storage Object
  res = list() 
  res$iter = iter
  res$rms = rep(NA,2)
  res$lrt = matrix(NA,nrow=1,ncol=length(analysis.days))
  res$curr.only = matrix(NA,nrow=4,ncol=length(analysis.days))
  res$curr.hist = matrix(NA,nrow=5,ncol=length(analysis.days))
  
  # Generate trial and treatment effect, calculate corresponding hazards
  trial.effect = 0 #Use to reproduce Table 1
  #trial.effect = runif(1)-0.5 #Use to reproduce Figures 2-4
  hazC = exp(trial.effect)*hazC0/(1-hazC0)/(1+exp(trial.effect)*hazC0/(1-hazC0))
  hazE = hazC 
  #Uncomment the next three lines for the alternative; leave commented out for null
  #rrr = c(rep(0.4,15),0.4-0.2*(1:15)/15,rep(0.2,30)) #Targeted Relative Risk Reduction
  #hazE[1] = (1-rrr[1])*hazC[1]
  #for(t in 2:60) hazE[t] = 1-(1-(1-rrr[t])*(1-prod(1-hazC[1:t])))/prod(1-hazE[1:(t-1)])
  res$rms[1] = round(-sum(diff(c(1,cumprod(1-hazE),0))*(0:60)),4) #True RMS in E arm
  res$rms[2] = round(-sum(diff(c(1,cumprod(1-hazC),0))*(0:60)),4) #True RMS in C arm

  # Generate Data on Concurrent Arms
  dat = get.current.data(nobs=nobs,hazE=hazE,hazC=hazC)
  
  # Posterior Statistics at each analysis
  for(analysis in 1:length(analysis.days)){
    # Current Data at Analysis <analysis>
    dat.curr = filter(dat,enroll<=analysis.days[analysis]) %>% mutate(y = ifelse(y.star==61,60,pmin(y.star,analysis.days[analysis] - enroll + 1)),delta = 1*(y.star<=analysis.days[analysis]-enroll+1 & y.star<=60)) %>% select(y,delta,arm)
      
    # Calculate log-rank test statistic
    res$lrt[1,analysis] = get.lrt.stat(dat=dat.curr)

    # Fit Curr Only
    res$curr.only[1:4,analysis] <- get.curr.only.stats(dat=dat.curr)

    # Fit Curr-Hist
    res$curr.hist[1:5,analysis] <- get.curr.hist.stats(dat=rbind(dat.curr,dat.hist))
  }

  # Write Out Statistics
  write(toJSON(res),file="Results/Null-Stats.json",append=TRUE) #Use for Table 1 Null Scenario
  #write(toJSON(res),file="Results/Alt-Stats.json",append=TRUE) #Use for Table 1 Alternative Scenario
  #write(toJSON(res),file="Results/Null-HB-Stats.json",append=TRUE) #Use for Figures 2-4 Null Scenario
  #write(toJSON(res),file="Results/Alt-HB-Stats.json",append=TRUE) #Use for Figures 2-4 Alternative Scenario
}








#######################
#######################
### Analyze Results ###
#######################
#######################
library(jsonlite); library(gsDesign)

# Determine Monitoring Boundaries
des = gsDesign(k=3,test.type=6,alpha=0.025,beta=0.1,sfu=sfHSD,sfupar=-4,sfl=sfHSD,sflpar=-2)
a = des$upper$bound; b = des$lower$bound #z-score boundaries 
p = pnorm(a); q = pnorm(b) #posterior probability boundaries

###############
### Table 1 ###
###############
#Null Scenario
stats <- lapply(readLines("Results/Null-Stats.json"),function(x) fromJSON(x))

#LRT
lrt.stat <- sapply(stats,function(x) x$lrt[1,])
lrt.recE <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]>=a[1],z[1]<a[1] & z[2]>=a[2],z[1]<a[1] & z[2]<a[2] & z[3]>=a[3])))
lrt.ess <- sum(956*c(305,610,855)/855*c(lrt.recE[1:2],1-sum(lrt.recE[1:2])))
# With stopping for futility
#lrt.recE <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]>=a[1],b[1]<z[1] & z[1]<a[1] & z[2]>=a[2],b[1]<z[1] & z[1]<a[1] & b[2]<z[2] & z[2]<a[2] & z[3]>=a[3])))
#lrt.recC <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]<=b[1],b[1]<z[1] & z[1]<a[1] & z[2]<=b[2],b[1]<z[1] & z[1]<a[1] & b[2]<z[2] & z[2]<a[2] & z[3]<=b[3])))
#lrt.ess <- sum(956*c(305,610,855)/855*(lrt.recE+lrt.recC))

#Curr Only
bco.stat <- sapply(stats,function(x) x$curr.only[1,])
bco.recE <- rowMeans(apply(bco.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
bco.ess <- sum(956*c(305,610,855)/855*c(bco.recE[1:2],1-sum(bco.recE[1:2])))
# With stopping for futility
#bco.recE <- rowMeans(apply(bco.stat,2,function(z) c(z[1]>=p[1],q[1]<z[1] & z[1]<p[1] & z[2]>=p[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]>=p[3])))
#bco.recC <- rowMeans(apply(bco.stat,2,function(z) c(z[1]<=q[1],q[1]<z[1] & z[1]<p[1] & z[2]<=q[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]<=q[3])))
#bco.ess <- sum(956*c(305,610,855)/855*(bco.recE+bco.recC))

#Curr Hist
bch.stat <- sapply(stats,function(x) x$curr.hist[1,])
bch.recE <- rowMeans(apply(bch.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
bch.ess <- sum(956*c(305,610,855)/855*c(bch.recE[1:2],1-sum(bch.recE[1:2])))
# With stopping for futility
#bch.recE <- rowMeans(apply(bch.stat,2,function(z) c(z[1]>=p[1],q[1]<z[1] & z[1]<p[1] & z[2]>=p[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]>=p[3])))
#bch.recC <- rowMeans(apply(bch.stat,2,function(z) c(z[1]<=q[1],q[1]<z[1] & z[1]<p[1] & z[2]<=q[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]<=q[3])))
#bch.ess <- sum(956*c(305,610,855)/855*(bch.recE+bch.recC))

# Null Results
tab1.null <- cbind(round(rbind(lrt.recE, bco.recE, bch.recE),3),round(rowSums(rbind(lrt.recE, bco.recE, bch.recE)),3),round(c(lrt.ess,bco.ess,bch.ess),1))


#Alternative Scenario
stats <- lapply(readLines("Results/Alt-Stats.json"),function(x) fromJSON(x))

#LRT
lrt.stat <- sapply(stats,function(x) x$lrt[1,])
lrt.recE <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]>=a[1],z[1]<a[1] & z[2]>=a[2],z[1]<a[1] & z[2]<a[2] & z[3]>=a[3])))
lrt.ess <- sum(956*c(305,610,855)/855*c(lrt.recE[1:2],1-sum(lrt.recE[1:2])))
# With stopping for futility
#lrt.recE <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]>=a[1],b[1]<z[1] & z[1]<a[1] & z[2]>=a[2],b[1]<z[1] & z[1]<a[1] & b[2]<z[2] & z[2]<a[2] & z[3]>=a[3])))
#lrt.recC <- rowMeans(apply(lrt.stat,2,function(z) c(z[1]<=b[1],b[1]<z[1] & z[1]<a[1] & z[2]<=b[2],b[1]<z[1] & z[1]<a[1] & b[2]<z[2] & z[2]<a[2] & z[3]<=b[3])))
#lrt.ess <- sum(956*c(305,610,855)/855*(lrt.recE+lrt.recC))

#Curr Only
bco.stat <- sapply(stats,function(x) x$curr.only[1,])
bco.recE <- rowMeans(apply(bco.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
bco.ess <- sum(956*c(305,610,855)/855*c(bco.recE[1:2],1-sum(bco.recE[1:2])))
# With stopping for futility
#bco.recE <- rowMeans(apply(bco.stat,2,function(z) c(z[1]>=p[1],q[1]<z[1] & z[1]<p[1] & z[2]>=p[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]>=p[3])))
#bco.recC <- rowMeans(apply(bco.stat,2,function(z) c(z[1]<=q[1],q[1]<z[1] & z[1]<p[1] & z[2]<=q[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]<=q[3])))
#bco.ess <- sum(956*c(305,610,855)/855*(bco.recE+bco.recC))

#Curr Hist
bch.stat <- sapply(stats,function(x) x$curr.hist[1,])
bch.recE <- rowMeans(apply(bch.stat,2,function(z) c(z[1]>=p[1],z[1]<p[1] & z[2]>=p[2],z[1]<p[1] & z[2]<p[2] & z[3]>=p[3])))
bch.ess <- sum(956*c(305,610,855)/855*c(bch.recE[1:2],1-sum(bch.recE[1:2])))
# With stopping for futility
#bch.recE <- rowMeans(apply(bch.stat,2,function(z) c(z[1]>=p[1],q[1]<z[1] & z[1]<p[1] & z[2]>=p[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]>=p[3])))
#bch.recC <- rowMeans(apply(bch.stat,2,function(z) c(z[1]<=q[1],q[1]<z[1] & z[1]<p[1] & z[2]<=q[2],q[1]<z[1] & z[1]<p[1] & q[2]<z[2] & z[2]<p[2] & z[3]<=q[3])))
#bch.ess <- sum(956*c(305,610,855)/855*(bch.recE+bch.recC))

# Alternative Results
tab1.alt <- cbind(round(rbind(lrt.recE, bco.recE, bch.recE),3),round(rowSums(rbind(lrt.recE, bco.recE, bch.recE)),3),round(c(lrt.ess,bco.ess,bch.ess),1))

# Table 1
tab1 = rbind(tab1.null,tab1.alt)
colnames(tab1) = c("1","2","3","Tot","mSS"); rownames(tab1) = rep(c("LRT","Curr Only","Curr+Hist"),2)
tab1









###################
### Figures 2-4 ###
###################
library(mgcv)

### Figure 2 ###
# Null Stats (For left panel of Figure 2)
stats <- c(lapply(readLines("Results/Null-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Null-HB-Stats.json"),function(x) fromJSON(x)))

# True RMS
rmsE = sapply(stats,function(x) x$rms[1])
rmsC = sapply(stats,function(x) x$rms[2])
rmsC0 = round(-sum(diff(c(1,cumprod(1-hazC0),0))*(0:60)),4)

# Prediction Points for \mu_C
pred.dat <- data.frame(rmsC=seq(min(rmsC),max(rmsC),length.out=100))

### Determine which Trials Recommend E and Estimate then Plot How Power Changes with \mu_C
lrt.recE <- apply(sapply(stats,function(x) x$lrt[1,]),2,function(z) z[1]>=a[1] | z[2]>=a[2] | z[3]>=a[3])
lrt.pred <- predict(gam(lrt.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

bco.recE <- apply(sapply(stats,function(x) x$curr.only[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
bco.pred <- predict(gam(bco.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

bch.recE <- apply(sapply(stats,function(x) x$curr.hist[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
bch.pred <- predict(gam(bch.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

# Left Panel of Figure 2
plot(pred.dat$rmsC,lrt.pred$fit,type="l",lwd=2,las=1,ylim=c(0,0.2),ylab="Probability of Recommending E",xlab=expression(mu[C]),main="Null")
lines(pred.dat$rmsC,bco.pred$fit,type="l",lwd=2,lty=2)
lines(pred.dat$rmsC,bch.pred$fit,type="l",lwd=2,lty=4)
abline(h=c(0.025,0.90),lty=3); abline(v=rmsC0,lty=3)
legend("topleft",c("LRT","Curr Only","Curr-Hist"),lty=c(1,2,4),bty='n',lwd=2,seg.len=2)

# Alternative Stats (For right panel of Figure 2)
stats <- c(lapply(readLines("Results/Alt-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Alt-HB-Stats.json"),function(x) fromJSON(x)))

# True RMS
rmsE = sapply(stats,function(x) x$rms[1])
rmsC = sapply(stats,function(x) x$rms[2])
rmsC0 = round(-sum(diff(c(1,cumprod(1-hazC0),0))*(0:60)),4)

# Prediction Points for \mu_C
pred.dat <- data.frame(rmsC=seq(min(rmsC),max(rmsC),length.out=100))

### Determine which Trials Recommend E and Estimate then Plot How Power Changes with \mu_C
lrt.recE <- apply(sapply(stats,function(x) x$lrt[1,]),2,function(z) z[1]>=a[1] | z[2]>=a[2] | z[3]>=a[3])
lrt.pred <- predict(gam(lrt.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

bco.recE <- apply(sapply(stats,function(x) x$curr.only[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
bco.pred <- predict(gam(bco.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

bch.recE <- apply(sapply(stats,function(x) x$curr.hist[1,]),2,function(z) z[1]>=p[1] | z[2]>=p[2] | z[3]>=p[3])
bch.pred <- predict(gam(bch.recE~s(rmsC),family=binomial),newdata=pred.dat,type="response",se.fit=TRUE)

# Right Panel of Figure 2
plot(pred.dat$rmsC,lrt.pred$fit,type="l",lwd=2,las=1,ylim=c(0,1),ylab="Probability of Recommending E",xlab=expression(mu[C]),main="Alternative")
lines(pred.dat$rmsC,bco.pred$fit,type="l",lwd=2,lty=2)
lines(pred.dat$rmsC,bch.pred$fit,type="l",lwd=2,lty=4)
abline(h=c(0.025,0.90),lty=3); abline(v=rmsC0,lty=3)








### Figure 3 ###
stats <- c(lapply(readLines("Results/Null-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Null-HB-Stats.json"),function(x) fromJSON(x)))
#stats <- c(lapply(readLines("Results/Alt-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Alt-HB-Stats.json"),function(x) fromJSON(x)))

# True RMS
rmsE = sapply(stats,function(x) x$rms[1])
rmsC = sapply(stats,function(x) x$rms[2])
rmsC0 = round(-sum(diff(c(1,cumprod(1-hazC0),0))*(0:60)),4)

# Prediction Points for \mu_C
pred.dat <- data.frame(rmsC=seq(min(rmsC),max(rmsC),length.out=100))

# Get E[\mu_C | Data] and estimate then plot how (E[\mu_C | D] - \mu_C) changes with \mu_C
bco.pmc <- sapply(stats,function(x) x$curr.only[3,3])
bco.pred <- predict(gam((bco.pmc-rmsC)~s(rmsC)),pred.dat,type="response",se.fit=TRUE)
#bco.pred <- predict(loess((bco.pmc-rmsC)~rmsC),pred.dat,se=TRUE)

bch.pmc <- sapply(stats,function(x) x$curr.hist[3,3])
bch.pred <- predict(gam((bch.pmc-rmsC)~s(rmsC)),pred.dat,type="response",se.fit=TRUE)
#bch.pred <- predict(loess((bch.pmc-rmsC)~rmsC),pred.dat,se=TRUE)

#Figure 3
par(mfrow=c(1,3),mar=c(5,4.3,2,1)+0.1)

plot(rmsC,bco.pmc-rmsC,pch=20,las=1,col='grey60',ylab=expression(hat(mu)[C]-mu[C]),xlab=expression(mu[C]),main="Curr Only",cex=0.2,ylim=c(min(c(bco.pmc-rmsC,bch.pmc-rmsC)),max(c(bco.pmc-rmsC,bch.pmc-rmsC))))
lines(pred.dat$rmsC,bco.pred$fit,lwd=3)
abline(h=0,lty=3);abline(v=rmsC0,lty=3)

plot(rmsC,bch.pmc-rmsC,pch=20,las=1,col='grey60',ylab=expression(hat(mu)[C]-mu[C]),xlab=expression(mu[C]),main="Curr-Hist",cex=0.2,ylim=c(min(c(bco.pmc-rmsC,bch.pmc-rmsC)),max(c(bco.pmc-rmsC,bch.pmc-rmsC))))
lines(pred.dat$rmsC,bch.pred$fit,lwd=3)
abline(h=0,lty=3);abline(v=rmsC0,lty=3)

plot(bco.pmc,bch.pmc,pch=20,las=1,col='grey60',cex=0.2,ylab=expression(paste(hat(mu)[C]," from Curr-Hist",sep="")),xlab=expression(paste(hat(mu)[C]," from Curr Only",sep="")),main="Curr-Hist vs Curr Only")
abline(h=rmsC0,lty=3); abline(v=rmsC0,lty=3); lines(c(0,55),c(0,55),lty=3)






### Figure 4 ###
stats <- c(lapply(readLines("Results/Null-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Null-HB-Stats.json"),function(x) fromJSON(x)))
#stats <- c(lapply(readLines("Results/Alt-Stats.json"),function(x) fromJSON(x)),lapply(readLines("Results/Alt-HB-Stats.json"),function(x) fromJSON(x)))

# Estimated RMS[C]'s
bco.pmc <- sapply(stats,function(x) x$curr.only[3,3])
bch.pmc <- sapply(stats,function(x) x$curr.hist[3,3])
rmsC0 = round(-sum(diff(c(1,cumprod(1-hazC0),0))*(0:60)),4)

# Calculate RIG and estimate then plot how it changes with \mu_C
rig <- sapply(stats,function(x) x$curr.only[4,3]/x$curr.hist[4,3]-1)

# Get Prob(nu = 1 | Data) and estimate then plot how it changes with \mu_C
ppe <- sapply(stats,function(x) x$curr.hist[5,3])

# Figure 4
par(mfrow=c(1,3),mar=c(5,4,1,1)+0.1)
plot(bco.pmc,ppe,pch=20,las=1,col='grey60',cex=0.2,ylab="PPE",xlab=expression(paste(hat(mu)[C]," from Curr Only",sep="")),main="")
abline(v=rmsC0,lty=3)

plot(bco.pmc,rig,pch=20,las=1,col='grey60',cex=0.2,ylab="RIG",xlab=expression(paste(hat(mu)[C]," from Curr Only",sep="")),main="")
abline(h=0,lty=3); abline(v=rmsC0,lty=3)

plot(ppe,rig,pch=20,las=1,col='grey60',cex=0.2,ylab="RIG",xlab="PPE")
abline(h=0,lty=3);abline(v=rmsC0,lty=3)
