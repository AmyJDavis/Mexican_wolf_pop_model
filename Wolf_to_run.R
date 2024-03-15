##################################################################################
###
### Mexican Wolf population model
### 
### September 5, 2019
### Amy J Davis
###
###################################################################################

### Libraries
library(viridis)
source('wolf_pop_MCMC.R')

### Data
mdat=read.csv("MexWolfDataFile.csv")

## Process data for analysis
mdat$Count=c(0,mdat$pop[-dim(mdat)[1]])
mdat$justpoach=mdat$poaching
mdat$poaching=mdat$poaching+mdat$unknown

mdat$mort=mdat$poaching+mdat$vehicle+mdat$natural+mdat$otherlegal
mdat$nadj=mdat$Count+mdat$Released+mdat$Translocations
mdat$MortNotPoach=mdat$vehicle+mdat$natural+mdat$otherlegal

### Set up tuning parameters and MCMC iteration size
betar.tune=0.1
betam.tune=0.1
betap.tune=0.1
n.mcmc=1000000
n.burn=n.mcmc/2
mortmis=c(0.8,0.2)


### Model we are using (commented out all of the model selection stuff)
xm=data.frame(Int=1)
xr=data.frame(Int=rep(1,dim(mdat)[1]))
xp=data.frame(Int=1,Transrate=(mdat$Translocations+mdat$Released)/mdat$pop,Legal=mdat$Removals/mdat$pop)


### Code to run MCMC model
wolf.results=wolf.Malth.Pois.sigT.2mort.mcmc(mdat,xm,xr,xp,mortmis,
                                             betar.tune,betam.tune,
                                             betap.tune,n.mcmc)

### Calculating posterior means and confidence intervals
npredmean=apply(wolf.results$npredsave[,n.burn:n.mcmc],1,function(x).Internal(mean(x)))
npredci=apply(wolf.results$npredsave[,n.burn:n.mcmc],1,function(x)quantile(x,c(0.025,0.975)))

Rmean=apply(wolf.results$Rtsave[,n.burn:n.mcmc],1,function(x).Internal(mean(x)))
Rci=apply(wolf.results$Rtsave[,n.burn:n.mcmc],1,function(x)quantile(x,c(0.025,0.975)))

Mmean=apply(wolf.results$Mtsave[,n.burn:n.mcmc],1,function(x).Internal(mean(x)))
Mci=apply(wolf.results$Mtsave[,n.burn:n.mcmc],1,function(x)quantile(x,c(0.025,0.975)))

pmmean=apply(wolf.results$pmsave[,n.burn:n.mcmc],1,function(x).Internal(mean(x)))
pmci=apply(wolf.results$pmsave[,n.burn:n.mcmc],1,function(x)quantile(x,c(0.025,0.975)))

if(dim(wolf.results$betamsave)[1]==1){
  betammean=mean(wolf.results$betamsave[n.burn:n.mcmc])
  betamci=quantile(wolf.results$betamsave[n.burn:n.mcmc],c(0.025,0.975))
}else{
  betammean=apply(wolf.results$betamsave[,n.burn:n.mcmc],1,function(x).Internal(mean(x)))
  betamci=apply(wolf.results$betamsave[,n.burn:n.mcmc],1,function(x)quantile(x,c(0.025,0.975)))
}

poachmean=apply(wolf.results$poacht[,n.burn:n.mcmc],1,function(x).Internal(mean(x)))
poachCI=apply(wolf.results$poacht[,n.burn:n.mcmc],1,function(x)quantile(x,c(0.025,0.975)))

poachratemean=apply(wolf.results$poachratet[,n.burn:n.mcmc],1,function(x).Internal(mean(x)))
poachrateCI=apply(wolf.results$poachratet[,n.burn:n.mcmc],1,function(x)quantile(x,c(0.025,0.975)))

combdat=cbind(mdat,npredmean,t(npredci),Rmean,t(Rci),Mmean,t(Mci),poachmean,t(poachCI),poachratemean,t(poachrateCI))
combdat


### Calculating rates of interest
combdat$Poachrate=combdat$poachmean/combdat$nadj
combdat$Removrate=combdat$Removals/combdat$nadj
combdat$TRrate=(combdat$Translocations+combdat$Released)/combdat$nadj
combdat$Releaserate=combdat$Released/combdat$nadj
combdat$logPoachrate=log(combdat$Poachrate+0.1)
combdat$logRemrate=log(combdat$Removrate+0.1)
combdat$Policy=rep(c(0,1),each=11)

### Regression looking at relationship between poaching rate and removal rate 
### and translocations and removals rate. Relationships shown in Figure 5.
prtr1=lm(data=combdat,logPoachrate~Removrate+TRrate)

### Repeat of the above lm without the first year. The population is small this year
###  and the rates are thus larger than in other years. We wanted to see if this
###  had an undue influence on the results
prtr1=lm(data=combdat[-1,],logPoachrate~Removrate+TRrate)


### Model results for Table 2.
remlm_null=lm(data=combdat,logRemrate~1)
remlm_policy=lm(data=combdat,logRemrate~Policy)
remlm_Release=lm(data=combdat,logRemrate~Releaserate)
remlm_PolicyRelease=lm(data=combdat,logRemrate~Policy+Releaserate)
remlm_PolicyByRelease=lm(data=combdat,logRemrate~Policy*Releaserate)