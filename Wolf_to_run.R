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
mdat=read.csv("MexWolfData2020.csv")


### Set up tuning parameters and MCMC iteration size
betar.tune=0.1
betam.tune=0.1
betap.tune=0.1
n.mcmc=100000
n.burn=n.mcmc/2
mortmis=c(0.8,0.2)


### Set up design matrices for mortality, reproduction, and proportion mortality
xm=data.frame(Int=1)
xr=data.frame(Int=rep(1,dim(mdat)[1]))
xp=data.frame(Int=1,Depredations=mdat$Depperwolf,Transrate=(mdat$Translocations+mdat$Released)/mdat$pop,Legal=mdat$Removals/mdat$pop)


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