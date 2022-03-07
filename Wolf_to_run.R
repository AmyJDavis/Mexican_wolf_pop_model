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
wolf.res.sigT.2mortk=wolf.Malth.Pois.sigT.2mort.mcmc(mdat,xm,xr,xp,mortmis,
                                                     betar.tune,betam.tune,
                                                     betap.tune,n.mcmc)