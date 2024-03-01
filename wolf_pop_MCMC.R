##################################################################################
###
### Mexican Wolf hierarchical population model 
###    - malthusian growth on counts
###    - covariates on reproduction and mortalities
###  *** New formulation with Nt-1=nt-1+St+Rt  ****
###    - Poisson connection from mort and recruit vals
### October 21, 2018
### Amy J Davis
###
###################################################################################

###
wolf.Malth.Pois.sigT.2mort.mcmc<-function(mdat,xm,xr,xp,mortmis=c(0.8,0.2),
                                          betar.tune=0.1,betam.tune=0.1,
                                          betap.tune=0.1,n.mcmc){
  
  ### mdat = data frame with the following columns
  ###      1.  year = Year of data
  ###      2.  pop = Minimum population size by year (count of wolf population)
  ###      3.  Reproduction = Number of reproduction events by year
  ###      4.  PupRecruitment = Minimum pup recruitment number
  ###      5.  Removals = Number of individuals removed by management by year
  ###      6.  Releases = Number of management released individuals
  ###      7.  Translocations = Number of management translocations by year
  ###      8.  mort = Number of deaths from known or legal causes
  ###      9.  poaching = Number of deaths from illegal causes and cryptic
  
  ###
  ### xm = design matrix for covariates on mortality
  ### xr = design matrix for covariates on reproduction
  ### xp = design matrix for covariates on mortality proportions
  ###
  ### mortmis = the a priori likelihoods of misclassifying mortality causes
  ### betar.tune = tuning parameter for covariates on reproduction
  ### betam.tune = tuning parameter for covariates on mortality 
  ### betap.tune = tuning parameter for covariates on mortality proportions
  ### n.mcmc = the number of MCMC iteations to run 
  
  
  
  ### libraries
  library(MCMCpack)
  library(boot)
  t_col <- function(color, percent = 30, name = NULL) {
    #	  color = color name; percent = % transparency; name = an optional name for the color
    rgb.val <- col2rgb(color) ## Get RGB values for named color
    ## Make new color using input color as base and alpha set by transparency
    t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255,alpha = (100-percent)*255/100,names = name)
    invisible(t.col) ## Save the color
  }
  
  ### 
  ### Data set up
  ###
  yrs=dim(mdat)[1]
  
  xms=dim(xm)[2]
  xrs=dim(xr)[2]
  p.levs=dim(xp)[2]
  
  ###
  ### Starting values
  ###
  ### Hyper priors
  beta.p=rnorm(p.levs,0,1)
  p=inv.logit(as.matrix(xp)%*%beta.p)
  
  ###
  ### Parameter starts  
  ###
  xsig=cbind(1,(seq(1,yrs,1)-mean(seq(1,yrs,1)))/sqrt(var(seq(1,yrs,1))))
  ### Reproduction and betas for reproduction
  beta.r=rnorm(xrs,0,0.1)
  Rt=rpois(yrs,mdat$Reproduction)
  
  ### Mortality rate and betas for mortality rate, and converted to total mortality
  beta.m=rnorm(xms,0,0.1)
  Mt=rpois(yrs,mdat$mort)
  
  ###
  npred=rpois(yrs,mdat$pop)
  
  ###
  ### Create save matricies
  ###
  n.burn=round(n.mcmc/2)           ### Create a burn in portion      
  betapsave=matrix(0,p.levs,n.mcmc)
  betarsave=matrix(0,xrs,n.mcmc)
  betamsave=matrix(0,xms,n.mcmc)
  
  pmsave=matrix(0,yrs,n.mcmc)
  Mtsave=matrix(0,yrs,n.mcmc)
  Rtsave=matrix(0,yrs,n.mcmc)
  
  mse.n=rep(0,n.mcmc)
  mse.npred=rep(0,n.mcmc)
  msediffsave=rep(0,n.mcmc)
  msesave=rep(0,n.mcmc)
  npredsave=matrix(0,yrs,n.mcmc)
  
  poacht=matrix(0,yrs,n.mcmc)
  poachratet=matrix(0,yrs,n.mcmc)
  
  
  pd2save=matrix(0,yrs,n.mcmc)
  
  
  ### Run MCMC
  for(k in 1:n.mcmc){
    if(k%%1000==0)cat(k," ");flush.console()
    
    ###########################
    ### 
    ### Sample from npred
    ###
    ###########################
    nchange=1+(Rt-Mt-mdat$Removals)/(npred+mdat$Released+mdat$Translocations)
    
    ### t middles
    for(t in 1:(yrs-1)){
      nt.star=rpois(1,npred[t])
      nchange.star=1+(Rt-Mt-mdat$Removals)/(nt.star+mdat$Released+mdat$Translocations)
      mh1=dpois(mdat$pop[t+1],(nt.star+mdat$Translocations[t]+mdat$Released[t])*nchange.star[t],log=TRUE)+
        dpois(nt.star,(mdat$Count[t]+mdat$Translocations[t]+mdat$Released[t])*nchange.star[t],log=TRUE)+dpois(npred[t],nt.star,log=TRUE)
      mh2=dpois(mdat$pop[t+1],(npred[t]+mdat$Translocations[t]+mdat$Released[t])*nchange[t],log=TRUE)+
        dpois(npred[t],(mdat$Count[t]+mdat$Translocations[t]+mdat$Released[t])*nchange[t],log=TRUE)+dpois(nt.star,npred[t],log=TRUE)
      mh=exp(mh1-mh2)
      tmp.keep=mh>runif(1)
      npred[t][tmp.keep]=nt.star[tmp.keep]
    }
    npred
    
    ### t=T
    nTstar=rpois(1,npred[yrs])
    mh1=dpois(mdat$pop[yrs],(nTstar+mdat$Translocations[t]+mdat$Released[t])*nchange.star[t],log=TRUE)+dpois(npred[yrs],nTstar,log=TRUE)
    mh2=dpois(mdat$pop[yrs],(npred[yrs]+mdat$Translocations[t]+mdat$Released[t])*nchange[t],log=TRUE)+dpois(nTstar,npred[yrs],log=TRUE)
    mh=exp(mh1-mh2)
    tmp.keep=mh>runif(1)
    npred[yrs][tmp.keep]=nTstar[tmp.keep]
    npred
    
    ###########################
    ### 
    ### Sample from Rt  
    ###
    ########################### 
    Rstar=rpois(yrs,Rt)
    nchangestar=1+(Rstar-Mt-mdat$Removals)/(npred+mdat$Released+mdat$Translocations)
    nchangestar
    
    npredt.1=c(mdat$pop[1],npred[-1])
    
    rtmhration=exp((dpois(npredt.1,(npred+mdat$Translocations+mdat$Released)*nchangestar,log=TRUE)+
                      dpois(Rstar,exp(as.matrix(xr)%*%beta.r),log=TRUE)+
                      dpois(mdat$PupRecruitment,Rstar,log=TRUE)+dpois(Rt,Rstar,log=TRUE))-
                     (dpois(npredt.1,(npred+mdat$Translocations+mdat$Released)*nchange,log=TRUE)+
                        dpois(Rt,exp(as.matrix(xr)%*%beta.r),log=TRUE)+
                        dpois(mdat$PupRecruitment,Rt,log=TRUE)+dpois(Rstar,Rt,log=TRUE)))
    
    tmp.keep=rtmhration>runif(yrs)
    tmp.keep[is.na(tmp.keep)]=FALSE
    Rt[tmp.keep]=Rstar[tmp.keep]
    nchange=1+(Rt-Mt-mdat$Removals)/(npred+mdat$Released+mdat$Translocations)
    
    
    
    ###########################
    ### 
    ### Sample from betas for r
    ###
    ########################### 
    beta.rstar=rnorm(xrs,beta.r,betar.tune)
    
    brration=exp((sum(dpois(Rt,exp(as.matrix(xr)%*%beta.rstar),log=TRUE))+dnorm(beta.rstar,0,1,log=TRUE))-
                   (sum(dpois(Rt,exp(as.matrix(xr)%*%beta.r),log=TRUE))+dnorm(beta.r,0,1,log=TRUE)))
    
    tmp.keep=brration>runif(xrs)
    beta.r[tmp.keep]=beta.rstar[tmp.keep]
    beta.r
    
    
    ###########################
    ### 
    ### Sample from Mt  
    ###
    ########################### 
    
    Mstar=rpois(yrs,Mt)
    
    nchangestar=1+(Rt-Mstar-mdat$Removals)/(npred+mdat$Released+mdat$Translocations)
    nchangestar
    
    mtmhrtio=exp((dpois(npredt.1,(npred+mdat$Translocations+mdat$Released)*nchangestar,log=TRUE)+
                    dpois(Mstar,exp(as.matrix(xm)%*%beta.m),log=TRUE)+dpois(mdat$mort,Mstar,log=TRUE)+dpois(Mt,Mstar,log=TRUE))-
                   (dpois(npredt.1,(npred+mdat$Translocations+mdat$Released)*nchange,log=TRUE)+
                      dpois(Mt,exp(as.matrix(xm)%*%beta.m),log=TRUE)+dpois(mdat$mort,Mt,log=TRUE)+dpois(Mstar,Mt,log=TRUE)))
    mtmhrtio[is.na(mtmhrtio)]=0
    tmp.keep=mtmhrtio>runif(yrs)
    Mt[tmp.keep]=Mstar[tmp.keep]
    Mt
    
    
    ###########################
    ### 
    ### Sample from betas for m
    ###
    ########################### 
    beta.mstar=rnorm(xms,beta.m,betam.tune)
    
    bmration=exp((sum(dpois(Mt,exp(as.matrix(xm)%*%beta.mstar),log=TRUE))+dnorm(beta.mstar,0,1,log=TRUE))-
                   (sum(dpois(Mt,exp(as.matrix(xm)%*%beta.m),log=TRUE))+dnorm(beta.m,0,1,log=TRUE)))
    
    tmp.keep=bmration>runif(xms)
    beta.m[tmp.keep]=beta.mstar[tmp.keep]
    beta.m
    
    
    ###########################
    ### 
    ### Sample from beta for the proportion of mortalities that are poaching
    ###
    ########################### 
    beta.pstar=rnorm(p.levs,beta.p,betap.tune)
    p.star=inv.logit(as.matrix(xp)%*%beta.pstar)
    
    pmhratio=exp((sum(dbinom(mdat$poaching,mdat$mort,p.star,log=TRUE))+dnorm(beta.pstar,0,1))-
                   (sum(dbinom(mdat$poaching,mdat$mort,p,log=TRUE))+dnorm(beta.p,0,1)))
    
    tmp.keep=pmhratio>runif(p.levs)
    beta.p[tmp.keep]=beta.pstar[tmp.keep]
    beta.p
    p=inv.logit(as.matrix(xp)%*%beta.p)
    
    
    ################################################
    ###
    ### Save samples
    ###
    ################################################
    betapsave[,k]=beta.p
    betarsave[,k]=beta.r
    betamsave[,k]=beta.m
    
    pmsave[,k]=p
    Mtsave[,k]=Mt
    Rtsave[,k]=Rt
    pd2save[,k]=dpois(npred,mdat$pop)
    npredsave[,k]=npred
    
    ###
    ### Cacluate MSE
    ###
    mse.n[k]=mean((mdat$Count[-1]-(mdat$Count[-yrs]+nchange[-yrs]))^2)
    mse.npred[k]=mean((npred[-1]-(mdat$Count[-yrs]+nchange[-yrs]))^2)
    msediffsave[k]=mse.npred[k]-mse.n[k]
    msesave[k]=mean((npred-mdat$Count)^2)
    
    
    
    ###
    ### Calculate poaching estimate and poaching rate
    ###
    poachest=rbinom(yrs,Mt,mean(p,mortmis[1]))
    poacht[,k]=poachest
    
    poachratet[,k]=poachest/npred
    
  }
  cat("\n")
  
  ###
  ### Make plots
  ###
  par(family="serif",mar=c(3,4,1,1)+0.1)
  
  ###
  ### Plots for Beta R
  ###
  brmin=min(betarsave[,n.burn:n.mcmc])
  brmax=max(betarsave[,n.burn:n.mcmc])
  layout(matrix(c(seq(1,2*min(10,dim(betarsave)[1]),by=2),seq(2,2*min(10,dim(betarsave)[1]),by=2),seq(2,2*min(10,dim(betarsave)[1]),by=2)),min(10,dim(betarsave)[1]),3))
  for(i in 1:dim(betarsave)[1]){
    plot(betarsave[i,],type="l",main=paste("Trace: beta r ",names(xr)[i]),ylab=paste("N ",i),xlab="MCMC Iteration")
    plot(density(betarsave[i,n.burn:n.mcmc]),lwd=2,main=paste("Posterior and Prior: beta r ",names(xr)[i]),xlim=c(brmin,brmax))
    curve(dnorm(x,0,1),col=2,lwd=3,lty=2,add=TRUE)
    legend("topright",col=c(1,2),lwd=c(2,3),lty=c(1,2),legend=c("Posterior","Prior"),bty='n')
  }
  
  ###
  ### Plots for Beta M
  ###
  bmmin=min(betamsave[,n.burn:n.mcmc])
  bmmax=max(betamsave[,n.burn:n.mcmc])
  layout(matrix(c(seq(1,2*min(10,dim(betamsave)[1]),by=2),seq(2,2*min(10,dim(betamsave)[1]),by=2),seq(2,2*min(10,dim(betamsave)[1]),by=2)),min(10,dim(betamsave)[1]),3))
  for(i in 1:dim(betamsave)[1]){
    plot(betamsave[i,],type="l",main=paste("Trace: beta m ",names(xm)[i]),ylab=paste("N ",i),xlab="MCMC Iteration")
    plot(density(betamsave[i,n.burn:n.mcmc]),lwd=2,main=paste("Posterior and Prior: beta m ",names(xm)[i]),xlim=c(bmmin,bmmax))
    curve(dnorm(x,0,1),col=2,lwd=3,lty=2,add=TRUE)
    legend("topright",col=c(1,2),lwd=c(2,3),lty=c(1,2),legend=c("Posterior","Prior"),bty='n')
  }
  
  
  ###
  ### Plots for Beta p
  ###
  bpmin=min(betapsave[,n.burn:n.mcmc])
  bpmax=max(betapsave[,n.burn:n.mcmc])
  layout(matrix(c(seq(1,2*min(10,dim(betapsave)[1]),by=2),seq(2,2*min(10,dim(betapsave)[1]),by=2),seq(2,2*min(10,dim(betapsave)[1]),by=2)),min(10,dim(betapsave)[1]),3))
  for(i in 1:dim(betapsave)[1]){
    plot(betapsave[i,],type="l",main=paste("Trace: beta p ",names(xp)[i]),ylab=paste("N ",i),xlab="MCMC Iteration")
    plot(density(betapsave[i,n.burn:n.mcmc]),lwd=2,main=paste("Posterior and Prior: beta p ",names(xp)[i]),xlim=c(bpmin,bpmax))
    curve(dnorm(x,0,1),col=2,lwd=3,lty=2,add=TRUE)
    legend("topright",col=c(1,2),lwd=c(2,3),lty=c(1,2),legend=c("Posterior","Prior"),bty='n')
  }
  
  
  layout(matrix(c(seq(1,2*min(10,6),by=2),seq(2,2*min(10,6),by=2),seq(2,2*min(10,6),by=2)),min(10,6),3))
  npredmin=min(npredsave)
  npredmax=max(npredsave)
  for(i in 1:yrs){
    plot(npredsave[i,],type="l",main=paste("Trace: npred ",i),ylab="pm ",xlab="MCMC Iteration")
    plot(density(npredsave[i,n.burn:n.mcmc]),lwd=2,main="Posterior and Prior: npred ",xlim=c(npredmin,npredmax))
    legend("topright",col=c(1,2),lwd=c(2,3),lty=c(1,2),legend=c("Posterior","Prior"),bty='n')
  }

  
  list(npredsave=npredsave,Mtsave=Mtsave,Rtsave=Rtsave,betarsave=betarsave,pmsave=pmsave,poacht=poacht,
       poachratet=poachratet,betamsave=betamsave,betapsave=betapsave)
}