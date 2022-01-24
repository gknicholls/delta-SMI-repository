#Here is the code to produce figure 4 in section 5.2 of the paper
#"Valid belief updates for prequentially additive loss functions arising in Semi-Modular Inference"
#Geoff K. Nicholls, Jeong Eun Lee, Chieh-Hsi Wu and Chris U. Carmona
#(2022)
#
###################################################################
setwd("~/collab - Kate/delta-SMI-repository")

rm(list=ls())
set.seed(0) 

#true generative model for Y
yim <- function(phi,theta,x,s,n,k=2){
  return(rnorm(n,theta*x^k+phi,s))
}

#ture generative model for Z
zim <- function(phi,s,n){
  return(rnorm(n,phi,s))
}

#simulate values for the covariates in the regression
xim <- function(n, centre, spread) {
  #return(rnorm(n,mean = centre, sd=spread))
  return(runif(n,min=centre-spread,max=centre+spread))
  #n=50; m=n/2; sy=0.5; sz=2#
  #phi_T=0; theta_T=1; X.centre=1; X.spread=1
}

#priors for theta and phi are ~1 improper

#observation model setup and data - covariate values are U(0,2)
n=50; m=n; sy=0.25; sz=3; X.centre=1; X.spread=1#
phi_T=0; theta_T=1; 

#the number of K-values and the list of K-values we consider - this is
#the power of the covariate X in the true generative model - the fitted model
#takes k=1 so linear regression
K.len=11; 
K=seq(from=1,to=2,length.out=K.len)

#the number of reps per k-val - the number of samples in each bar in the boxplot
#each rep is an independent data set and a complete refit of Bayes/Cut/delta-SMI
N.trials=100

#for each k,trial pair we get the MSE and ELPD_z values and store them 
MSE.phi=MSE.theta=MSE.phi.B=MSE.theta.B=MSE.phi.C=MSE.theta.C=delta.opt.loocv=delta.opt.exact=matrix(NA,K.len,N.trials,dimnames = list(K))
ELPDz.excact.chosen=ELPDz.excact.Bayes=ELPDz.excact.Cut=matrix(NA,K.len,N.trials,dimnames = list(K))
  
T.sim=1e3 #the number of samples of each posterior used for example in estimating the posterior mean square errors used to form fig 4 top
v.gam <- 5^seq(from=-2,to=1,length.out=100); #these are the delta values considered - when we find delta^* we chose the one with the largest ELPD_z
n.gam=length(v.gam) #in this setting v.gam[1] should be a decent approx to Bayes and v.gam[n.gam] should give cut (to a decent approx)

#LOOCV or exact ELPD - choose delta-star using LOOCV (feasible in real applications) or compute exact ELPD. Paper uses
#LOOCV to estimate delta and then evaluates success by seeing how well we did on exact elpd (bottom graph fig 4)
USE.LOOCV=TRUE

for(trial in 1:N.trials){
  for (k.ind in 1:K.len) {
    #simulate some new synthetic data
    xobs <- xim(n=n,centre=X.centre, spread=X.spread)
    yobs <- yim(phi=phi_T,theta=theta_T,x=xobs,s=sy,n=n,k=K[k.ind]); 
    zobs <- zim(phi=phi_T,s=sz,m); 
    
    ELPDz=rep(0,n.gam)
    
    #we have the exact posteriors which are normal
    mx=mean(xobs); mxy=mean(xobs*yobs); mxx=mean(xobs^2)
    my=mean(yobs)
    
    rho=(sy^2+v.gam^2)*(m-1)/(sz^2*n)
    #LOOCV on the entries in Z
    for (ii in 1:m) {
      mz=mean(zobs[-ii])
      s.phi=sqrt(rho*sz^2/(m-1))/sqrt(rho+1-mx^2/mxx)
      mu.phi=(rho*mz+my*(1-mx*mxy/(mxx*my)))/(rho+1-mx^2/mxx)
      for (j in 1:n.gam) {ELPDz[j] = ELPDz[j]+dnorm(zobs[ii],mean=mu.phi[j],sd=sqrt(s.phi[j]^2+sz^2),log=T) }
    }
    ELPDz=ELPDz/m
    
    which.delta.loocv=which.max(ELPDz) #if we choose delta based on ELPDz
    delta.star.loocv=v.gam[which.delta.loocv]
    delta.opt.loocv[k.ind,trial]=delta.star.loocv
    
    #OK LOOCV is done so get mean of all data
    mz=mean(zobs)
    
    rho=(sy^2+v.gam^2)*m/(sz^2*n)
    s.thGphY=sy/sqrt(n*mxx)
    s.phi=sqrt(rho*sz^2/m)/sqrt(rho+1-mx^2/mxx)
    mu.phi=(rho*mz+my*(1-mx*mxy/(mxx*my)))/(rho+1-mx^2/mxx)
    s22=s.phi^2+sz^2
    
    ELPDz.exact = -0.5*log(2*pi)-0.5*log(s22)-(0.5/s22)*(sz^2+phi_T^2-2*phi_T*mu.phi+mu.phi^2)
    which.delta.exact=which.max(ELPDz.exact) #if we choose delta based on ELPDz
    delta.star.exact=v.gam[which.delta.exact]
    delta.opt.exact[k.ind,trial]=delta.star.exact
    ELPDz.excact.chosen[k.ind,trial]=ELPDz.exact[which.delta.exact]
    ELPDz.excact.Bayes[k.ind,trial]=ELPDz.exact[1]
    ELPDz.excact.Cut[k.ind,trial]=ELPDz.exact[n.gam]
    
    which.delta.use=ifelse(USE.LOOCV,which.delta.loocv,which.delta.exact)
    
    #SMI simulation using nested monte carlo at delta.star.loocv
    PHI=rnorm(T.sim,mean=mu.phi[which.delta.use],sd=s.phi[which.delta.use]) #we need sampled phi's to sample theta - could calculate theta marginal
    mu.thGphY=(mxy-PHI*mx)/mxx
    THETA=rnorm(T.sim,mean=mu.thGphY,sd=s.thGphY)
    
    #just trusting the limits v.gam[1] and v.gam[n.gam] when delta~0 and delta~inf
    PHI.BAYES=rnorm(T.sim,mean=mu.phi[1],sd=s.phi[1])
    mu.thGphY=(mxy-PHI.BAYES*mx)/mxx
    THETA.BAYES=rnorm(T.sim,mean=mu.thGphY,sd=s.thGphY)
    
    PHI.CUT=rnorm(T.sim,mean=mu.phi[n.gam],sd=s.phi[n.gam])
    mu.thGphY=(mxy-PHI.CUT*mx)/mxx
    THETA.CUT=rnorm(T.sim,mean=mu.thGphY,sd=s.thGphY)
    
    MSE.phi[k.ind,trial]=mean( (PHI-phi_T)^2 )
    MSE.theta[k.ind,trial]=mean( (THETA-theta_T)^2 )
    
    MSE.phi.B[k.ind,trial]=mean( (PHI.BAYES-phi_T)^2 )
    MSE.theta.B[k.ind,trial]=mean( (THETA.BAYES-theta_T)^2 )
    
    MSE.phi.C[k.ind,trial]=mean( (PHI.CUT-phi_T)^2 )
    MSE.theta.C[k.ind,trial]=mean( (THETA.CUT-theta_T)^2 )
    
  }
  print(trial)
}

#Fig 4 top
#pdf(file="PMSE-LOOCV-regression-example.pdf",width=6,height=5)
boxplot(t(MSE.phi),boxwex=0.15,border = par("bg"),ylim=c(0,0.6),col=1,outcol=0,names=as.character(round(K,1)),
        xlab="covariate power k",ylab="posterior mean squared error (phi)")
boxplot(t(MSE.phi.B),boxwex=0.15,border = par("bg"),add=TRUE,names=NULL,col=2,outcol=0,pch=2,axes=FALSE,at=1:K.len-0.2)
boxplot(t(MSE.phi.C),boxwex=0.15,border = par("bg"),add=TRUE,names=NULL,col=3,outcol=0,pch=3,axes=FALSE,at=1:K.len+0.2)
legend("topleft", legend = c("Bayes","delta-SMI","Cut") , 
       col = c(2,1,3) , bty = "n", lwd=c(3,3,3))
#dev.off()

#Fig 4 bottom
#pdf(file="ELPD-EXACT-regression-example.pdf",width=6,height=5)
boxplot(t(ELPDz.excact.chosen),boxwex=0.15,border = par("bg"),col=1,outcol=0,names=as.character(round(K,1)),
        xlab="covariate power k",ylab="ELPD_z")
boxplot(t(ELPDz.excact.Bayes),boxwex=0.15,border = par("bg"),add=TRUE,names=NULL,col=2,outcol=0,pch=2,axes=FALSE,at=1:K.len-0.2)
boxplot(t(ELPDz.excact.Cut),boxwex=0.15,border = par("bg"),add=TRUE,names=NULL,col=3,outcol=0,pch=3,axes=FALSE,at=1:K.len+0.2)
legend("bottomleft", legend = c("Bayes","delta-SMI","Cut") , 
       col = c(2,1,3) , bty = "n", lwd=c(3,3,3))
#dev.off()

#Like fig 4 top but for theta - less interesting than phi as the Y-module is misspecified so not in paper
boxplot(t(MSE.theta),boxwex=0.15,ylim=c(0,1),names=as.character(round(K,1)),
        xlab="covariate power k",ylab="posterior mean squared error (theta)")
boxplot(t(MSE.theta.B),boxwex=0.15,add=TRUE,names=NULL,col=2,outcol=2,pch=2,axes=FALSE,at=1:K.len-0.2)
boxplot(t(MSE.theta.C),boxwex=0.15,add=TRUE,names=NULL,col=3,outcol=3,pch=3,axes=FALSE,at=1:K.len+0.2)
legend("topleft", legend = c("Bayes","SMI at delta*","Cut") , 
       col = c(2,1,3) , bty = "n", pch=c(2,1,3))

#compare the LOOCV estimates of delta* with exact delta-star (well, pretty accurate) - pretty rough
boxplot(t(delta.opt.loocv),boxwex=0.2,pch=1,names=as.character(round(K,1)),main="Delta-opt values - LOOCV and Exact",
        xlab="covariate power k",ylab="optimal delta")
boxplot(t(delta.opt.exact),boxwex=0.2,add=TRUE,names=NULL,col=2,outcol=2,pch=2,axes=FALSE,at=1:K.len+0.3)
legend("topleft", legend = c("LOOCV","Exact") , 
       col = c(1,2) , bty = "n", pch=c(1,2))
#dev.off()

 
