###################################################################
setwd("~/collab - Kate/delta-SMI-repository")

rm(list=ls())
set.seed(0) 

yim <- function(phi,theta,x,s,n,k=2){
  return(rnorm(n,theta*x^k+phi,s))
}

zim <- function(phi,s,n){
  return(rnorm(n,phi,s))
}

xim <- function(n, centre, spread) {
  #return(rnorm(n,mean = centre, sd=spread))
  return(runif(n,min=centre-spread,max=centre+spread))
  #n=50; m=n/2; sy=0.5; sz=2#
  #phi_T=0; theta_T=1; X.centre=1; X.spread=1
}

#priors for theta and phi are ~1 improper

#observation model setup and data
n=30; m=n/3; sy=0.25; sz=1; X.centre=1; X.spread=1#
phi_T=0; theta_T=1; 
#X.centre=2/3; X.spread=2/3 #these are nice values for illustration - artificial because the cut limit is exactly correct

#the number of K-values and the list of K-values we consider
K.len=30; 
K=seq(from=0.1,to=2.3,length.out=K.len)

#the number of reps per k-val
N.trials=100

#for each k,trial pair we get the MSE and store it - we stor cut and bayes as well
MSE.phi=MSE.theta=MSE.phi.B=MSE.theta.B=MSE.phi.C=MSE.theta.C=delta.opt.loocv=delta.opt.exact=matrix(NA,K.len,N.trials,dimnames = list(K))

T.sim=1e3
v.gam <- 5^seq(from=-2,to=1,length.out=100); #my delta is sqrt(gamma) from Kate's 
n.gam=length(v.gam)

#trial=1; k.ind=1

for(trial in 1:N.trials){
  for (k.ind in 1:K.len) {
    #simulate some new synthetic data
    xobs <- xim(n=n,centre=X.centre, spread=X.spread)
    yobs <- yim(phi=phi_T,theta=theta_T,x=xobs,s=sy,n=n,k=K[k.ind]); 
    zobs <- zim(phi=phi_T,s=sz,m); 
    
    ELPDz=rep(0,n.gam)
    
    mx=mean(xobs); mxy=mean(xobs*yobs); mxx=mean(xobs^2)
    my=mean(yobs)
    
    rho=(sy^2+v.gam^2)*(m-1)/(sz^2*n)
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
    
    #now delta is chosen sample post
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
    
    #SMI sim at delta.star.loocv
    PHI=rnorm(T.sim,mean=mu.phi[which.delta.loocv],sd=s.phi[which.delta.loocv]) #we need sampled phi's to sample theta - could calculate theta marginal
    mu.thGphY=(mxy-PHI*mx)/mxx
    THETA=rnorm(T.sim,mean=mu.thGphY,sd=s.thGphY)
    
    #just trusting the limits when delta~0 and delta~inf
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

#pdf(file="PMSE-regression-example.pdf")
# boxplot(t(MSE.phi),boxwex=0.15,ylim=c(0,0.7),names=as.character(round(K,1)),
#         xlab="covariate X power k",ylab="posterior mean squared error (phi)")
# boxplot(t(MSE.phi.B),boxwex=0.15,add=TRUE,names=NULL,col=2,outcol=2,pch=2,axes=FALSE,at=1:K.len-0.2)
# boxplot(t(MSE.phi.C),boxwex=0.15,add=TRUE,names=NULL,col=3,outcol=3,pch=3,axes=FALSE,at=1:K.len+0.2)
boxplot(t(MSE.phi),boxwex=0.15,ylim=c(0,0.7),col=1,outcol=0,names=as.character(round(K,1)),
        xlab="covariate X power k",ylab="posterior mean squared error (phi)")
boxplot(t(MSE.phi.B),boxwex=0.15,add=TRUE,names=NULL,col=2,outcol=0,pch=2,axes=FALSE,at=1:K.len-0.2)
boxplot(t(MSE.phi.C),boxwex=0.15,add=TRUE,names=NULL,col=3,outcol=0,pch=3,axes=FALSE,at=1:K.len+0.2)
legend("bottomright", legend = c("Bayes","delta-SMI","Cut") , 
       col = c(2,1,3) , bty = "n", lwd=c(3,3,3))
#dev.off()

boxplot(t(MSE.theta),boxwex=0.15,ylim=c(0,1),names=as.character(round(K,1)),
        xlab="covariate power k",ylab="posterior mean squared error (theta)")
boxplot(t(MSE.theta.B),boxwex=0.15,add=TRUE,names=NULL,col=2,outcol=2,pch=2,axes=FALSE,at=1:K.len-0.2)
boxplot(t(MSE.theta.C),boxwex=0.15,add=TRUE,names=NULL,col=3,outcol=3,pch=3,axes=FALSE,at=1:K.len+0.2)
legend("topleft", legend = c("Bayes","SMI at delta*","Cut") , 
       col = c(2,1,3) , bty = "n", pch=c(2,1,3))

boxplot(t(delta.opt.loocv),boxwex=0.2,pch=1,names=as.character(round(K,1)),main="Delta-opt values - LOOCV and Exact",
        xlab="covariate power k",ylab="optimal delta")
boxplot(t(delta.opt.exact),boxwex=0.2,add=TRUE,names=NULL,col=2,outcol=2,pch=2,axes=FALSE,at=1:K.len+0.3)
legend("topleft", legend = c("LOOCV","Exact") , 
       col = c(1,2) , bty = "n", pch=c(1,2))
#dev.off()
