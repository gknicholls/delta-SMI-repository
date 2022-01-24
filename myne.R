#Here is the code to produce figures 2 and 3 in section 5.1 of the paper
#"Valid belief updates for prequentially additive loss functions arising in Semi-Modular Inference"
#Geoff K. Nicholls, Jeong Eun Lee, Chieh-Hsi Wu and Chris U. Carmona
#(2022)
#
###################################################################
setwd("~/collab - Kate/SMI-ABC")
#pdf(file="simple-normal-example.pdf",paper="a4")

rm(list=ls())
set.seed(6) # tried 1-6 wanted good separation of the cut and Bayes - slightly dodgy 
#but we were explicit about this in the paper


#true observation model for Y (is same as fitted)
yim <- function(phi,theta,s,n){
  return(rnorm(n,theta+phi,s))
}

#ditto Z
zim <- function(phi,s,n){
  return(rnorm(n,phi,s))
}

#sd for normal prior for theta
s.theta=0.33 #quite extreme

#observation model setup and data
n=50; m=n/2; sy=1; sz=3 
phi_T=0; theta_T=1
yobs <- yim(phi=phi_T,theta=theta_T,s=sy,n); 
zobs <- zim(phi=phi_T,s=sz,m); 

 
T.sim=1e5 #the number of samples of each posterior in scatter plot fig 2 and in estimating the posterior mean square error fig 3 top right 
v.gam <- 0.1*1000^seq(from=0,to=1,length.out=40); #these are the delta values considered - when we find delta^* we chose the one with the largest ELPD_z
n.gam=length(v.gam) #in this setting v.gam[1] should be a decent approx to Bayes and v.gam[n.gam] should give cut (to a decent approx)
PHI=THETA=matrix(NA,T.sim,n.gam) #save posterior samples

#calculate the posteriors analyticaly as simple normal - seee qns bottom page 16 top 17 of arXiv version
my=mean(yobs)
mz=mean(zobs)
rho=(s.theta^2/(s.theta^2+sy^2/n))
s.thGphY=1/sqrt(n/sy^2+1/s.theta^2)
s.delta=mu.delta=matrix(NA,n.gam,1)

#smi-abc posterior simulation - go through each of the delta values (confusingly called v.gam)
for (j in 1:n.gam) {
  delta=v.gam[j]
  lam=(m/sz^2)/((m/sz^2)+ (n/(sy^2+delta^2+n*s.theta^2)))
  s.delta[j]=sz*sqrt(lam/m)
  mu.delta[j]=lam*mz+(1-lam)*my
  PHI[,j]=rnorm(T.sim,mean=mu.delta[j],sd=s.delta[j]) #we need sampled phi's to sample theta - could calculate theta marginal
  mu.thGphY=rho*(my-PHI[,j])
  THETA[,j]=rnorm(T.sim,mean=mu.thGphY,sd=s.thGphY) #nested monte carlo
}
 
#taking the limit delta->0 and delta->infty respectively
#full bayes  
lam.bayes=(m/sz^2)/((m/sz^2)+ (n/(sy^2+n*s.theta^2)))
s.bayes.phi=sz*sqrt(lam.bayes/m)
mu.bayes.phi=lam.bayes*mz+(1-lam.bayes)*my
PHI.BAYES=rnorm(T.sim,mean=mu.bayes.phi,sd=s.bayes.phi) #we need sampled phi's to sample theta - could calculate theta marginal
mu.thGphY=rho*(my-PHI.BAYES)
THETA.BAYES=rnorm(T.sim,mean=mu.thGphY,sd=s.thGphY) #this stuff doesnt change

#full cut  
s.cut.phi=sz/sqrt(m)
mu.cut.phi=mz
PHI.CUT=rnorm(T.sim,mean=mu.cut.phi,sd=s.cut.phi) #we need sampled phi's to sample theta - could calculate theta marginal
mu.thGphY=rho*(my-PHI.CUT)
THETA.CUT=rnorm(T.sim,mean=mu.thGphY,sd=s.thGphY) #this stuff doesnt change


#means of phi and theta as delta varies
#plot(log(v.gam),apply(PHI,2,mean)); #we dont need this as we have the formula
par(mfrow=c(2,1));
par(mar=c(1,4,4,2))
plot(log(v.gam),mu.delta,type='l',xlab='',ylab="SMI-Post mean phi",main="SMI-Posterior Means with log(delta)"); 
points(log(v.gam[1]),mu.bayes.phi,pch=2,col=2);points(log(v.gam[n.gam]),mu.cut.phi,pch=6,col=3)
abline(h=phi_T,lty=2)
legend('bottomleft',y.intersp=0.9,cex=0.9,bty="n",col=c(1,1,2,3),lty=c(1,2,NA,NA),pch=c(NA,NA,2,6),legend=c("Post mean phi","True phi","Bayes mean","Cut mean"))
par(mar=c(5,4,2,2))
plot(log(v.gam),apply(THETA,2,mean),type='l',xlab="log(delta)",ylab="SMI-Post mean theta"); 
points(log(v.gam[1]),mean(THETA.BAYES),pch=2,col=2); points(log(v.gam[n.gam]),mean(THETA.CUT),pch=6,col=3)
abline(h=theta_T,lty=2)
legend('bottomright',y.intersp=0.9,cex=0.9,bty="n",col=c(1,1,2,3),lty=c(1,2,NA,NA),pch=c(NA,NA,2,6),legend=c("Post mean phi","True theta","Bayes mean","Cut mean"))

#MSE's 
par(mfrow=c(2,1));
par(mar=c(3,4,2,2))
plot(log(v.gam),(mu.delta-phi_T)^2+s.delta^2,type='l',xlab='',ylab="MSE phi",main="MSE phi (top) and theta (bottom)")
par(mar=c(5,4,0,2))
plot(log(v.gam),(apply(THETA,2,mean)-theta_T)^2+apply(THETA,2,sd)^2,type='l',xlab="log(delta)",ylab="MSE theta")


####
#important bits now - calculate the exact ELPD

#exact ELPD - just y - formula item (10) in notes
par(mfrow=c(3,1));
par(mar=c(3,4,2,2))
s1.sq=(1-rho)^2*s.delta^2+s.thGphY^2+sy^2
mu1=mu.delta*(1-rho)+rho*my
ELPDy=-(1/2)*log(2*pi)-(1/2)*log(s1.sq)-sy^2/(2*s1.sq)-(phi_T+theta_T-mu1)^2/(2*s1.sq)
plot(log(v.gam),ELPDy,type='l',xlab="",main="ELPD's just y (top), just z (mid) and (y,z) (bottom)")

#exact ELPD - just z - forula item (11) in notes
s2.sq=s.delta^2+sz^2
mu2=mu.delta
ELPDz=-(1/2)*log(2*pi)-(1/2)*log(s2.sq)-sz^2/(2*s2.sq)-(phi_T-mu2)^2/(2*s2.sq)
plot(log(v.gam),ELPDz,type='l',xlab="")

#exact ELPD - full elpd y and z - formula end of item (12) in notes
a=s1.sq; b=(1-rho)*s.delta^2; d=s2.sq; e=mu1; f=mu2
Delta=a*d-b^2
expr1=Delta^(-1)*(d*sy^2+d*(phi_T+theta_T-e)^2-2*b*(phi_T+theta_T-e)*(phi_T-f)+a*sz^2+a*(phi_T-f)^2)
ELPDyz=-(1/2)*expr1-(1/2)*log(Delta)-log(2*pi)
par(mar=c(5,4,0,2))
plot(log(v.gam),ELPDyz,type='l',xlab="log(delta)")
points(log(v.gam),ELPDz+ELPDy,col=2) #after all that calculation the ELPD joint for y and z is basicly equal ELPDy+ELPDz
which.delta=which.max(ELPDyz)
delta.star=v.gam[which.delta]
eta.star=sy^2/(sy^2+delta.star^2) #for comparison with eta-SMI
abline(v=log(delta.star),lty=2)
legend('topright',bty="n",col=c(1,2,1),lty=c(1,NA,2),pch=c(NA,1,NA),legend=c("ELPD(y,z)","ELPD(y)+ELPD(z)","log(delta*)"))

#density at cut and post
par(mfrow=c(2,1));
par(mar=c(4,4,2,2))
plot(density(PHI.BAYES),xlab="phi",ylab="Post density phi",main="SMI-Post density phi (top) and theta (bottom)",xlim=c(-2.5,2)); 
lines(density(PHI.CUT),col=3); abline(v=phi_T,col=2)
lines(density(PHI[,which.delta]),col=4)
legend('topleft',bty="n",col=c(1,4,3,2),lty=c(1,1,1,1),legend=c("Bayes","SMI-optimal","Cut","True phi"))
par(mar=c(5,4,1,2))
plot(density(THETA.BAYES),xlab="theta",ylab="Post density theta",main="",xlim=c(-1,3)); 
lines(density(THETA.CUT),col=3); abline(v=theta_T,col=2)
lines(density(THETA[,which.delta]),col=4)
legend('topleft',bty="n",col=c(1,4,3,2),lty=c(1,1,1),legend=c("Bayes","SMI-optimal","Cut","True theta"))

#joint dbns
par(mfrow=c(1,1))
#plot(PHI[1:(T.sim/10),which.delta],THETA[1:(T.sim/10),which.delta],col=2,xlab="phi",ylab="theta",pch=1)
#points(PHI.BAYES[1:(T.sim/10)],THETA.BAYES[1:(T.sim/10)],col=1,pch=1)
plot(PHI.BAYES[1:(T.sim/10)],THETA.BAYES[1:(T.sim/10)],col=1,pch='.',xlab="phi",ylab="theta",xlim=c(-2,2),ylim=c(-0.5,2.5))
points(PHI.CUT[1:(T.sim/10)],THETA.CUT[1:(T.sim/10)],col=3,pch='.')
abline(v=phi_T,h=theta_T,lty=2)

pi.smi<-function(delta,theta,phi){
  #SMI-ABC joint posterior
  rho=(s.theta^2/(s.theta^2+sy^2/n))
  s.thGphY=1/sqrt(n/sy^2+1/s.theta^2)
  lam=(m/sz^2)/((m/sz^2)+ (n/(sy^2+delta^2+n*s.theta^2)))
  s.delta=sz*sqrt(lam/m)
  mu.delta=lam*mz+(1-lam)*my
    PHI.d=dnorm(phi,mean=mu.delta,sd=s.delta) #we need sampled phi's to sample theta - could calculate theta marginal
    mu.thGphY=rho*(my-phi)
    THETA.d=dnorm(theta,mean=mu.thGphY,sd=s.thGphY)
    return(THETA.d*PHI.d)
}

M=300; u=seq(-2,2,length.out=M); v=seq(-0.5,2.5,length.out=M); z=matrix(NA,M,M); 
for (i in 1:M) { 
  for (j in 1:M) { 
    z[i,j]=pi.smi(delta.star,v[j],u[i])
  } 
}
contour(u,v,z,col=4,nlevels=7,add=T,lwd=2)
legend('bottomleft',bty='n',lty = c(NA,NA,1,2),pch=c(16,16,NA,NA),col = c(1,3,4,1),legend=c("Bayes","Cut","SMI-optimal","True"))

#dev.off()

##############
# pdfs for paper

#Fig 3 left column - marginal posterior density at cut, selected delta*-SMI and post
pdf('NormEx-MarginalDenisties.pdf')
par(mfrow=c(2,1));
par(mar=c(3.5,4,1,2))
plot(density(PHI.BAYES),xlab="",ylab="Post density phi",main="",xlim=c(-2.5,2)); 
title(xlab="phi", line=2)
lines(density(PHI.CUT),col=3); abline(v=phi_T,col=2)
lines(density(PHI[,which.delta]),col=4)
legend('topleft',bty="n",col=c(1,4,3,2),lty=c(1,1,1,1),legend=c("Bayes","SMI-optimal","Cut","True phi"))
par(mar=c(3.5,4,1,2))
plot(density(THETA.BAYES),xlab="",ylab="Post density theta",main="",xlim=c(-1,3)); 
title(xlab="theta", line=2)
lines(density(THETA.CUT),col=3); abline(v=theta_T,col=2)
lines(density(THETA[,which.delta]),col=4)
legend('topleft',bty="n",col=c(1,4,3,2),lty=c(1,1,1),legend=c("Bayes","SMI-optimal","Cut","True theta"))
dev.off()

#Fig 3 right column - MSE & exact ELPD_{y,z}
#MSE
pdf("NormEx-MSE-ELPD.pdf")
par(mfrow=c(2,1));
par(mar=c(3,4,1,2))
plot(log(v.gam),(mu.delta-phi_T)^2+s.delta^2,type='l',xlab='',ylab="PMSE",ylim=c(0.1,0.8))
lines(log(v.gam),(apply(THETA,2,mean)-theta_T)^2+apply(THETA,2,sd)^2,lty=2)
legend('bottomright',bty="n",col=c(1,1),lty=c(1,2),pch=c(NA,NA),legend=c("PMSE phi","PMSE theta"))
#exact ELPD - full elpd y and z - formula end of item (12) in notes
a=s1.sq; b=(1-rho)*s.delta^2; d=s2.sq; e=mu1; f=mu2
Delta=a*d-b^2
expr1=Delta^(-1)*(d*sy^2+d*(phi_T+theta_T-e)^2-2*b*(phi_T+theta_T-e)*(phi_T-f)+a*sz^2+a*(phi_T-f)^2)
ELPDyz=-(1/2)*expr1-(1/2)*log(Delta)-log(2*pi)
par(mar=c(3.5,4,0.5,2))
plot(log(v.gam),ELPDyz,type='l',xlab="")
title(xlab="log(delta)", line=2)
#points(log(v.gam),ELPDz+ELPDy,col=2) #after all that calculation the ELPD joint for y and z is basicly equal ELPDy+ELPDz
which.delta=which.max(ELPDyz)
delta.star=v.gam[which.delta]
abline(v=log(delta.star),lty=2)
legend('topright',bty="n",col=c(1,1),lty=c(1,2),pch=c(NA,NA),legend=c("ELPD(y,z)","log(delta*)"))
dev.off()

#Fig 1 - scatter plot of theta v. phi
#joint dbns
pdf("NormEx-Scatter.pdf")
par(mfrow=c(1,1))
par(mar=c(4,4,2,2))
#plot(PHI[1:(T.sim/10),which.delta],THETA[1:(T.sim/10),which.delta],col=2,xlab="phi",ylab="theta",pch=1)
#points(PHI.BAYES[1:(T.sim/10)],THETA.BAYES[1:(T.sim/10)],col=1,pch=1)
plot(PHI.BAYES[1:(T.sim/10)],THETA.BAYES[1:(T.sim/10)],col=1,pch='.',xlab="",ylab="",xlim=c(-2.5,2),ylim=c(-0.5,3.0))
title(xlab="phi", ylab="theta",line=2.5)

points(PHI.CUT[1:(T.sim/10)],THETA.CUT[1:(T.sim/10)],col=3,pch='.')
abline(v=phi_T,h=theta_T,lty=2)

M=300; u=seq(-2,2,length.out=M); v=seq(-0.5,2.5,length.out=M); z=matrix(NA,M,M); 
for (i in 1:M) { 
  for (j in 1:M) { 
    z[i,j]=pi.smi(delta.star,v[j],u[i])
  } 
}
contour(u,v,z,col=4,nlevels=7,add=T,lwd=2)
legend('bottomleft',bty='n',lty = c(NA,NA,1,2),pch=c(16,16,NA,NA),col = c(1,3,4,1),legend=c("Bayes","Cut","SMI-optimal","True"))
dev.off()
