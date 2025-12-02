single.rd.tald.mcmc <- function(y,X,b,n.mcmc){
	
###
###  Subroutines
###
	
library(statmod)	
	
pALD<-function(q,mu,sigma,tau){
  prob=ifelse(test=q<mu,yes=tau*exp((1-tau)*(q-mu)/sigma),no=1-(1-tau)*exp(-(tau)*(q-mu)/sigma))
  prob 
}	

pTALD<-function(q,mu,sigma,tau,b){
	if(!is.matrix(b)){
    F.q=pALD(q,mu,sigma,tau)
    F.u=pALD(b,mu,sigma,tau) 
    F.l=pALD(b-100,mu,sigma,tau) 
    prob=(F.q-F.l)/(F.u-F.l) 
    prob[q<(b-100)]=0
    prob[q>b]=1
	}
	if(is.matrix(b)){
		p.fcn<-function(q,mu,sig,tau,b){
      F.q=pALD(q,mu,sigma,tau)
      F.u=pALD(b,mu,sigma,tau) 
      F.l=pALD(b-100,mu,sigma,tau) 
      prob=(F.q-F.l)/(F.u-F.l) 
      prob[q<(b-100)]=0
      prob[q>b]=1
      prob
		}
		prob=apply(b,2,p.fcn,q=q,mu=mu,sig=sigma,tau=tau)
	}
	prob
}

pATALD <- function(q,mu,sig,tau,b,theta.big,w.big) { # vectorized
	B=outer(b,theta.big/2+.5,"+")
	prob=pTALD(q,mu,sig,tau,B)%*%w.big/2
  prob	
}
	
ldpATALD <- function(d,mu,sig,tau,b,theta,w,theta.big,w.big) { # vectorized
  thresh=1e-20
  d=as.vector(d)
  mu=as.vector(mu)
  n.theta=length(theta)
  D.1=outer(d,(theta+1)/2,"+")
  D.2=outer(d,(theta-1)/2,"+")
  F.1=matrix(pATALD(as.vector(D.1),mu,sig,tau,b,theta.big,w.big),,n.theta)
  F.2=matrix(pATALD(as.vector(D.2),mu,sig,tau,b,theta.big,w.big),,n.theta)
  dens=(F.1-F.2)%*%w/2
  dens[dens<thresh]=thresh
  sum(log(dens))
}

ldIG <- function(x,q,r){
  -(q+1)*log(x)-1/r/x-q*log(r)-lgamma(q)
}

###
###  Setup Variables 
###

n.burn=round(.2*n.mcmc)
pp=dim(X)[2]	
m=length(y)

beta.save=matrix(0,pp,n.mcmc)
sig.save=rep(0,n.mcmc)
s2.save=rep(0,n.mcmc)
tau.save=rep(0,n.mcmc)
Dbar.save=rep(0,n.mcmc)

###
###  Priors and Starting Values 
###

mu.beta=rep(0,pp)
sig.beta=1000000
Sig.beta.inv=diag(pp)/(sig.beta^2)
q=.001
r=1000

beta=as.vector(coef(lm(y ~ 0+X)))
Xbeta=X%*%beta

sig=sd(y)/20
s2=sig^2
tau=.5

tau.tune=.05
beta.tune=.1
s2.tune=.5

###
###  Get Large Set of G-L Quadrature Weights (for pATALD)
###

n.quad.big <- 1001 # needs to be big for approximation to work well in this case
quad.big <- gauss.quad(n.quad.big,kind="legendre")
theta.big <- quad.big$nodes
w.big <- quad.big$weights

###
###  Get Small Set of G-L Quadrature Weights (for ldpATALD)
###

n.quad <- 7 
quad <- gauss.quad(n.quad,kind="legendre")
theta <- quad$nodes
w <- quad$weights

###
###  Begin MCMC Loop
###
	
for(k in 1:n.mcmc){
  if(k%%10==0){cat(k," ")}	
	
  ###
  ### Sample s2
  ###

  s2.star=rnorm(1,s2,s2.tune)
  if(s2.star > 0){
    sig.star=sqrt(s2.star)
    mh.1=ldpATALD(y,Xbeta,sig.star,tau,b,theta,w,theta.big,w.big)+ldIG(s2.star,q,r) 
    mh.2=ldpATALD(y,Xbeta,sig,tau,b,theta,w,theta.big,w.big)+ldIG(s2,q,r)
    mh=exp(mh.1-mh.2)
    if(mh>runif(1)){
      s2=s2.star	
      sig=sig.star
    }
  }
	
  ###
  ### Sample beta 
  ###

  beta.star=rnorm(pp,beta,beta.tune)
  Xbeta.star=X%*%beta.star
  mh.1=ldpATALD(y,Xbeta.star,sig,tau,b,theta,w,theta.big,w.big)+sum(dnorm(beta.star,mu.beta,sig.beta,log=TRUE))
  mh.2=ldpATALD(y,Xbeta,sig,tau,b,theta,w,theta.big,w.big)+sum(dnorm(beta,mu.beta,sig.beta,log=TRUE))
  mh=exp(mh.1-mh.2)
  if(mh>runif(1)){
    beta=beta.star	
    Xbeta=Xbeta.star
  }
	
  ###
  ### Sample tau
  ###

  tau.star=rnorm(1,tau,tau.tune)
  if(tau.star > 0 & tau.star < 1){
    mh.1=ldpATALD(y,Xbeta,sig,tau.star,b,theta,w,theta.big,w.big)
    mh.2=ldpATALD(y,Xbeta,sig,tau,b,theta,w,theta.big,w.big)
    mh=exp(mh.1-mh.2)
    if(mh>runif(1)){
      tau=tau.star	
    }
  }
  	
  ###
  ### Save Sample 
  ###

  Dbar.save[k]=-2*ldpATALD(y,Xbeta,sig,tau,b,theta,w,theta.big,w.big)
  
  beta.save[,k]=beta
  tau.save[k]=tau
  s2.save[k]=s2
  sig.save[k]=sig
		
};cat("\n")

###
###  DIC Calc  
###

Dbar=mean(Dbar.save[n.burn:n.mcmc])	
beta.mn=apply(beta.save[,n.burn:n.mcmc],1,mean)
s2.mn=mean(s2.save[n.burn:n.mcmc])
sig.mn=mean(sig.save[n.burn:n.mcmc])
tau.mn=mean(tau.save[n.burn:n.mcmc])

Dhat=-2*ldpATALD(y,X%*%beta.mn,sig.mn,tau.mn,b,theta,w,theta.big,w.big)
pD=Dbar-Dhat
DIC=Dhat+2*pD

cat("DIC:",DIC,"pD:",pD,"Dbar:",Dbar,"Dhat:",Dhat,"\n")

###
###  Write Output
###
		
list(y=y,X=X,b=b,n.mcmc=n.mcmc,beta.save=beta.save,sig.save=sig.save,s2.save=s2.save,tau.save=tau.save,DIC=DIC,pD=pD,Dhat=Dhat,Dbar=Dbar,n.burn=n.burn)	
	
}
