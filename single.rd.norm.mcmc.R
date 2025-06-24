single.rd.norm.mcmc <- function(y,X,n.mcmc){
	
###
###  Subroutines
###

library(statmod)	

loglik <- function(d,mu,sig,theta,w) { # vectorized
	thresh=.000000001
	d=as.vector(d)
	mu=as.vector(mu)
	D.1=outer(d,(theta+1)/2,"+")
	D.2=outer(d,(theta-1)/2,"+")
	dens=(pnorm(D.1,mu,sig)-pnorm(D.2,mu,sig))%*%w/2
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
s2.save=rep(0,n.mcmc)
sig.save=rep(0,n.mcmc)
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

sig=sd(y)/10
s2=sig^2

beta.tune=.1
s2.tune=.5

###
###  Get G-L Quadrature Weights 
###

n.quad <- 7  # or more for higher accuracy
quad <- gauss.quad(n.quad,kind="legendre")
theta <- quad$nodes
w <- quad$weights

###
###  Begin MCMC Loop
###
	
for(k in 1:n.mcmc){
  if(k%%1000==0){cat(k," ")}	
	
  ###
  ### Sample s2
  ###

  s2.star=rnorm(1,s2,s2.tune)
  if(s2.star > 0){
  sig.star=sqrt(s2.star)
    mh.1=loglik(y,Xbeta,sig.star,theta,w)+ldIG(s2.star,q,r)
    mh.2=loglik(y,Xbeta,sig,theta,w)+ldIG(s2,q,r)
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
  mh.1=loglik(y,Xbeta.star,sig,theta,w)+sum(dnorm(beta.star,mu.beta,sig.beta,log=TRUE))
  mh.2=loglik(y,Xbeta,sig,theta,w)+sum(dnorm(beta,mu.beta,sig.beta,log=TRUE))
  mh=exp(mh.1-mh.2)
  if(mh>runif(1)){
    beta=beta.star	
    Xbeta=Xbeta.star
  }
	
  ###
  ### Save Sample 
  ###

  Dbar.save[k]=-2*loglik(y,Xbeta,sig,theta,w)
  
  beta.save[,k]=beta
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

Dhat=-2*loglik(y,X%*%beta.mn,sig.mn,theta,w)
pD=Dbar-Dhat
DIC=Dhat+2*pD

cat("DIC:",DIC,"pD:",pD,"Dbar:",Dbar,"Dhat:",Dhat,"\n")

###
###  Write Output
###
		
list(y=y,X=X,n.mcmc=n.mcmc,beta.save=beta.save,s2.save=s2.save,sig.save=sig.save,DIC=DIC,pD=pD,Dbar=Dbar,Dhat=Dhat,n.burn=n.burn)	
	
}
