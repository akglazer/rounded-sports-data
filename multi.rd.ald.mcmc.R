multi.rd.ald.mcmc <- function(y,X,player,n.mcmc){
	
###
###  Subroutines
###
	
library(statmod)
	
pALD<-function(q,mu,sigma,tau){
  prob=ifelse(test=q<mu,yes=tau*exp((1-tau)*(q-mu)/sigma),no=1-(1-tau)*exp(-(tau)*(q-mu)/sigma))
  prob 
}	
	
lddALD <- function(y,mu,sig,tau,theta,w) { # vectorized
	thresh=.000000001
	y=as.vector(y)
	mu=as.vector(mu)
	Y.1=outer(y,(theta+1)/2,"+")
	Y.2=outer(y,(theta-1)/2,"+")
	dens=(pALD(Y.1,mu,sig,tau)-pALD(Y.2,mu,sig,tau))%*%w/2
	dens[dens<thresh]=thresh
	sum(log(dens))
}

lddALD.nosum <- function(y,mu,sig,tau,theta,w) { # vectorized for y and mu
	thresh=.000000001
	y=as.vector(y)
	mu=as.vector(mu)
	Y.1=outer(y,(theta+1)/2,"+")
	Y.2=outer(y,(theta-1)/2,"+")
	dens=(pALD(Y.1,mu,sig,tau)-pALD(Y.2,mu,sig,tau))%*%w/2
	dens[dens<thresh]=thresh
	log(dens)
}

ldIG <- function(x,q,r){
  -(q+1)*log(x)-1/r/x-q*log(r)-lgamma(q)
}

intloglik <- function(y,mu0,sig0,mu,sig,tau,theta,w,z.gh,w.gh,grp){  # G-H Quad 
  n=max(grp)
  tmp.sum=0
  for(i in 1:n){
    beta0.mat=matrix(mu0+sqrt(2)*sig0*z.gh,length(y[grp==i]),length(z.gh),byrow=TRUE)
    ll.vec=apply(matrix(lddALD.nosum(rep(y[grp==i],length(z.gh)),c(beta0.mat+mu[grp==i]),sig,tau,theta,w),length(y[grp==i]),length(z.gh)),2,sum)
    ll.max=max(ll.vec)
    tmp.sum=tmp.sum+ll.max+log(sum(w.gh*exp(ll.vec-ll.max)))-log(pi)/2
  }
  tmp.sum
}

###
###  Get G-H Quadrature Weights 
###

L.gh=11
gh=gauss.quad(L.gh,kind="hermite")
z.gh=gh$nodes
w.gh=gh$weights

###
###  Get G-L Quadrature Weights 
###

n.quad <- 7  
quad <- gauss.quad(n.quad,kind="legendre")
theta <- quad$nodes
w <- quad$weights
		
###
###  Setup Variables 
###

n.burn=round(.2*n.mcmc)
pp=dim(X)[2]	
m=length(y)
n=length(unique(player))

beta.0.save=matrix(0,n,n.mcmc)
beta.save=matrix(0,pp,n.mcmc)
sig.save=rep(0,n.mcmc)
s2.save=rep(0,n.mcmc)
tau.save=rep(0,n.mcmc)
mu.0.save=rep(0,n.mcmc)
s2.0.save=rep(0,n.mcmc)
Dbar.save=rep(0,n.mcmc)

###
###  Priors and Starting Values 
###

mu.beta=rep(0,pp)
sig.beta=1000000
Sig.beta.inv=diag(pp)/(sig.beta^2)
q=.001
r=1000

q.0=.001
r.0=1000

mu.0=0
mu.00=0
s2.00=1000
sig.00=sqrt(s2.00)

beta=as.vector(coef(lm(y ~ 1+X)))[-1]
beta.0=rep(0,n)
for(i in 1:n){
  beta.0[i]=mean(y[player==i]-X[player==i,]%*%beta)
}
s2.0=var(beta.0)
sig.0=sqrt(s2.0)
beta.0.big=beta.0[player]

beta=as.vector(coef(lm(y-beta.0.big ~ 0+X)))
#beta=rep(0,dim(X)[2])

Xbeta=X%*%beta

sig=sd(y)/20
s2=sig^2
tau=.5

beta.0.tune=1
tau.tune=.01
beta.tune=.02
s2.tune=.1

###
###  Begin MCMC Loop
###
	
for(k in 1:n.mcmc){
  if(k%%100==0){cat(k," ")}	
	
  ###
  ### Sample s2
  ###

  s2.star=rnorm(1,s2,s2.tune)
  if(s2.star > 0){
    sig.star=sqrt(s2.star)
    mh.1=sum(lddALD(y,beta.0.big+Xbeta,sig.star,tau,theta,w))+ldIG(s2.star,q,r)
    mh.2=sum(lddALD(y,beta.0.big+Xbeta,sig,tau,theta,w))+ldIG(s2,q,r)
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
  mh.1=sum(lddALD(y,beta.0.big+Xbeta.star,sig,tau,theta,w))+sum(dnorm(beta.star,mu.beta,sig.beta,log=TRUE))
  mh.2=sum(lddALD(y,beta.0.big+Xbeta,sig,tau,theta,w))+sum(dnorm(beta,mu.beta,sig.beta,log=TRUE))
  mh=exp(mh.1-mh.2)
  if(mh>runif(1)){
    beta=beta.star	
    Xbeta=Xbeta.star
  }
	
  ###
  ### Sample beta.0
  ###

  for(i in 1:n){
    beta.0.star=rnorm(1,beta.0[i],beta.0.tune)
    mh.1=sum(lddALD(y[player==i],beta.0.star+Xbeta[player==i],sig,tau,theta,w))+dnorm(beta.0.star,mu.0,sig.0,log=TRUE)
    mh.2=sum(lddALD(y[player==i],beta.0[i]+Xbeta[player==i],sig,tau,theta,w))+dnorm(beta.0[i],mu.0,sig.0,log=TRUE)
    mh=exp(mh.1-mh.2)
    if(mh>runif(1)){
      beta.0[i]=beta.0.star	
    }
  }
  beta.0.big=beta.0[player]
  
  ###
  ### Sample mu.0 
  ###
  
  tmp.var=1/(n/s2.0+1/s2.00)
  tmp.mn=tmp.var*(sum(beta.0)/s2.0+mu.00/s2.00)	
  mu.0=rnorm(1,tmp.mn,sqrt(tmp.var))
  
  ###
  ### Sample s2.0 
  ###
  
  q.tmp=n/2+q.0 
  r.tmp=1/(sum((beta.0-mu.0)^2)/2+1/r.0)
  s2.0=1/rgamma(1,q.tmp,,r.tmp)
  sig.0=sqrt(s2.0)
  
  ###
  ### Sample tau
  ###

  tau.star=rnorm(1,tau,tau.tune)
  if(tau.star > 0 & tau.star < 1){
    mh.1=sum(lddALD(y,beta.0.big+Xbeta,sig,tau.star,theta,w))
    mh.2=sum(lddALD(y,beta.0.big+Xbeta,sig,tau,theta,w))
    mh=exp(mh.1-mh.2)
    if(mh>runif(1)){
      tau=tau.star	
    }
  }
  	
  ###
  ### Save Sample 
  ###

  Dbar.save[k]=-2*intloglik(y,mu.0,sig.0,Xbeta,sig,tau,theta,w,z.gh,w.gh,player)
  
  beta.save[,k]=beta
  beta.0.save[,k]=beta.0
  tau.save[k]=tau
  mu.0.save[k]=mu.0
  s2.0.save[k]=s2.0
  s2.save[k]=s2
  sig.save[k]=sig
		
};cat("\n")

###
###  DIC Calc  
###

Dbar=mean(Dbar.save[n.burn:n.mcmc])	
Dhat=0
beta.mn=apply(beta.save[,n.burn:n.mcmc],1,mean)
mu.0.mn=mean(mu.0.save[n.burn:n.mcmc])
sig.0.mn=mean(sqrt(s2.0.save[n.burn:n.mcmc]))
sig.mn=mean(sig.save[n.burn:n.mcmc])
tau.mn=mean(tau.save[n.burn:n.mcmc])

Dhat=-2*intloglik(y,mu.0.mn,sig.0.mn,X%*%beta.mn,sig.mn,tau.mn,theta,w,z.gh,w.gh,player)
pD=Dbar-Dhat
DIC=Dhat+2*pD

cat("DIC:",DIC,"pD:",pD,"Dbar:",Dbar,"Dhat:",Dhat,"\n")

###
###  Write Output
###
		
list(y=y,X=X,n.mcmc=n.mcmc,beta.save=beta.save,sig.save=sig.save,s2.save=s2.save,beta.0.save=beta.0.save,mu.0.save=mu.0.save,s2.0.save=s2.0.save,tau.save=tau.save,DIC=DIC,pD=pD,Dhat=Dhat,Dbar=Dbar)	
	
}
