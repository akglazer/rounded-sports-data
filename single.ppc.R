single.ppc <- function(out){
	
###
### Subroutines 
###

library(parameters)
	
rALD <- function(n,mu,sig,tau){
  z=rbinom(n,1,tau)	
  y=rep(0,n)
  if(length(mu)==1){
    y[z==1]=mu-rexp(sum(z),(1-tau)/sig)
    y[z==0]=mu+rexp(sum(1-z),tau/sig)
  }
  if(length(mu)==n){
    y[z==1]=mu[z==1]-rexp(sum(z),(1-tau)/sig)
    y[z==0]=mu[z==0]+rexp(sum(1-z),tau/sig)
  }
  y
}	
	
###
### Make variables
###

n.mcmc=out$n.mcmc	
n.burn=out$n.burn
beta.save=out$beta.save
sig.save=out$sig.save
if(!is.null(out$tau.save)){
	  tau.save=out$tau.save
}
y=out$y
X=out$X
m=length(y)

p.val.skew.save=rep(0,n.mcmc)
p.val.kurt.save=rep(0,n.mcmc)

###
### PP Loop
###

for(k in 1:n.mcmc){
	if(k%%10000==0){cat(k," ")}

	Xbeta=X%*%beta.save[,k]
	
  if(is.null(out$tau.save)){
  	y.tilda=rnorm(m,Xbeta,sig.save[k])
  }
	
  if(!is.null(out$tau.save)){
  	y.tilda=rALD(m,Xbeta,sig.save[k],tau.save[k])
  }
	
  y.pred=round(y.tilda+runif(m,-.5,.5))
  skew.data=skewness(y-Xbeta,type=2)$Skewness
  skew.pred=skewness(y.pred-Xbeta,type=2)$Skewness
  kurt.data=kurtosis(y-Xbeta)$Kurtosis
  kurt.pred=kurtosis(y.pred-Xbeta)$Kurtosis
	p.val.skew.save[k]=skew.pred>skew.data
	p.val.kurt.save[k]=kurt.pred>kurt.data
		
};cat("\n")
	
###	
### Write Output	
###	

p.val.skew=mean(p.val.skew.save)
p.val.kurt=mean(p.val.kurt.save)
cat("skewness p-value:",p.val.skew,"kurtosis p-value:",p.val.kurt)
	
list(p.val.skew=p.val.skew,p.val.kurt=p.val.kurt)	
	
}
