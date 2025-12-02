###
###  Prepare Data
###

rb.df=read.csv("nfl-rb-pbp-rushing.csv")
rb.sm.df=rb.df[rb.df$week<18 & !is.na(rb.df$opponent_rush_epa_allowed) & !is.na(rb.df$rushing_yards),]
m=dim(rb.sm.df)[1]
player.uniq=unique(rb.sm.df$rusher_player_id)
n=length(player.uniq)
rb.sm.df$player=rep(0,m)
for(i in 1:n){
  rb.sm.df$player[rb.sm.df$rusher_player_id==player.uniq[i]]=i	
}

player.name=rep(0,n)
n.i.vec=rep(0,n)
for(i in 1:n){
  player.name[i]=rb.sm.df$rusher_player_name[rb.sm.df$rusher_player_id==player.uniq[i]][1]	
  n.i.vec[i]=sum(rb.sm.df$player==i)
}

###
###  Prepare Single Player Data
###

rb.player.df=rb.sm.df[rb.sm.df$rusher_player_id==player.uniq[16],] # Bijan Robinson 
rb.player.df$down_cat=as.factor(rb.player.df$down)

y=rb.player.df$rushing_yards
X=model.matrix(~down+I(down^2)+game_seconds_remaining+score_differential+ydstogo+yardline_100+rest+opponent_rush_epa_allowed+is_home+travel_distance,data=rb.player.df)
X.sc=X
X.sc[,-c(1,10)]=scale(X[,-c(1,10)])
X.sc[,3]=X.sc[,2]^2
b=rb.player.df$yardline_100

datawide.df=data.frame(y,X[,-1])
colnames(datawide.df)[1]="rushing_yards"

library(ggplot2)
library(tidyr)
library(dplyr)

datalong.df <- pivot_longer(datawide.df, cols = -rushing_yards, names_to = "Variable", values_to = "covariate") %>%
  mutate(Variable = factor(Variable, levels = colnames(datawide.df)[-1]))  

datalong.df <- datalong.df %>%
  mutate(Variable = factor(Variable,
                           levels = colnames(datawide.df)[-1],
                           labels = c("down","down^2","game seconds remaining","score differential","yards to go","line of scrimmage","days of rest","opp. EPA","homefield status","travel distance")))

ggplot(datalong.df, aes(x = covariate, y = rushing_yards)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ Variable, scales = "free_x") +
  theme_minimal()+
	theme(panel.border = element_rect(color = "black", fill = NA),axis.title.x = element_blank())+
  labs(y = "rushing yards") 

###
###  Fit Rounded Norm Model
###

library(tictoc)
n.mcmc=100000
n.burn=round(.2*n.mcmc)

source("single.rd.norm.mcmc.R")
tic()
rd.norm.out=single.rd.norm.mcmc(y,X.sc,n.mcmc)
rd.norm.toc=toc()

out=rd.norm.out
matplot(t(out$beta.save),type="l",lty=1)
plot(out$s2.save,type="l")

library(vioplot)
vioplot(data.frame(t(out$beta.save)),names=colnames(X.sc))
abline(h=0,col=rgb(0,0,0,.25),lty=2)

###
###  Fit rounded ALD Model
###

source("single.rd.ald.mcmc.R")
tic()
rd.ald.out=single.rd.ald.mcmc(y,X.sc,n.mcmc)
rd.ald.toc=toc()

out=rd.ald.out
matplot(t(out$beta.save),type="l",lty=1)
plot(out$s2.save,type="l")
plot(out$tau.save,type="l")

library(vioplot)
vioplot(data.frame(t(out$beta.save)),names=colnames(X.sc))
abline(h=0,col=rgb(0,0,0,.25),lty=2)

###
###  Fit rounded TALD Model
###

source("single.rd.tald.mcmc.R")
tic()
rd.tald.out=single.rd.tald.mcmc(y,X.sc,b,n.mcmc)
rd.tald.toc=toc()

out=rd.tald.out
matplot(t(out$beta.save),type="l",lty=1)
plot(out$s2.save,type="l")
plot(out$tau.save,type="l")

library(vioplot)
vioplot(data.frame(t(out$beta.save)),names=colnames(X.sc))
abline(h=0,col=rgb(0,0,0,.25),lty=2)

###
###  Compare rd-ALD and rd-TALD Inference 
###

plot(apply(rd.ald.out$beta.save,1,mean),apply(rd.tald.out$beta.save,1,mean),asp=TRUE)
abline(0,1,col=2)

coef.idx=0:(dim(X.sc)[2]-1)
p=length(coef.idx)
ald.mn=apply(rd.ald.out$beta.save,1,mean)
ald.l=apply(rd.ald.out$beta.save,1,quantile,.025)
ald.u=apply(rd.ald.out$beta.save,1,quantile,.975)
tald.mn=apply(rd.tald.out$beta.save,1,mean)
tald.l=apply(rd.tald.out$beta.save,1,quantile,.025)
tald.u=apply(rd.tald.out$beta.save,1,quantile,.975)

matplot(coef.idx,cbind(ald.l,ald.u),type="n",xlab="",ylab="",xaxt="n")
segments(coef.idx-.1,ald.l,coef.idx-.1,ald.u,lwd=2)
points(coef.idx-.1,ald.mn,pch=16,lwd=2)
segments(coef.idx+.1,tald.l,coef.idx+.1,tald.u,lwd=2,col=2)
points(coef.idx+.1,tald.mn,pch=16,lwd=2,col=2)
axis(1,coef.idx,labels=as.expression(sapply(0:(p-1),function(j){bquote(beta[.(j)])})))
abline(h=0,lty=2,col=rgb(0,0,0,.4))
legend("top",bty="n",col=1:2,lwd=2,legend=c("ALD","trunc. ALD"))

###
###  Compare rd-Norm and rd-ALD Inference 
###
	
cbind(apply(rd.norm.out$beta.save,1,mean),apply(rd.ald.out$beta.save,1,mean))
plot(apply(rd.norm.out$beta.save,1,mean),apply(rd.ald.out$beta.save,1,mean),asp=TRUE)
abline(0,1,col=2)

cbind(apply(rd.norm.out$beta.save,1,sd),apply(rd.ald.out$beta.save,1,sd))
plot(apply(rd.norm.out$beta.save,1,sd),apply(rd.ald.out$beta.save,1,sd),asp=TRUE)
abline(0,1,col=2)

###
###  Posterior Predictive Check for Skewness (b/w rd-norm and rd-ald) 
###

source("single.ppc.R")
p.val.norm=single.ppc(rd.norm.out)
p.val.ald=single.ppc(rd.ald.out)
c(p.val.norm,p.val.ald)

###
###  Summarize all DIC
###

DIC.mat=matrix(0,3,4)
DIC.mat[1,]=c(rd.norm.out$DIC,rd.norm.out$pD,rd.norm.out$Dbar,rd.norm.out$Dhat)
DIC.mat[2,]=c(rd.ald.out$DIC,rd.ald.out$pD,rd.ald.out$Dbar,rd.ald.out$Dhat)
DIC.mat[3,]=c(rd.tald.out$DIC,rd.tald.out$pD,rd.tald.out$Dbar,rd.tald.out$Dhat)
DIC.mat

###
###  Summarize Posteriors 
###

round(cbind(apply(rd.norm.out$beta.save,1,mean), apply(rd.norm.out$beta.save,1,sd), t(apply(rd.norm.out$beta.save,1,quantile,c(0.025,.975)))),2)
round(cbind(apply(rd.ald.out$beta.save,1,mean), apply(rd.ald.out$beta.save,1,sd), t(apply(rd.ald.out$beta.save,1,quantile,c(0.025,.975)))),2)

round(cbind(apply(rd.norm.out$beta.save,1,mean), apply(rd.norm.out$beta.save,1,sd), t(apply(rd.norm.out$beta.save,1,quantile,c(0.05,.95)))),2)
round(cbind(apply(rd.ald.out$beta.save,1,mean), apply(rd.ald.out$beta.save,1,sd), t(apply(rd.ald.out$beta.save,1,quantile,c(0.05,.95)))),2)

round(c(mean(rd.ald.out$sig.save), sd(rd.ald.out$sig.save), quantile(rd.ald.out$sig.save,c(0.05,.95))),2)
round(c(mean(rd.ald.out$tau.save), sd(rd.ald.out$tau.save), quantile(rd.ald.out$tau.save,c(0.05,.95))),2)


