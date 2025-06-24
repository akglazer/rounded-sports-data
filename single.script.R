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
y=rb.player.df$rushing_yards
X=model.matrix(~game_seconds_remaining+score_differential+ydstogo+yardline_100+rest+opponent_rush_epa_allowed+is_home+travel_distance,data=rb.player.df)
X.sc=X
X.sc[,-c(1,8)]=scale(X[,-c(1,8)])

plot(table(c(y,min(y):max(y)))-1,type="h",lwd=2,xlab="yards gained",ylab="frequency",main="")
abline(v=0,col=rgb(0,0,1,.75),lty=3,lwd=2)

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
                           labels = c("game seconds remaining","score differential","yards to go","line of scrimmage","days of rest","opp. EPA","homefield status","travel distance")))

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
###  Compare rd-Norm and rd-ALD Inference 
###

cbind(apply(rd.norm.out$beta.save,1,mean),apply(rd.ald.out$beta.save,1,mean))
plot(apply(rd.norm.out$beta.save,1,mean),apply(rd.ald.out$beta.save,1,mean),asp=TRUE)
abline(0,1,col=2)

cbind(apply(rd.norm.out$beta.save,1,sd),apply(rd.ald.out$beta.save,1,sd))
plot(apply(rd.norm.out$beta.save,1,sd),apply(rd.ald.out$beta.save,1,sd),asp=TRUE)
abline(0,1,col=2)

###
###  Posterior Predictive Check for Skewness 
###

source("single.ppc.R")
p.val.norm=single.ppc(rd.norm.out)
p.val.ald=single.ppc(rd.ald.out)

###
###  Summarize all DIC
###

DIC.mat=matrix(0,2,4)
DIC.mat[1,]=c(rd.norm.out$DIC,rd.norm.out$pD,rd.norm.out$Dbar,rd.norm.out$Dhat)
DIC.mat[2,]=c(rd.ald.out$DIC,rd.ald.out$pD,rd.ald.out$Dbar,rd.ald.out$Dhat)
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

