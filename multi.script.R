###
###  Read in Data
###

rb.df=read.csv("nfl-rb-pbp-rushing.csv")

###
###  Prepare Data
###

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

y=rb.sm.df$rushing_yards
X=model.matrix(~0+age+I(age^2)+down+I(down^2)+game_seconds_remaining+score_differential+ydstogo+yardline_100+rest+opponent_rush_epa_allowed+is_home+travel_distance,data=rb.sm.df)
X.sc=scale(X)
X.sc[,2]=X.sc[,1]^2
X.sc[,4]=X.sc[,3]^2
X.sc[,11]=X[,11]

###
###  Fit Rounded Norm Model
###

library(tictoc)

n.mcmc=200000
n.burn=round(.2*n.mcmc)

source("multi.rd.norm.mcmc.R")
tic()
rd.norm.out=multi.rd.norm.mcmc(y,X.sc,c(rb.sm.df$player),n.mcmc)
rd.norm.toc=toc()

###
###  Fit Rounded ALD Model
###

source("multi.rd.ald.mcmc.R")
tic()
rd.ald.out=multi.rd.ald.mcmc(y,X.sc,c(rb.sm.df$player),n.mcmc)
rd.ald.toc=toc()

###
###  Summarize all DIC
###

DIC.mat=matrix(0,2,4)
DIC.mat[1,]=c(rd.norm.out$DIC,rd.norm.out$pD,rd.norm.out$Dbar,rd.norm.out$Dhat)
DIC.mat[2,]=c(rd.ald.out$DIC,rd.ald.out$pD,rd.ald.out$Dbar,rd.ald.out$Dhat)

###
###  Save Output
###

save.image(file="multi.RData")

###
### Check Convergence
###
### Note:  These draw slowly for large n.mcmc
###

load("multi.RData")

#matplot(t(rd.norm.out$beta.0.save[1:4,]),type="l",lty=1)
#matplot(t(rd.norm.out$beta.save),type="l",lty=1)
#plot(rd.norm.out$mu.0.save,type="l")
#plot(rd.norm.out$s2.0.save,type="l")

#matplot(t(rd.ald.out$beta.0.save[1:4,]),type="l",lty=1)
#matplot(t(rd.ald.out$beta.save),type="l",lty=1)
#plot(rd.ald.out$s2.save,type="l")
#plot(rd.ald.out$mu.0.save,type="l")
#plot(rd.ald.out$s2.0.save,type="l")
#plot(rd.ald.out$tau.save,type="l")

###
###  DIC Summary 
###

c(rd.norm.out$DIC,rd.norm.out$pD)
c(rd.ald.out$DIC,rd.ald.out$pD)

###
###  Beta_0 Plots
###

out=rd.norm.out
beta.0.mn=apply(out$beta.0.save[,-(1:n.burn)],1,mean)
beta.0.lu=apply(out$beta.0.save[,-(1:n.burn)],1,quantile,c(0.025,.975))
xlim=range(c(beta.0.lu))
png(file="beta0intervals.png",height=16,width=10,units="in",res=400)
par(mar=c(5,6,4,2))
layout(matrix(1:2,1,2))
plot(0,type="n",yaxt="n",xlab=bquote(beta[0]),xlim=xlim,ylim=c(1,n),ylab="",axes=FALSE,main="rounded-Norm")
axis(1)
segments(beta.0.lu[1,order(beta.0.mn)],1:n,beta.0.lu[2,order(beta.0.mn)],1:n,col=rgb(0,0,0,.6),lwd=2)
points(sort(beta.0.mn),1:n,col=rgb(0,0,0,.8),pch=16,cex=.5)
mtext(player.name[order(beta.0.mn)],side=2,at=1:n,las=2,line=.1,cex=.6)

out=rd.ald.out
beta.0.mn=apply(out$beta.0.save[,-(1:n.burn)],1,mean)
beta.0.lu=apply(out$beta.0.save[,-(1:n.burn)],1,quantile,c(0.025,.975))
xlim=range(c(beta.0.lu))
plot(0,type="n",yaxt="n",xlab=bquote(beta[0]),xlim=xlim,ylim=c(1,n),ylab="",axes=FALSE,main="rounded-ALD")
axis(1)
segments(beta.0.lu[1,order(beta.0.mn)],1:n,beta.0.lu[2,order(beta.0.mn)],1:n,col=rgb(0,0,0,.6),lwd=2)
points(sort(beta.0.mn),1:n,col=rgb(0,0,0,.8),pch=16,cex=.5)
mtext(player.name[order(beta.0.mn)],side=2,at=1:n,las=2,line=.1,cex=.6)
dev.off()

plot(apply(rd.norm.out$beta.0.save[,-(1:n.burn)],1,mean),apply(rd.ald.out$beta.0.save[,-(1:n.burn)],1,mean),pch=16,col=rgb(0,0,0,.4),asp=TRUE,type="n")
text(apply(rd.norm.out$beta.0.save[,-(1:n.burn)],1,mean),apply(rd.ald.out$beta.0.save[,-(1:n.burn)],1,mean),labels=player.name,cex=.5)
abline(0,1,col=rgb(1,0,0,.4))

###
###  Beta Posterior Summary 
###

round(cbind(apply(rd.norm.out$beta.save,1,mean),apply(rd.norm.out$beta.save,1,sd),t(apply(rd.norm.out$beta.save,1,quantile,c(0.025,.975)))),2)

round(cbind(apply(rd.ald.out$beta.save,1,mean),apply(rd.ald.out$beta.save,1,sd),t(apply(rd.ald.out$beta.save,1,quantile,c(0.025,.975)))),2)

library(vioplot)
vioplot(data.frame(t(rd.ald.out$beta.save)),names=c("age","age^2","down","down^2","game seconds remaining","score differential","yards to go","line of scrimmage","days of rest","opp. EPA","homefield status","travel distance"),las=2,mar=c(10,4,4,2))
abline(h=0,col=2,lty=2)

###
###  tau Posterior Summary 
###

quantile(rd.ald.out$tau.save,c(.025,.975))


