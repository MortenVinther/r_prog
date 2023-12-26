#HCR icons
library(tidyverse)

#fixed F

Biomass<-0:100
FF<-Biomass
FF[]<-1


png(file=file.path(www,"fixed_F.png"),width=500,height=500)
par(mar=c(5,6,1,1)) #c(bottom, left, top, right
plot(Biomass,FF,type='l',col='red',lwd=5,xlab=list('Biomass',cex=4),ylab=list('F',cex=4),xaxt='n',yaxt='n',ylim=c(0.0,1.3))
dev.off()

T1<-20
T2<-60
FF[Biomass<T1]<-0
a<-(1-0)/(T2-T1)
b<-0-a*T1
FF[Biomass>=T1 & Biomass<T2]<-a*Biomass[Biomass>=T1 & Biomass<T2]+b


png(file=file.path(www,"AR_F_SSB.png"),width=500,height=500)
par(mar=c(5,6,1,1)) #c(bottom, left, top, right
plot(Biomass,FF,type='l',col='red',lwd=5,xlab=list('Spawning Biomass',cex=4),ylab=list('F',cex=4),xaxt='n',yaxt='n',ylim=c(0.0,1.3))

lines(x=c(T1,T1),y=c(0,1.0),lwd=4,col='blue')
text(x=T1,y=1.2,'T1',cex=4,col='blue')

lines(x=c(T2,T2),y=c(0,1.0),lwd=4,col='green')
text(x=T2,y=1.2,'T2',cex=4,col='green')

dev.off()


png(file=file.path(www,"AR_F_TSB.png"),width=500,height=500)
par(mar=c(5,6,1,1)) #c(bottom, left, top, right
plot(Biomass,FF,type='l',col='red',lwd=5,xlab=list('Total Biomass',cex=4),ylab=list('F',cex=4),xaxt='n',yaxt='n',ylim=c(0.0,1.3))

lines(x=c(T1,T1),y=c(0,1.0),lwd=4,col='blue')
text(x=T1,y=1.2,'T1',cex=4,col='blue')

lines(x=c(T2,T2),y=c(0,1.0),lwd=4,col='green')
text(x=T2,y=1.2,'T2',cex=4,col='green')

dev.off()

#### SSB rec plot


T1<-0
T2<-60
FF[Biomass<T1]<-0
a<-(1-0)/(T2-T1)
b<-0-a*T1
FF[Biomass>=T1 & Biomass<T2]<-a*Biomass[Biomass>=T1 & Biomass<T2]+b


png(file=file.path(www,"rec_deter.png"),width=500,height=500)

par(mar=c(5,6,1,1)) #c(bottom, left, top, right
plot(Biomass,FF,type='l',col='red',lwd=5,xlab=list('Spawning Biomass',cex=4),ylab=list('Recruitment',cex=4),xaxt='n',yaxt='n',ylim=c(0.0,1.7))
n<-15
recB<-runif(n,15,100)

rec<-recB*a+b
rec[recB>T2]<-1
points(x=recB,y=rec,cex=2,pch=22,col='blue')

dev.off()


png(file=file.path(www,"rec_stoch.png"),width=500,height=500)

par(mar=c(5,6,1,1)) #c(bottom, left, top, right
plot(Biomass,FF,type='l',col='red',lwd=5,xlab=list('Spawning Biomass',cex=4),ylab=list('Recruitment',cex=4),xaxt='n',yaxt='n',ylim=c(0.0,1.7))
rec<-rec*rlnorm(n, meanlog = 0, sdlog = 0.4)
points(x=recB,y=rec,cex=2,pch=22,col='blue')

dev.off()



