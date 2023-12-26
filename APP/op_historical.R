dat<-Read.summary.data(extend=FALSE,read.init.function=F)
dat<-subset(dat,Species.n>=first.VPA)
summary(dat)

dat<-data.frame(dat,deadM1w=dat$M1*dat$N.bar*dat$west, deadM1w.core=dat$M1*dat$N.bar*dat$west*dat$prop.in,deadM1=dat$M1*dat$N.bar,
                    deadM2w=dat$M2*dat$N.bar*dat$west, deadM2w.core=dat$M2*dat$N.bar*dat$west*dat$prop.in,deadM2=dat$M2*dat$N.bar,
                    deadZ=dat$Z*dat$N.bar)


dat<-aggregate(cbind(deadM1,deadM1w,deadM1w.core, deadM2,deadM2w,deadM2w.core, deadZ,Z)~Species.n+Species+Year+Age,data=dat,sum)

anno.M<-data.frame(Year=dat$Year,Species.n=dat$Species.n,Age=dat$Age,M2=dat$deadM2/dat$deadZ*dat$Z)
write.csv(anno.M,file=file.path(out_op,'hist_anno_M.csv'),row.names=FALSE) 

a<-Read.summary.table()
a<-subset(a,Year<=SMS.control@last.year.model)

b<-aggregate(cbind(DeadM1=deadM1w,deadM2=deadM2w, DeadM1.core=deadM1w.core,deadM2.core=deadM2w.core)~Species.n+Year,data=dat,sum)

b<-merge(a,b,by=c('Species.n','Year'))

b<-data.frame(Year=b$Year,Species.n=b$Species.n,Yield=b$Yield,yield.core=b$Yield.core,DeadM1=b$DeadM1,DeadM2=b$deadM2,DeadM1.core=b$DeadM1.core,DeadM2.core=b$deadM2.core,Fbar=b$mean.F,SSB=b$SSB,TSB=b$TSB,Recruits=b$Rec)
write.csv(b,file=file.path(out_op,'hist_condensed.csv'),row.names=FALSE)           
