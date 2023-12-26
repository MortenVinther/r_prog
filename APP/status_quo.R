a<-subset(Read.summary.table(),Year==SMS.control@last.year.model,select=c(-SOP.hat, -Yield.hat))

#sum of quarterly F
b<-subset(Read.summary.data(),Year==SMS.control@last.year.model,select=c(Species,Age,F))

ff<-as.data.frame(SMS.control@avg.F.ages)
ff$Species<-rownames(ff)

b<-merge(b,ff)

colnames(b)<-c("Species","Age","sum.q.F","first","last")
b<-subset(b,Age>=first & Age<=last)

b<-aggregate(sum.q.F~Species+Age,data=b,sum)
b<-aggregate(sum.q.F~Species,data=b,mean)

a<-merge(a,b)

a<-a[order(a$Species.n),]
a

write.csv(a,file=file.path(out_op,'status_quo.csv'),row.names=FALSE)
a

ff<-a$mean.F
cat(c(1,'\n',ff,'\n'),file=file.path(out_op,"op_multargetf.in"))

