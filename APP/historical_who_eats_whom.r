
dat<-Read.summary.data(read.init.function=F)
dat<-subset(dat,Year<=SMS.control@last.year.model )

dat<-data.frame(dat, deadM2=ifelse(dat$M2>=0,dat$M2*dat$N.bar*dat$west*dat$prop.in,0), deadM1=ifelse(dat$M1>=0,dat$M1*dat$N.bar*dat$west*dat$prop.in,0))
dat<-subset(dat,select=c(Species, Year, Quarter, Species.n, Age, M2,deadM1,deadM2))

M2<-Read.part.M2.data()

a<-merge(x=dat,y=M2, by.x = c("Year","Quarter","Species","Age"), by.y = c("Year","Quarter","Prey","Prey.age"))
a$eatenW<- a$deadM2*a$Part.M2/a$M2

a$Prey<-a$Species
a$Prey.age<-a$Age
a$tot.M2.prey<-a$M2

b<-subset(a,select=c( Year, Quarter, Predator,Predator.age, Prey, Prey.age,Prey.no, eatenW, Part.M2,tot.M2.prey))
#bb<-droplevels(aggregate(list(eatenW=b$eatenW),list(Year=b$Year, Quarter=b$Quarter, Predator=b$Predator,Prey=b$Prey),sum))
bbb<-droplevels(aggregate(list(eatenW=b$eatenW),list(Year=b$Year, Predator=b$Predator,Prey=b$Prey,Prey.no=b$Prey.no),sum))

pred_format<-read.csv(file.path(data.path,'pred_format.csv'),header=TRUE)

s<-merge(x=bbb,y=pred_format,by.x='Prey',by.y='old',all.x=TRUE)
s$Prey<-s$new; s$new<-NULL
s$Prey.no<-s$new_no; s$new_no<-NULL

s<-merge(x=s,y=pred_format,by.x='Predator',by.y='old',all.x=TRUE)
a<-aggregate(s$eatenW,list(s$new,s$Year,s$new_no,s$Prey,s$Prey.no),sum)
names(a)<-c("Predator","Year","Predator.no","Prey","Prey.no","eatenW")
###

sort(unique(a$Prey))

dat1<-Read.summary.table(read.init.function=F)
dat1<-subset(dat1,Year<=SMS.control@last.year.model )
human<-data.frame(Predator='Humans',Year=dat1$Year,Predator.no=0,Prey=dat1$Species,Prey.no=dat1$Species.n,eatenW=dat1$SOP.core)
human<-merge(x=human,y=pred_format,by.x='Prey',by.y='old',all.x=TRUE)
human$Prey<-human$new; 
human$Prey.no<-human$new_no; 
human<-subset(human,select=c(Predator, Year, Predator.no, Prey, Prey.no,   eatenW))
head(human)
a<-bind_rows(a,human)

sort(unique(a$Prey))

#####


dat<-subset(dat,deadM1>=0)
b<-aggregate(dat$deadM1,list(dat$Species,dat$Species.n,dat$Year),sum)
names(b)<-c("Species","Species.no","Year","deadM1")
b<-merge(x=b,y=pred_format,by.x='Species',by.y='old',all.x=TRUE)
b$Prey<-b$new
b$Prey.no<-b$new_no; 
b<-data.frame(Predator='Residual mortality',Year=b$Year,Predator.no=-1,Prey=b$Prey,Prey.no=b$Prey.no,eatenW=b$deadM1)
sort(unique(b$Prey))

a<-bind_rows(a,b)


pp<-unique(subset(pred_format,select=c(new,new_no)))
pp<-pp[order(pp$new_no),]
pformat<-pp$new

  
write.table(a,file=file.path(out_op,'who_eats_whom_historical.csv'),sep=',',row.names = F)
