# user options

#Legend position  Valid options are "bottomright", "bottom", "bottomleft", 
#"left", "topleft", "top", "topright", "right" and "center". Use NULL for no legend
legend.position <-  "topleft"
legend.ncol     <-  2

#weigting factor for survey (or other obs) observations, start with default value to get a solid line for default
w.fac<-c(1.0,3,2,0.33,0.25)

 #w.fac<-c(1991,1993,1995,1997,1999,2001,2002)  # used for Blue whiting
dev<-'screen';portrait <- TRUE
nox<-1; noy<-1;
cleanup()

######   end user options  ###################



setwd(data.path)
bio.interact<-FALSE

# retrospective runs are made in a separate dirictory
retro.dir<-"retro"

dir.create(retro.dir,showWarnings = FALSE)

SMS.files.single<-c("natmor.in","canum.in","west.in","weca.in","propmat.in","fleet_catch.in","fleet_names.in","fleet_info.dat","just_one.in")

for (from.file in SMS.files.single) {
  to.file<-file.path(data.path,retro.dir,from.file)
  file.copy(from.file, to.file, overwrite = TRUE)
}

# read data and options into FLR objects
control<-read.FLSMS.control()

control@read.HCR<-0

setwd(retro.dir)

for (fac in (w.fac)){
    control@obj.func.weight[1,"survey"]<-fac
   # control@obj.func.weight[1,"catch"]<-fac

  # control@catch.sep.year[[1]][2]<-fac


   if (F) {
     fl.dat<-fleet.info<-scan('fleet_info.dat',comment.char='#')
     fl.dat[1]<-fac
     write(fl.dat, file = "fleet_info.dat")
   }

    
    write.FLSMS.control(control,write.multi=bio.interact)

    shell(paste( file.path(root,"program","sms.exe"),"-nox -nohess",sep=" "), invisible = TRUE)
    if (!file.copy("summary.out", paste("summary",fac,".out",sep=""), overwrite = TRUE)) stop("Retro stopped: something went wrong in copying summary.dat")
    file.remove("summary.out")

    if (!file.copy("sms.par", paste("sms",fac,".par",sep=""), overwrite = TRUE)) stop("Retro stopped: something went wrong in copying sms.par")
    file.remove("sms.par")

    if (!file.copy("objective_function.out", paste("objective_function",fac,".out",sep=""), overwrite = TRUE)) stop("Retro stopped: something went wrong in copying objective_function.out")
    file.remove("objective_function.out")

}

av.F.age<-control@avg.F.ages 
  
lab<-''
i<-1
for (fac in(w.fac)) {
  file<-paste('summary',fac,'.out',sep='')
  s<-read.table(file,header=TRUE)
  s<-data.frame(Species=control@species.names[s$Species.n],s)
  rec<-subset(s,Quarter==control@rec.season & Age==control@first.age & Z>-1)
  rec<-tapply(rec$N,list(rec$Species.n,rec$Year),sum)/1000

  s1<-subset(s,Quarter==1)
  ssb<-tapply(s1$SSB,list(s1$Species.n,s1$Year),sum)/1000

  if (i==1) {
   SSB<-array(NA,dim=c(nsp,dim(ssb)[2],length(w.fac)),
         dimnames=list(control@species.names,seq(control@first.year,control@last.year.model+1),(w.fac)))
   FI<-array(NA,dim=c(nsp,dim(rec)[2],length(w.fac)),
         dimnames=list(control@species.names,seq(control@first.year,control@last.year.model),(w.fac)))
   REC<-FI
  }
  SSB[,1:dim(ssb)[2],i]<-ssb
  REC[,1:dim(rec)[2],i]<-rec
  
  s<-subset(s,Z>-1)
  F.at.age<-tapply(s$F,list(s$Year,s$Age,s$Species.n),sum,na.rm=T)
  for (sp in (1:nsp)){
    FI[,1:(dim(rec)[2]),i]<-apply(F.at.age[ ,(av.F.age[1]-fa+1):(av.F.age[2]-fa+1),sp],c(1),mean)  ### kun med en art-- fejl
  }
 i<-i+1
}  


#Plot results using the plot compare function
noxy<-nox*noy
newplot(dev,nox,noy,portrait)

#source(file.path(prog.path,"plot_triple.R"))
plot.triple(legend.position,legend.ncol)


if (F) {
  i<-1
  for (fac in(w.fac)) {
    file<-paste('sms',fac,'.par',sep='')
    s<-scan(file,comment.char = "#" )
    if (i==1) sms.par<-s else sms.par<-cbind(sms.par,s)
   i<-i+1
  }

  npar<- dim(sms.par)[1]
  first.par<-npar-18-7
  last.par<-npar-18
  bias.fac<-(sms.par[first.par:last.par,])
  matplot(bias.fac,type='b')
  #abline(h=1)
}


if (T) {
  i<-1
  for (fac in(w.fac)) {
    file<-paste('objective_function',fac,'.out',sep='')
     s<-read.table(file,header=TRUE)
     s<-data.frame(factor=fac,s)
    if (i==1) obj<-s else obj<-rbind(obj,s)
   i<-i+1
  }


}

setwd(data.path)

