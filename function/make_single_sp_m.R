#  this script calculates the mean M (=M1+M2) over all years and write it to a file
m1m2_calcs<-function(calcDone=c("averageByAge","averageBySpeciesQuarter","averageBySpecies","simpleSum","smoothed")[1]) {
  stopifnot(calcDone %in% c("averageByAge","averageBySpeciesQuarter","averageBySpecies","simpleSum","smoothed"))
  dat<-Read.summary.data(read.init.function=T)
  oth<-Read.other.predator()
  oth$M<- -77
  oth<-subset(oth,select=c(Species.n,Year,Quarter,Age,M))
  
  dat$M<-dat$M1+dat$M2
  dat<-subset(dat,select=c(Species.n,Year,Quarter,Age,M))
  dat<-rbind(dat,oth)
  
  #dat<-droplevels(subset(dat,Species.n>=first.VPA))
  #maxAge<-SMS.control@max.age.all
  
  if (calcDone=="averageByAge") {  # average M over all years
    M<-tapply(dat$M,list(dat$Quarter,dat$Age,dat$Species.n),mean)
    
    M[is.na(M)]<- -1
     outf<-file.path(data.path,'natmor_avg.in')
    cat('# single species M from average multispecies M1+M2 over the full period\n',file=outf)
    
    for (sp in (first.VPA:nsp)) {
      for (y in (SMS.control@first.year:SMS.control@last.year )) {
      
      cat(paste('# ',sp.names[sp],y,'\n'),file=outf,append=T)
      write.table(round(M[,,as.character(sp)],3),append=T,file=outf,row.names = F,col.names = F)
    }
    }
  } 
  if (calcDone=="averageBySpeciesQuarter") { # average M 
    nage<-SMS.control@max.age.all-SMS.control@first.age+1
    
    M<-tapply(dat$M,list(dat$Quarter,dat$Species.n),mean)

    outf<-file.path(data.path,'natmor_by_species_quarter.in')
    cat('# single species M from average multispecies M1+M2 by species and quarter\n',file=outf)
    
    for (sp in (first.VPA:nsp)) {
      MM<-M[,as.character(sp)]
      MM<-matrix(rep(MM,nage),SMS.control@last.season)
      for (y in (SMS.control@first.year:SMS.control@last.year )) {
        
        cat(paste('# ',sp.names[sp],y,'\n'),file=outf,append=T)
        write.table(round(MM,3),append=T,file=outf,row.names = F,col.names = F)
      }}    
  }
  
  if (calcDone=="averageBySpecies") { # average M 
    nage<-SMS.control@max.age.all-SMS.control@first.age+1
    
    M<-tapply(dat$M,list(dat$Species.n),mean)
     
    outf<-file.path(data.path,'natmor_by_species.in')
    cat('# single species M from average multispecies M1+M2 by species \n',file=outf)
    
    for (sp in (first.VPA:nsp)) {
      MM<-M[as.character(sp)]
      MM<- matrix(MM,ncol=nage,nrow=SMS.control@last.season)
      
      for (y in (SMS.control@first.year:SMS.control@last.year )) {
        
        cat(paste('# ',sp.names[sp],y,'\n'),file=outf,append=T)
        write.table(round(MM,3),append=T,file=outf,row.names = F,col.names = F)
      }}    
  }
    
  if (calcDone=="simpleSum") { # Sum of M1 and M2 
    M<-tapply(dat$M,list(dat$Year,dat$Quarter,dat$Age,dat$Species.n),sum)
    M[is.na(M)]<- -1
    
    outf<-file.path(data.path,'natmor_sum.in')
    cat('# single species from sum of  multispecies M1+M2\n',file=outf)
    
    for (sp in (first.VPA:nsp)) {
      for (y in (SMS.control@first.year:SMS.control@last.year )) {
        
        cat(paste('# ',sp.names[sp],y,'\n'),file=outf,append=T)
        write.table(round(M[as.character(y),,,as.character(sp)],3),append=T,file=outf,row.names = F,col.names = F)
      }}    
   }


  if (calcDone=="smoothed") { # smoothed M 
    library(mgcv)
    outf<-file.path(data.path,'natmor_smoothed.in')
    cat('# single species M from year-smoothed multispecies M1+M2 by species and quarter\n',file=outf)
    
    for (sp in (first.VPA:nsp)) {

      b<-subset(dat,Species.n==sp)
      bb<- by(b,list(b$Species.n,b$Quarter,b$Age),function(x) {
        a<-gam(M ~s(Year,bs='cs'),data=x)
        return(data.frame(Species.n=x$Species.n,Quarter=x$Quarter,Year=x$Year,Age=x$Age,M=predict(a)))
      })
      
      dat2<-do.call(rbind,bb)
      dat2<-rbind(dat2,filter(dat,Species.n!=sp))  # to get the full age range
      
      M<-tapply(dat2$M,list(dat2$Year,dat2$Quarter,dat2$Age,dat2$Species.n),mean)
      M[is.na(M)]<- -1
      for (y in (SMS.control@first.year:SMS.control@last.year )) {
        cat(paste('# ',sp.names[sp],y,'\n'),file=outf,append=T)
        write.table(round(M[as.character(y),,,as.character(sp)],3),append=T,file=outf,row.names = F,col.names = F)
      }}    
  }
  
 
  cat(-999," # checksum",file=outf,append=T)
  cat("file: ",outf,' has been written\n')
}
