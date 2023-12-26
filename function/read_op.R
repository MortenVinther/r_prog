
#Function to read rep file
Read.OP.rep<-function(dir=data.path,infile='op.rep')
{
  file<-file.path(dir,infile)
  s<-readLines(file)
  s<-suppressWarnings(as.numeric(unlist(strsplit(s,' '))))
  s<-s[!is.na(s)]
  s<-list(objFunc=s[1],npar=s[2],maxGrad=s[3],Penalty=s[4],Ftarget=s[5:(5+first.VPA-1)],Penalty.sp=tail(s,first.VPA-1))
  return(s)
}


#Function to OP read size data
Read.OP.size<-function(dir=data.path,infile='op_size.in')
{
  
  
  n.VPA<-sum(SMS.control@species.info[,'predator'] %in% c(0,1))
  n.other.pred<-sum(SMS.control@species.info[,'predator'] %in% c(2))
  n.pred<-sum(SMS.control@species.info[,'predator'] %in% c(1,1)) 
  n.age<-SMS.control@max.age.all-SMS.control@first.age+1
  n.sp<- SMS.control@no.species
  
  OP<-read.FLOP.control(file="op.dat",path=data.path,
                        n.VPA,
                        n.other.pred,
                        n.pred)
  
  file<-file.path(dir,infile)
  s<-scan(file,comment.char = "#",quiet=TRUE,nmax=n.sp*n.age*SMS.control@last.season)
  wpred<-array(data=s,dim=c(n.age,  n.sp,SMS.control@last.season))
  
  dimnames(wpred)<-list(Age=SMS.control@first.age:SMS.control@max.age.all,
                        Species=sp.names,
                        Quarter=1:SMS.control@last.season)
  wpred<-arr2dfny(wpred,"size")
  wpred$type<-"predator W"
  
  
  ###
  s<-scan(file,comment.char = "#",quiet=TRUE)
  s<-s[(n.sp*n.age*SMS.control@last.season+1):length(s)]
  s<-head(s,nmax=n.VPA*n.age*SMS.control@last.season)
  wprey<-array(s,dim=c(n.age,n.VPA,SMS.control@last.season))
  
  dim(wprey)
  dimnames(wprey)<-list(Age=SMS.control@first.age:SMS.control@max.age.all,
                        Species=sp.names[ (n.sp-n.VPA+1):nsp],
                        Quarter=1:SMS.control@last.season)
  wprey<-arr2dfny(wprey,"size")
  wprey$type<-"prey W"
  
  
  ###
  s<-scan(file,comment.char = "#" ,quiet=TRUE)
  s<-tail(s,n.age*n.sp*SMS.control@last.season)
  len<-array(data=s,dim=c(n.age,  n.sp,SMS.control@last.season))
  dimnames(len)<-list(Age=SMS.control@first.age:SMS.control@max.age.all,
                      Species=sp.names,
                      Quarter=1:SMS.control@last.season)
  len<-arr2dfny(len,"size")
  len$type<-"length"
  
  
  out<-rbind(rbind(wpred,wprey),len) %>%  mutate_if(is.factor, as.character)
  return(out)

} 


#Function to read par file
Read.OP.par<-function(dir=data.path,infile='op_optim_par.out')
{
  file<-file.path(dir,infile)
  s<-read.table(file,header=TRUE)
  data.frame(Species=sp.names[s$Species.n],s)
}

#Function to read summary data
Read.OP.condensed<-function(dir=data.path,infile='op_condensed.out')
{
  file<-file.path(dir,infile)
  s<-read.table(file,header=TRUE)
  data.frame(Species=sp.names[s$Species.n],s)
}

#Function to read variable size and growth other pred
Read.op.other.sp.var<-function(dir=data.path,infile='op_other_sp_var.out')
{
  file<-file.path(dir,infile)
  s<-read.table(file,header=TRUE)
  data.frame(Predator=sp.names[s$Species.n],s)
}


#Function to read summary data
Read.OP.community.indicator<-function(dir=data.path,infile='op_indicator_system.out')
{
  file<-file.path(dir,infile)
  read.table(file,header=TRUE)
}

Read.OP.other<-function(dir=data.path,infile="op_other_sp.out")
{
  file<-file.path(dir,infile)
  s<-read.table(file,header=TRUE)
  data.frame(Predator=sp.names[s$Species.n],s)
}

#Read refence points from file reference_points.in 

Read.reference.points.OP<-function(dir=data.path){
    a<-scan(file.path(dir,"op_reference_points.in"),comment.char = "#",quiet = TRUE)
    b<-matrix(a,ncol=4,byrow=TRUE)
    colnames(b)<-c("Flim","Fpa","Blim","Bpa")   
    rownames(b)<-sp.names[first.VPA:nsp]
    b
}
