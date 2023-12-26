control<-read.FLSMS.control()
indices<-SMS2FLIndices(control) 

nc<-trunc.FLSMS.control(control,inclVPA=c(16,22)) 
nc 
  

write.FLSMS.control(control=nc,file="sms2.dat",path=data.path,write.multi=FALSE,nice=TRUE,writeSpNames=F,expand=F) 
  

