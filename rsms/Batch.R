#
### init: Run: init.R in the SMS directory and _init_rsms.R
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"quarter2annual.R"))
#source(file.path(rsms.root.prog,"rsms_function.R"))
source(file.path(rsms.root.prog,"map_param.R"))
source(file.path(rsms.root.prog,"lowerUpper.R"))

### Extract data from SMS

rm(data,parameters)


doIt<-function(sms,runno){
  write.RSMS.control(sms,file=file.path(data.path,"rsms.dat"),write.multi=T,nice=T)
  inp<-make_rsms_data(dir=data.path,sms.dat='rsms.dat',adj_san0=TRUE,seasFfrom=c('F','catch')[2])

  #### prepare to run
  data<-inp[['data']]
  
  parameters<-inp[['parameters']]
  data$Debug<-1L
  my.map<-map_param(data,parameters)
  #source(file.path(rsms.root.prog,"map_param.R"))
  cat('map done\n')
  random=c("Un","Uf")
  source(file.path(rsms.root.prog,"rsms_function.R"),local=TRUE)
  obj <- MakeADFun(func,parameters, random,silent=T,map=my.map)
  cat("MakeADFun done\n")
  lu<-lowerUpper(obj,data,parameters )
  cat('lowerUpper is done\n')
  rm(opt,sdrep)
  system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1000,eval.max=1000)))
  announce(opt)
  
  sdrep <- sdreport(obj); 
  cat('Hesssian:',sdrep$pdHess,'\n')
  out<-as.list(sdrep, "Est", report=TRUE)
  sdPar<-calcCV(sdrep,data)
  plts<-summary_plot(obj,data,sdrep,out,incl_ICES_plot=FALSE,pLabel=paste(' Run:',runno)) 
  print(plts)
  fplts<-plotF(obj,sdrep,data,combineAges=FALSE) 
  print(fplts)
  a<-list(run=runno,opt=opt,Hessian=sdrep$pdHess,sdrep=summary(sdrep, "fixed", p.value = FALSE),obj=obj,sdPar=sdPar,plts=plts,fplts=fplts)
  return(a)
}

sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat")
x<- doIt(sms=sms,runno=1) # just cheking

run<-1
allruns<-list(list(run=run,sms=sms))

##
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-1
allruns<-c(allruns,list(list(run=run,sms=sms)))

##
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-2
allruns<-c(allruns,list(list(run=run,sms=sms)))

##
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-3
allruns<-c(allruns,list(list(run=run,sms=sms)))

cleanup()
results<-lapply(allruns,function(x) doIt(sms=x$sms,runno=x$run))



x<-lapply(results,function(x) cat(' run:',x$run,announce(x$opt),'Hessian:',x$Hessian))

x<-lapply(results,function(x) {
  a<-x$sdPar %>% filter(is.na(sd)) %>% mutate(run=x$run)
 if (dim(a)[[1]]<4)  print(a) else cat(dim(a)[[1]],'parameters without std\n')
})

x<-lapply(results,function(x) print(x$plts))
x<-lapply(results,function(x) print(x$fplts))

