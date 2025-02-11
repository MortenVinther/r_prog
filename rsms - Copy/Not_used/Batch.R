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

printPlot=T
doIt<-function(sms,runno){
  write.RSMS.control(sms,file=file.path(data.path,"rsms.dat"),write.multi=T,nice=T)
  inp<-make_rsms_data(dir=data.path,sms.dat='rsms.dat',seasFfrom=c('F','catch')[2])

  if (exists('opt')) rm(opt);if (exists('sdrep')) rm(sdrep);if (exists('obj')) rm(obj)

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
  system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1000,eval.max=1000)))
  announce(opt)
  
  sdrep <- sdreport(obj); 
  cat('Hesssian:',sdrep$pdHess,'\n')
  out<-as.list(sdrep, "Est", report=TRUE)
  sdPar<-calcCV(sdrep,data)
  plts<-summary_plot(obj,data,sdrep,out,incl_ICES_plot=FALSE,pLabel=paste(' Run:',runno)) 
  if (printPlot) print(plts)
  fplts<-plotF(obj,sdrep,data,combineAges=T,pLabel=paste(' Run:',runno)) 
  if (printPlot) print(fplts)
  a<-list(run=runno,opt=opt,Hessian=sdrep$pdHess,sdrep=sdrep,sdreps=summary(sdrep, "fixed", p.value = FALSE),ssb=summary(sdrep, "report"),obj=obj,sdPar=sdPar,plts=plts,fplts=fplts)
  return(a)
}

sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat")
x<- doIt(sms=sms,runno=1) # just cheking

run<-1
allruns<-list(list(run=run,sms=sms))

if (FALSE) {
##
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-1
allruns<-c(allruns,list(list(run=run,sms=sms)))

##
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-2
allruns<-c(allruns,list(list(run=run,sms=sms)))
}

sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-3
allruns<-c(allruns,list(list(run=run,sms=sms)))


## no F at age zero, one F group 
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-1
sms@species.info[,2]<-1L
sms@keyLogFsta[[1]]<-c(1L)
sms@catch.s2.group[[1]]<-c(1L,3L)
allruns<-c(allruns,list(list(run=run,sms=sms)))


## no F at age zero, two F groups 
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-1
sms@species.info[,2]<-1L
sms@keyLogFsta[[1]]<-c(1L,2L)
sms@catch.s2.group[[1]]<-c(1L,3L)
allruns<-c(allruns,list(list(run=run,sms=sms)))


## no F at age zero, two F groups 
sms<-read.RSMS.control(dir=data.path,file="rsms_source.dat"); run<-run+1
sms@SSB.R<-1
sms@species.info[,2]<-1L
sms@keyLogFsta[[1]]<-c(1L,3L)
sms@catch.s2.group[[1]]<-c(1L,3L)
allruns<-c(allruns,list(list(run=run,sms=sms)))


# x<- doIt(sms=sms,runno=1) # just checking


cleanup()
results<-lapply(allruns,function(x) doIt(sms=x$sms,runno=x$run))



x<-lapply(results,function(x) cat(' run:',x$run,announce(x$opt),'Hessian:',x$Hessian))

x<-lapply(results,function(x) {
  a<-x$sdPar %>% filter(is.na(sd)) %>% mutate(run=x$run)
 if (dim(a)[[1]]<4)  print(a) else cat(dim(a)[[1]],'parameters without std\n')
})

x<-lapply(results,function(x) aa<-summary(x$sdrep, "fixed", p.value = FALSE) )

# CV on SSB
yssb<-5
x<-lapply(results,function(x) {
  a<-tail(as.data.frame(x$ssb),yssb) ;
  ssbCV<-a[,2]/a[,1]
})
xx<-do.call(rbind,x);rownames(xx)<-paste('run',1:run); colnames(xx)<- paste((sms@last.year.model-yssb+1):sms@last.year.model  )
xx

x<-lapply(results,function(x) print(x$plts))
x<-lapply(results,function(x) print(x$fplts))

