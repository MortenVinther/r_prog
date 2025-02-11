rsms<-read.RSMS.control(dir=data.path,file='rsms.dat')
rsms@incl.process.noise[]<-pn
write.RSMS.control(rsms,file=file.path(data.path,"rsms.dat"),write.multi=T,nice=T)

doMultiExtract<-FALSE
if (TRUE) {  # transform  SMS data into RSMS format 
  switch(my.stock.dir,
         "rsms_input"       = inp_all<-make_rsms_data(dir=data.path,sms.dat='rsms.dat',seasFfrom=c('F','catch')[2],multi=doMultiExtract),
         stop(paste0("Not valid stock dir: ",my.stock.dir))
  ) #end switch
  save(inp_all,file=file.path(data.path,"rsms_input_all.Rdata")) 
} else load(file=file.path(data.path,"rsms_input_all.Rdata"),verbose=TRUE)

smsConf<-0L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated

# select a combination of species from the (full) data set, also including multi species information
runName<-paste0('Single_',pn,'_',paste(my.ps,collapse='_'))
print(runName)
my.pso<-c(0L)

inp<-pick_species(ps=my.ps,pso=my.pso, inp=inp_all,smsConf) 
inp$data$sms.mode<-smsConf

#### prepare to run
data<-inp[['data']]
parameters<-inp[['parameters']]
cleanrun(silent=TRUE)
myMap<-map_param(data,parameters)
random=c("Un","Uf")
print(system.time(obj <- MakeADFun(func, parameters, random,silent=T,map=myMap)))
lu<-lowerUpper(obj,data,parameters )
system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1000,eval.max=1000)))
print(announce(opt))
system.time(sdrep <- sdreport(obj))
print(paste('Hesssian:',sdrep$pdHess))
myRep<-obj$report()
sms<-saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
runs<-c(runs,runName)
