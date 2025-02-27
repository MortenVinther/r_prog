write.RSMS.control(sms,file=file.path(data.path,"rsms.dat"),write.multi=T,nice=T)
inp<-make_rsms_data(dir=data.path,sms.dat='rsms.dat',adj_san0=FALSE,seasFfrom=c('F','catch')[2])

#### prepare to run
source(file.path(rsms.root.prog,"rsms_function.R"))
data<-inp[['data']]

parameters<-inp[['parameters']]
data$Debug<-1L
my.map<-map_param(data,parameters)
#source(file.path(rsms.root.prog,"map_param.R"))
cat('map done\n')
random=c("Un","Uf")
obj <- MakeADFun(func,parameters, random,silent=T,map=my.map)
cat("MakeADFun done\n")
source(file.path(rsms.root.prog,"lowerUpper.R"))

rm(opt,sdrep)
system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper,control=list(iter.max=1000,eval.max=1000)))
announce(opt)

sdrep <- sdreport(obj); 
cat('Hesssian:',sdrep$pdHess,'\n')
#sdrep

plts<-summary_plot(obj,data,sdrep,incl_ICES_plot=FALSE) 
a<-list(run=runno,opt=opt,Hessian=sdrep$pdHess,sdrep=summary(sdrep, "fixed", p.value = FALSE),obj=obj,plts=plts)
