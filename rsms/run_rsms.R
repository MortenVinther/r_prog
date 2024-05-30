#
### init: Run: init.R in the SMS directory and _init_rsms.R
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"quarter2annual.R"))
source(file.path(rsms.root.prog,"rsms_function.R"))
source(file.path(rsms.root.prog,"summary_plot.R"))
source(file.path(rsms.root.prog,"map_param.R"))
source(file.path(rsms.root.prog,"lowerUpper.R"))

### Extract data from SMS

if (TRUE) {  # transform  SMS data into RSMS format 
  if (my.stock.dir=="rsms_input") inp_all<-make_rsms_data(dir=data.path,sms.dat='rsms.dat',seasFfrom=c('F','catch')[2])
  
  # 
  if (my.stock.dir=="rsms_SAN-area-1r")  inp_all<-make_rsms_data(dir=data.path,sms.dat='rsms.dat',seasFfrom=c('F','catch')[2])
  
  # 
  # if (my.stock.dir=="rsms_SAN-area-3r") inp_all<-make_rsms_data(dir=my.stock.dir,sms.dat='rsms.dat',
  #                                                               seasFfrom=c('F','catch')[2], effort_fl=c("Effort season 1","Effort season 2") )
  # 
  # 
  # if (my.stock.dir=="rsms_Sprat-div-4_plus_IIIa") inp_all<-make_rsms_data(dir=my.stock.dir,sms.dat='rsms.dat',adj_san0=FALSE,
  #                                                               seasFfrom=c('F','catch')[2] )
  # 
 
  save(inp_all,file=file.path(data.path,"rsms_input_all.Rdata"))
}

load(file=file.path(data.path,"rsms_input_all.Rdata"),verbose=TRUE)

annualData<-FALSE

# select a combination of species from the (full) data set
#inp<-pick_species(ps=c(1L,3L,4L,6L), inp=inp_all) # example with more species, convergence and Hessian
#inp<-pick_species(ps=c(1L,2L,3L,4L,5L,6L,7L,9L), inp=inp_all) # 8,10 no Hessian

#inp<-pick_species(ps=c(1L,2L,3L,4L,5L,6L), inp=inp_all)  #ok

#inp<-pick_species(ps=c(1L,6L,7L,8L), inp=inp_all)  
#inp<-pick_species(ps=c(7L,8L), inp=inp_all) 
inp=inp_all

#  transform quarterly data into to annual data (testing)
if (annualData) inp<-into_annual(inp)

inp$data$spNames
inp$data$fleetNames

#### prepare to run
data<-inp[['data']]
parameters<-inp[['parameters']]

data$Debug<-1L

#### Run rsms
my.map<-map_param(data,parameters)
random=c("Un","Uf")

obj <- MakeADFun(func, parameters, random,silent=T,map=my.map)

#obj$simulate()
# checkConsistency(obj);
               
lu<-lowerUpper(obj,data,parameters )

# round(ftable(data$seasFprop[[2]][,2,]),3)

# 
# CC<-as.data.frame(data$keyCatch)
# CC$canum<-exp(data$logCatchObs)
# CC$canumR<-round(exp(data$logCatchObs))
# round(ftable(xtabs(canum~s+year+age,data=CC)))
# round(ftable(tapply(CC$canum,list(CC$s,CC$year,CC$age),FUN=sum)))
# filter(CC,year==2012)
# round(ftable(data$seasFprop[[1]][45,,]),5)
# arrange(CC,canum)


rm(opt,sdrep)
system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1000,eval.max=1000)))
announce(opt)

sdrep <- sdreport(obj); 
cat('Hesssian:',sdrep$pdHess,'\n')
sdrep

FF<-FToDF(obj,sdrep,data) %>% mutate(species=data$spNames[s])
by(FF,FF$species,function(x) round(ftable(xtabs(FF~year+age,data=x)),6))

sF<-seasonalF(obj,sdrep,data)
by(sF,sF$species,function(x) round(ftable(xtabs(SF~year+q+age,data=x)),6))

filter(sF,year >2021 & age <2)
data.frame(name=attr(sdrep$par.fixed,'names'),value=sdrep$par.fixed,gradient=sdrep$gradient.fixed)

print(calcCV(sdrep,data) %>%arrange(sp.no,desc(sd)),n=40)
print(calcCV(sdrep,data) %>%arrange(desc(abs(gradient))),n=20)

attr(sdrep$par.fixed,'names')
####  Re-run with estimated parameters
if (FALSE) {
  # using opt$par(with the estimated parameters) instead of  obj$par
  system.time(opt1 <- nlminb(opt$par, obj$fn, obj$gr, lower=lower, upper=upper,control=list(iter.max=2000,eval.max=2000)))
  announce(opt1)
  sdrep <- sdreport(obj); cat('Hesssian:',sdrep$pdHess,'\n')

  if (FALSE) { # another way of doing the same, requires call to MakeADFun, but faster optimization
     # no need to repeat if sdreport done already
    system.time(sdrep<-sdreport(obj,ignore.parm.uncertainty=TRUE))
    
    newPar<-as.list(sdrep, what="Est")  #parameters and random effects
    system.time(obj2 <- MakeADFun(func, newPar, random,silent=T,map=my.map))
     #obj2$simulate()
    system.time(opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, lower=lower, upper=upper,control=list(iter.max=1300,eval.max=1300)))
    announce(opt2)
    sdrep <- sdreport(obj2); cat('Hesssian:',sdrep$pdHess,'\n')
    data.frame(name=names(obj$par),ini=obj$par,opt=opt$par,exp_opt=exp(opt$par),opt1=opt1$par,opt2=opt2$par)
  }
}  

