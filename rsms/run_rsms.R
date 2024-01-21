#
### init: Run: init.R in the SMS directory and _init_rsms.R
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"quarter2annual.R"))
source(file.path(rsms.root.prog,"rsms_function.R"))

### Extract data from SMS

if (FALSE) {  # transform  SMS data into RSMS format 
  inp_all<-make_rsms_data(dir="rsms_input",outDir=rsms.root,sms.dat='rsms.dat')
  save(inp_all,file=file.path(rsms.root,"rsms_input_all.Rdata"))
  load(file=file.path(rsms.root,"rsms_input_all.Rdata"),verbose=TRUE)
}
load(file=file.path(rsms.root,"rsms_input_all.Rdata"),verbose=TRUE)

# read a data set on SMS format
# inp<-make_rsms_data(dir="S21",outDir=rsms.root)
annualData<-F

# select a combination of species from the (full) data set
#inp<-pick_species(ps=c(1L,3L,4L,6L), inp=inp_all) # example with more species, convergence and Hessian
inp<-pick_species(ps=c(1L,2L,3L,4L,6L), inp=inp_all)
inp<-pick_species(ps=c(5L), inp=inp_all)  
#inp=inp_all

#  transform quarterly data into to annual data (testing)
if (annualData) inp<-into_annual(inp)


inp$data$spNames

#### prepare to run

data<-inp[['data']]
parameters<-inp[['parameters']]

#data$stockRecruitmentModelCode[]<-1L
data$Debug<-1L


#### Run rsms


#Going back to our objective function f, first step is to check that you can evaluate the function as a normal R function:
# func(parameters)   # KALDET VIL PÅVIRKE KALDET TIL MakeAdFun !!!?
#An error at this point is obviously not due to RTMB.

# adjust if the are species/year combination with zero catches (or assumed F is very low and highly uncertain)
if (data$zeroCatchYearExists==1) {
  UfMap<-matrix(1L:(dim(parameters$Uf)[[1]]*dim(parameters$Uf)[[2]]),nrow=sum(data$nlogF),byrow=TRUE)
  for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0 ) {
    zy<-data$zeroCatchYear[[s]]
    fromTo<-data$nlogFfromTo[s,]
    UfMap[fromTo[1]:fromTo[2],zy]<-NA
    parameters$Uf[fromTo[1]:fromTo[2],zy]<-log(0.001)
  }
  UfMap<-factor(UfMap)
}

if (data$zeroCatchYearExists==1) my.map<-list(Uf=UfMap) else my.map=list()

if (any(data$stockRecruitmentModelCode==0)) { #random walk recruitment, no need for recruitment parameters
  aMap<-1L:length(parameters$rec_loga)
  bMap<-1L:length(parameters$rec_logb)
  aMap[data$stockRecruitmentModelCode==0]<-NA
  bMap[data$stockRecruitmentModelCode==0]<-NA
  aMap<-factor(aMap)
  bMap<-factor(bMap)
  my.map<-c(my.map,list(rec_loga=aMap),list(rec_logb=bMap))
}
  #logSdLogObsSurvey=factor(rep(NA,length(parameters$logSdLogObsSurvey))),
  # logSdLogN =factor(rep(NA,length(parameters$logSdLogN)))
  #rho =factor(rep(NA,length(parameters$rho)))

#obj <- MakeADFun(func, parameters,silent=FALSE,map=my.map); obj$simulate()
#obj <- MakeADFun(func, parameters, random=c("Uf"),silent=FALSE,map=my.map); obj$simulate()
#obj <- MakeADFun(func, parameters, random=c("Un"),silent=FALSE,map=my.map); obj$simulate()  # problems ?
random=c("Un","Uf")
obj <- MakeADFun(func, parameters, random,silent=T,map=my.map)

#obj$simulate()
# checkConsistency(obj);
               
source(file.path(rsms.root.prog,"lowerUpper.R"))
#t(rbind(lower,upper)) 

opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper,control=list(iter.max=300,eval.max=500))

announce(opt)

sdrep <- sdreport(obj); 
cat('Hesssian:',sdrep$pdHess,'\n')

####  Re-run with estimated parameters
if (FALSE) {
  # using opt$par(with the estimated parameters) instead of  obj$par
  system.time(opt1 <- nlminb(opt$par, obj$fn, obj$gr, lower=lower, upper=upper,control=list(iter.max=300,eval.max=300)))
  announce(opt1)
  sdrep <- sdreport(obj); cat('Hesssian:',sdrep$pdHess,'\n')

  if (FALSE) { # another way of doing the same, requires call to MakeADFun, but faster optimization
     # no need to repeat if sdreport done already
    system.time(sdrep<-sdreport(obj,ignore.parm.uncertainty=TRUE))
    
    newPar<-as.list(sdrep, what="Est")  #parameters and random effects
    system.time(obj2 <- MakeADFun(func, newPar, random,silent=T,map=my.map))
     #obj2$simulate()
    system.time(opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, lower=lower, upper=upper,control=list(iter.max=300,eval.max=300)))
    announce(opt2)
    sdrep <- sdreport(obj2); cat('Hesssian:',sdrep$pdHess,'\n')
    data.frame(name=names(obj$par),ini=obj$par,opt=opt$par,exp_opt=exp(opt$par),opt1=opt1$par,opt2=opt2$par)
    
  }
}  

