#
### init: Run: init.R in the SMS directory and _init_rsms.R
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"quarter2annual.R"))
source(file.path(rsms.root.prog,"rsms_function.R"))

### Extract data from SMS

if (FALSE) {  # transform  SMS data into RSMS format 
  if (my.stock.dir=="rsms_input") inp_all<-make_rsms_data(dir=my.stock.dir,outDir=rsms.root,sms.dat='rsms.dat',adj_san0=TRUE,seasFfrom=c('F','catch')[1],
                          effort_fl=c("NSA Commercial 1983-1998","NSA Commercial 1999-2022", "SSA Commercial 1983-2002","SSA Commercial 2003-2022") )
                      #    effort_fl="none") 
  
  if (my.stock.dir=="rsms_SAN-area-1r") inp_all<-make_rsms_data(dir=my.stock.dir,outDir=rsms.root,sms.dat='rsms.dat',adj_san0=FALSE,
                        seasFfrom=c('F','catch')[2], effort_fl=c("Effort season 1","Effort season 2") )
  
  save(inp_all,file=file.path(rsms.root,"rsms_input_all.Rdata"))
}

load(file=file.path(rsms.root,"rsms_input_all.Rdata"),verbose=TRUE)

annualData<-FALSE

# select a combination of species from the (full) data set
#inp<-pick_species(ps=c(1L,3L,4L,6L), inp=inp_all) # example with more species, convergence and Hessian
#inp<-pick_species(ps=c(1L,2L,3L,4L,5L,6L,7L,9L), inp=inp_all) # 8,10 no Hessian
inp<-pick_species(ps=c(7L), inp=inp_all)  
#inp=inp_all

#  transform quarterly data into to annual data (testing)
if (annualData) inp<-into_annual(inp)


inp$data$spNames

#### prepare to run

data<-inp[['data']]
parameters<-inp[['parameters']]

#data$stockRecruitmentModelCode[]<-1L
data$Debug<-1L


#data$seasFprop[[1]]

#### Run rsms


#Going back to our objective function f, first step is to check that you can evaluate the function as a normal R function:
# func(parameters)   # KALDET VIL PÅVIRKE KALDET TIL MakeAdFun !!!?
#An error at this point is obviously not due to RTMB.


  #logSdLogObsSurvey=factor(rep(NA,length(parameters$logSdLogObsSurvey))),
  # logSdLogN =factor(rep(NA,length(parameters$logSdLogN)))
  #rho =factor(rep(NA,length(parameters$rho)))

#obj <- MakeADFun(func, parameters,silent=FALSE,map=my.map); obj$simulate()
#obj <- MakeADFun(func, parameters, random=c("Uf"),silent=FALSE,map=my.map); obj$simulate()
#obj <- MakeADFun(func, parameters, random=c("Un"),silent=FALSE,map=my.map); obj$simulate()  # problems ?

source(file.path(rsms.root.prog,"map_param.R"))
random=c("Un","Uf")

obj <- MakeADFun(func, parameters, random,silent=T,map=my.map)

#obj$simulate()
# checkConsistency(obj);
               
source(file.path(rsms.root.prog,"lowerUpper.R"))
#t(rbind(lower,upper)) 

system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper,control=list(iter.max=1000,eval.max=1000)))
announce(opt)

sdrep <- sdreport(obj); 
cat('Hesssian:',sdrep$pdHess,'\n')

sdrep


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

