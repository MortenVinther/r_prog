extractDataAndRunSingle<-function(runName='Single',
                            my.ps=1,my.pso=0,
                            rsmsControl='rsms.dat',
                            doMultiExtract=FALSE,
                            dir=data.path,
                            silent=TRUE,
                            fleetIn="new_fleet_info.dat",
                            smsConf=0L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
) 
{
  
  inp_all<-make_rsms_data(dir=dir,rsmsControl,fleetIn=fleetIn,multi=doMultiExtract) 
  inp<-pick_species(ps=my.ps,pso=my.pso, inp=inp_all,smsConf) 
  
  inp$data$sms.mode<-smsConf
  
  #### prepare to run
  data<-inp[['data']]
  data$silent<-FALSE
  parameters<-inp[['parameters']]
  
  myMap<-map_param(data,parameters)
  random=c("Un","Uf")
  environment(func) <- environment() # see https://github.com/kaskr/RTMB/issues/26
  system.time(obj <- MakeADFun(func, parameters, random,silent=silent,map=myMap))
  
  lu<-lowerUpper(obj,data,parameters )
  opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1000,eval.max=1000))
  announce(opt)
  
  sdrep <- sdreport(obj)
  cat('Hesssian:',sdrep$pdHess,'\n')
  
  myRep<-obj$report()
  saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
  
  return(TRUE)
}

extractDataAndRunMulti<-function(runName='Multi',
                            my.ps=1,my.pso=0,
                            rsmsControl,
                            doMultiExtract=TRUE,
                            dir=data.path,
                            silent=TRUE,
                            fleetIn="new_fleet_info.dat",
                            smsConf=1L, # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
                            parametersFrom="Single",
                            doLock=TRUE,
                            maxIter=3000
) 
{
  load(file=file.path(data.path,paste0(singleRun,".Rdata")),verbose=TRUE)
  cleanrun(silent)
  sdrep<-sms$sdrep # from single species run
  cat('Hesssian from Single run:',sdrep$pdHess,'\n')
  newPar<-as.list(sms$sdrep, what="Est")  #parameters and random effects
  
  inp_all<-make_rsms_data(dir,rsmsControl=rsmsControl,fleetIn=fleetIn,multi=doMultiExtract) 
  inp<-pick_species(ps=my.ps,pso=my.pso, inp=inp_all,smsConf) 
  
  inp$data$sms.mode<-smsConf
  
  #### prepare to run
  data<-inp[['data']]
  data$silent<-FALSE
  
  # add MS parameters
  newPar$vulnera<-inp$parameters$vulnera
  newPar$overlapP<-inp$parameters$overlapP
  newPar$logStomObsVar<-inp$parameters$logStomObsVar
  parameters<-newPar
  
  myMap<-map_param(data,parameters)
  
  #lock single species parameters,if necessary
  if (doLock) {
    lockP<-c('logCatchability','logSdLogN','logSdLogFsta','logSdLogObsSurvey','logSdLogObsCatch','rho','rec_loga','rec_logb','logSeparF')  #,'Uf','Un'
    myMap<-lock_param(data,parameters,myMap,lockP) 
  }

  random=c("Un","Uf")
  environment(func) <- environment() # see https://github.com/kaskr/RTMB/issues/26
  system.time(obj <- MakeADFun(func, parameters, random,silent=silent,map=myMap))
  
  lu<-lowerUpper(obj,data,parameters )
  opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=maxIter,eval.max=maxIter))
  announce(opt)
  
  sdrep <- sdreport(obj)
  cat('Hesssian:',sdrep$pdHess,'\n')
  
  myRep<-obj$report()
  saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
  
  return(TRUE)
}

# make default rsms.dat (control) file
batch_default_configuration <-function(outfile='rsms.dat',dir=data.path,writeConrol=TRUE) {
  
  a<-RSMS.control(first.year=1974,last.year=2022,no.species=27,last.season=4,max.age.all=10,
                  no.other.predators=15,no.VPA.predators=5,
                  species.names=c("FUL", "GLT", "HEG", "KTW", "GBG", "GNT", "PUF", "RAZ", "RAJ", "GUR", "WHM", "NHM", "GSE", "HBP", "HKE", "COD", "WHG", "HAD", "POK", "MAC", "HER", "NSA", "SSA", "NOP", "SPR", "PLE", "SOL"),
                  species.names.long=c("Fulmar","Guillemot","Her.Gull","Kittiwake","GBB.Gull","Gannet","Puffin","Razorbill","A.radiata","G.gurnards","W.horse.mac","N.horse.mac","Grey.seal","H.porpoise","Hake","Cod","Whiting","Haddock","Saithe","Mackerel","Herring","N.sandeel","S.sandeel","Nor.pout","Sprat","Plaice","Sole")
     )
  
  nonPrey<-c('POK','MAC','PLE','SOL')
  a@species.info[nonPrey,'prey']<-0L
  
  firstAgeF1<-c('COD','MAC','SSA','SPR','PLE','SOL')
  a@species.info[firstAgeF1,'first-age F>0']<-1L
  firstAgeF3<-('POK')
  a@species.info[firstAgeF3,'first-age F>0']<-3L
  a@species.info[,'last-age'] <-as.integer(
    # FUL GLT HEG KTW GBG GNT PUF RAZ RAJ GUR WHM NHM GSE HBP HKE COD WHG HAD POK MAC HER NSA SSA,NOP SPR PLE SOL
    c(1,    1,  1,  1,  1,  1,  1,  1,  3,  4,  3,  6,  1,  1,  9, 10,  6, 10, 10, 10,  8,  4,  4,  3,  3, 10, 10))
  
  no.other.predators<-sum(a@species.info[,'predator']==2)
  
  
                 #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@SSB.R<-as.integer(c(0,    0,    0,    1,    0,    1,    1,    1,    2,    2,    0,    0) )
  
  a@rec.season<-3L
  
  bySpAge<-list(
    c(1,2,4), #COD   
    c(0,1,3),   #WHG   
    c(0,1,2), #HAD   
    c(3,4), #POK   
    c(1,2), #MAC   
    c(0,1, 4), #HER   
    c(0,1), #NSA
    c(1,3), #SSA   
    c(0,1), #NOP
    c(0,1), #SPR
    c(0,1,5), #PLE
    c(0,1,5)  #SOL 
  )
  a@catch.s2.group<-bySpAge
 
  
                             #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@incl.process.noise<-as.integer(c(1,    1,    1,    1,    1,    1,    2,    1,    1,    1,    1,    1) )
  
  
  bySpAge<-list(
    c(0,1), #COD   
    c(0,1),   #WHG   
    c(0,1), #HAD   
    c(0,3), #POK   
    c(0,1), #MAC   
    c(0,1), #HER   
    c(0), #NSA
    c(0,1), #SSA   
    c(0,1), #NOP
    c(0,1), #SPR
    c(0,1), #PLE
    c(0,1)  #SOL 
  )
  a@keyVarLogN<-bySpAge
  
  bySpAge<-list(
    c(1,2,3), #COD   
    c(0,1,2),     #WHG   
    c(0,1,2),   #HAD   
    c(3,4,5),   #POK   
    c(1),   #MAC   
    c(0,1,2), #HER   
    c(0,1),   #NSA
    c(1,3),   #SSA   
    c(0,1,2), #NOP
    c(1),     #SPR
    c(1,2,3), #PLE
    c(1,2,3)  #SOL 
  );
  a@keyLogFsta<-bySpAge
  
  
  #                      COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@use.rho<-as.integer(c(1,    1,    0,    0,    0,    1,    0,    0,    0,    0,    1,    1) )
  
  avgF<-matrix(as.integer(c(
    ## first and last age in calculation of average F by species (option avg.F.ages)
  2, 4,  # COD 
  2, 5,  # WHG 
  2, 4,  # HAD 
  4, 7,  # POK 
  4, 8,  # MAC 
  2, 6,  # HER 
  1, 2,  # NSA 
  1, 2,  # SSA 
  1, 2,  # NOP 
  1, 2,  # SPR 
  2, 6,  # PLE 
  2, 6  # SOL 
  )),ncol=2, byrow = TRUE)
  dimnames(avgF)<-dimnames(a@avg.F.ages)
  a@avg.F.ages<-avgF
  
  if (writeConrol) write.RSMS.control(a,file=file.path(data.path,outfile))
  invisible(a)
}



# make  rsms.dat (control) file for estimation of proportion of F by season for species with seasonal catches
batch_seasonal_catch_configuration <-function(outfile='rsms.dat',dir=data.path,writeConrol=TRUE) {
  
  # use the default control file as starting point
  a<-batch_default_configuration(outfile,dir=data.path,writeConrol=FALSE)
  
  
                           #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@combined.catches<-as.integer(c(1,    1,    1,    1,    1,    3,    3,    3,    3,    3,    1,    1) )
  
  
                   #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@fModel  <-as.integer(c(1,    1,    1,    1,    1,    2,    2,    2,    2,    2,    1,    1) )
  
  bySpAge<-list(
    c(1,2,3), #COD   
    c(0,1,2),     #WHG   
    c(0,1,2),   #HAD   
    c(3,4,5),   #POK   
    c(1),   #MAC   
    c(0), #HER   
    c(0),   #NSA
    c(1),   #SSA   
    c(0), #NOP
    c(1),     #SPR
    c(1,2,3), #PLE
    c(1,2,3)  #SOL 
  );
  a@keyLogFsta<-bySpAge
  
  
                              #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@incl.process.noise<-as.integer(c(1,    1,    1,    1,    1,    1,    2,    1,    1,    1,    1,    1) )
 
  
  bySpAge<-list(
    c(0,1), #COD   
    c(0,1),   #WHG   
    c(0,1), #HAD   
    c(0,3), #POK   
    c(0,1), #MAC   
    c(0,1), #HER   
    c(0), #NSA
    c(0,1), #SSA   
    c(0,1), #NOP
    c(0,1), #SPR
    c(0,1), #PLE
    c(0,1)  #SOL 
  )
  a@keyVarLogN<-bySpAge
  
  bySpAge<-list(
    c(1,2,3), #COD   
    c(0,1,2),     #WHG   
    c(0,1,2),   #HAD   
    c(3,4,5),   #POK   
    c(1),   #MAC   
    c(0,1,2), #HER   
    c(0,1),   #NSA
    c(1,2),   #SSA   
    c(0,1,2), #NOP
    c(1),     #SPR
    c(1,2,3), #PLE
    c(1,2,3)  #SOL 
  );
  a@keyLogFsta<-bySpAge
 
 #                                COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
 a@firstAgeYearEffect<-as.integer(c(0,    0,    0,    0,    0,    0,    0,    1,    0,    1,    0,    0) )
 
   
  bySpAge<-list(
    c(0), #COD   
    c(0),     #WHG   
    c(0),   #HAD   
    c(0),   #POK   
    c(0),   #MAC   
    c(0,1,2), #HER   
    c(0,1,3),   #NSA
    c(1,2),   #SSA   
    c(0,1,3), #NOP
    c(1,3),     #SPR
    c(0), #PLE
    c(0)  #SOL 
  );
 a@catch.sep.age<-  bySpAge
  
  bySpYear<-list(
    c(1974), #COD   
    c(1974),     #WHG   
    c(1974),   #HAD   
    c(1974),   #POK   
    c(1974),   #MAC   
    c(1974,1983,1998), #HER   
    c(1974,2005),   #NSA
    c(1974,2005),   #SSA   
    c(1974,2003), #NOP
    c(1974),     #SPR
    c(1974), #PLE
    c(1974)  #SOL 
  );
  a@catch.sep.year<-bySpYear
  
  #                      COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@use.rho<-as.integer(c(1,    1,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1) )
  
  
  if (writeConrol) write.RSMS.control(a,file=file.path(data.path,outfile))
  invisible(a)
}




# final (fine tuned) configuration
batch_final_single_configuration <-function(outfile='rsms.dat',dir=data.path,writeConrol=TRUE) {
  
  # use the default control file as starting point
  a<-batch_default_configuration(outfile,dir=data.path,writeConrol=FALSE)
  
  
                           #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@combined.catches<-as.integer(c(1,    1,    1,    1,    1,    2,    2,    2,    2,   2,     1,    1) )
  
                              #     COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@incl.process.noise<-as.integer(c(1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1,    1) )
  
  bySpAge<-list(
    c(0,1), #COD   
    c(0,1),   #WHG   
    c(0,1), #HAD   
    c(0,3), #POK   
    c(0,1), #MAC   
    c(0,1), #HER   
    c(0,1), #NSA
    c(0,1), #SSA   
    c(0,1), #NOP
    c(0,1), #SPR
    c(0,1), #PLE
    c(0,1)  #SOL 
  )
  a@keyVarLogN<-bySpAge
  
  #                                COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@fModel  <-          as.integer(c(2,    1,    1,    1,    1,    2,    1,    1,    1,    1,    2,    2) )
  
  #                                COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@firstAgeYearEffect<-as.integer(c(3,    0,    0,    0,    0,    2,    0,    0,    0,    0,    3,    3) )
  
  bySpAge<-list(
    c(3,6), #COD   
    c(0),     #WHG   
    c(0),   #HAD   
    c(0),   #POK   
    c(0),   #MAC   
    c(2,5), #HER   
    c(0),   #NSA
    c(0),   #SSA   
    c(0), #NOP
    c(0),     #SPR
    c(3,6), #PLE
    c(3,6)  #SOL 
  );
  a@catch.sep.age<-  bySpAge
  
  
  if (writeConrol) write.RSMS.control(a,file=file.path(data.path,outfile))
  invisible(a)
}


# final (fine tuned) configuration
batch_final_single_configuration2 <-function(outfile='rsms.dat',dir=data.path,writeConrol=TRUE) {
  
  # use the default control file as starting point
  a<-batch_final_single_configuration(outfile,dir=data.path,writeConrol=FALSE)

  #                       COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@fModel  <-as.integer(c(2,    1,    1,    1,    1,    2,    1,    2,    1,    1,    2,    2) )
  
  #                                COD   WHG   HAD   POK   MAC   HER   NSA   SSA   NOP   SPR   PLE   SOL 
  a@firstAgeYearEffect<-as.integer(c(3,    0,    0,    0,    0,    2,    0,    1,    0,    0,    3,    3) )
  
  bySpAge<-list(
    c(3,6), #COD   
    c(0),     #WHG   
    c(0),   #HAD   
    c(0),   #POK   
    c(0),   #MAC   
    c(2,5), #HER   
    c(0),   #NSA
    c(1,3),   #SSA   
    c(0), #NOP
    c(0),     #SPR
    c(3,6), #PLE
    c(3,6)  #SOL 
  );
  a@catch.sep.age<-  bySpAge
  
  bySpYear<-list(
    c(1974), #COD   
    c(1974),     #WHG   
    c(1974),   #HAD   
    c(1974),   #POK   
    c(1974),   #MAC   
    c(1974), #HER   
    c(1974),   #NSA
    c(1974,2005),   #SSA   
    c(1974), #NOP
    c(1974),     #SPR
    c(1974), #PLE
    c(1974)  #SOL 
  );
 # a@catch.sep.year<-bySpYear
  
  if (writeConrol) write.RSMS.control(a,file=file.path(data.path,outfile))
  invisible(a)
}

batch_final_multi_configuration <-function(outfile='rsms.dat',dir=data.path,writeConrol=TRUE) {
  
  # use the default control file as starting point
  a<-batch_final_single_configuration(outfile,dir=data.path,writeConrol=FALSE)

  a@incl.stom.all <- 1L
    a@consum<- 1L
  
  #                         FUL   GLT   HEG   KTW   GBG   GNT   PUF   RAZ   RAJ   GUR   WHM   NHM   GSE   HBP   HKE   COD   WHG   HAD   POK   MAC 
  a@stomach.variance<-     c(4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    4,    5,    5,    5,    5,    5,    5,    5,    5,    5)
  a@size.selection<-       c(0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0)
  a@size.range.lower<-     c(0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    1,    1,    1,    1,    1,    1,    1)
  a@size.range.upper<-  c(100,   100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,   99,   99,   99,   99,   99,   99)
  a@stom.max.sumP   <-     c(5,    5,    5,    5,    5,    5,    5,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100,  100)
  a@size.other.food.suit<- c(0,    0,    0,    0,    0,    0,    0,    0,    1,    0,    0,    0,    0,    0,    0,    0,    1,    0,    1,    0)
  a@prey.pred.size.fac<-   c(5,    5,    5,    5,    5,    5,    5,    5,  0.5,  0.5,  0.5,  0.5,   50,   50,  0.9,  0.5,  0.9,  0.5,   0.5,  0.5)
  a@stom.type.include<-    c(2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2,    2)
  
  if (writeConrol) write.RSMS.control(a,file=file.path(data.path,outfile))
  invisible(a)
}
