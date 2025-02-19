

if (TRUE) {
  
  rsmsControl<-'rsms.dat'
  
  doMultiExtract<-FALSE
  
  runName<-'old_SMS_like'
  rsms<-batch_SMS_old_like_configuration(outfile=rsmsControl,writeConrol=T)
}

### Extract data from SMS
doMultiExtract<-FALSE

if (TRUE) {  # transform  SMS data into RSMS format 
  switch(my.stock.dir,
   "rsms_input"       = inp_all<-make_rsms_data(dir=data.path,rsmsControl='rsms.dat',multi=doMultiExtract),
   "rsms_SAN-area-1r"       = inp_all<-make_rsms_data(dir=data.path,rsmsControl='rsms.dat',multi=FALSE),
    stop(paste0("Not valid stock dir: ",my.stock.dir))
  ) #end switch
  save(inp_all,file=file.path(data.path,"rsms_input_all.Rdata")) 
} else load(file=file.path(data.path,"rsms_input_all.Rdata"),verbose=TRUE)


smsConf<-0L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated


#runName<-'Single'
# select a combination of species from the (full) data set, also including multi species information
my.ps=c(1,2,6,7,8)
my.ps=c(1,7,8)
my.pso<-c(0L)
#my.pso<-13L:27L


inp<-pick_species(ps=my.ps,pso=my.pso, inp=inp_all,smsConf) 
#inp=inp_all

inp$data$sms.mode<-smsConf

inp$data$spNames

#### prepare to run
data<-inp[['data']]
data$silent<-FALSE
parameters<-inp[['parameters']]

if (FALSE) {
 load(paste0(runName,'.Rdata'),verbose=TRUE)
  parameters<-as.list(sms$sdrep, what="Est")  #parameters and random effects
  exp(parameters$logSeparF)
  exp(parameters$logFSeasonal)
}

#### Run rsms
cleanrun(silent=TRUE)

myMap<-map_param(data,parameters)

if (length(myMap$Un)>0) random<-c("Un") else random<-NULL
if (length(myMap$Uf)>0) random<-c(random,"Uf")


system.time(obj <- MakeADFun(func, parameters, random,silent=T,map=myMap))

# checkConsistency(obj);
lu<-lowerUpper(obj,data,parameters )

system.time(opt <-nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1000,eval.max=1000)))
announce(opt)

system.time(sdrep <- sdreport(obj))
cat('Hesssian:',sdrep$pdHess,'\n')

a<-extractParameters(sdrep,myMap,data)[[2]]
print(a,n=600)

myRep<-obj$report()
a<-myRep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a,1)

sms<-saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)

plotCompareRunSummary(Type=c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat")[2],showSpecies=1:12,
                      inpRdata=list(runName,"SMS_old"),
                      labels=c(runName,"SMS_old"),
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:10, ncol=3,allInOne=T,
                      longSpNames=FALSE, fileLabel='single')
-

plotSeasonalData(inp=runName,Type="FiProp", #Type="FiProp",
                 outFormat=c('screen','pdf','png')[1],
                 showAges=0:8,
                 multN=0.001,
                 ncols=3,
                 cummulate=T,
                 fileLabel='pl')


 
plotCompareRunSummary(Type=c(c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat"))[2],showSpecies=1:12,
                                 inpRdata=list("Single","SMS_old"),
                                 labels=c("single","SMS_old"),
                                 outFormat=c('screen','pdf','png')[1],
                                 multN=0.000001,
                                 longSpNames=FALSE, fileLabel='single')

 
plotCompareRunSummary(Type=c(c("SummaryConf","Summary","SSBrecConf","SSBrec","M2","F","N","ExpPat"))[1],showSpecies=1:12,
                       inpRdata=list("Single"),
                       labels=c("single"),
                       outFormat=c('screen','pdf','png')[1],
                       multN=0.000001,
                       longSpNames=FALSE, fileLabel='single')
 
   
   AICCompare( 
     inpRdata=list("Single3","Single4"),
     labels=c("Single3","Single4")
   )
   
   # read results from ICES stock assessment into a suitable format,saved as "outSet" .Rdata
   transExternalSummary(inp='summary_table_raw_ICES.out',outSet='ICES_single_sp',spNames=data$spNames) #, exSpeciesNames=data$spNames)
   transExternalSummary(inp='summary_table_raw.out',outSet='SMS_old',spNames=data$spNames) #, exSpeciesNames=data$spNames)
   transExternalData(inp='summary.out',outSet='SMS_old_detailed',spNames=data$spNames)
   transExternalSummary(inp='ICES_SAN_R1_2025.out',outSet='SandeelR1',spNames=data$spNames,exSpeciesNames=c('SA1'))
   transExternalData(inp='ICES_SAN_R1_2025_detailed.out',outSet='SandeelR1_detailed',spNames=data$spNames, exSpeciesNames=data$spNames)
   transExternalData(inp='ICES_SAN_R1_2025_detailed_twoBlocks.out',outSet='SandeelR1_detailed_twoBlocks',spNames=data$spNames, exSpeciesNames=data$spNames)
   
   

if (FALSE) plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[2],showSpecies=1:12,
                                 inpRdata=list("noTechCreep","TechCreep"),
                                 labels=c("No TechCreep","TechCreep"),
                                 outFormat=c('screen','pdf','png')[1],
                                 multN=0.000001,
                                 longSpNames=FALSE, fileLabel='techCreep_')


# plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[5],
#                       showSpecies=1:12,
#                       inpRdata=list('SandeelR1_detailed','SandeelR1_detailed_twoBlocks'),
#                       labels=c('ICES 2025','ICES 2025 two blocks'),
#                       outFormat=c('screen','pdf','png')[1],
#                       showAges=0:4,
#                       longSpNames=FALSE, fileLabel='san1R_')


# F at age
plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[4],showSpecies=1:12,
                      inpRdata=list('Single','SandeelR1_detailed'),
                      labels=c('RSMS','ICES 2025'),
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      longSpNames=FALSE, fileLabel='san1R_')

plotCompareRunSummary(Type=c("compSummaryConf","compSummary")[2],showSpecies=1:12,
                      inpRdata=list('Single','SandeelR1'),
                      labels=c('RSMS','ICES 2025'),
                      outFormat=c('screen','pdf','png')[1],
                      longSpNames=FALSE, fileLabel='san1R_')

if (FALSE) plotCompareRunSummary(Type=c("compSummaryConf","compSummary")[2],showSpecies=1:12,
                                 inpRdata=list('Single','SMS_old'),
                                 labels=c('Single sp','SMS_old'),
                                 outFormat=c('screen','pdf','png')[1],
                                 longSpNames=FALSE, fileLabel='SS_oldSMS_')


if (FALSE) plotCompareRunSummary(Type=c("compSummaryConf","compSummary")[2],showSpecies=1:12,
                      inpRdata=list('Single','ICES_single_sp'),
                      labels=c('Single sp','ICES'),
                      outFormat=c('screen','pdf','png')[1],
                      longSpNames=TRUE, fileLabel='SS_ICES_')


if (FALSE) plotCompareRunSummary(Type=c("compSummaryConf","compSummary")[2],showSpecies=1:12,
                                 inpRdata=list('Single'),
                                 labels=c('Single sp'),
                                 outFormat=c('screen','pdf','png')[3],
                                 longSpNames=FALSE, fileLabel='techCreepNO_')



###################  multi species mode 1
runName<-'Single'; load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)

runName<-'Multi'
smsConf<-1L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated

cleanrun(silent=FALSE)

sdrep<-sms$sdrep # from single species run
cat('Hesssian,single:',sdrep$pdHess,'\n')

newPar<-as.list(sms$sdrep, what="Est")  #parameters and random effects

my.pso<-c(13L:27L)
my.ps; my.pso
inp<-pick_species(ps=my.ps,pso=my.pso,inp=inp_all,smsConf) # to extract multi species data and parameters

inp$data$sms.mode<-smsConf
data<-inp[['data']]

# add MS parameters
newPar$vulnera<-inp$parameters$vulnera
newPar$overlapP<-inp$parameters$overlapP
newPar$logStomObsVar<-inp$parameters$logStomObsVar
parameters<-newPar

myMap<-map_param(data,parameters)

#lock single species parameters,if necessary
doLock<-TRUE
if (doLock) {
  lockP<-c('logCatchability','logSdLogN','logSdLogFsta','logSdLogObsSurvey','logSdLogObsCatch','rho','rec_loga','rec_logb','logSeparF')  #,'Uf','Un'
   myMap<-lock_param(data,parameters,myMap,lockP) 
}

system.time(obj <- MakeADFun(func, parameters, random,silent=T,map=myMap))

lu<-lowerUpper(obj,data,parameters)
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1300,eval.max=1300)))
announce(opt)

system.time(sdrep <- sdreport(obj)); 
cat('Hesssian:',sdrep$pdHess,'\n')
#sdrep

a<-extractParameters(sdrep,myMap,data)[[2]] 
arrange(a,desc(abs(gradient)))
arrange(a,desc(is.na(estimate.sd)))
#rint(a,n=500)


myRep<-obj$report()
a<-myRep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a)

sms<-saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
#load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)


if (FALSE) {
  plotCompareRunSummary(Type=c("compSummaryConf","compSummary")[2],showSpecies=1:10,
                      inpRdata=list('Single','Multi'),
                      labels=c('Single sp','Multi sp'),
                      outFormat=c('screen','pdf','png')[1],
                      longSpNames=TRUE,fileLabel='SS_MS1_')


}


###################  multi species mode 2
runName<-'Multi'; load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)

runName<-'Multi2'
smsConf<-2L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated

cleanrun(silent=FALSE)

sdrep<-sms$sdrep # from single species run
cat('Hesssian:',sdrep$pdHess,'\n')
newPar<-as.list(sms$sdrep, what="Est")  #parameters and random effects

inp$data$sms.mode<-smsConf
data<-inp[['data']]
parameters<-newPar
myMap<-map_param(data,parameters)

system.time(obj <- MakeADFun(func, parameters, random,silent=T,map=myMap))

lu<-lowerUpper(obj,data,parameters )
system.time(opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1300,eval.max=1300)))
announce(opt)

system.time(sdrep <- sdreport(obj)) 
cat('Hesssian:',sdrep$pdHess,'\n')
#sdrep


a<-extractParameters(sdrep,myMap,data)[[2]] %>% mutate(expParm=round(exp(estimate),2))
#arrange(a,desc(abs(gradient)))
#arrange(a,desc(is.na(estimate.sd)))
#print(a,n=500)

myRep<-obj$report()
a<-myRep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a)

sms<-saveResults(runName=runName,data=data,parameters=parameters,obj=obj,opt=opt,lu=lu,map=myMap,random=random,rep=myRep,sdrep=sdrep)
#load(file=file.path(data.path,paste0(runName,".Rdata")),verbose=TRUE)

a<-sms$rep$nlls; a<-rbind(a,all=colSums(a)); a<-cbind(a,all=rowSums(a));round(a)

plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2")[3],
                      showSpecies=1:10,M2ages=0:4,ncol=3,
                      inpRdata=list('Multi2','SMS_old_detailed'),
                      labels=c('Multi sp2','old SMS'),
                      outFormat=c('screen','pdf','png')[3],
                      longSpNames=TRUE,fileLabel='Multi2_SMSold')

plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2")[2],
                      showSpecies=1:10,ncol=3,
                      inpRdata=list('Multi2','SMS_old'),
                      labels=c('Multi sp2','old SMS'),
                      outFormat=c('screen','pdf','png')[3],
                      longSpNames=TRUE,fileLabel='Multi2_SMSold')

plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2")[1],showSpecies=1:12,
                      inpRdata=list('Multi2'),
                      labels='Multi sp',
                      outFormat=c('screen','pdf','png')[3],
                      longSpNames=TRUE,fileLabel='Multi2')


hertil




inp$data$stom[1,'data'][[1]][[1]]
rep1$stom[1,'data'][[1]][[1]]

inp$data$suitIdx
inp$data$suitIdx[1,'data'][[1]][[1]]
inp$data$suitIdx[1,'data'][[1]][[1]][['data']][[2]]


inp$data$stom[1,'data'][[1]][[1]][['data']][[2]]
rep1$stom[1,'data'][[1]][[1]][['data']][[2]]

inp$data$partM2Idx
inp$data$partM2Idx[1,'data'][[1]][[1]]
inp$data$partM2Idx[1,'data'][[1]][[1]][['data']][[2]]

inp$data$partM2Idx[1,'data'][[1]][[1]][['data']][[2]]
rep1$partM2Idx[1,'data'][[1]][[1]][['data']][[2]]

stom2<-unnest(rep1$stom,cols = c(data))
stom2<-unnest(stom2,cols = c(data))

stom2 %>% select(pred,prey,vulneraIdx) %>% unique() %>% arrange(pred,prey)

rep1$res
M2 for 0-gruppe i Q1 ????
  
  
xx<-rep1$res %>% group_by(species,age) %>% mutate(sumM=sum(M2,na.rm=TRUE)) %>% filter(sumM>0) %>% 
  group_by(s,species,year,age) %>% summarize(M2=sum(M2,na.rm=TRUE))  %>% ungroup() %>% mutate(age=factor(age))
xx

  
  
ggplot(xx,aes(x=year,y=M2,shape=age,col=age))+
  geom_point(size=2)  + geom_line()+labs(ylab='M2',xlab='Year')+
  facet_wrap(vars(data$spNames[s]),ncol=2,scales="free_y")

rep1$nlls
rep1$availFood
# p<-extractParameters(sdrep1)
#view(p[[2]])

save(parameters,data,obj1,my.map1,newPar1,lu,random,file=file.path(data.path,"multi_sp1.Rdata"))
load(file=file.path(data.path,"multi_sp1.Rdata"),verbose=T)

##########################
smsConf<-2L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
load(file=file.path(data.path,"multi_sp1.Rdata"),verbose=T)

inp$data$sms.mode<-smsConf
data<-inp[['data']]

parameters<-newPar1

my.map<-map_param(data,parameters)
my.map2<-my.map
lockParameters<-FALSE
if (lockParameters) {
  #lockP<-c('logCatchability','logSdLogN','logSdLogFsta','logSdLogObsSurvey','logSdLogObsCatch','rho','rec_loga','rec_logb')
  lockP<-c('Uf','Un')
  my.map2<-lock_param(data,parameters,my.map,lockP)
}

system.time(obj2 <- MakeADFun(func, parameters, random,silent=F,map=my.map2))

lu<-lowerUpper(obj2,data,parameters)

if (exists('opt2')) rm(opt2);if (exists('sdrep2')) rm(sdrep2);
system.time(opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, lower=lu[['lower']], upper=lu[['upper']],control=list(iter.max=1300,eval.max=1300)))
announce(opt2)

sdrep2 <- sdreport(obj2,ignore.parm.uncertainty=FALSE); # ignore.parm.uncertainty=TRUE to speed up
cat('Hesssian:',sdrep2$pdHess,'\n')
sdrep2

a<-extractParameters(sdrep2)[[2]]
arrange(a,desc(abs(gradient)))


vulnerabilityTab<-function(sdrep,expIt=FALSE,asDF=FALSE) {
  vul<-sdrep[['par.fixed']]
  vul<-vul[names(vul)=='vulnera']
  if (expIt) vul<-exp(vul)
  vulnerability<-data$vulneraIdx
  vulnerability[match(1:length(vul),vulnerability)]<-vul
  vulnerability<-vulnerability[,data$preds]

  if (asDF) {
    vul.se<-sqrt(diag(sdrep$cov.fixed))
    vul.se<-vul.se[names(vul.se)=='vulnera']
    #if (expIt) vul.se<-exp(vul.se)
    vulnerability.se<-data$vulneraIdx
    vulnerability.se[match(1:length(vul.se),vulnerability.se)]<-vul.se
    vulnerability.se<-vulnerability.se[,data$preds]
    vulnerability.se<- array2DF(vulnerability.se)
    colnames(vulnerability.se)<-c('prey','pred','vulnerability.se')
    
    
    vulnerability<- array2DF(vulnerability)
    colnames(vulnerability)<-c('prey','pred','vulnerability')
    
    vulnerability<-left_join(vulnerability,vulnerability.se,by = join_by(prey, pred))
    vulnerability<-filter( vulnerability,!(vulnerability==0))
  }
  return(vulnerability)
}

vulnerabilityTab(sdrep1,expIt=TRUE,asDF=TRUE)
round(vulnerabilityTab(sdrep1,expIt=TRUE,asDF=FALSE),3) 

round(vulnerabilityTab(sdrep2,expIt=TRUE,asDF=FALSE),3) 

overlapTab<-function(sdrep,expIt=FALSE) {
  ov<-sdrep[['par.fixed']]
  ov<-ov[names(ov)=='overlapP']
  if (expIt) ov<-exp(ov)

  overlap<-data$overlapIdx
  overlap$overlap<-ov[overlap$overlap]

  return(overlap)
}

overlapTab(sdrep=sdrep1,expIt=TRUE) 
overlapTab(sdrep=sdrep2,expIt=TRUE) 


if (FALSE) {
   res<-obj2$report();  
   resul<-res$res
   resulAn<-res$resAnnual
   filter(resul,species=='WHG' & quarter==3)
   res$nlls; 
 
   M2Graph<-function(preys='COD',inp=resul,annual=TRUE,maxAgePlot=4) {
     res<-filter(inp,species %in% preys)
     res<-res %>% filter(age<=maxAgePlot) %>% group_by(species,year,age) %>% summarize(M2=sum(M2))  %>% mutate(age=factor(age))
    
     p<-ggplot(data=res, aes(x=year, y=M2, group=age)) +
       geom_line(aes(linetype=age,col=age))+
       geom_point(aes(shape=age,col=age))+
       facet_grid(species ~ ., scales="free_y")
       #ggtitle(unlist(data[1,'species']))
     print(p)
    ftable(round(xtabs(M2~species+year+age,data=filter(res,year %in% c(1974:1979,2020:2022))),6))
   }
   
   M2Graph(preys='COD') 
   M2Graph(preys=c('COD','WHG')) 
   M2Graph(preys=c('HER','NOP')) 
   M2Graph(preys=c('NSA','SSA')) 
 
   M2Graph(preys=c('HER')) 
   M2Graph(preys=c('NOP'))  # age 3 is zero ???????
}
