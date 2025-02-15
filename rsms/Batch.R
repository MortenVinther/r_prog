#### prepare the "baseline run"
my.ps=c(1,2,3,4,5,6,7,8,9,10,11,12)
my.pso<-c(0L)


rsmsControl<-'rsms.dat'

doMultiExtract<-FALSE

################################
# make default control file and run
runName<-'RUN1'
rsms<-batch_default_configuration(outfile=rsmsControl,writeConrol=T)

extractDataAndRunSingle(
  runName,
  my.ps,my.ps0,
  rsmsControl,
  doMultiExtract,
  dir=data.path,
  silent=TRUE,
  fleetIn="new_fleet_info.dat",
  smsConf=0L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
) 

###########################
# Estimate seasonal F proportions from species with seasonal catches
runName<-'RUN2'
rsmsRun2<-batch_seasonal_catch_configuration(outfile=rsmsControl,writeConrol=TRUE)
extractDataAndRunSingle(runName, my.ps, my.ps0, rsmsControl)

writeSeasonalF(inp=runName,outfile="proportion_of_annual_f.in",dir=data.path) 
#############################

# make default control file and run again with the estimated seasonal F proportions
runName<-'RUN3'
rsms<-batch_default_configuration(outfile=rsmsControl,writeConrol=T)
extractDataAndRunSingle(runName, my.ps, my.ps0, rsmsControl)

########
# final fine tuning
runName<-'Single'
batch_final_single_configuration(outfile=rsmsControl,dir=data.path,writeConrol=TRUE) 
extractDataAndRunSingle(runName, my.ps, my.ps0, rsmsControl)

# getParameters(runName)

########
# final fine tuning testing
runName<-'Single_test'
batch_final_single_configuration2(outfile=rsmsControl,dir=data.path,writeConrol=TRUE) 
extractDataAndRunSingle(runName, my.ps, my.ps0, rsmsControl)


inpRdata<- list('Single','Single_test')
labels<-      c('Single',"Single_test")

AICCompare(inpRdata,labels) 


plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[2],showSpecies=8,
                      inpRdata,
                      labels,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      longSpNames=FALSE, fileLabel='san')

###############  Multi sp 

runName<-'Multi1'
my.ps=c(1,2,3,4,5,6,7,8,9,10,11,12)
my.pso<-c(0L)
#my.pso<-13L:27L

extractMulti<-TRUE
rsms<-batch_final_single_configuration(outfile=rsmsControl,writeConrol=T)

extractDataAndRunMulti(
  runName,
  my.ps,my.pso,
  rsmsControl,
  doMultiExtract=extractMulti,
  parametersFrom="Single",
  doLock=TRUE,
  dir=data.path,
  silent=TRUE,
  fleetIn="new_fleet_info.dat",
  smsConf=1L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
) 


runName<-'Multi2'
extractDataAndRunMulti(
  runName,
  my.ps,my.pso,
  rsmsControl,
  doMultiExtract=extractMulti,
  parametersFrom="Multi1",
  doLock=FALSE,
  dir=data.path,
  silent=TRUE,
  fleetIn="new_fleet_info.dat",
  maxIter=3000,
  smsConf=2L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
) 




getParameters(runName)


inpRdata<- list('Single','Multi2')
labels<-      c('single','Multi2')


#AICCompare(inpRdata,labels) 


plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[4],showSpecies=8,
                      inpRdata,
                      labels,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      longSpNames=FALSE, fileLabel='san')


plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[4],showSpecies=8,
                      inpRdata=list('Multi2'),
                      labels=c('Multi2'),
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      longSpNames=FALSE, fileLabel='san')



### Extract data from SMS
doMultiExtract<-FALSE


fleetIn<-"new_fleet_info_template.dat"

system.time(res<-extractDataAndRun(runName,my.ps=ps,my.pso=pso,doMultiExtract=FALSE,sms.dat=sms.dat.in,fleetIn=fleetIn,silent=TRUE) )


load(file=file.path(data.path,paste0(runName,'.Rdata')),verbose=T)
a<-extractParameters(sms$sdrep,sms$map,sms$data)[[2]]
print(a,n=300)

#Used later on
off.age<-sms$data$off.age
off.year<-sms$data$off.year
nSpecies<-sms$data$nSpecies
speciesNames<-sms$data$spNames

###





fleetOut<-paste0("new_fleet_info_",runName,".dat")
key<-readNewFleetInfo(dir=data.path,of=fleetIn,off.age=off.age,off.year=off.year,nSpecies=1,verbose=FALSE) 

# change fleet
key[['k']][,'useFleet']<-0L
key[['k']][c(1,2,3),'useFleet']<-1L
writeNewFleetInfo(dir=dir,of=fleetOut,key,speciesNames)



inpRdata<- list('RUN1','RUN2','RUN3','RUN6')
labels<- c('RUN1','RUN2','RUN3','RUN6')


AICCompare(inpRdata,labels) 
  
plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[4],showSpecies=1:12,
                      inpRdata=list(runName,'SandeelR1_detailed'),
                      labels=c(runName,'ICES 2025'),
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      longSpNames=FALSE, fileLabel='san1R_')

plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[5],showSpecies=1:12,
                      inpRdata=list(runName,'SandeelR1_detailed'),
                      labels=c(runName,'ICES 2025'),
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:4,
                      multN=0.000001,
                      longSpNames=FALSE, fileLabel='san1R_')

