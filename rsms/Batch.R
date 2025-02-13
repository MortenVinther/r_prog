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

###############  Multi sp 

runName<-'Multi'

my.ps=c(1,2,3,4,5,6,7,8,9,10,11,12)
my.pso<-c(0L)
#my.pso<-13L:27L

extractMulti<-TRUE

rsms<-batch_default_configuration(outfile=smsControlFile,writeConrol=T)

extractDataAndRun(
  runName,
  my.ps,my.ps0,
  doMultiExtract=extractMulti,
  dir=data.path,
  silent=TRUE,
  smsControlFile=my.smsControlFile,
  fleetIn="new_fleet_info.dat",
  smsConf=1 # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
) 







extractDataAndRun<-function(runName=runName,
                            my.ps=1,my.pso=0,
                            doMultiExtract=FALSE,
                            dir=data.path,
                            silent=TRUE,
                            sms.dat='rsms.dat',
                            fleetIn="new_fleet_info.dat",
                            smsConf=0L # 0=single species, 1=multi species, but fixed single species parameters, 2=multi species, all parameters are estimated
) 
  
 


  
inpRdata<- list('RUN3','Single')
labels<- c('RUN3','Single')


AICCompare(inpRdata,labels) 


plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[2],showSpecies=1:12,
                      inpRdata,
                      labels,
                      outFormat=c('screen','pdf','png')[1],
                      showAges=0:7,
                      longSpNames=FALSE, fileLabel='san')


plotCompareRunSummary(Type=c("compSummaryConf","compSummary","compM2","compF","compN")[1],showSpecies=1:12,
                      inpRdata=list('RUN3'),
                      labels=c('RUN3'),
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

