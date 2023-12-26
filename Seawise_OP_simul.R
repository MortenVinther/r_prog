# this script does
# 0. some administration
# 1. Used NS 2020 keyrun (last year 2019) to produce Stock numbers in 2020 (done by SMS automatically)
# 2. Update N to the beginning of 2021 from input F(2020) taken from ?


################ 0 ´###################
# you first have to run an SMS to produce OP data

if (FALSE) {
  do.a.full.SMS.run(label="run_",                   # label for output
                    cleanup=T,                      # delete files in the deleteFiles variable?
                    do.single=T,                    # run SMS in single species mode
                    do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                    do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                    do.multi.2.redo=T,              # Run the full model, with simultaneously estimation of all parameters
                    do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                    do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                    shake.ms2.par=F,
                    SSB.R.seperate=F,               # Estimate S/R parameters in a separate step  
                    do.MCMC=F,                      # Prepare for MCMC analysis
                    mcmc=0,mcsave=0,                # Options for MCMC analysis
                    do.prediction=F,                # Make a prediction
                    pause=F,                        # Make a confirm between each stage
                    Screen.show=F,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                    do.run=T,                       # Make the run immediately, or just make the batch file for the run
                    deleteFiles=deleteFiles,        # clean up in files before the run is made
                    HPC=F)                          # run it as batch program on the UNIX High Performance Computer 
  
}

library(tidyverse)
scenF<-readxl::read_excel(file.path("H:","Seawise","scenarierTilAlex","Summary.table_NS_fishingMort_stock.xlsx"), sheet = "F per stock")
sort(unique(scenF$stock))
ss<-c("COD-NS","WHG-NS","HAD","POK","PLE-NS","SOL-NS" )
sn<-c(16,17,18,19,26,27)
names(sn)<-ss
sn
scenF<-scenF %>% mutate(Species.n=sn[stock]) %>% rename(FF=`F`) %>% dplyr::select(scenario,year,Species.n,FF)

scens<-readxl::read_excel(file.path("H:","Seawise","scenarierTilAlex","scenarios_AK.xlsx"), sheet = "Metadata") %>%
   mutate(description=NULL) %>% filter(UsedInRep!="Deleted")

scenF<-left_join(scens,scenF,by = "scenario") %>% mutate( scenario=UsedInRep,UsedInRep=NULL )

ar<-expand.grid(Species.n=unique(scenF$Species.n),year=unique(scenF$year),scenario='ICES-AR')
targetF<-scan(file=file.path(data.path,"op_mulTargetf_FMSY.in"),skip=1L)
targetF<-data.frame(targetF=targetF,Species.n=first.VPA:nsp)
ar<-left_join(ar,targetF,by = "Species.n") %>% rename(FF=targetF) %>% tibble()
ar
scenF<-rbind(scenF,ar) %>% mutate(stock=sp.names[Species.n])
# ftable(round(tapply(scenF$FF,list(scenF$scenario,scenF$year,scenF$stock),FUN=sum),3))

scenarios<-unique(scenF$scenario)
scenario.dirs<-paste("NS_2022_Seawise",scenarios,sep='_')


# reset skill input files
#skill_fleet_make(control=SMS.control,make_excl=TRUE,write_data=TRUE,exclFile='skill_exclude.in') 
#makes the files skill_cpue.in and skill_exclude.in  

if (FALSE) {
  labels<-c("2020 Keyrun","2023 Seawise")  
  dirs<-c("Seawise_NorthSeaKeyRun_2020","NS_2022_seawise")
  
  compare_runs(
    dirs=dirs,
    labels=labels,
    nox=2, noy=3,
    paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
    run.ID='ICES_com',         # file id used for paper output
    doGrid=TRUE,
    extent.SSB=FALSE,  # plot SSB for the year after last assessment year
    first.year.on.plot=1974,
    last.year.on.plot=2020,
    plot.MCMC=FALSE,                        # plot values from MCMC scenarios. FALSE=plot hindcast values from "summary_table_raw.out"
    single.species=FALSE,                   # single species mode or multi species mode
    include.assess.forcast.line=FALSE,      # vertical line at last assessment year
    include.F.reference.points=FALSE,
    include.SSB.reference.points=FALSE,
    include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
    include.2.std=FALSE,
    #incl.sp=c('Herring'),                      # species number to be included. Numbers or "all"
    #incl.sp="all",
    first.pch=0,    # first pch symbol
    first.color=1,   # first color
    palette="default"               # good for colour full plots
    #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
   ) 
  
  
  compare_runs_stock_rec(dirs,labels,first.year.on.plot=1975,last.year.on.plot=2020,incl.sp="all",
                                   include.CV=TRUE,include.CV2=TRUE,include.mean=TRUE,
                                   palette="R3", makeAllGraphs=FALSE,nox=1, noy=2, w8=8,w11=8,
                                   include.year.labels=TRUE,incl_not_used=TRUE,run.ID='SBB_rec', 
                                   paper=TRUE,facSSB=1000,facRec=1000000,
                                   compare.dir=data.path,verbose=FALSE,writeData=TRUE) 
  
  compare_runs_stock_rec(dirs[2],labels[2],first.year.on.plot=1975,last.year.on.plot=2020,incl.sp="all",
                         include.CV=TRUE,include.CV2=TRUE,include.mean=TRUE,
                         palette="R3", makeAllGraphs=FALSE,nox=1, noy=1, w8=5,w11=7,
                         include.year.labels=TRUE,incl_not_used=TRUE,run.ID='SBB_rec_1', 
                         paper=TRUE,facSSB=1000,facRec=1000000,
                         compare.dir=data.path,verbose=FALSE,writeData=TRUE) 
  
  compare_runs_stock_rec(dirs[2],labels[2],first.year.on.plot=1975,last.year.on.plot=2020,incl.sp="all",
                         include.CV=TRUE,include.CV2=TRUE,include.mean=TRUE,
                         palette="R3", makeAllGraphs=FALSE,nox=3, noy=4, w8=15,w11=11.5,
                         include.year.labels=TRUE,incl_not_used=TRUE,run.ID='SBB_rec_2', 
                         paper=TRUE,facSSB=1000,facRec=1000000,
                         compare.dir=data.path,verbose=FALSE,writeData=TRUE,combinePlot=TRUE) 
}


hind<-Read.summary.table() %>% filter(Year<=SMS.control@last.year.model)


#################    2.  ##############
do_scenario<-function(scenarioNo) {
 # test scenarioNo<-2  
  
  data.path<<-file.path(root,my.stock.dir)
  setwd(data.path)
  scenario<-scenarios[scenarioNo]
  scenario.dir.full<-file.path(root,scenario.dirs[scenarioNo])
  cat(scenario.dir.full,'\n')
  if (file.exists(scenario.dir.full)) unlink(scenario.dir.full,recursive = T)
  dir.create(scenario.dir.full,showWarnings = FALSE)
  
  OP.files<-c("sms.dat","op.exe","area_names.in","species_names.in","op.dat","op_trigger.dat","op_config.dat","op_msfd.dat","just_one.in",
              "op_consum.in","op_f.in","op_m1.in","op_m.in","op_n.in","op_propmat.in","op_prop_landed.in","op_size.in","op_wcatch.in","op_wsea.in",
              "op_growth_type1.in","op_consum_ab.in","op_other_n.in","op_exploitation.in","op_reference_points.in","covariance_rec.in","op_price.in",
              "op_ssb_rec_residuals.in","op_length_weight_relations.in",'op_eqsim.in','op_eqsim_stoch.in','op_n_proportion_m2.in','op_seed.in',
              "reference_points.in","summary.out")
  
  
  
  for (from.file in OP.files) {
    to.file<-file.path(scenario.dir.full,from.file)
    if (!file.copy(from.file, to.file, overwrite = TRUE)) cat('problems with file: ',from.file,'\n')
  }
  
   
  # 
  
  source(file=file.path(prog.path.func,'hcr_op_batch_common.R'))
  
  res<-make.OP.dat(my.area='North Sea',my.last.year=2020,first.year.output=2020,do.indicators=F,stochastic.recruitment=0,recruit.adjust.CV=0, recruit.adjust.factor=1,    
                        y.end=2019,y.first=2017,y.first.other=2019)
  
  
  SMS<-res[["SMS"]]
  
  nsp<-SMS@no.species
  n.other.pred<-sum(SMS@species.info[,'predator']==2)
  n.pred<-n.other.pred+sum(SMS@species.info[,'predator']==1)
  n.vpa<-nsp-n.other.pred
  n.vpa.pred<-sum(SMS@species.info[,'predator']==1)
 
  N<-read.table(file=file.path(scenario.dir.full,"op_n.in"))
  N<-as.matrix(N,byrow=TRUE)
  colnames(N)<-SMS.control@first.age:SMS.control@max.age.all
  rownames(N)<-SMS.control@species.names[(n.other.pred+1):nsp]
  N
  N['Sole',"2"]<- N['Sole',"2"]*0.5
  N['Plaice',"2"]<- N['Plaice',"2"]*0.75
  write.table(N,file=file.path(scenario.dir.full,"op_n.in"),col.names = FALSE,row.names = FALSE)
  
  OP<-res[["OP"]]
  OP@F.or.C[1,]<-31  #Update N by target F values from file op_mulTargetf.in
  OP@output<-15
  OP@recruit.min[]<-10
  OP@recruit.adjust[,'Saithe']<-0.9
  OP@recruit.adjust.CV[,'Haddock']<-2
  OP@recruit.adjust.CV[,'N.sandeel']<-2
  OP@recruit.adjust.CV[,'S.sandeel']<-2
  OP@recruit.adjust.CV[,'Nor.pout']<-2
  OP@recruit.adjust.CV[,'Sprat']<-2
  
  
  write.FLOP.control(control=OP,path=scenario.dir.full)
  
  file.copy("op_mulTargetf_2020.in",file.path(scenario.dir.full,"op_mulTargetf.in"),overwrite = TRUE)
  
  opt<-read.FLOPtrigger.control(file="op_trigger.dat",path=scenario.dir.full,n.VPA=n.vpa,n.other.pred)

  opt@HCR[1,]<-1 # all constant F from input
  write.FLOPtrigger.control(opt,path=scenario.dir.full,nice=FALSE,writeSpNames=FALSE)
  
  
  #run OP command
  op.command<-'op.exe'
  comm.sms<-paste(op.command,"-maxfn 0 -nohess ",sep=" ")
  setwd(scenario.dir.full)
  shell(comm.sms, invisible = TRUE)
  setwd(data.path)
  
  
  N<-read.table(file.path(scenario.dir.full,"op_n.out"))
  file.copy(from=file.path(scenario.dir.full,"op_n.out"),to=file.path(scenario.dir.full,"op_n.in"),overwrite = TRUE)
  
  
  
  targetF<-scan(file=file.path(data.path,"op_mulTargetf_FMSY.in"),skip=1L)
  targetF<-data.frame(targetF=targetF,Species.n=first.VPA:nsp)
  
  master<-read.table(file=file.path(scenario.dir.full,"op_condensed.out"),header=TRUE)
  master2<-read.table(file=file.path(scenario.dir.full,"op_summary.out"),header=TRUE)
  
  # change HCR for 2021 onwards, Fixed F for demersal (FLBIEA) stocks and AR or "escapement proxy" for the rest 
  #               Cod Whiting Haddock Saithe Mackerel Herring N.sandeel S.sandeel Nor.pout Sprat Plaice Sole
  if (scenario!='ICES-AR') opt@HCR[1,]<-c(   1,     1,       1,    1,       2,       2,      22,        22,      22,   22,      1,    1)
  if (scenario=='ICES-AR') opt@HCR[1,]<-c(   2,     2,       2,    2,       2,       2,      22,        22,      22,   22,      2,    2) #ICES AR
  
  
  do_one_year<-function(sYear=2022,sScenario="Baseline_CaseStudy") {
  
   ff<-filter(scenF,scenario==sScenario & year==sYear) %>% dplyr::select(Species.n,targetF=FF) 
   ff<- rbind(filter(targetF,!(Species.n %in% sn)),ff) %>% arrange(Species.n)
   ff$targetF
   cat(1,'\n', ff$targetF,"\n",file=file.path(scenario.dir.full,"op_mulTargetf.in"))
   
   opt@first.year<-sYear
   opt@last.year<-sYear
   write.FLOPtrigger.control(opt,path=scenario.dir.full,nice=FALSE,writeSpNames=FALSE)
    
   
   OP@first.year <-sYear
   OP@first.year.out <-sYear 
   OP@last.year <-sYear
   write.FLOP.control(OP,path=scenario.dir.full,nice=F,writeSpNames=F)
   N<-read.table(file.path(scenario.dir.full,"op_n.out"))
  
   file.copy(from=file.path(scenario.dir.full,"op_n.out"),to=file.path(scenario.dir.full,"op_n.in"),,overwrite = TRUE)
   setwd(scenario.dir.full)
   shell(comm.sms, invisible = TRUE)
   setwd(data.path)
   
   master<<-rbind(master,read.table(file=file.path(scenario.dir.full,"op_condensed.out"),header=TRUE))
   master2<<-rbind(master2,read.table(file=file.path(scenario.dir.full,"op_summary.out"),header=TRUE))
  }
  
  
  for (y in 2022:2060) do_one_year(sYear=y,sScenario=scenarios[scenarioNo])
  
  mm<-master %>% transmute(Species=sp.names[Species.n],Species.n=Species.n,Year=Year,Rec=recruit,SSB=SSB,
                       TSB=TSB,SOP=CWsum,SOP.hat=CWsum,SOP.core=CWsum.core,Yield=CWsum,Yield.hat=CWsum,
                       Yield.core=CWsum,mean.F=Fbar,Eaten=eaten)
  
  mm<-rbind(hind,mm)
  
  write.table(mm, file=file.path(scenario.dir.full,'summary_table_raw.out'),col.names = TRUE,row.names = FALSE)
  
  write.table(master2, file=file.path(scenario.dir.full,'op_summary.out'),col.names = TRUE,row.names = FALSE)
  
  data.path<<-scenario.dir.full
  plot_summary_ices_multi(
    Portrait=T,                 # graphical output orientation
    include.terminal.year= FALSE,          # plot terminal year (last assessment year +1) as well?
    include.last.assess.year.recruit=FALSE,          # plot recruits terminal year as well?
    
    first.year= -1974,                #first year on plot, negative value means value defined by data
    last.year= 2050,             #last year on plot
    incl.M2.plot=TRUE,
    incl.reference.points=TRUE,
    incl.TSB=FALSE,
    splitLine=TRUE,
    OperatingModel=TRUE,
    redefine.scenario.manually=FALSE,
    output.dir=scenario.dir.full,
    op.dir=scenario.dir.full,,
    my.dev=c('screen','wmf', 'png', 'pdf')[3]
  )
  
  data.path<<-file.path(root,my.stock.dir)
}
# do_scenario(scenarioNo=1)
for (i in (1:length(scenarios))) do_scenario(scenarioNo=i)

 getwd()   
 
if (TRUE) {
  compare_runs(
      dirs=c(scenario.dirs),  # directory files to be compared
      labels=c(scenarios),  # legends
    
      nox=2, noy=3,
      portrait=TRUE,
      paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
      run.ID='all',         # file id used for paper output
      doGrid=TRUE,
      extent.SSB=FALSE,  # plot SSB for the year after last assessment year
      first.year.on.plot=1974,
      last.year.on.plot=2060,
      plot.MCMC=FALSE,                        # plot values from MCMC scenarios. FALSE=plot hindcast values from "summary_table_raw.out"
      single.species=FALSE,                   # single species mode or multispecies mode
      plot.yield=TRUE,
      include.assess.forcast.line=TRUE,      # vertical line at last assessment year
      include.F.reference.points=FALSE,
      include.SSB.reference.points=FALSE,
      include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
      include.2.std=FALSE,
      std.first.only=FALSE,
      w8=9,w11=11,
      makeAllGraphs=FALSE, # make plots for HTML output
      compare.dir=data.path,
      incl.sp="all"
   )
  
  a<-read_csv(file.path(data.path,"R_F_SSB.csv"))
  ac<-filter(a,Species=='Cod' & variable=='F')
  round(tapply(ac$value,list(ac$scenario,ac$Year),sum),3)
  } 
 
