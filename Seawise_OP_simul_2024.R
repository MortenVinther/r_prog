# this script does
# 0. some administration
# 1. Used NS 2020 keyrun (last year 2019) to produce Stock numbers in 2020 (done by SMS automatically)
# 2. Update N to the beginning of 2021 from input F(2020) taken from ?


output.path<-file.path('H:',"Seawise","scenarierTilAlex_2024")

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

## recruitments
sp<-c('Cod','Herring','Plaice','Saithe','Whiting')

a<-lapply(sp,function(x) {
  print(x)
  rec<-readxl::read_excel(file.path("H:","Seawise","scenarierTilAlex_2024","EMSRR_scaled_recruitment_noCC.xlsx"), sheet =x)
})
rec<-do.call(rbind,a)

rec<-rec %>% rename(rcp45=rcp45_rec_percent,rcp85=rcp85_rec_percent) 
rec<-rec %>% pivot_longer(
  cols = starts_with("rcp"),
  names_to = "scenario",
  values_to = "rec"
)  
  
  
ggplot(data = filter(rec,stock=='Cod'), aes(x = year, y = rec,col=scenario)) + 
  geom_point(size = 1) +
  facet_wrap(vars(ssb))

ggplot(data = filter(rec,stock=='Cod' & rec>0), aes(x = year, y = rec,col=scenario)) + 
  geom_point(size = 1) +
  facet_wrap(vars(scenario))

#rec2<-filter(rec,ssb>0 &  stock=='Herring')
rec2<-filter(rec,ssb>0 &  stock=='Plaice')
library(mgcv)
a<-gam(rec~year+ssb+scenario,data=rec2)
plot(rec2$year,predict(a))
summary(a)

a1<-gam(rec~s(year)+ssb+scenario,data=rec2)
plot(rec2$year,predict(a1))
summary(a1)

a2<-gam(rec~s(year)+scenario,data=rec2)
plot(rec2$year,predict(a2))
summary(a2)

data.frame(stock=rec2$stock,scenario=rec2$scenario,year=rec2$year,recpct=predict(a2)) %>% unique() %>%arrange(stock,scenario,year)

b<-lapply(sp,function(x){
  rec2<-filter(rec,ssb>0 &  stock==x)
  a2<-gam(rec~s(year)+scenario,data=rec2)
  data.frame(stock=rec2$stock,scenario=rec2$scenario,year=rec2$year,recpct=predict(a2)) %>% unique() %>%arrange(stock,scenario,year)
})
recruit<-do.call(rbind,b) %>% as_tibble()
recruit

pp<-ggplot(data = filter(recruit,year>=2020), aes(x = year, y = recpct,col=scenario)) + 
  geom_point(size = 1) +
  ylab('Percentage of unchanged climate recruitment')+
  facet_wrap(vars(stock))


png(file.path(output.path,'climate_rec_fac.png'),width = 600, height = 800, units = "px", pointsize = 12)
print(pp)
dev.off()


unique(recruit$scenario)
recruits<-rbind(recruit,filter(recruit,scenario=='rcp45') %>% mutate(scenario='noCC',recpct=100)) %>% mutate(recFac=recpct/100,recpct=NULL)

summary(recruits)
dd<-sort(unique(recruits$stock))
vpas<-SMS.control@species.names[16:27]
miss<- vpas[!(vpas %in% dd)]
recruits<-rbind(recruits,expand.grid(stock=miss,scenario=unique(recruits$scenario),year=unique(recruits$year),recFac=1))
recruits<-pivot_wider(recruits,names_from=stock,values_from=recFac) %>%
  select( scenario,year,Cod, Whiting, Haddock, Saithe, Mackerel, Herring, N.sandeel, S.sandeel, Nor.pout, Sprat,Plaice,Sole)

recruits<-recruits %>% filter(year>=2019)

aa<-recruits %>% filter(year==2029 & scenario=='rcp45')
unlist(aa[1,3:14])


# scenario F
scenF<-readxl::read_excel(file.path("H:","Seawise","scenarierTilAlex_2024","Summary.table_NS_fishingMort_stock_DeliverableMonth36.xlsx"), sheet = "F per stock")
sort(unique(scenF$stock))
ss<-c("COD-NS","WHG-NS","HAD","POK","PLE-NS","SOL-NS" )
sn<-c(16,17,18,19,26,27)
names(sn)<-ss
sn
scenF<-scenF %>% mutate(scenario=paste(HCR,ClimScen,sep='_'),Species.n=sn[stock],variable=NULL) %>% rename(FF=median) %>% dplyr::select(scenario,year,Species.n,stock,FF) %>% as_tibble()
sort(unique(scenF$scenario))

ar<-expand.grid(stock=unique(scenF$stock),year=unique(scenF$year),scenario=c('ICES-AR_noCC','ICES-AR_rcp45','ICES-AR_rcp85')) %>% mutate(Species.n=sn[stock])  %>% as_tibble()
targetF<-scan(file=file.path(data.path,"op_mulTargetf_FMSY.in"),skip=1L)
targetF<-data.frame(targetF=targetF,Species.n=first.VPA:nsp)
ar<-left_join(ar,targetF,by = "Species.n") %>% rename(FF=targetF) %>% tibble()
ar
scenF<-rbind(scenF,ar) %>% mutate(stock=sp.names[Species.n])




sort(unique(scenF$scenario))

if (FALSE) {  #not used in 2024
  ara<-gam(rec~s(year)+ssb+scenario,data=rec2)
summary(a)
plot(rec2$year,predict(a))
plot(a)


}
ftable(round(tapply(scenF$FF,list(scenF$scenario,scenF$year,scenF$stock),FUN=sum),3))
scenF

scenarios<-unique(scenF$scenario)
scenario.dirs<-paste("NS_2022_seawise_ver_2024",scenarios,sep='_')


# reset skill input files
#skill_fleet_make(control=SMS.control,make_excl=TRUE,write_data=TRUE,exclFile='skill_exclude.in') 
#makes the files skill_cpue.in and skill_exclude.in  


labels<-c("2020 Keyrun","2024 Seawise")  
dirs<-c("Seawise_NorthSeaKeyRun_2020","NS_2022_seawise_ver_2024")


if (FALSE) {
  
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
}

  if (FALSE) {
    compare_runs_stock_rec(dirs,labels,first.year.on.plot=1975,last.year.on.plot=2020,incl.sp="all",
                                     include.CV=TRUE,include.CV2=TRUE,include.mean=TRUE,
                                     palette="R3", makeAllGraphs=FALSE,nox=1, noy=2, w8=8,w11=8,
                                     include.year.labels=TRUE,incl_not_used=TRUE,run.ID='SBB_rec', 
                                     paper=TRUE,facSSB=1000,facRec=1000000,
                                     compare.dir=data.path,verbose=FALSE,writeData=TRUE) 
  } 
    compare_runs_stock_rec(dirs=dirs[2],labels=labels[2],first.year.on.plot=1975,last.year.on.plot=2020,incl.sp="all",
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



hind<-Read.summary.table() %>% filter(Year<=SMS.control@last.year.model)


unique(recruit$scenario)


#################    2.  ##############
do_scenario<-function(scenarioNo) {
 # test scenarioNo<-2  
  
  data.path<<-file.path(root,my.stock.dir)
  setwd(data.path)
  scenario<-scenarios[scenarioNo]

  scRec<-strsplit(scenario,'_')[[1]][2]
  #filter(recruits,year==2025)
         
  usedRec<-filter(recruits,scenario==scRec) %>% mutate(scenario=NULL) 
  
  scenario.dir.full<-file.path(root,scenario.dirs[scenarioNo])
  cat(scenario.dir.full,'\n')
  if (file.exists(scenario.dir.full)) unlink(scenario.dir.full,recursive = T)
  dir.create(scenario.dir.full,showWarnings = FALSE)
  
  OP.files<-c("sms.dat","op.exe","area_names.in","species_names.in","op.dat","op_trigger.dat","op_config.dat","op_msfd.dat","just_one.in",
              "op_consum.in","op_f.in","op_m1.in","op_m.in","op_n.in","op_propmat.in","op_prop_landed.in","op_size.in","op_wcatch.in","op_wsea.in",
              "op_growth_type1.in","op_consum_ab.in","op_other_n.in","op_exploitation.in","op_reference_points.in","covariance_rec.in","op_price.in",
              "op_ssb_rec_residuals.in","op_length_weight_relations.in",'op_eqsim.in','op_eqsim_stoch.in','op_n_proportion_m2.in','op_seed.in',
              "reference_points.in","summary.out")
  
  OP.files<-c("sms.dat","op.exe","area_names.in","species_names.in","op.dat","op_trigger.dat","op_config.dat","op_msfd.dat","just_one.in",
              "op_consum.in","op_f.in","op_m1.in","op_m.in","op_n.in","op_propmat.in","op_prop_landed.in","op_size.in","op_wcatch.in","op_wsea.in",
              "op_other_n.in","op_exploitation.in","op_reference_points.in","op_price.in",
              "op_ssb_rec_residuals.in","op_length_weight_relations.in",'op_n_proportion_m2.in','op_seed.in',
              "reference_points.in","summary.out")
  
  
  for (from.file in OP.files) {
    to.file<-file.path(scenario.dir.full,from.file)
    if (!file.copy(from.file, to.file, overwrite = TRUE)) cat('problems with file: ',from.file,'\n')
  }
  
   
  # 
  
  source(file=file.path(prog.path.func,'hcr_op_batch_common.R'))
  
  recfactor<-unlist(filter(usedRec,year==2020))[2:13]
  recfactor
  
  res<-make.OP.dat(my.area='North Sea',my.last.year=2020,first.year.output=2020,do.indicators=F,stochastic.recruitment=0,recruit.adjust.CV=0, recruit.adjust.factor=recfactor,    
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
  #                                               Cod Whiting Haddock Saithe Mackerel Herring N.sandeel S.sandeel Nor.pout Sprat Plaice Sole
  if (!grepl('ICES-AR',scenario)) opt@HCR[1,]<-c(   1,     1,       1,    1,       2,       2,      22,        22,      22,   22,      1,    1)
  if ( grepl('ICES-AR',scenario)) opt@HCR[1,]<-c(   2,     2,       2,    2,       2,       2,      22,        22,      22,   22,      2,    2) #ICES AR
  
  
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
   OP@recruit.adjust[1,]<-unlist(filter(usedRec,year==sYear))[2:13]
   
   cat(scenario,scRec,OP@recruit.adjust[1,],'\n')
   
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
#for (i in (1:2)) do_scenario(scenarioNo=i)
#for (i in (13:14)) do_scenario(scenarioNo=i)


 getwd()   
 
if (TRUE) {
 
 plot_it<-function(sce,sces,ncolLeg=1,fypl=1974) { 
  compare_runs(
   dirs=sces$dirs,  # directory files to be compared
   labels=sces$labels,  # legends
   ncollegd=ncolLeg,
   nox=2, noy=3,
   portrait=TRUE,
   paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
   run.ID=sce,         # file id used for paper output
   doGrid=TRUE,
   extent.SSB=FALSE,  # plot SSB for the year after last assessment year
   do_op_summary.out=TRUE,
   first.year.on.plot=fypl,
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
   a<-read_csv(file.path(data.path,"R_F_SSB.csv"),show_col_types = FALSE)
   write_csv(a,file.path(data.path,paste("R_F_SSB",sce,".csv",sep='_')))
 }
 
 
 allScen<- data.frame(dirs=scenario.dirs,labels=scenarios)

  #
  sce<-'All'
  sces<-allScen
  plot_it(sce,sces,ncolLeg=2)
  
  allScen<-filter(allScen, labels!= "PGY-Min_noCC")
 
  sce<-'All_short'
  sel<-grepl('ICES-AR',allScen$labels)
  sces<-allScen[!sel,] 
  plot_it(sce,sces,ncolLeg=2,fypl=2020)
  
  sce<-'All_short_all'
  sces<-allScen
  plot_it(sce,sces,ncolLeg=2,fypl=2020)
  
  #
  sce<-'ICES-AR'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020)

  sce<-'ICES-AR'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,]
  sce<-'ICES-AR-long'
  plot_it(sce,sces,fypl=1975)
  
  sce<-'PGY-Min'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020)
 
  sce<-'Status-quo'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020) 
  
  sce<-'FMSY-Min'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020) 
  
  
  sce<-'Case-study'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020) 
  
  sce<-'noCC'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020) 
  
  sce<-'rcp45'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020) 
  
  sce<-'rcp85'
  sel<-grepl(sce,allScen$labels)
  sces<-allScen[sel,] 
  plot_it(sce,sces,fypl=2020) 
}

 
 # a<-read_csv(file=file.path(data.path,'bernhard',"R_F_SSB_All_.csv"))
# a<-read_csv(file=file.path(data.path,"R_F_SSB.csv"))
 a<-read_csv(file=file.path(data.path,"R_F_SSB_All_.csv"))
 
 a<-a[a$variable=='Rec' & a$Species=='Cod',] %>%  mutate(scGroup=sapply(strsplit(scenario,'_'),function(x) x[1]),scClimate=sapply(strsplit(scenario,'_'),function(x) x[2]))
 
 arrange(unique(select(a,scenario,scGroup,scClimate)))
 
 ggplot(data = filter(a,value>0), aes(x = Year, y = value,col=scClimate)) + 
   geom_point(size = 1) +
   facet_wrap(vars(scGroup),scales='free_y')
 
##################### one third for the birds

 a<- Read.summary.data() %>% filter(Quarter==3) %>% 
   filter((Species=='Herring' & Age<=15) | (Species %in% c('S.sandeel','N.sandeel','Sprat'))) %>% 
   group_by(Species.n,Species,Year) %>% summarize(TSB=sum(BIO)) %>% ungroup()

 SSB.R.year.first<-SMS.control@SSB.R.year.first
 SSB.R.year.last <-SMS.control@SSB.R.year.last
 SSB.R.year.first[SSB.R.year.first==-1]<-SMS.control@first.year.model
 SSB.R.year.last[SSB.R.year.last==-1]<-SMS.control@last.year.model
 
 fl<-data.frame(Species.n=first.VPA:nsp,first=SSB.R.year.first,last=SSB.R.year.last)
 
 a<-left_join(a,fl,by = join_by(Species.n)) %>% filter(Year>=first & Year<=last) %>% mutate(TSB=TSB/1000) 

 oneThird<-a %>% group_by(Species.n,Species) %>% summarize(maxTSB=max(TSB),pct95=quantile(TSB, probs = c(.95))) %>%
   mutate(oneThird=maxTSB/3,OTpct95=pct95/3)
 a<-left_join(a,oneThird,by = join_by(Species, Species.n))
 
hist<-select(a, Species.n, Species , Year,   TSB)
 
save(a,hist,file=file.path('one_third.Rdata'))
p<- ggplot(data = a, aes(x=Year, y=TSB)) +
   geom_hline(aes(yintercept = oneThird), col="red",size=1.5) + 
   geom_hline(aes(yintercept = OTpct95), col="blue",size=1.5) + 
   geom_bar(stat = 'identity') +
   facet_wrap(.~ Species, ncol=2,scales = "free_y") + 
   ylab('Total biomas, 1. July (1000 tonnes)')+
   theme(axis.text.x=element_text(angle=90,
                                  hjust=1,
                                  vjust=.5,
                                  colour='gray50'))


png(file.path(output.path,'one_third_bird.png'),width = 900, height = 800, units = "px", pointsize = 12)
print(p)
dev.off()
###

## scenario

a<-read_csv(file=file.path(data.path,"R_F_SSB_All_.csv")) %>% 
 filter(Species %in% c('S.sandeel','N.sandeel','Sprat','Herring' ) & variable=='TSBq3' ) %>%
 mutate(scGroup=sapply(strsplit(scenario,'_'),function(x) x[1]),scClimate=sapply(strsplit(scenario,'_'),function(x) x[2]))


ah<-a %>%select(Species,scenario,scClimate,scGroup) %>% unique()
ah<-left_join(hist,ah,by = join_by(Species),relationship = "many-to-many")
ah<-ah %>%rename(value=TSB) %>% mutate(Species.n=NULL,variable="TSBq3")
a<-rbind(a,ah)
a
a<-left_join(a,oneThird,by = join_by(Species))

arrange(unique(select(a,scenario,scGroup,scClimate)))

aa<-a %>% mutate(scGroup=paste(Species,scGroup))
by(aa,list(aa$Species),function(x){
  p<-ggplot(data = x, aes(x = Year, y = value,col=scClimate)) + 
    geom_point(size = 1) +
    ylab('Total biomas, 1. July (1000 tonnes)')+
    geom_hline(aes(yintercept = OTpct95), col="blue",size=1.0) + 
    facet_wrap(vars(scGroup))
  
  
  png(file.path(output.path,paste0('one_third_bird_',x[1,"Species"],'.png')),width = 900, height = 800, units = "px", pointsize = 12)
  print(p)
  dev.off()

})



a<-read_csv(file=file.path(data.path,"R_F_SSB_All_short_all_.csv")) %>% 
  mutate(scGroup=sapply(strsplit(scenario,'_'),function(x) x[1]),scClimate=sapply(strsplit(scenario,'_'),function(x) x[2]))


b<-a %>% filter(Year>= 2051 & Year <=2060) %>% mutate(mGroup=if_else (grepl('ICES',scGroup),'Dynamic F','FLBEIA')) %>%
  group_by(scenario, Species, variable, scGroup , scClimate, mGroup) %>% summarize(value=mean(value),meanValue=mean(value)) %>% ungroup()


b2<-b %>%  group_by(Species, variable, scGroup) %>% summarize(minv=min(value),maxv=max(value),mean=mean(value)) %>% 
   ungroup() %>% mutate(range=paste(round(minv),round(maxv),sep=':')) %>% filter(Species!='Mackerel') %>%
  rename(Combination=scGroup)

b3<-b %>%  group_by(Species, variable, scClimate,) %>% summarize(minv=min(value),maxv=max(value),mean=mean(value)) %>% 
  ungroup() %>% mutate(range=paste(round(minv),round(maxv),sep=':')) %>% filter(Species!='Mackerel') %>%
  rename(Combination=scClimate)

b33<-b %>%  group_by(Species, variable,  mGroup,) %>% summarize(minv=min(value),maxv=max(value),mean=mean(value)) %>% 
  ungroup() %>% mutate(range=paste(round(minv),round(maxv),sep=':')) %>% filter(Species!='Mackerel') %>%
  rename(Combination=mGroup)


b4<-rbind(b2,b3,b33) 
sort(unique(b4$Combination))
summary(b4)

unique(b4$variable)
labs<-data.frame(scen=c("Eaten","Rec" ,  "SSB", "TSBq3","Yield",'F'),
           label=c("Eaten Biomass (1000 tonnes)",'Recruitment (billions)','SSB (1000 tonnes)', 'Total biomass (1000 tonnes)','Yield (1000 tonnes)','F')
)


for (i in (1:dim(labs)[[1]])) {
 b5<-filter(b4,variable==labs[i,'scen'])
p<-ggplot(b5,aes(factor(Combination,levels=c("ICES-AR" ,"Case-study", "FMSY-Min", "PGY-Min", "Status-quo","noCC","rcp45","rcp85","Dynamic F","FLBEIA" )), mean)) + 
  geom_linerange(aes(ymin=minv, ymax=maxv), size=1, color="blue") +
  geom_point(aes(x=factor(Combination,levels=c("ICES-AR" ,"Case-study", "FMSY-Min", "PGY-Min", "Status-quo","noCC","rcp45","rcp85","Dynamic F","FLBEIA" )), y=mean), size=4, shape=21, fill="white") +
  facet_wrap(vars(Species),scales = "free_y")+
  xlab('Scenario combination')+
  ylab(labs[i,'label'])+
  theme(axis.text.x=element_text(angle=90,
                                 hjust=1,
                                 vjust=.5,
                                 colour='gray50'))

png(file.path(output.path,paste0('comb_',labs[i,'scen'],'.png')),width = 700, height = 800, units = "px", pointsize = 20)
print(p)
dev.off()
}

