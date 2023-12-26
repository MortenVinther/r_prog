update_environment<-function(my.area) {
  if (my.area=='North Sea') {local_data_dir <- "Data_northsea" }
  if (my.area=='Baltic Sea'){local_data_dir<- "Data_baltic"}
  
 
  
  # control objects for predictions
  #data.path<-data_dir # used by the control objects
  
  SMS<-read.FLSMS.control(file='sms.dat',dir=out_op)
  #get information from SMS run
  n.species.tot <- SMS@no.species  # number of species including "other predators"
  n.pred.other<-sum(SMS@species.info[,'predator']==2) #number of "other predators"
  n.VPA<-n.species.tot-n.pred.other #number of species with analytical assessment
  n.pred<-n.pred.other+sum(SMS@species.info[,'predator']==1) # number of predators
  n.fleet <- n.VPA   # in this case it is just the number of species, one "fleet" per species
  stq_year<-SMS@last.year.model
  fy_year_hist<-SMS@first.year.model
  n.seasons<-SMS@last.season
  first.age<-SMS@first.age
  spNames<-SMS@species.names
  firstVPA<-n.species.tot-n.VPA+1   #first species with analytical assessment
  VPA.spNames<-spNames[firstVPA:n.species.tot]
  last.age<-SMS@species.info[VPA.spNames,'last-age']
  max.age<-SMS@max.age.all
  n.age<-max.age-first.age+1
  Fages<-SMS@avg.F.ages
  
  spOtherNames<-c('Other.food',spNames)
  other.spNames<-spNames[1:n.pred.other]
  
  
  # group predator names for eaten biomass
  pred_format<-read.csv(file.path(out_op,'pred_format.csv'),header=TRUE, stringsAsFactors = FALSE)
  #pp<-unique(subset(pred_format,new_no>=0,select=c(new,new_no,group,is.predator,is.prey)))
  pp<-unique(subset(pred_format,select=c(new,new_no,group)))
  pp<-pp[order(pp$new_no),]
  predPreyFormat<-pp$new
  pp_short<-filter(pp,new_no >=0) 
  
  recruitMode<-c('Determenistic','Stochastic')[1]
  
  # options for predictions, reset from master version
  OP<-read.FLOP.control(file="op_master.dat",path=out_op,n.VPA=n.VPA,n.other.pred=n.pred.other,n.pred=n.pred)
  OP.trigger<-read.FLOPtrigger.control(file="op_trigger_master.dat",path=out_op,n.VPA=n.VPA,n.other.pred=n.pred.other)
  
  
  # read units and label used for plots etc.
  a<-read_csv(file=file.path(out_op,'units.csv'),col_types = cols()) 
  plotUnits<-a$plotUnits; names(plotUnits)<-a$type
  plotLabels<-a$plotLabels; names(plotLabels)<-a$type
  plotLabels[is.na(plotLabels)]<-" "              
  roundUnits<-a$roundUnits; names(roundUnits)<-a$type
  
  
  # read status quo values (= values in the terminal hindcast SMS model year)  
  status_quo<-read_csv(file=file.path(out_op,'status_quo.csv'),col_types = cols()) %>%
    mutate(Rec=Rec*plotUnits['Recruits'],SSB=SSB*plotUnits['SSB'],TSB=TSB*plotUnits['TSB'],SOP=SOP*plotUnits['Yield'],
           Yield=Yield*plotUnits['Yield'],mean.F=mean.F*plotUnits['Fbar'],Eaten=Eaten*plotUnits['DeadM'])
  
  # read hindcast SMS values  
  histAnnoM<-read_csv(file=file.path(out_op,'hist_anno_M.csv'),col_types = cols()) %>%mutate(Species=spNames[Species.n])
  
  histCondensed<-read_csv(file=file.path(out_op,'hist_condensed.csv'),col_types = cols())%>%mutate(Species=spNames[Species.n]) %>%
    mutate(Recruits=Recruits*plotUnits['Recruits'],SSB=SSB*plotUnits['SSB'],TSB=TSB*plotUnits['TSB'],
           Yield=Yield*plotUnits['Yield'],yield.core=yield.core*plotUnits['Yield'],
           Fbar=Fbar*plotUnits['Fbar'],
           DeadM1=DeadM1*plotUnits['DeadM'],DeadM2=DeadM2*plotUnits['DeadM'],
           DeadM1.core=DeadM1.core*plotUnits['DeadM'],DeadM2.core=DeadM2.core*plotUnits['DeadM'])
  
  histEaten <- read_csv(file=file.path(out_op,'who_eats_whom_historical.csv'),col_types = cols()) %>%
    mutate(Predator=parse_factor(Predator,levels=predPreyFormat),Prey=parse_factor(Prey,levels=predPreyFormat),eatenW=eatenW*plotUnits['DeadM'])
  
  refPoints<-matrix(scan(file.path(out_op,"op_reference_points.in"),quiet = TRUE, comment.char = "#"),ncol=4,byrow=TRUE)
  rownames(refPoints)<-VPA.spNames
  colnames(refPoints)<-c('Flim','Fpa','Blim','Bpa')
  refPoints[,3:4] <- refPoints[,3:4]* plotUnits['SSB']
  
  explPat_in<-scan(file.path(out_op,"op_exploitation.in"),quiet = TRUE, comment.char = "#")
  explPat<-array(explPat_in,dim=c(n.age,n.VPA,n.seasons),dimnames=list(paste('age',as.character(first.age:max.age)),VPA.spNames,paste0('Q_',1:n.seasons)))
  annExplPat<-apply(explPat,c(1,2),sum) # annual exploitation pattern
  
  # SMS output values in the last year
  base_SSB<-stqSSB<-status_quo$SSB
  base_F<-stqF<-status_quo$mean.F
  #stqF<-status_quo$sum.q.F
  base_Yield<-stqYield<-status_quo$Yield
  base_Rec<-stqRec<-status_quo$Rec
  
  # write status quo F
  #cat("1\n",base_F,"\n",file=file.path(local_data_dir,"op_multargetf.in")) # write F values
  
  # read various setting for options files
  hcr_ini<-read.csv(file.path(out_op,'HCR_ini.csv'),header=TRUE,stringsAsFactors = FALSE)
  
  ### change option values (could have been done in the master files!)
  OP@rec.noise['lower',]<-hcr_ini$noise.low
  OP@rec.noise['upper',]<-hcr_ini$noise.high
  OP@recruit.adjust.CV[1,]<-hcr_ini$rec.adjust.CV.single
  OP@recruit.adjust[1,]<-hcr_ini$rec.adjust.single
  #
  OP.trigger@Ftarget['init',]<-base_F
  OP.trigger@trigger['T1',]<-hcr_ini$T1
  OP.trigger@trigger['T2',]<-hcr_ini$T2
  OP.trigger@HCR[1,]<-1
  
  # write option files to be used
  write.FLOP.control(OP,file=file.path(out_op,"op.dat"),nice=TRUE)
  write.FLOPtrigger.control(OP.trigger,file="op_trigger.dat",path=out_op, nice=TRUE)
  
  
  fleetNames<-paste0('fl_',VPA.spNames) # just this special case where we have no fleets
  
  #recruitment parameters from op_config.dat
  rec<-readLines(file.path(out_op,"op_config.dat"))
  found<-grep("#model alfa  beta std info1 info2",rec)
  rec<-scan(file.path(out_op,"op_config.dat"),skip=found,comment.char = "#",nlines=n.VPA,quiet = TRUE)
  rec<-matrix(rec,nrow=n.VPA,byrow=TRUE)
  colnames(rec)<-c('model','a','b','s','o1','o2')
  rownames(rec)<-VPA.spNames
  
  
  # maximum recruits
  max_rec<-rep(0,n.VPA);names(max_rec)<-VPA.spNames
  i<-rec[,'model']==100;max_rec[i]<-exp(rec[i,'a'])*rec[i,'b']  # Hockey stick
  i<-rec[,'model']==1;  max_rec[i]<-rec[i,'a']/(rec[i,'b']*exp(1))  # Ricker
  i<-rec[,'model']==2;  max_rec[i]<-rec[i,'a']/rec[i,'b']  # B & H
  i<-rec[,'model']==3;  max_rec[i]<-exp(rec[i,'a'])  # GM
  i<-OP@recruit.adjust.CV==2; max_rec[i]<-max_rec[i]*exp((rec[i,'s']^2)/2)
  max_rec<-max_rec*plotUnits['Recruits']
  
  # values for baselines option lists
  bsF<-list(Names=list('No change',paste0('F(',stq_year,')'),'Most recent results'),Values=list(0,1,2))
  bsSSB<-list(Names=list('No change',paste0('SSB(',stq_year,')'),'Most recent results'),Values=list(0,1,2))
  bsYield<-list(Names=list('No change',paste0('Yield(',stq_year,')'),'Most recent results'),Values=list(0,1,2))
  bsRec<-list(Names=list('No change',paste0('Recruitment(',stq_year,')'),'Most recent results','Maximum recrutiment'),Values=list(0,1,2,3))
  
  
  F_mult<-1.0  
  
  Fvalues<-stqF*F_mult
  
  oldwd<-getwd()
  
   get_terminal_year<-function(OP){
    return(OP@last.year)
  }
  
  termYear<-get_terminal_year(OP)
  
  
  get_other_predators<-function(){
    First.year<-rep(1.0,length(other.spNames)) ;names(First.year)<-other.spNames
    Last.year<-Total.change<-First.year
    
    for (sp in other.spNames)  {
      if (OP@other.predator['first',sp]== -1)   First.year[sp]<-stq_year+1 else First.year[sp] <- OP@other.predator['first',sp]
      if (OP@other.predator['second',sp]== -1)  Last.year[sp]<-OP@other.predator['second',sp]<- termYear else  Last.year[sp] <- OP@other.predator['second',sp]
      if ((OP@other.predator['first',sp]== -1) || (OP@other.predator['second',sp]== -1)) Total.change[sp]<- OP@other.predator['factor',sp]** (OP@other.predator['second',sp]-OP@other.predator['first',sp]+1) else Total.change[sp]<-1
    }
    
    (data.frame(Predator=other.spNames,
                change=OP@other.predator['factor',], 
                First.year=as.integer(First.year), 
                Last.year=as.integer(Last.year),
                Total.change=Total.change , stringsAsFactors = FALSE))
  }
  other_predators<-get_other_predators()
  
  
  
  hcrlab <- c("Fixed F", "F from SSB", "F from TSB")
  hcr <- data.frame(val = hcrlab)
  hcrval<-c(1,2,22);names(hcrval)<-hcrlab 
  
  
  get_op_Fmodel<-function(){
    HCR<-OP.trigger@HCR
    trigger<-OP.trigger@trigger*plotUnits['SSB']
    Ftarget<-OP.trigger@Ftarget['init',]
    return(data.frame(Species=VPA.spNames,target.F=Ftarget, HCR=names(hcrval[match(HCR[1,],hcrval)]),T1=trigger[1,],T2=trigger[2,],stringsAsFactors = FALSE))
  }
  
  
  Foption_tab<-get_op_Fmodel()
  
  
  doRunModel<-TRUE  # flag for re-running the prediction model
  doWriteOptions<-TRUE  # flag for writing option files for the prediction model
  doWriteExplPattern<- FALSE # flag for writing exploitation pattern file (op_exploitation.in)
  
  # icons for HCR options
  hcr$img <-c(
    sprintf("<img src='fixed_F.png' width=100px><div class='jhr'>%s</div></img>", hcr$val[1]),
    sprintf("<img src='AR_F_SSB.png' width=100px><div class='jhr'>%s</div></img>", hcr$val[2]),
    sprintf("<img src='AR_F_TSB.png' width=100px><div class='jhr'>%s</div></img>", hcr$val[3])
  )
  oldFvals<-rep(1.0,n.fleet)
  
  
  library(shiny)
  doSlider<-function(fleets,ars){
     slider <- div()
    for (i in (1:length(fleets))) {
      slider <- tagAppendChild(slider, sliderInput(inputId =paste0(ars,'_',fleets[i]),
                                                   label =  fleets[i],
                                                   min = 0.25, max = 4.0, value = 1.0, step = 0.05))
    }
    return(slider)
  }

  fl_baltic<-c("fl_Herring", "fl_Sprat")
  fl_northSea<-c("fl_Cod","fl_Whiting","fl_Haddock","fl_Saithe","fl_Mackerel","fl_Herring","fl_N.sandeel","fl_S.sandeel","fl_Nor.pout","fl_Sprat","fl_Plaice","fl_Sole") 
  
  cat("fleetNames:  ",fleetNames,'\n')
  if (my.area=='North Sea')  stopifnot(all(fl_northSea==fleetNames))
  if (my.area=='Baltic Sea')  stopifnot(all(fl_baltic==fleetNames))
  
  slidersNS<-doSlider(fl_northSea,ars='NS')
  slidersBS<-doSlider(fl_baltic,ars='BS')
  
  
  selectedAreas<-c('North Sea','Baltic Sea')
  selectedArea<-selectedAreas[1]
  

  

  
  save(
    local_data_dir,
    SMS,
   # n.species.tot,
   # n.pred.other,
    n.VPA,
    n.pred,
    n.fleet, 
    stq_year,
    fy_year_hist,
    n.seasons,
    first.age,
    spNames,
    firstVPA,
    VPA.spNames,
    last.age,
    max.age,
    n.age,
    Fages,
    spOtherNames,
    other.spNames,
    pred_format,
  # pp,
    predPreyFormat,
    pp_short,
    recruitMode,
    OP.trigger,
    plotUnits,
    plotLabels,
    roundUnits,
    status_quo,
    histAnnoM,
    histCondensed,
    histEaten, 
    refPoints,
    explPat,
    annExplPat,
    base_SSB,
    stqSSB,
    base_F,
    stqF,
    base_Yield,
    stqYield,
    base_Rec,
    stqRec,
    hcr_ini,
    fleetNames,
    rec,
    max_rec,
    bsF,
    bsSSB,
    bsYield,
    bsRec,
    F_mult, 
    Fvalues,
   # oldwd,
    termYear,
    other_predators,
    hcrlab, 
    hcr, 
    hcrval,
    Foption_tab,
    doRunModel,
    doWriteOptions,
    doWriteExplPattern,
    hcr,
    oldFvals,
    OP.trigger,
    slidersNS,
    slidersBS,
    selectedAreas,
    OP,file=file.path(out_op,paste0('environment.Rdata'))
  )    
  
  
  
  return('done')
}  #end update_environment
