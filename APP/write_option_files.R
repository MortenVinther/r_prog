
write_options<-function(rec.mode=c(0,1),APP_final=TRUE,change_op_config=FALSE,make_AMOEBA_data=FALSE) {
  
  if ((APP_final | make_AMOEBA_data) & rec.mode==1) cat('Final data for the App will be produced (with determenistic recruitment)\n')
  if (APP_final | make_AMOEBA_data) rec.mode<-0
    
  SMS<-SMS.control
  HCR<-read.csv(file.path(data.path,'HCR_ini.csv'))
  
 
   yLong<-100
  yShort<-40
  
  ### OP.dat file
  
  nsp<-SMS@no.species
  n.other.pred<-sum(SMS@species.info[,'predator']==2)
  n.pred<-n.other.pred+sum(SMS@species.info[,'predator']==1)
  n.vpa<-nsp-n.other.pred
  n.vpa.pred<-sum(SMS@species.info[,'predator']==1)
  
  OP<-read.FLOP.control(file="op.dat",n.VPA=n.vpa,n.other.pred=n.other.pred,n.pred=n.pred) 
  
  OP@first.year <-SMS@last.year.model+1
  OP@first.year.out<-SMS@last.year.model+1
  OP@indicator<-0
  OP@stochastic.recruitment[1,]<-rec.mode 
  OP@F.or.C[1,]<-31
  
  if (rec.mode==1){ # stochastic
    OP@output<-26
    OP@last.year<-SMS@last.year.model+ yLong
    OP@rec.noise['lower',]<-HCR$noise.low
    OP@rec.noise['upper',]<-HCR$noise.high
    OP@recruit.adjust[1,]<-HCR$rec.adjust
    OP@recruit.adjust.CV[1,] <-0
    OP@recruit.min[1,]<-10
  } else if (rec.mode==0) {  #deterministic
    if (APP_final) OP@output<- 20 else  OP@output<- 26
    OP@last.year<-SMS@last.year.model+ yShort
    OP@rec.noise['lower',]<- -2
    OP@rec.noise['upper',]<-  2
    OP@recruit.adjust[1,]<-   HCR$rec.adjust.single
    OP@recruit.adjust.CV[1,]<-HCR$rec.adjust.CV.single
    OP@recruit.min[1,]<-10
  }
  my.last.year<<-OP@last.year
  
  if (make_AMOEBA_data & APP_final) {
    OP@first.year.out<-OP@last.year-1
    OP@output<-14
  }
  
  write.FLOP.control(OP,file="op_master.dat",path=out_op) 
  write.FLOP.control(OP,file="op.dat") 
  
  #### trigger
  trigger<-read.FLOPtrigger.control(file="op_trigger.dat",n.VPA=n.vpa,n.other.pred=n.other.pred) 
  trigger@first.run.no<-1
  trigger@first.iter.no<-1
  trigger@no.iter<-1
 
  trigger@first.year<-SMS@last.year.model+1 
  trigger@Ftarget['lower',]<-0
  trigger@Ftarget['higher',]<-3
  trigger@Ftarget['phase',]<-1
  trigger@Ftarget['init',]<-HCR$Ftarget    # trig@Ftarget is really not used !!!! (OP uses F from  op_multargetf.in as target F in simulation mode)
  
  trigger@at.age.weighting<-0
  
  trigger@trigger['T1',]<-HCR$T1
  trigger@trigger['T2',]<-HCR$T2
  
  
  if (rec.mode==1){ # stochastic
    trigger@last.year<-SMS@last.year.model+ yLong
    trigger@HCR[1,]<-HCR$HCR 

    # Ftarget is really not used !!!! (OP uses F from  op_multargetf.in as target F in simulation mode)
  } else if (rec.mode==0) {  #deterministic
       trigger@last.year<-SMS@last.year.model+ yShort
       trigger@first.year<-SMS@last.year.model+1
       trigger@HCR[1,]<-1
  }
  write.FLOPtrigger.control(trigger,file="op_trigger.dat")
  write.FLOPtrigger.control(trigger,file="op_trigger_master.dat",path=out_op)
  
  
  stq<-read.csv(file=file.path(out_op,'status_quo.csv'))
  
  
  
  if (rec.mode==0 & make_AMOEBA_data==FALSE ) cat('1\n',stq$mean.F,'\n',file='op_multargetf.in')
  if (rec.mode==0 & make_AMOEBA_data==TRUE ) cat('1\n',HCR$AMOEBA_F,'\n',file='op_multargetf.in')
  if (rec.mode==1) cat('1\n',HCR$Ftarget,'\n',file='op_multargetf.in')
  
  
  
  #recruitment parameters from op_config.dat
  if (change_op_config) {
    if (file.exists(file.path(data.path,"op_config - Copy.dat")) ) {
      rec<-readLines(file.path(data.path,"op_config - Copy.dat"))
      found<-grep("#model alfa  beta std info1 info2",rec)
      n<-found+9 # Norway pout
      np<-rec[n]
      np<-scan(text=np,quiet = TRUE, comment.char = "#")
      np
      np[2]<-np[2]*1.5 #peak at lower SSB
      np[3]<-np[3]*1.5
      rec[n]<-paste(np,collapse = " ")
      writeLines(rec,con=file.path(data.path,"op_config.dat"))
    } else stop("File 'op_config - Copy.dat' does not exist")
   }

}


