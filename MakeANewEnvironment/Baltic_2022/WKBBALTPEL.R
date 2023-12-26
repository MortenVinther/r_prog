
library(tidyverse)
library(reshape2)
exchangeDir<-"Data_baltic/2022-data"   # directory with data  files on spread sheet format

source(file.path(prog.path,"compare_WKBBALTPE.R"))

runs<-list(lim_10=0.1,lim_25=0.25,key_run=0.5,lim_75=0.75,lim_90=0.90)

do_data<-FALSE
do_runs<-FALSE


id<-"WKBBALTPEL"

if (do_data) {
  cod<-read.csv(file.path(root,exchangeDir,'Cod',"cod_length.csv"))
  cod<-as_tibble(cod)
  tst<-filter(cod,year %in% c(1974,1990,2000,2021))
  cat("Cod length\n");ftable(round(tapply(tst$mean_l,list(tst$year,tst$quarter,tst$age),sum),3))
  
  
  
  #consumption
  conf<-read.csv(file.path(root,exchangeDir,'Cod',"Average QUARTERLY consum param.csv"))
  conf$species<-'COD'
  conf<-as_tibble(conf)
  
  
  do_ab<-function(conf,lim=0.05) {
   out<-conf %>% mutate(error_a= qnorm(lim)*conf$a_std, error_b= qnorm(lim)*conf$b_std) %>%
        mutate(a= a+error_a,b=b+error_b) 
    #return(dplyr::select(out,-error_a,-error_b))
   return(out)
  }
  
  
  confs<-lapply(runs,function(x)do_ab(conf,lim=x))
  names(confs)
  confs
  
  if (TRUE) { # just testing
    calc_consum<-function(l=10:60,r=1,conf){
      a<-as.numeric(conf[r,'a']); b<-as.numeric(conf[r,'b'])
      consum<-a*l**b
      names(consum)<-as.character(l)
      return(consum)
    }
    
    cons<-lapply(confs,function(x) calc_consum(conf=x))
    cons
    png(file.path(data.path,'gem','food_cons.png'))
    plot(x=as.numeric(names(cons[['key_run']])) ,y=cons[['key_run']],t='l',col='black',ylim=c(0,1800),lwd=4,xlab='Length (cm)',ylab='Food ration (kg/quarter)')
    lapply(cons ,function(x) lines(as.numeric(names(x)),y=x,col='green',lwd=2))
    lines(x=as.numeric(names(cons[['key_run']])) ,y=cons[['key_run']],col='black',lwd=4)
    cleanup()
  }
  
  
  prop<-read.csv(file.path(root,exchangeDir,'Cod',"Quarterly_consumption_proportions.csv"))
  prop<-melt(prop,id=1:4)
  prop$quarter<-as.numeric(substr(prop$variable,2,2)); 
  prop<-mutate(prop,variable=NULL,species='COD')
  
  
  do_cons<-function(l=cod,conf) {
    a<-left_join(l,conf) %>% filter(year>=year.start & year<=year.stop) %>% mutate(CONSUM=a*(mean_l/10)**b /1000,SMS_area=1 )
    ap<-left_join(dplyr::select(a,c(-a,-b,-a_std, -b_std)),prop) %>% filter(mean_l >=l.start & mean_l<=l.stop) %>%
     mutate(l.start=NULL,l.stop=NULL, CONSUM=CONSUM*4*value)
    ap0<-subset(ap,age==1) %>% mutate(age=0L,CONSUM=0)
    ap<-rbind(ap0,ap)
   return(ap)
  }
  
  cons<-lapply(confs,function(x)do_cons(conf=x))
  
  
  names(cons)
  
  tst<-filter(cons[['key_run']],year %in% c(1974,1990,2000,2021))
  cat("Consumption\n");ftable(round(tapply(tst$CONSUM,list(tst$year,tst$quarter,tst$age),sum),3))
  
  write.cons<-function(x,fname) {
    apOut<-x %>% dplyr::select(year,species,quarter,age,CONSUM) %>% mutate(CONSUM=round(CONSUM,3)) %>%
      pivot_wider(names_from=age,values_from=CONSUM) %>% mutate(comment=paste(" #",year,quarter)) %>%
      dplyr::select(-year,-species,-quarter)
    ofile<-file.path(data.path,paste0('Cod_',fname,"_consum.in"))
    write.table(apOut,row.names=FALSE,col.names=FALSE,file=ofile,quote = FALSE)
    cat("-999\n",file=ofile,append=TRUE)
  }  
   
  for (i in (1:length(cons))) {
    write.cons(x=cons[[i]],fname=names(cons[i]))
  }
}



if (do_runs) {
  lapply(names(cons), function(x){
    # test x<-names(cons)[1]
    outputDir<-file.path(root,paste0('Baltic-2022_',id,"_",x))
    if (file.exists(outputDir)) unlink(outputDir,recursive = T)
    dir.create(outputDir,showWarnings = FALSE)
  
    for (from.file in list.files(path=data.path,pattern='*.in')) {
      to.file<-file.path(outputDir,from.file)
      print(paste(file.copy(from.file, to.file, overwrite = TRUE),from.file))
    }
    for (from.file in list.files(path=data.path,pattern='*.dat')) {
      to.file<-file.path(outputDir,from.file)
      print(paste(file.copy(from.file, to.file, overwrite = TRUE),from.file))
    }
    file.copy(file.path(root,'program','sms.exe'), file.path(outputDir,"sms.exe"), overwrite = TRUE)
    
    file.copy(file.path(data.path,paste0('Cod_',x,"_consum.in")), file.path(outputDir,"consum.in"), overwrite = TRUE)
    
  })
  
  
  
  deleteFiles<-c('*.png','Her_*.in','Cod_*.in')
  
  lapply(names(cons), function(x){
    # test x<-names(cons)[1]
   outputDir<-file.path(root,paste0('Baltic-2022_',id,"_",x))
   do.a.full.SMS.run(label="run_",                   # label for output
                     rundir=outputDir,outdir=outputDir,
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
                    Screen.show=T,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                    do.run=T,                       # Make the run immediately, or just make the batch file for the run
                    deleteFiles=deleteFiles,        # clean up in files before the run is made
                    HPC=F)                          # run it as batch program on the UNIX High Performance Computer 
  })
 
} 


labels<-names(runs)
dirs<-paste0('Baltic-2022_',id,"_",labels)

compare_WKBBALTPEL(labels,dirs,IDout='consum')


################################  M1 values ##################################

#runsM1<-list(M1_005=c(0.04,0.06)/4, M1_010=c(0.08,0.12)/4,M1_key_run=c(0.1,0.1)/4,M1_020=c(0.17,0.23)/4,M1_020_uniform=c(0.20,0.20)/4)

runsM1<-list(M1_005=c(0.04,0.06)/4, M1_010=c(0.08,0.12)/4,M1_key_run=c(0.1,0.1)/4,M1_020=c(0.17,0.23)/4)

if (do_data) {
  M1<-expand.grid(Species.n=2:3,Year=1974:2021,Quarter=1:4, Age=0:11,M1=0.1)
  M1[M1$Species.n==3,'M1']=0.2/4  #Sprat
  M1$first<-TRUE
  M1[M1$Year>=2000,'first']<-FALSE
  head(M1)
  tail(M1) 
  
  m1s<-lapply(runsM1,function(x){
    m1<-M1
    m1[m1$Species.n==2 &  m1$first,'M1'] <-x[1]
    m1[m1$Species.n==2 & !m1$first,'M1'] <-x[2]
    return(m1)
  })             
  
  
  write.m1s<-function(x,fname) {
    # test x<-m1s[[1]]
    mOut<-as_tibble(x) %>% dplyr::select(Year,Species.n,Quarter,Age,M1) %>% mutate(M1=round(M1,5)) %>%
      arrange(Species.n,Year,Quarter,Age) %>% 
      pivot_wider(names_from=Age,values_from=M1) %>% mutate(comment=paste(" #",sp.names[Species.n], Year,Quarter)) %>%
      dplyr::select(-Year,-Species.n,-Quarter)
    ofile<-file.path(data.path,paste0('Her_',fname,".in"))
    
    write.table(mOut,row.names=FALSE,col.names=FALSE,file=ofile,quote = FALSE)
    cat("-999\n",file=ofile,append=TRUE)
  }  
  
  for (i in (1:length(m1s))) {
    write.m1s(x=m1s[[i]],fname=names(m1s[i]))
  }
}


if (do_runs) {
  lapply(names(m1s), function(x){
    # test x<-names(cons)[1]
    outputDir<-file.path(root,paste0('Baltic-2022_',id,"_",x))
    if (file.exists(outputDir)) unlink(outputDir,recursive = T)
    dir.create(outputDir,showWarnings = FALSE)
    
    for (from.file in list.files(path=data.path,pattern='*.in')) {
      to.file<-file.path(outputDir,from.file)
      print(paste(file.copy(from.file, to.file, overwrite = TRUE),from.file))
    }
    for (from.file in list.files(path=data.path,pattern='*.dat')) {
      to.file<-file.path(outputDir,from.file)
      print(paste(file.copy(from.file, to.file, overwrite = TRUE),from.file))
    }
    file.copy(file.path(root,'program','sms.exe'), file.path(outputDir,"sms.exe"), overwrite = TRUE)
    
    file.copy(file.path(data.path, paste0('Her_',x,".in")), file.path(outputDir,"natmor1.in"), overwrite = TRUE)

  })
  
  
  
  deleteFiles<-c('*.png','Her_*.in','Cod_*.in')
  
  lapply(names(m1s), function(x){
    # test x<-names(cons)[1]
    outputDir<-file.path(root,paste0('Baltic-2022_',id,"_",x))
    do.a.full.SMS.run(label="run_",                   # label for output
                      rundir=outputDir,outdir=outputDir,
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
                      Screen.show=T,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                      do.run=T,                       # Make the run immediately, or just make the batch file for the run
                      deleteFiles=deleteFiles,        # clean up in files before the run is made
                      HPC=F)                          # run it as batch program on the UNIX High Performance Computer 
  })
  
} 

labels<-names(runsM1)
dirs<-paste0('Baltic-2022_',id,"_",labels)

compare_WKBBALTPEL(labels,dirs,IDout='M1')


###  all

labels<-c(names(runs),names(runsM1))
labels<-labels[labels!="M1_key_run"]
dirs<-paste0('Baltic-2022_',id,"_",labels)
compare_WKBBALTPEL(labels,dirs,IDout='all')


### output

source(file.path(prog.path,"r_prog_less_frequently_used",'plot_Z.R'))

for (i in (1:length(dirs))) {
  cat(dirs[i],labels[i],'\n')
  do_Z(data_in=file.path(root,dirs[i]),label=' ',  first_cohort=1973, last_cohort=2014,inclSPn=2:3,age_ranges=NULL) 
}

a<-lapply(dirs,function(x) {
  read_csv(file.path(root,x,'WKBBALTPEL.csv'))
})

b<-do.call(rbind,a)
write_csv(b,file='WKBBALTPEL_all.csv')
