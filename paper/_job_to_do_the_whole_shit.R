# first run init.r as usual

useNbar<-TRUE
UseLengthSMS<-FALSE


reRunMulti<-FALSE
reRun<-FALSE


if (UseLengthSMS) {
  if (useNbar) data.path.multi<-file.path(root,"NS_2020_multi_len")  else data.path.multi<-file.path(root,"NS_2020_multi_no_Nbar_len")
  
} else {
 if (useNbar) data.path.multi<-file.path(root,"NS_2020_multi")  else data.path.multi<-file.path(root,"NS_2020_multi_no_Nbar")
}
data.path.multi



# just two functions
transf<-function(a,file.name,dig) {
  out<-file.path(data.path,file.name)
  unlink(out)
  years<-  dimnames(a)[[2]]
  for (sp in dimnames(a)[[1]]) {
    cat(sp,'\n')
    out1<-a[sp,,,]
    for (y in years){
      out2<-out1[y,,]
      out2[is.na(out2)]<- -1
      cat(paste("#",sp.names[as.numeric(sp)] ,"year:",y,"\n"),file=out,append=TRUE)
      write.table(format(round(out2,dig),width=11),file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)    
    }
  }
  cat(-999,file=out,append=TRUE)
}

# gradient function, in case of poor fit
gr_check<-function(grfile="run_ms0.grd",dir=data.path) {
  g<-Read.gradient.dat(gfile=grfile,dir=dir)
  g$parNo<-1:dim(g)[[1]]
  g<-g[order(abs(g$Gradient)),]
  print(tail(g,10))
}  
library(mgcv)

smoothM<-function(x,doPlot=FALSE) {
  cat(x$Species[1],'Age:',x$Age[1],'Q:',x$Quarter[1],'\n')
  if(min(x$MM) != max(x$MM)){
    y<-gam(MM ~s(Year,bs='cs'),data=x)
    yy<-predict(y)
    if (doPlot) {
      plot(x$Year,x$MM,type='b',main=paste(x$Species[1],x$Age[1],x$Quarter[1]))
      lines(x=x$Year,y=yy,col='red')
    }
  } else {
    yy<-x$MM
  }
  return(data.frame(Species=x$Species,Species.n=x$Species.n,Age=x$Age, Year=x$Year,Quarter=x$Quarter, MM=x$MM,predM=yy))
}

#########
# Path to data directory
data.path<-data.path.multi
setwd(data.path)

# run to  get multi species M
if (reRunMulti) {
  
  # delete files before the first run is made, only in force if do.single=T
  
  deleteFiles<-c("*.?0?","*.out","*.mcm","*.bin","admodel.*","*.csv","*.std","*.bar","*.mc2","*.cor","*.psv","*.ecm",
                 "*.xls","*.html", "mcout*.all","_*.txt","OP_BATCH*.*","SMS.o*","gradient.dat","*.grd","*.Rdata",
                 "*.wmf","*.png","*.ps","*.lg","*.log","ud.dat","gradient.dat","op_*.out","iter*.*","stom_and_noise*.in","canum_and_noise*","survey_and_noise*","baseline*.*",
                 "HCR_prob.dat","HCR_yield.dat","HCR_SSB.dat","*.par","*.rep","*.hst","*.eva","*.tmp","amoeba*.*","covariance_*.*","forecast*.*")

  do.a.full.SMS.run(label="run_",                   # label for output
                    cleanup=T,                      # delete files in the deleteFiles variable?
                    do.single=T,                    # run SMS in single species mode
                    do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                    do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                    do.multi.2.redo=!useNbar,              # Run the full model, with simultaneously estimation of all parameters
                    do.multi.2.redo.Nbar=useNbar,    # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                    do.hessian=T,                   # Make the Hessian matrix and estimate uncertainties
                    shake.ms2.par=F,
                    SSB.R.seperate=F,               # Estimate S/R parameters in a separate step  
                    do.MCMC=F,                      # Prepare for MCMC analysis
                    mcmc=0,mcsave=0,                # Options for MCMC analysis
                    do.prediction=F,                # Make a prediction
                    pause=F,                        # Make a confirm between each stage
                    Screen.show=T,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                    do.run=T,                       # Make the run immediately, or just make the batch file for the run
                    deleteFiles=deleteFiles,        # clean up in files before the run is made
                    HPC=F)                          # run it as batch program on the UNIX High  Performance Computer 


  if (useNbar) gr_check(grfile="run_ms3Nbar.grd") else gr_check(grfile="run_ms3.grd")
  
  SMS<-read.FLSMS.control()
  
  #read M

  dat<-Read.summary.data(read.init.function=F)
  dat<-subset(dat,Species.n>=first.VPA)
  dat$MM<-dat$M1+dat$M2
  save(dat,file=file.path(data.path,'in_dat.Rdata'))
  
  #re-format
  a<-tapply(dat$MM,list(dat$Species.n,dat$Year,dat$Quarter,dat$Age),sum)
  dimnames(a)
  a['16','1974',,]
  Mmulti<-a
  save(Mmulti,dat,file=file.path(data.path,'dat_M_multi.Rdata'))
  
  
  # smoothed version
  bb<-by(dat,list(dat$Age,dat$Quarter,dat$Species.n,dat$Species),smoothM,simplify=FALSE )
  dat<-do.call(rbind,bb)
  a<-tapply(dat$predM,list(dat$Species.n,dat$Year,dat$Quarter,dat$Age),sum) #re-format
  a['16','1974',,]
  Mmulti<-a
  save(Mmulti,dat,file=file.path(data.path,'dat_M_smooth.Rdata'))
  
  
  
  # Average M
  load(file=file.path(data.path,'dat_M_multi.Rdata'),verbose=T)
  a<-apply(Mmulti,c(1,3,4),mean,na.rm=TRUE)
  for ( s in dimnames(Mmulti)[[1]]) for (y in dimnames(Mmulti)[[2]]) Mmulti[s,y,,]<-a[s,,]
  Mmulti['16','1974',,]
  save(Mmulti,file=file.path(data.path,'dat_M_average.Rdata'))
  
  # M=0.2 for multispecies sp
  load(file=file.path(data.path,'dat_M_multi.Rdata'),verbose=T)
  Mmulti['16','1974',,]
  dimnames(Mmulti)
  for (s in c('16','17','18','21','22','23','24','25')) Mmulti[s,,,]<-0.05
  save(Mmulti,file=file.path(data.path,'dat_M_02.Rdata'))

}

###### establish  a new directory for single species SMS, Copy MS M 
scenario.dirs<-  c("NS_2020_s_01_MS_M","NS_2020_s_02_smooth_MS_M","NS_2020_s_03_avg_MS_M","NS_2020_s_04_M02")
scenario.labels<-c("s_01_MS_M","s_02_smooth_MS_M","s_03_avg_MS_M","s_04_M=02")

multi.dir<-c(data.path.multi)
setwd(multi.dir)
SMS<-read.FLSMS.control()
              
if (reRun) for (i in (1:length(scenario.dirs))) {
    setwd(multi.dir)
    copySMSfiles(control=SMS,file.path("..",scenario.dirs[i]),doSingle=TRUE,doMulti=FALSE,doArea=FALSE) 
      
    my.stock.dir<-scenario.dirs[i]
    
    # Path to data directory
    data.path<-file.path(root,scenario.dirs[i])
    setwd(data.path)
    
    if (scenario.dirs[i]=="NS_2020_s_01_MS_M")        {load(file=file.path(multi.dir,'dat_M_multi.Rdata')) ;transf(Mmulti,file.name="natmor.in",dig=5) }
    if (scenario.dirs[i]=="NS_2020_s_02_smooth_MS_M") {load(file=file.path(multi.dir,'dat_M_smooth.Rdata')) ;transf(Mmulti,file.name="natmor.in",dig=5) }
    if (scenario.dirs[i]=="NS_2020_s_03_avg_MS_M")    {load(file=file.path(multi.dir,'dat_M_average.Rdata')) ;transf(Mmulti,file.name="natmor.in",dig=5) }
    if (scenario.dirs[i]=="NS_2020_s_04_M02")         {load(file=file.path(multi.dir,'dat_M_02.Rdata')) ;transf(Mmulti,file.name="natmor.in",dig=5) }
    
    # run in single species mode (with MS M)
    do.a.full.SMS.run(label="run_",                   # label for output
                                 cleanup=F,                      # delete files in the deleteFiles variable?
                                 do.single=T,                    # run SMS in single species mode
                                 SSB.R.seperate=T,               # Estimate S/R parameters in a separate step 2 
                                 do.hessian=T,                   # Make the Hessian matrix and estimate uncertainties
                                 do.run=T)                       # Make the run immediately, or just make the batch file for the run
    
  # gr_check(dir=data.path,grfile="run_ms0.grd")
}


##  compare runs

scenario.dir<-data.path.multi
my.stock.dir<- scenario.dir

# Path to data directory
data.path<-file.path(data.path.multi)
setwd(data.path)

getwd()

# check, some runs might have a too large gradient, re-run with estimated parameters
for (dd in scenario.dirs) { print(dd);gr_check(grfile="run_ms0.grd",dir=file.path(root,dd))}


dirs<-basic_dirs<-c("NorthSeaKeyRun_2020","NS_2020_multi_no_Nbar","NS_2020_multi","NS_2020_multi_no_Nbar_len",scenario.dirs);
labels<-basic_labels<-labels<-c('1. key run',"2. MS N","3. MS Nbar","4. MS Length",scenario.labels)
data.frame(dirs,labels)


source(file.path(prog.path,"compare_runs_objective_function.R"))


a<-read_csv(file="compare_like_species.csv") %>% select(label,Species,Species.n,neg_like) %>% filter(Species.n>=16 & Species.n<=25)
a
a<-spread(a, key =label, value = neg_like) %>% arrange(Species.n)
write_csv(a,file='species_like.csv')
##

Read.objective.function.fleet()

mySpecies<-'Cod'; myAge<-1; myQuarter<-'1'


load(file=file.path(multi.dir,'dat_M_smooth.Rdata'),verbose=TRUE) 
dat<-subset(dat,Species==mySpecies & Age==myAge & Quarter==myQuarter)

png(file=file.path(data.path.multi,paste0('compare',mySpecies,'age_',myAge, 'Q_',myQuarter,'.png')), width = 1000, height = 1000, units = "px", pointsize = 20)
plot(dat$Year,dat$MM,ylab='M1+M2',xlab='',type='b',ylim=c(0,max(dat$MM)),main=paste(mySpecies,'age:',myAge, 'Q:',myQuarter),lwd=2)
lines(dat$Year,dat$predM,col='red',lwd=2)
abline(a=mean(dat$MM),b=0,col='green',lwd=2)
abline(a=0.05,b=0,col='blue',lwd=2)
cleanup()


# multi species
dirs<-basic_dirs[1:3]
labels<-basic_labels[1:3]
data.frame(dirs,labels)

compare_runs(
  dirs=dirs,
  labels=labels,
  nox=2, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID='multi_sp',         # file id used for paper output
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
  incl.sp="all",
  first.pch=0,    # first pch symbol
  first.color=1,   # first color
  palette="default"               # good for colour full plots
  #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
)  


# single species
dirs<-basic_dirs[4:7]
labels<-basic_labels[4:7]
data.frame(dirs,labels)

compare_runs(
  dirs=dirs,
  labels=labels,
  nox=2, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID='single_sp',         # file id used for paper output
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
  incl.sp="all",
  first.pch=0,    # first pch symbol
  first.color=1,   # first color
  palette="default"               # good for colour full plots
  #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
)  




# multi species
dirs<-basic_dirs[1:4]
labels<-basic_labels[1:4]
data.frame(dirs,labels)

compare_runs_stock_rec(dirs,labels,first.year.on.plot=1975,last.year.on.plot=2020,incl.sp="all",include.CV=TRUE,include.CV2=TRUE,include.mean=TRUE,
                       nox=2, noy=2, run.ID='SBB_rec_multi', include.year.labels=TRUE,incl_not_used=TRUE,paper=FALSE) 



# single species
dirs<-basic_dirs[5:8]
labels<-basic_labels[5:8]
data.frame(dirs,labels)


compare_runs_stock_rec(dirs,labels,first.year.on.plot=1975,last.year.on.plot=2020,incl.sp="all",include.CV=TRUE,include.CV2=TRUE,include.mean=TRUE,
                       nox=2, noy=2, run.ID='SBB_rec_single', include.year.labels=TRUE,incl_not_used=T,paper=TRUE) 




#residuals


for (dir in dirs) {
  if ( file.access(file.path(root,dir,"sms.dat"), mode = 0)!=0)  stop(paste('Directory',dir,'does not exist')) #else cat(dir,'\n')
} 


read.survey.residuals<-function(data.path,standardize=F) {
  
  fleet.names<-Read.fleet.names()
  
  file<-file.path(data.path,'fleet_info.dat') 
  finfo<-scan(file,comment.char = "#") 
  
  file<-file.path(data.path,'catch_survey_residuals.out')
  res<-read.table(file,comment.char = "#",header=T)
  res<-subset(res,data=='survey')
  if (standardize) res$residual<- res$stand.residual
  res[res$residual==-99.9 ,'residual']<-NA
  res$Species<-sp.names[res$Species.n]
  res$Fleet<-' '
  for (i in (1:dim(res)[[1]])) res[i,'Fleet']<-fleet.names[res[i,"Species"],res[i,"fleet"]]
  return(res)
}


Init.function() # get SMS.contol object  including sp.names
for (dir in dirs) {  
  a<-read.survey.residuals(data.path=file.path(root,dir),standardize=F)
  a$dir<-dir
  a$label<-labels[ which(dirs==dir)]
  
  if (dir==dirs[1]) aa<-a else aa<-rbind(aa,a)
  
}

head(aa)

#a<-subset(aa,Species==sp & fleet==fl & Age==age)
aa<-subset(aa,!(Species %in% c('Saithe','Mackerel','Plaice','Sole')))
l<-unique(aa$label)
length(l)

if (TRUE) { #individual plots
  dev='png'
  nox=2
  noy=2
  by(aa,list(aa$Species.n,aa$fleet,aa$Age), function(a){
    
    sp<-a[1,"Species"]
    fl<-a[1,'fleet']
    age<-a[1,'Age']
    
    bb<-tapply(a$residual,list(a$Year,a$label),sum,na.rm=TRUE)
 
    newplot(dev,nox,noy,Portrait=TRUE,filename=paste('Survey',sp,fl,age,sep='_'))
    by(a,list(a$label) ,function(x) {
      b<-ts(x$residual, start=min(x$Year), end=max(x$Year), frequency=1) # Yearly Data
      plot(b,main=x[1,'label'],ylab='log(observed/predicted)',xlab='Year');abline(h=0,col='red')
      
     # hw<-HoltWinters(b, beta=FALSE, gamma=FALSE)
    #  hw
    #  plot(hw)
    })
    
    #matplot(x=dimnames(bb)[[1]],y=bb,type='b',ylab='log(observed/predicted)',xlab='Year',main=paste(sp,'fleet:',fl,'Age:',age));abline(h=0,col='red')
 
    newplot(dev,nox,noy,Portrait=TRUE,filename=paste('Survey_pred',sp,fl,age,sep='_'))
    by(a,list(a$label) ,function(x) {
       plot(x=x$model,y=x$residual,main=x[1,'label'],ylab='log(observed/predicted)',xlab='Predicted');abline(h=0,col='red')
    })

    cleanup()
  })
}

#test a<-subset(aa,Species=='Cod' & fleet==1 & Age==1)
if (FALSE) { #combined plots
  dev='png'
  nox=2
  noy=2
  by(aa,list(aa$Species.n,aa$fleet,aa$Age), function(a){
    
    sp<-a[1,"Species"]
    fl<-a[1,'fleet']
    age<-a[1,'Age']
    
    newplot(dev,nox,noy,Portrait=TRUE,filename=paste('Survey',sp,fl,age,sep='_'))
    bb<-tapply(a$residual,list(a$Year,a$label),sum,na.rm=TRUE)
    matplot(x=dimnames(bb)[[1]],y=bb,type='b',ylab='log(observed/predicted)',xlab='Year');abline(h=0,col='red')
    by(a,list(a$label) ,function(x) {
      b<-ts(x$residual, start=min(x$Year), end=max(x$Year), frequency=1) # Yearly Data
      plot(b,main=x[1,'label'],ylab='log(observed/predicted)',xlab='Year');abline(h=0,col='red')
      
      # hw<-HoltWinters(b, beta=FALSE, gamma=FALSE)
      #  hw
      #  plot(hw)
    })
    cleanup()
  })
}



plotResidDist <- function(x)
{
  # make a histogram of the forecast errors:
  mybinsize <- IQR(x)/4
  mysd   <- sd(x)
  mymin  <- min(x) - mysd*5
  mymax  <- max(x) + mysd*3
  # generate normally distributed data with mean 0 and standard deviation mysd
  mynorm <- rnorm(10000, mean=0, sd=mysd)
  mymin2 <- min(mynorm)
  mymax2 <- max(mynorm)
  if (mymin2 < mymin) { mymin <- mymin2 }
  if (mymax2 > mymax) { mymax <- mymax2 }
  # make a red histogram of the forecast errors, with the normal distributed data overlaid:
  mybins <- seq(mymin, mymax, mybinsize)
  hist(x, col="red", freq=FALSE, breaks=mybins,main='')
  # freq=FALSE ensures the area under the histogram = 1
  # generate normally distributed data with mean 0 and standard deviation mysd
  myhist <- hist(mynorm, plot=FALSE, breaks=mybins)
  # plot the normal curve as a blue line on top of the histogram of forecast errors:
  points(myhist$mids, myhist$density, type="l", col="blue", lwd=2)
}

#plotResidDist(b)


#acf(b, lag.max=20)
#Box.test(b, lag=20, type="Ljung-Box")
