# Begin by running init.R and define my.stock.dir

my.stock.dir

skill_single<- file.path(root,"Baltic-2022-skill_01_single")
skill_multi <- file.path(root,"Baltic-2022-skill_02_multi")


# reset skill input files
skill_fleet_make(control=SMS.control,make_excl=TRUE,write_data=TRUE,exclFile='skill_exclude.in') 
#makes the files skill_cpue.in and skill_exclude.in  
  
###### make single species M from the average of the sum of M1 and M2

# make first a multispecies key-run to reset everything

do.a.full.SMS.run(label="run_",                   # label for output
                  cleanup=F,                      # delete files in the deleteFiles variable?
                  do.single=T,                    # run SMS in single species mode
                  do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                  do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                  do.multi.2.redo=F,              # Run the full model, with simultaneously estimation of all parameters
                  do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                  do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                  SSB.R.seperate=F,               # Estimate S/R parameters in a separate step  
                   do.prediction=F,                # Make a prediction
                  pause=F,                        # Make a confirm between each stage
                  Screen.show=T,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                  do.run=T)                       # Make the run immediately, or just make the batch file for the run


# make average M
m1m2_average(averageByYear=FALSE)
file.copy(from=file.path(root,my.stock.dir,"natmor_avg.in"),to=file.path(root,my.stock.dir,"natmor.in"),overwrite=TRUE)



# make single species run directory
copySMSfiles(control=SMS.control,scenario.dir=skill_single,doSingle=TRUE,doMulti=FALSE,verbose=TRUE)
file.copy(from=file.path(root,my.stock.dir,"skill_exclude.in"),to=file.path(skill_single,"skill_exclude.in"),overwrite=TRUE)

# make a command file (single_do.bat) to run single (and run it to check it works) 
do.a.full.SMS.run(label="single_",                   # label for output
                  rundir=skill_single,
                  outdir=skill_single,
                  cleanup=F,                      # delete files in the deleteFiles variable?
                  do.single=T,                    # run SMS in single species mode
                  do.multi.1=F,                   # Make preliminary estimate of "predation parameters"
                  do.multi.2=F,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                  do.multi.2.redo=F,              # Run the full model, with simultaneously estimation of all parameters
                  do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                  do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                  SSB.R.seperate=F,               # Estimate S/R parameters in a separate step  
                  do.prediction=F,                # Make a prediction
                  pause=F,                        # Make a confirm between each stage
                  Screen.show=T,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                  do.run=F)                       # Make the run immediately, or just make the batch file for the run


shell(cmd=file.path(skill_single,"single_do.bat"), invisible = FALSE) 



# make multi species species run directory

copySMSfiles(control=SMS.control,scenario.dir=skill_multi,doSingle=TRUE,doMulti=TRUE,verbose=TRUE)
file.copy(from=file.path(root,my.stock.dir,"skill_exclude.in"),to=file.path(skill_multi,"skill_exclude.in"),overwrite=TRUE)
# make a command file (multi_do.bat) to run multi (and run it to check it works) 
do.a.full.SMS.run(label="multi_",                   # label for output
                  rundir=skill_multi,
                  outdir=skill_multi,
                  cleanup=F,                      # delete files in the deleteFiles variable?
                  do.single=T,                    # run SMS in single species mode
                  do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                  do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                  do.multi.2.redo=T,              # Run the full model, with simultaneously estimation of all parameters
                  do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                  do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                  SSB.R.seperate=F,               # Estimate S/R parameters in a separate step  
                  do.prediction=F,                # Make a prediction
                  pause=F,                        # Make a confirm between each stage
                  Screen.show=T,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                  do.run=F)                       # Make the run immediately, or just make the batch file for the run


shell(cmd=file.path(skill_multi,"multi_do.bat"), invisible = FALSE) 


#######  We are now ready for the skill assessment


excl<-read.csv(file.path(skill_single,"skill_exclude.in"))
ftable(xtabs(~sp2+year2+fl2,data=excl))

#define scenario (what to leave out)
excl<-subset(excl,sp2==1 & fl2==2 & year2 %in% 2013:2017) %>% mutate(excl=0)
excl

# single
skill_fleet_make(SMS.control,path=skill_single,write_data=TRUE,adj=excl) 
shell(cmd=file.path(skill_single,"single_do.bat"), invisible = FALSE) 

#multi'
skill_fleet_make(SMS.control,path=skill_multi,write_data=TRUE,adj=excl) 
shell(cmd=file.path(skill_multi,"multi_do.bat"), invisible = FALSE) 

  


useLogValues<-FALSE

get_residuals<-function(path,excl) {
  a<-Read.catch.survey.residuals(path=path)
  excl$sp2<-excl$sp2+(nsp-first.VPA)
  left_join(x=excl,y=a,by=c("sp2"="Species.n","year2"="Year","fl2"="fleet" )) %>% mutate(excl=NULL,col=NULL,data=NULL)
}

ss<-get_residuals(path=skill_single,excl) %>% mutate(model='single')
ms<-get_residuals(path=skill_multi,excl) %>% mutate(model='multi')

b<-rbind(ss,ms) %>% mutate(f=paste('Age',Age," Fleet",fl2))

ggplot(b, aes(x=model, y=stand.residual)) + 
  geom_boxplot(outlier.colour="red", outlier.shape=8,
               outlier.size=4)+
  facet_wrap(vars(f))



a$logObs=log(a$observed)
a$logPred=log(a$model)
head(a)

b<-aggregate(cbind(meanO=observed,meanLO=logObs,meanLP=logPred)~data+Species+Quarter+fleet+Age,mean,data=a)
head(b)

a<-merge(x=a,y=b,all.x=TRUE)
head(a)
tail(a)


aa<-by(a,list(a$data,a$Species.n,a$Quarter,a$fleet,a$Age),simplify=FALSE,function(x) {
  r<-cor(x=x$logObs, y=x$logPred, method = c("pearson", "kendall", "spearman")[1])
  return(data.frame(Species=x[1,"Species"],Species.n=x[1,"Species.n"],Age=x[1,"Age"],Quarter=x[1,"Quarter"], fleet=x[1,"fleet"],r=r))
})
aa<-do.call(rbind,aa)
aa



aa<-by(a,list(a$data,a$Species.n,a$fleet),simplify=FALSE,function(x) {
  r<-cor(x=x$logObs, y=x$logPred, method = c("pearson", "kendall", "spearman")[1])
  return(data.frame(Species=x[1,"Species"],Species.n=x[1,"Species.n"], data=x[1,'data'],fleet=x[1,"fleet"],r=r))
})
aa<-do.call(rbind,aa)



