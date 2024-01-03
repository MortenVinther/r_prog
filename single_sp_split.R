# user options

# species combinations
spComb<-list(s16=c(16),s17=c(17),s18=c(18),s19=c(19),s20=c(20),s21=c(21),s22=c(22),s23=c(23),s24=c(24),s25=c(25),s26=c(26),s27=c(27),s16_21=c(16,21),s16_17_21=c(16,17,21))
#spComb<-list(s16=c(16),s21=c(21),s22=c(22),s23=c(23),s16_21=c(16,21))
#spComb<-list(s16=c(16),s16_21=c(16,21))
#spComb<-list(s16=c(16))
codes<-c('COD','WHG','HAD','POK','MAC','HER','NSA','SSA','NOP','SPR','PLE','SOL')

make_runs<-TRUE

dev<-'screen'
nox<-1; noy<-3;

######   end user options  ###################
cleanup()

setwd(data.path)
oldDir<-data.path

if  (make_runs) {

  source(file.path(prog.path,"trunc_input.R"))
  # read data and options into FLR objects
  control<-read.FLSMS.control()
  indices<-SMS2FLIndices(control)
  flinfo<-readFleetInfo() %>% select(species.n=species,fleet,fleetName) %>% unique() %>% arrange(species.n,fleet)
  flinfo$fl.n<-1:dim(flinfo)[[1]]



  # write all the files in a number of directories
  for (sps in spComb) {
     # sps<-spComb[[1]]

     # runs are made in a separate directory
     retro.dir<-file.path(oldDir,paste0('S',paste0(sps,collapse = '_S')))
     copySMSfiles(control,scenario.dir=retro.dir,doSingle=TRUE,doMulti=FALSE,doArea=FALSE,verbose=T)

     nc<-trunc.FLSMS.control(control,inclVPA=sps)
     nc@avg.F.ages
     write.FLSMS.control(control=nc,file="sms.dat",path=retro.dir,write.multi=FALSE,nice=TRUE,writeSpNames=T,expand=F)

     trunc_input(outDir=retro.dir,inclVPA=sps)

     fll<-filter(flinfo,species.n %in% sps)
     fl.n<-as.vector(fll$fl.n)
     outIndi<-indices[fl.n]

     FLIndices2SMS(out.path=retro.dir,indices=outIndi,control=nc)

     fname<-'f_q_ini.in'
     a<-readLines(file.path(data.path,fname))
     b<-lapply(sps,function(x) {print(x);lapply(x,function(xx) {print(xx);a[grep(codes[xx-first.VPA+1],a)]})})
     b<-c(unlist(do.call(rbind,b)),' -999 #check')
     writeLines(text=b,con=file.path(retro.dir,fname))

    # setwd(retro.dir)

    cat("\nDoing run ",sps,"\n")

      do.a.full.SMS.run(outdir=retro.dir, rundir=retro.dir,
                    label="run_",                   # label for output
                    cleanup=T,                      # delete files in the deleteFiles variable?
                    do.single=T,                    # run SMS in single species mode
                    do.multi.1=F,        # Make preliminary estimate of "predation parameters"
                    do.multi.2=F,        # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                    do.multi.2.redo=F,   # Run the full model, with simultaneously estimation of all parameters
                    do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters, Use mean stock numbers (Nbar) for predation
                    do.hessian=T,                   # Make the Hessian matrix and estimate uncertainties
                    do.MCMC=F,                      # Prepare for MCMC analysis
                    mcmc=1000,mcsave=100,                # Options for MCMS analysis
                    do.prediction=F,                # Make a prediction
                    pause=F,                        # Make a pause between each stage
                    Screen.show=F,                  # show the output on screen, or save it in file
                    do.run=T,                       # Make the run immediately, or just make the batch file for the run
                    deleteFiles=NA      )       # clean up in files before the run is made
  }
}


dirs<-do.call(c,lapply(spComb,function(x) paste0('S',paste0(x,collapse = '_S'))))
label<-names(dirs)
setwd(oldDir)
oldRoot<-root; root<-oldDir

for (dir in dirs) {
  if ( file.access(file.path(root,dir,"sms.dat"), mode = 0)!=0)  stop(paste('Directory',dir,'does not exist')) #else cat(dir,'\n')
}

read.fit<-function(dir=data.path){
  # Function to read a basic fit
  k<-2
  parfile<-as.numeric(scan(file.path(dir,"sms.par"),
                           what='', n=16, quiet=TRUE)[c(6,11,16)])
  data.frame(n.par=as.integer(parfile[1]),
             neg.log.like=parfile[2],
             max.grad=parfile[3],
             AIC=round(k*parfile[1]+2*parfile[2],2)
  )
}

a<-lapply(dirs,function(x) read.fit(dir=file.path(root,x)))
a<-do.call(rbind,a); a$run<-rownames(a)
a


b<-lapply(dirs,function(x) {
  a<-Read.SMS.std(dir=file.path(root,x)) %>% filter(name %in% c('hist_SSB','avg_F')) %>%
  mutate(label=label[ which(dirs==x)]) %>% select(sp=species,value,std,year,label,CV=CV.round,name)
})
b<-do.call(rbind,b);
b$sp<-factor(b$sp)

head(b)

ggplot(filter(b,name=='hist_SSB'), aes(x=year, y=CV,  shape=sp, color=sp)) + geom_point()+
  facet_wrap(~ label, ncol=3, scales='free')


ggplot(filter(b,name=='avg_F'), aes(x=year, y=CV,  shape=sp, color=sp)) + geom_point()+
  facet_wrap(~ label, ncol=3, scales='free')

