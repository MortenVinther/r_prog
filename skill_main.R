area<-c('NorthSea','Baltic')[1]
make_dirs<-TRUE
use_simul_data<-F

if (area=='Baltic') {
 main_skill_dir<-"Baltic-2022-skill_00"  #should be used in init.r
 conf1<-matrix(
       c("Baltic-2022-skill_1a","1. MS M","yes",
        "Baltic-2022-skill_1b","2. MS M, independent fit","yes",
        "Baltic-2022-skill_1c","3. Smoothed MS M","yes",
        "Baltic-2022-skill_1d","4. Avg, age dep. MS M","yes",
        "Baltic-2022-skill_1e","5. Avg. MS M","yes"
        ),nrow=3,byrow=FALSE)
}
if (area=='NorthSea') {
  if (use_simul_data)  main_skill_dir<-"NS_2020_skill_00_simul"  #should be used in init.r 
  if (!use_simul_data) main_skill_dir<-"NS_2020_skill_00"  #should be used in init.r
 
  conf1<-matrix(
    c("NS_2020_skill_1a","1. MS model","yes",
      "NS_2020_skill_1b","2. M sum","yes",
      "NS_2020_skill_1c","3. M smooth","yes",
      "NS_2020_skill_1d","4. M age level","yes",
      "NS_2020_skill_1e","5. M sp. level","yes"
    ),nrow=3,byrow=FALSE)
}


stopifnot(data.path==file.path(root,main_skill_dir))

t(conf1)
dirs<-conf1[1,]
labels<-conf1[2,]
if (use_simul_data) {
 dirs<-paste(dirs,'simul',sep='_') 
 labels<-paste(labels,'simul',sep='_') 
}

dirsFull<-file.path(root,dirs)

doRun1<-rep(TRUE,length(dirs))
doRun1[conf1[2,]=='no']<-FALSE
doRun1

makeMyDir<-function(my.dir,doSingle=TRUE,doMulti=TRUE){
  myDir<-file.path(root,my.dir)
  if (make_dirs) {
    copySMSfiles(control=SMS.control,scenario.dir=myDir,doSingle,doMulti,verbose=TRUE)
  } 
}

deleteFiles<-c("*.?0?","*.out","*.mcm","*.bin","admodel.*","*.csv","*.std","*.bar","*.mc2","*.cor","*.psv","*.ecm",
               "*.xls","*.html", "mcout*.all","_*.txt","OP_BATCH*.*","SMS.o*","gradient.dat","*.grd","*.Rdata",
               "*.wmf","*.png","*.ps","*.lg","*.log","ud.dat","gradient.dat","op_*.out","iter*.*","stom_and_noise*.in","canum_and_noise*","survey_and_noise*","baseline*.*",
               "HCR_prob.dat","HCR_yield.dat","HCR_SSB.dat","*.par","*.rep","*.hst","*.eva","*.tmp","amoeba*.*","covariance_*.*","forecast*.*")

#deleteFiles<-NA
###################### to be absolutely sure that the newest version is used


# reset skill input files
skill_fleet_make(control=SMS.control,make_excl=TRUE,write_data=TRUE,exclFile='skill_exclude.in') 
#makes the files skill_cpue.in and skill_exclude.in  

if (FALSE) do.a.full.SMS.run(label="run_",                   # label for output
                  cleanup=T,                      # delete files in the deleteFiles variable?
                  do.single=T,                    # run SMS in single species mode
                  do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                  do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                  do.multi.2.redo=T,              # Run the full model, with simultaneously estimation of all parameters
                  do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                  do.hessian=T,                   # Make the Hessian matrix and estimate uncertainties
                  Screen.show=F,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                  deleteFiles=deleteFiles,        # clean up in files before the run is made
                  do.run=T)                       # Make the run immediately, or just make the batch file for the run
  


m1m2_calcs(calcDone=c("averageByAge","averageBySpeciesQuarter","averageBySpecies","simpleSum","smoothed")[1]) # natmor_avg.in
m1m2_calcs(calcDone=c("averageByAge","averageBySpeciesQuarter","averageBySpecies","simpleSum","smoothed")[2]) # natmor_by_species_quarter.in
m1m2_calcs(calcDone=c("averageByAge","averageBySpeciesQuarter","averageBySpecies","simpleSum","smoothed")[3]) # natmor_by_species.in
m1m2_calcs(calcDone=c("averageByAge","averageBySpeciesQuarter","averageBySpecies","simpleSum","smoothed")[4]) # natmor_sum.in
m1m2_calcs(calcDone=c("averageByAge","averageBySpeciesQuarter","averageBySpecies","simpleSum","smoothed")[5]) # natmor_smoothed.in 
labels
M_files<-c("natmor.in","natmor_sum.in","natmor_smoothed.in","natmor_avg.in","natmor_by_species_quarter.in")

################# Make dir 1
dirNo<-1
makeMyDir(dirs[dirNo])

# Multi species run
if (doRun1[dirNo]) do.a.full.SMS.run(label="run_",                   # label for output
                             rundir=dirsFull[dirNo],
                             outdir=dirsFull[dirNo],
                             cleanup=T,                      # delete files in the deleteFiles variable?
                             do.single=T,                    # run SMS in single species mode
                             do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                             do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                             do.multi.2.redo=T,              # Run the full model, with simultaneously estimation of all parameters
                             do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                             do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                             Screen.show=F,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                             deleteFiles=deleteFiles,        # clean up in files before the run is made
                             do.run=T)                       # Make the run immediately, or just make the batch file for the run
  
# single species run
do_run_M<-function(dirNo) {
  makeMyDir(dirs[dirNo])
  file.copy(from=file.path(data.path,M_files[dirNo]),to=file.path(root,dirs[dirNo],"natmor.in"),overwrite=TRUE)
  
  if (doRun1[dirNo]) do.a.full.SMS.run(label="run_",                   # label for output
                                       rundir=dirsFull[dirNo],
                                       outdir=dirsFull[dirNo],
                                       cleanup=F,                      # delete files in the deleteFiles variable?
                                       do.single=T,                    # run SMS in single species mode
                                       do.multi.1=F,                   # Make preliminary estimate of "predation parameters"
                                       do.multi.2=F,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                                       do.multi.2.redo=F,              # Run the full model, with simultaneously estimation of all parameters
                                       do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                                       do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                                       Screen.show=F,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                                       do.run=T)                       # Make the run immediately, or just make the batch file for the run
}

do_run_M(2)
do_run_M(3)
do_run_M(4)
do_run_M(5)

#####################
IDout<-'Skill_1'
if (use_simul_data) IDout<-'Smoothed'

compare_runs_M(
  dirs=dirs,
  labels=labels,
  sumQuarterly=FALSE,  # calc M as sum of quarterly M2
  nox=3, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID=paste0('M_',IDout),         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1974,
  last.year.on.plot=2020,
  include.assess.forcast.line=FALSE,      # vertical line at last assessment year
  include.F.reference.points=FALSE,
  include.SSB.reference.points=FALSE,
  include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
  include.2.std=TRUE,
  # incl.sp=c('Herring'),                      # species number to be included. Numbers or "all"
  #incl.sp="all",
  first.pch=0,    # first pch symbol
  first.color=1,   # first color
  palette="default",               # good for colourful plots
  verbose=TRUE
) 



compare_runs(
  dirs=dirs,
  labels=labels,
  nox=2, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID=paste0('ICES_',IDout) ,         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1974,
  last.year.on.plot=2021,
  plot.MCMC=FALSE,                        # plot values from MCMC scenarios. FALSE=plot hindcast values from "summary_table_raw.out"
  single.species=TRUE,                   # single species mode or multi species mode
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


compare_runs_Z(
  dirs=dirs,
  labels=labels,
  nox=3, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID=paste0("Z_",IDout),         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1975,
  last.year.on.plot=2020,
  include.assess.forcast.line=FALSE,      # vertical line at last assessment year
  include.F.reference.points=FALSE,
  include.SSB.reference.points=FALSE,
  include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
  include.2.std=FALSE,
  #incl.sp=c("Cod","Whiting","Haddock","Saithe",'Herring',"Sprat",'Nor. pout'),                      # species number to be included. Numbers or "all"
  incl.sp="all",
  first.pch=0,    # first pch symbol
  first.color=1,   # first colour
  palette="default"               # good for colour full plots
  #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
) 



source(file.path(prog.path,"compare_runs_objective_function.r"))


compare_runs_stock_rec(dirs,labels,first.year.on.plot=1974,last.year.on.plot=2020,incl.sp="all",
                                 include.CV=TRUE,include.CV2=FALSE,include.mean=TRUE,
                                 palette="R3", makeAllGraphs=FALSE,nox=2, noy=3, w8=8,w11=11,
                                 include.year.labels=TRUE,incl_not_used=TRUE,run.ID='SBB_rec', 
                                 paper=TRUE,facSSB=1000,facRec=1000000,
                                 compare.dir=data.path,verbose=TRUE) 


b<-read_csv(file.path(data.path,"compare_like_Skill_1.csv"))
b<-b %>%dplyr::select( Run=label,catch ,survey=CPUE ,`SSB/R`=SSB.Rec, stomach=stomachs, `no. parameters` =n.par, neg.log.like )
kbl(b, booktabs = F,digits=0,caption='Likelihood contribution from catch, survey , stock-recruitment relation (SSB/R), stomach contents, number of parameters  and weighted sum of liklihood contributions (neg.log.like).')


r<-read_csv(file=file.path(data.path,paste0("compare_like_species_Skill_1.csv"))) %>%  
  filter(Species.n>=first.VPA & !(Species %in% c('Sole','Plaice','Saithe','Mackerel'))) %>%
  mutate(lab2=as.numeric(substr(label,1,1)),label=factor(label),labShort=str_sub(label,1,unlist(gregexpr('m=', label))-2)) 
b<-r%>%dplyr::select( Run=label,Species,catch ,survey=CPUE ,`SSB/R`=SSB.Rec, stomach=stomachs,  neg.log.like=neg_like )

kbl(b, booktabs = T,digits=0,caption='Likelihood contribution from catch, survey , stock-recruitment relation (SSB/R), stomach contents, number of parameters  and weighted sum of liklihood contributions (neg.log.like).')


##################################################################################### 



# from simulated data

sims<-c(-9,0.5,0.75,1,1.25,1.5) # multiplier on cv
#sims<-c(-1.0) # multiplier on cv

dirs<-conf1[1,]
labels<-conf1[2,]

a<-expand.grid(dirs=dirs,sims=sims)
a$dirNo<-match(a$dirs,dirs)
a$labels<-labels
a$labels<-paste0(a$labels,'_','m=',a$sims)
a$des_dir<-file.path(root,paste(a$dirs,a$sims,sep='_'))
a$dirs<-paste(a$dirs,a$sims,sep='_')
a

a<-filter(a,dirNo>1)  # only single species run


do_run_M_sim<-function(dirNo,dest_dir,dirs,smult) {
  makeMyDir(dirs)
  make_survey_data(smult=smult)
  make_catch_data(smult=smult)

  file.copy(from=file.path(data.path,M_files[dirNo]),to=file.path(dest_dir,"natmor.in"),overwrite=TRUE)
  file.copy(from=file.path(data.path,paste0('survey_mult_',smult,'.in')),to=file.path(dest_dir,"fleet_catch.in"),overwrite=TRUE)
  file.copy(from=file.path(data.path,paste0('catch_mult_',smult,'.in')),to=file.path(dest_dir,"canum.in"),overwrite=TRUE)
  

   do.a.full.SMS.run(label="run_",                   # label for output
                     rundir=dest_dir,
                     outdir=dest_dir,
                     cleanup=F,                      # delete files in the deleteFiles variable?
                     do.single=T,                    # run SMS in single species mode
                     do.multi.1=F,                   # Make preliminary estimate of "predation parameters"
                     do.multi.2=F,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                     do.multi.2.redo=F,              # Run the full model, with simultaneously estimation of all parameters
                     do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                     do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                     Screen.show=F,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                     do.run=T)                       # Make the run immediately, or just make the batch file for the run
}
a
i<-14 #test
for (i in (1:dim(a)[[1]])){
  print(a[i,])
  do_run_M_sim(dirNo=a[i,"dirNo"],dest_dir=a[i,"des_dir"],dirs=a[i,"dirs"],smult=a[i,"sims"]) 
}

a


IDout<-'simul'

dirs<-a$dirs
labels<-a$labels


source(file.path(prog.path,"compare_runs_objective_function.r"))

add_mult<-function(infi,outfi) {
  cat(infi,'\n',outfi,"\n")
  r<-read_csv(file=file.path(data.path,infi))
  r$mult<-str_sub(r$label,unlist(gregexpr('m=', r$label)),50)
  r$multn<-as.numeric(str_sub(r$label,unlist(gregexpr('m=', r$label))+2,50))
  nr<-names(r)
  r<-arrange(r,mult,label)
  if ("Species.n" %in% nr) r<-arrange(r,Species.n,mult,label)
  if ("Species.n" %in% nr & "sp_fl" %in% nr) r<-arrange(r,Species.n,sp_fl,mult,label)
  r<-r %>%mutate(dif.like=NULL, dif.n.par=NULL, prob.likelihood.ratio=NULL,meanCPUE=NULL,dir=NULL,stomachs=NULL,penalty=NULL,catch_w=NULL,ssb_rec_w=NULL, stom_w=NULL,	survey_w=NULL) 
  write_csv(r,file=file.path(data.path,outfi) )
}


add_mult(infi=paste0("compare_like_",IDout,".csv"),outfi=paste0("compare_like_",IDout,"2.csv"))   
add_mult(infi=paste0("compare_like_species_",IDout,".csv"),outfi=paste0("compare_like_species_",IDout,"2.csv"))  
add_mult(infi=paste0("compare_like_fleet_",IDout,".csv"),outfi=paste0("compare_like_fleet_",IDout,"2.csv"))  


r<-read_csv(file=file.path(data.path,paste0("compare_like_species_",IDout,"2.csv"))) %>%  
    filter(Species.n>=first.VPA & !(Species %in% c('Sole','Plaice','Saithe','Mackerel'))) %>%
  mutate(lab2=as.numeric(substr(label,1,1)),label=factor(label),labShort=str_sub(label,1,unlist(gregexpr('m=', label))-2)) %>%
  filter(multn %in% c(-9,0.75,1,1.25))



ggplot(r, aes(x=lab2, y=CPUE))   +geom_smooth(method = "lm", se = FALSE)+  geom_point()+facet_grid( Species~mult,scales='free_y')

ggplot(r, aes(x=lab2, y=CPUE))   +geom_smooth(method = "lm", se = FALSE)+  geom_point()+facet_grid(mult~ Species,scales='free_y')

#########################

dirs<-c("NS_2020_skill_00","NS_2020_skill_00_simul")
labels<-c('default','smoothed input')


IDout<-'Simulated'

compare_runs_M(
  dirs=dirs,
  labels=labels,
  sumQuarterly=FALSE,  # calc M as sum of quarterly M2
  nox=3, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID=paste0('M_',IDout),         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1974,
  last.year.on.plot=2020,
  include.assess.forcast.line=FALSE,      # vertical line at last assessment year
  include.F.reference.points=FALSE,
  include.SSB.reference.points=FALSE,
  include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
  include.2.std=TRUE,
  # incl.sp=c('Herring'),                      # species number to be included. Numbers or "all"
  #incl.sp="all",
  first.pch=0,    # first pch symbol
  first.color=1,   # first color
  palette="default",               # good for colourful plots
  verbose=TRUE
) 



compare_runs(
  dirs=dirs,
  labels=labels,
  nox=2, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID=paste0('ICES_',IDout) ,         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1974,
  last.year.on.plot=2021,
  plot.MCMC=FALSE,                        # plot values from MCMC scenarios. FALSE=plot hindcast values from "summary_table_raw.out"
  single.species=TRUE,                   # single species mode or multi species mode
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


compare_runs_Z(
  dirs=dirs,
  labels=labels,
  nox=3, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID=paste0("Z_",IDout),         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1975,
  last.year.on.plot=2020,
  include.assess.forcast.line=FALSE,      # vertical line at last assessment year
  include.F.reference.points=FALSE,
  include.SSB.reference.points=FALSE,
  include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
  include.2.std=FALSE,
  #incl.sp=c("Cod","Whiting","Haddock","Saithe",'Herring',"Sprat",'Nor. pout'),                      # species number to be included. Numbers or "all"
  incl.sp="all",
  first.pch=0,    # first pch symbol
  first.color=1,   # first colour
  palette="default"               # good for colour full plots
  #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
) 



source(file.path(prog.path,"compare_runs_objective_function.r"))


compare_runs_stock_rec(dirs,labels,first.year.on.plot=1974,last.year.on.plot=2020,incl.sp="all",
                       include.CV=TRUE,include.CV2=FALSE,include.mean=TRUE,
                       palette="R3", makeAllGraphs=FALSE,nox=2, noy=3, w8=8,w11=11,
                       include.year.labels=TRUE,incl_not_used=TRUE,run.ID='SBB_rec', 
                       paper=TRUE,facSSB=1000,facRec=1000000,
                       compare.dir=data.path,verbose=TRUE) 
###### 


