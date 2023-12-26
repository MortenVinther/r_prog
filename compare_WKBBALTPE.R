


compare_WKBBALTPEL<-function(labels,dirs,IDout='') {
  
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
  
  
  
  compare_runs_M2(
    dirs=dirs,
    labels=labels,
    sumQuarterly=FALSE,  # calc M2 as sum of quarterly M2
    nox=3, noy=2,
    paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
    run.ID=paste0('M2_',IDout),         # file id used for paper output
    doGrid=TRUE,
    extent.SSB=FALSE,  # plot SSB for the year after last assessment year
    first.year.on.plot=1974,
    last.year.on.plot=2021,
    include.assess.forcast.line=FALSE,      # vertical line at last assessment year
    include.F.reference.points=FALSE,
    include.SSB.reference.points=FALSE,
    include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
    include.2.std=FALSE,
    #incl.sp=c('Herring'),                      # species number to be included. Numbers or "all"
    incl.sp="all",
    first.pch=0,    # first pch symbol
    first.color=1,   # first color
    palette="default"               # good for clolorfull plots
    #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
  ) 
  
  
  compare_runs_M(
    dirs=dirs,
    labels=labels,
    sumQuarterly=FALSE,  # calc M as sum of quarterly M2
    nox=3, noy=2,
    paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
    run.ID=paste0('M',IDout),         # file id used for paper output
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
    palette="default"               # good for clolorfull plots
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
    first.color=1,   # first color
    palette="default"               # good for clolor full plots
    #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
  ) 
  
  
  
  ###  all
  
  
  source(file.path(prog.path,"compare_runs_objective_function.R"))
}
