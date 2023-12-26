
if (FALSE) {
 
  
  

  
  ################### 
  
 
   dirs<-c("NS_2023_01_2020_key_run","NS_2023_04__Simple_0001_haul_as_observed","NS_2023_04__Boots_0500_haul_as_observed","NS_2023_04__Boots_0500_haul_mu")
   labels<-c("2020 key run","simple","boot as obs","boot mu")
 
   dirs<-c("NS_2023_04__Boots_0500_haul_as_observed","NS_2023_04__Boots_0500_haul_mu")
   labels<-c("boot as obs","boot mu")
  
   
   dirs<-c("NS_2023_04__Simple_0001_haul_as_observed","NS_2023_04__Boots_0500_haul_as_observed","NS_2023_04__Boots_0500_haul_mu")
   labels<-c("simple","boot as obs","boot mu")
   
   
   dirs<-c("NS_2023_04__Simple_0001_haul_as_observed","NS_2023_04__Boots_0500_haul_mu","NS_2023_04__Boots_0500_haul_as_observed","NS_2023_04__Boots_0500_haul_mu_max")
   labels<-c("default","alpha prey","alpha 0","alpha max")
   
   
    
  }

data.frame(dirs,labels)

compare_runs_M2(
  dirs=dirs,
  labels=labels,
  sumQuarterly=FALSE,  # calc M2 as sum of quarterly M2
  nox=3, noy=2,
  paper=T,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID='SMS',         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1974,
  last.year.on.plot=2022,
  include.assess.forcast.line=FALSE,      # vertical line at last assessment year
  include.F.reference.points=FALSE,
  include.SSB.reference.points=FALSE,
  include.1.std=FALSE,                   # Include values plus/minus 1 times the standard deviation
  include.2.std=FALSE,
  #incl.sp=c('Herring'),                      # species number to be included. Numbers or "all"
  incl.sp="all",
  first.pch=0,    # first pch symbol
  first.color=1,   # first color
  palette="default"               # good for colourful plots
  #palette(gray(seq(0,.9,len=10)))  # gray scale for papers, use len =500 to get black only
) 
  
  

