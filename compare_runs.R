
if (FALSE) {

 
dirs<-c("NS_2023_1_noChange","NS_2023_2_new_stom_extract")  
labels<-c("2020 key run","new extract Stomachs")

dirs<-c("NS_2023_02__FishStomach_Boots_haul",
        "NS_2023_02__FishStomach_Boots_sample_id",
        "NS_2023_02__FishStomach_Simple_haul",
        "NS_2023_02__FishStomach_Simple_sample_id")
labels<-c("boots_haul","boots_sample","simple_haul","simpe_sample")


dirs<-c("NorthSeaKeyRun_2020","NS_2023_02__Simple_0001_haul_as_observed")

labels<-c("key_2020","Recalc")


}


compare_runs(
  dirs=dirs,
  labels=labels,
  nox=2, noy=2,
  paper=TRUE,      # graphics on paper=file (TRUE) or on screen (FALSE)
  run.ID='ICES_com',         # file id used for paper output
  doGrid=TRUE,
  extent.SSB=FALSE,  # plot SSB for the year after last assessment year
  first.year.on.plot=1974,
  last.year.on.plot=2022,
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
  
