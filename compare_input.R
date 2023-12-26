if (FALSE) {
  
  dirs.keyrun<-c("NS_2023_01_2020_key_run","NS_2023_04")
  labels.keyrun<-c("2020","2023")
  
  my.first.year.on.plot<-1974
  my.last.year.on.plot<-2024
}

for (v in c("M1", "c.obs","propmat", "west", "weca","ration")) {
#  for (v in c( "c.obs", "weca")){
  compare_runs_various(
    paper=TRUE,   # graphs on file (paper=TRUE) or on screen (paper=FALSE)
    do.log=v %in% c ("c.obs", "west", "weca","ration"),
    first.year.on.plot=my.first.year.on.plot,
    last.year.on.plot=my.last.year.on.plot,
    vari=v,
    maxAge=5,     # age > maxage are put in a separate plot 
    nonFish=c('Fulmar','Gannet','GBB. Gull','Grey seal','Guillemot','H. porpoise','Her. Gull','Kittiwake','Puffin','Razorbill'),
    makeAllGraphs=FALSE, # make plots for HTML output
    dirs=dirs.keyrun,
    labels=labels.keyrun,
    compare.dir=data.path
    #my.dir=my.dir  # output directory
  ) 
}
