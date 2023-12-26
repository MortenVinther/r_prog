# do a series of runs and output for input directories.

indir<-c("NS_2023_key_run","NS_2023_04__Boots_0500_haul_mu","NS_2023_04__Simple_0001_haul_as_observed","NS_2023_04__Boots_0500_haul_as_observed","NS_2023_04__Boots_0500_haul_as_observed_nomagic","NS_2023_04__Boots_0500_haul_mu_max","NS_2023_04__Boots_0500_haul_mu_plaiceM1","NS_2023_04__Boots_0500_haul_mu_nomagic")

indir<-c("NS_2023_key_run","NS_2023_04__Boots_0500_haul_mu","NS_2023_04__Simple_0001_haul_as_observed","NS_2023_04__Boots_0500_haul_as_observed","NS_2023_04__Boots_0500_haul_meanBoots")

data.frame(n=1:length(indir),indir)


#indir<-indir[c(1,2,3)]
#indir<-indir[c(4,5)]
indir<-indir[c(3)]
indir

old.data.path<-data.path
dorun<-T
do.hessian<-T
doretro<-T
dofig<-T
doknitr<-T


#just checking file existence
for (i in (indir)) {
  f<-file.path(root,i,'SMS.dat')
  cat(i,f,' exists:',file.exists(f),'\n')
}

# test my.stock.dir<-indir[1]
for (my.stock.dir in indir) {
  cat(my.stock.dir,'\n')
  data.path<-data.path<-file.path(root,my.stock.dir)
  setwd(data.path)
  
  if (dorun) do.a.full.SMS.run(label="run_",                   # label for output
                    cleanup=T,                      # delete files in the deleteFiles variable?
                    do.single=T,                    # run SMS in single species mode
                    do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                    do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                    do.multi.2.redo=T,              # Run the full model, with simultaneously estimation of all parameters
                    do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter, Use mean stock numbers (Nbar) for predation
                    do.hessian=do.hessian,                   # Make the Hessian matrix and estimate uncertainties
                    shake.ms2.par=F,
                    SSB.R.seperate=F,               # Estimate S/R parameters in a separate step  
                    do.MCMC=F,                      # Prepare for MCMC analysis
                    mcmc=0,mcsave=0,                # Options for MCMC analysis
                    do.prediction=F,                # Make a prediction
                    pause=F,                        # Make a confirm between each stage
                    Screen.show=F,                  # show the output on screen (TRUE), or save it in files "*.lg" (FALSE)
                    do.run=T,                       # Make the run immediately, or just make the batch file for the run
                    deleteFiles=deleteFiles,        # clean up in files before the run is made
                    HPC=F)                          # run it as batch program on the UNIX High Performance Computer 
  
  
  if (doretro) source(file.path(prog.path,"retrospectiv_multi_sp.R"))
  
  globalStockDir<<- my.stock.dir
  if (dofig) source(file.path(prog.path,"makeAllGraphs_NorthSea.R"))

  if (doknitr) {
    save(list = ls(all.names = TRUE), file = file.path(data.path,"SMS.RData"), envir = .GlobalEnv)
    inpf<-file.path(data.path,'HTML',paste0(my.stock.dir,'.Rmd'))
    outdir<-file.path(data.path,'HTML')
     rmarkdown::render(input=inpf,output_dir=outdir,output_format="html_document")
  }

}

data.path<-old.data.path

