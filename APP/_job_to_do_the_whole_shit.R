APP_r<-file.path(prog.path,'APP')
SMSApp0<-file.path("C:","MV","GitHub","SMSApp")


my.area<-c('North Sea','Baltic Sea')[2]
if (my.area=='North Sea') {
  y.end<-2019;y.first=2017; APP_root<-file.path(root,'APP_northSea'); SMSApp<-file.path(SMSApp0,"Data_northsea")
}
if (my.area=='Baltic Sea') {
  y.end<-2021;y.first=2019;APP_root<-file.path(root,'APP_baltic'); SMSApp<-file.path(SMSApp0,"Data_baltic") 
}




out_op<-file.path(APP_root,'op')
out_amoeba<-file.path(APP_root,'amoeba')

SMS<-SMS.control  # just a shorter name

##################################################################
# 1. Make a SMS and files for OP forcast 

# Make an OP control object. File='OP.dat' used for data extraction from SMS (which you have to run afterwards)
source(file=file.path(prog.path.func,'hcr_op_batch_common.r'))

# first  year for calculation of mean stock numbers of other predators  
#                                      Fulmar   Guillemot    Her.Gull   Kittiwake    GBB.Gull      Gannet      Puffin   Razorbill   A.radiata  G.gurnards W.horse.mac N.horse.mac   Grey.seal  H.porpoise      Hake 
if (my.area=='North Sea') y.first.other<-c(2019,        2019,       2019,       2019,      2019,        2019,      2019,         2019,       2019,      2018,        2019,      2019,        2019,      2019,      2015) 

if (my.area=='Baltic Sea') y.first.other<-c(2021) 
        
OP<-make.OP.dat(my.area,my.last.year=SMS@last.year.model+1,first.year.output=SMS@last.year.model+10,
                 stochastic.recruitment=0,  y.end=y.end,y.first=y.first,y.first.other=y.first.other)[["OP"]] 
write.FLOP.control(OP,file="op.dat")

# run to get updated OP_config.dat file
if (FALSE) do.a.full.SMS.run(label="run_",                   # label for output
                  cleanup=F,                      # delete files in the delete Files variable?
                  do.single=F,                    # run SMS in single species mode
                  do.multi.1=T,                   # Make preliminary estimate of "predation parameters"
                  do.multi.2=T,                   # Run the full model, with simultaneously estimation of all parameters except the stomach variance parameter
                  do.multi.2.redo=T,              # Run the full model, with simultaneously estimation of all parameters
                  do.multi.2.redo.Nbar=F,         # Run the full model, with simultaneously estimation of all parameters, Use mean stock numbers (Nbar) for predation
                  do.hessian=F,                   # Make the Hessian matrix and estimate uncertainties
                  do.run=T)                       # Make the run immediately, or just make the batch file for the run



#################################  files for APP

www<-file.path(out_op,'www')

source(file.path(APP_r,'HCR_icons.R'))
source(file.path(APP_r,'make_pred_format.R'))
source(file.path(APP_r,"op_historical.R"))
source(file.path(APP_r,"historical_who_eats_whom.R"))
source(file.path(APP_r,"status_quo.R"))


source(file.path(APP_r,"write_option_files.R"))




if (FALSE) {
 do_plot<-function() {
  plot_summary_ices_multi(
    Portrait=T,                 # graphical output orientation
    include.terminal.year= FALSE,          # plot terminal year (last assessment year +1) as well?
    include.last.assess.year.recruit=FALSE,          # plot recruits terminal year as well?
    
    first.year= -1974,                #first year on plot, negative value means value defined by data
    last.year= 2050,             #last year on plot
    incl.M2.plot=TRUE,
    incl.reference.points=TRUE,
    incl.TSB=FALSE,
    splitLine=FALSE,
    OperatingModel=TRUE,
    redefine.scenario.manually=FALSE,
    output.dir=data.path,
    op.dir=data.path,
    my.dev=c('screen','wmf', 'png', 'pdf')[1]
  )
  }  
  
  # do some op test runs
  rec.mode<-0   # 0=deterministic, 1=Stochastic
  write_options(rec.mode,APP_final=FALSE,change_op_config=FALSE,make_AMOEBA_data=FALSE)
  shell('OP -maxfn 0 -nohess', invisible = TRUE) 
  do_plot()
  
  # do some test runs
  rec.mode<-1   # 0=deterministic, 1=Stochastic
  write_options(rec.mode,APP_final=FALSE,change_op_config=FALSE,make_AMOEBA_data=FALSE)
  shell('OP -maxfn 0 -nohess', invisible = TRUE) 
  do_plot()
  
  # amoeba test, to see if chosen F status quo values are "reasonable"
  rec.mode<-0   # 0=deterministic, 1=Stochastic
  write_options(rec.mode,APP_final=FALSE,change_op_config=FALSE,make_AMOEBA_data=TRUE)
  shell('OP -maxfn 0 -nohess', invisible = TRUE) 
  do_plot()
}


## make the op.dat and op-trigger.dat for the App
rec.mode<-0   # 0=deterministic, 1=Stochastic
write_options(rec.mode,APP_final=TRUE,change_op_config=FALSE,make_AMOEBA_data=FALSE)

### copy OP files
op.files<-c("area_names.in","just_one.in","op_c.in","op_consum.in","op_exploitation.in","op_f.in","op_length_weight_relations.in",
            "op_m.in","op_m1.in","op_m1m2.in","op_multargetf.in","op_n.in","op_n_proportion_m2.in","op_other_n.in","op_price.in","op_prop_landed.in","op_propmat.in",
            "op_reference_points.in","op_seed.in","op_size.in","op_ssb_rec_residuals.in","op_wcatch.in","op_wsea.in","species_names.in",
            "op_msfd.dat","op.dat","op_config.dat","op_trigger.dat","sms.dat")

for (from.file in op.files) {
  to.file<-file.path(out_op,from.file)
  print(paste(file.copy(from.file, to.file, overwrite = TRUE),from.file))
}

write.FLSMS.control(SMS.control,path=out_op,writeSpNames = FALSE)


file.copy("op_config.dat", file.path(out_op,"op_config_master.dat"), overwrite = TRUE)
file.copy("op_exploitation.in", file.path(out_op,"op_exploitation_master.in"), overwrite = TRUE)
file.copy(file.path(root,'program','op.exe'), file.path(out_op,"op.exe"), overwrite = TRUE)
file.copy(file.path(root,'program','op.tpl'), file.path(out_op,"op.tpl"), overwrite = TRUE)


file.copy(file.path(APP_root,"gemmes","units.csv"),file.path(out_op,"units.csv"), overwrite = TRUE)
file.copy(file.path(APP_root,"gemmes","HCR_ini.csv"),file.path(out_op,"HCR_ini.csv"), overwrite = TRUE)


source(file.path(APP_r,"update_environment.R"))
update_environment(my.area)




### copy OP files to SMSApp local Data (Data_northsea or Data_baltic dir)

op.files2<-c(op.files,'environment.Rdata')
for (from.f in op.files2) {
  from.file<-file.path(out_op,from.f)
  to.file<-file.path(SMSApp,from.f)
  print(paste(file.copy(from.file, to.file, overwrite = TRUE),from.file,to.file))
}





### copy OP files to SMSApp Data dir

op.files2<-c("op.exe","op.tpl")
sort(op.files2)
for (from.f in op.files2) {
  from.file<-file.path(out_op,from.f)
  to.file<-file.path(SMSApp0,'Data',from.f)
  print(paste(file.copy(from.file, to.file, overwrite = TRUE),from.file,to.file))
}



#################### make  amoeba files
if (FALSE) {
  rec.mode<-0   # 0=deterministic, 1=Stochastic
  write_options(rec.mode,APP_final=TRUE,change_op_config=FALSE,make_AMOEBA_data=TRUE)
  
  make_data<-TRUE
  source(file.path(APP_r,"hcr_op_batch_amoeba.r"))
  
  ### copy AMOEBA files
  amoeba.files<-c("yield110.csv","status_q_2020.csv","fleetNames_2020.csv","ssb110_2020.csv")
  for (from.file in amoeba.files) {
    to.file<-file.path(out_amoeba,"new_data",from.file)
    file.copy(from.file, to.file, overwrite = TRUE)
  }
  
}  
