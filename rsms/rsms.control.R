### class ######################################################################


setClass("RSMS.control",
    representation(
        test.output         ="numeric",
        VPA.mode            ="numeric",
        first.year          ="numeric",
        first.year.model    ="numeric",
        last.year           ="numeric",
        last.year.model     ="numeric",
        last.season         ="numeric",
        no.species          ="numeric",
        species.names       ="vector",
        first.age           ="numeric",
        rec.season          ="numeric",
        max.age.all         ="numeric",
        species.info        ="matrix",
        SSB.R                ="vector",
        min.catch.CV        ="numeric",
        combined.catches    ="vector",
        
      #  nkeyVarLogN          ="vector",
        keyVarLogN          ="vector",
        
        fModel              ="vector",
       # nkeyLogFsta          ="vector",
        keyLogFsta          ="vector",
        catch.s2.group      ="vector",
  
        firstAgeYearEffect  ="vector",
        catch.sep.age       ="vector",
        catch.sep.year      ="vector",
      
        use_rho             ="vector",
        avg.F.ages          ="matrix",
        discard             ="vector",
        zero.catch.year.season="numeric",
        zero.catch.season.age="numeric",
      checkValueSingle="numeric",
      
      incl.stom.all       ="numeric",
      use.Nbar            ="numeric",
      M2.iterations       ="numeric",
      max.M2.sum2         ="numeric",
      stomach.variance    ="numeric",
      simple.ALK          ="numeric",
      consum              ="numeric",
      size.select.model   ="numeric",
      size.selection      ="vector",
      sum.stom.like       ="vector",
      #stom.obs.var        ="vector",
      stom.max.sumP        ="vector",
      #var.scale.stom      ="vector",
      size.other.food.suit="vector",
      min.stom.cont       ="vector",
      max.stom.sampl      ="vector",
      prey.pred.size.fac  ="vector",
      stom.type.include   ="vector",
      use.overlap         ="numeric",
      checkValueMulti="numeric"
  )
  ,
  prototype=prototype(
        test.output     =0,
        VPA.mode        =0,
        first.year      =1900, 
        first.year.model=1901,                                            
        last.year       =1901,                                            
        last.year.model =1901,                                            
        last.season     =1,                                            
        no.species      =1,  
        species.names   =as.vector("sp1",mode="character"),                                         
        first.age       =0,                                            
        rec.season      =1,                                            
        max.age.all     =0,                                            
        species.info    =matrix(0,ncol=9,nrow=1,dimnames=list(c("sp1"),c("last-age","first-age F>0","+group","predator","prey","SpawningQ","add1","add2","add3"))),
        SSB.R           =as.vector(0,mode="numeric"),
        min.catch.CV    =0.2,                                                
        discard         =as.vector(0,mode="list"),
        combined.catches=as.vector(0,mode="list"),                                                 
        keyVarLogN      =as.vector(0,mode="list"),
        fModel          =as.vector(1,mode="numeric"),
        keyLogFsta      =as.vector(0,mode="list"),
        catch.s2.group  =as.vector(0,mode="list"),  
        firstAgeYearEffect=as.vector(0,mode="numeric"),
        catch.sep.age   =as.vector(0,mode="list"), 
        catch.sep.year  =as.vector(0,mode="list"), 
        use_rho=         as.vector(0,mode="numeric"),
        avg.F.ages      =matrix(0,ncol=2,nrow=1,dimnames=list(c("sp1"),c("first-age","last-age"))),                                        
        discard         =as.vector(0,mode="numeric"),
        zero.catch.year.season =0,
        zero.catch.season.age =0,
        checkValueSingle= -999L,
        incl.stom.all   =0, 
        use.Nbar        =0,
        M2.iterations   =3,
        max.M2.sum2     =0.0, 
        stomach.variance    =1,                                              
        simple.ALK          =0, 
        consum              =0,                                            
        size.select.model   =2, 
        size.selection      =as.vector(0,mode="numeric"),  
        sum.stom.like       =as.vector(0,mode="numeric"),
          stom.max.sumP        =as.vector(0,mode="numeric"),
        size.other.food.suit=as.vector(0,mode="numeric"),
        min.stom.cont       =as.vector(0,mode="numeric"),
        max.stom.sampl      =as.vector(0,mode="numeric"),
        prey.pred.size.fac  =as.vector(0,mode="numeric"), 
        stom.type.include   =as.vector(1,mode="numeric"),                               
        use.overlap         =0,                                              
        checkValueMulti= -999L
        
        )                                                           
)



# in final version remove(validSMS.control)  # We do not need this function any more


### End class ###########################################################


### Methods #############################################################


write.RSMS.control<-function(control,file="rsms.dat",path=NULL,write.multi=TRUE,nice=TRUE,writeSpNames=T,expand=F,w=5) {
  
  wr.matrix<-function(m,text){
    cat("# ",text,"\n",file=file,append=TRUE)
    for (j in (1:dim(m)[1])) cat(m[j,],"\n",file=file,append=TRUE)
  }
  
  wr.matrix.nice<-function(m,sp){
    for (j in (1:dim(m)[1])) cat(m[j,]," #",sp[j],"\n",file=file,append=TRUE)
  }
  
  wr.vector.nice<-function(m,sp){
    cat("# ",formatC(sp,width=w),"\n  ",formatC(m,width=w),"\n",file=file,append=TRUE)
  }
  
  wr.vector.expand<-function(m,sp){
    for (j in 1:length(m)) cat(formatC(m[j],width=w)," # ",formatC(sp[j],width=w),"\n",file=file,append=TRUE)
  }
  
  
  wr.list<-function(l,text1,text2){
    for (j in 1:length(l)) cat(length(l[[j]])," ",file=file,append=TRUE) 
    cat("\t#",text1,"\n# ",text2,"\n",file=file,append=TRUE)  
    for (j in 1:length(l)) cat(l[[j]],"\n",file=file,append=TRUE)    
  }
  
  wr.list.nice<-function(l,text1,text2,sp){
    cat("#",text1,"\n#",formatC(sp,width=w),"\n ",file=file,append=TRUE)
    for (j in 1:length(l)) cat(formatC(length(l[[j]]),width=w+1),file=file,append=TRUE) 
    cat("\n# ",text2,"\n",file=file,append=TRUE)  
    for (j in 1:length(l)) cat(l[[j]],"\t# ",sp[j],"\n",file=file,append=TRUE)    
  }
  
  wr.list.expand<-function(l,text1,text2,sp){
    #cat("#",text1,"\n#",formatC(sp,width=11),"\n",file=file,append=TRUE)
    cat("#",text1,"\n",file=file,append=TRUE)
    for (j in 1:length(l)) cat(formatC(length(l[[j]]),width=w), "  #",formatC(sp[j],width=w),'\n', file=file,append=TRUE) 
    cat("\n# ",text2,"\n",file=file,append=TRUE)  
    for (j in 1:length(l)) cat(l[[j]],"\t# ",sp[j],"\n",file=file,append=TRUE)    
  }
  
  wr.list2<-function(l,text1,text2){
    for (j in 1:length(l)) {
      ifelse(l[[j]]==0, out<-0,out<-length(l[[j]]))
      cat(out," ",file=file,append=TRUE) 
    }
    cat("\t#",text1,"\n# ",text2,"\n",file=file,append=TRUE)  
    for (j in 1:length(l)) cat(l[[j]],"\n",file=file,append=TRUE)    
  }
  wr.list2.nice<-function(l,text1,text2,sp){
    cat("#",text1,"\n#",formatC(sp,width=w),"\n",file=file,append=TRUE)
    for (j in 1:length(l)) {
      ifelse(l[[j]]==0, out<-0,out<-length(l[[j]]))
      cat(formatC(out,width=w),file=file,append=TRUE) 
    }
    cat("\n# ",text2,"\n",file=file,append=TRUE)  
    for (j in 1:length(l)) cat(l[[j]],"\t# ",sp[j],"\n",file=file,append=TRUE)    
  }
  
  if (!inherits(control, "RSMS.control"))
    stop(paste("control" ,"must be an 'RSMS.contol' object!")) 
  
  old.path<-getwd()
  if (!is.null(path)) setwd(path) else path<-old.path
  
  sepLine<-"########################################\n"
  inpOnly<-"(for input transformation only)"
  
  last.pred<-1
  sp.names<-slot(control,"species.names")
  nsp<<-control@no.species
  for (ii in (1:nsp)) if (control@species.info[ii,'predator']!=2) {first.VPA<-ii; break;} #first VPA  species number
  VPA.species<-sp.names[first.VPA:length(sp.names)]
  for (ii in (1:nsp)) if (control@species.info[ii,'predator']==0) {last.pred<-ii-1; break;} #first VPA  species number
  pred.species<-sp.names[1:last.pred]
  
  
  cat("# sms.dat option file\n",file=file)
  cat('# the character "#" is used as comment character, such that all text and numbers\n# after # are skipped by the SMS program\n#\n',file=file, append=TRUE)
  n.<-slotNames(control)
  for (x in n.) {
    switch(x,
           "test.output"        ={ if (nice)  {
             cat(sepLine,file=file,append=TRUE) 
             cat("# Produce test output (option test.output)\n",
                 "#  0 no test output\n",
                   slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },
           "VPA.mode"          ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("# Single/Multispecies mode (option VPA.mode)\n",
                 "# 0=single species mode\n",
                 "# 1=multi species mode, but Z=F+M (used for initial food suitability parm. est.)\n",
                 "# 2=multi species mode, Z=F+M1+M2\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "no.areas"          ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("# Number of areas for multispecies run (default=1)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           
           "first.year"        ={if (nice) {
             cat("#\n#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n#\n",
                 "# single species parameters\n#\n",
                 "#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n#\n",
                 file=file,append=T,sep='')
             cat("## first year of input data (option first.year)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "first.year.model"   ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## first year used in the model (option first.year.model)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           
           
           "last.year"         ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## last year of input data (option last.year)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "last.year.model"   ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## last year used in the model (option last.year.model)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "last.season"       ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("##  number of seasons (option last.season). Use 1 for annual data\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "last.season.last.year"={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## last season last year (option last.season.last.year). Use 1 for annual data\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "no.species"        ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## number of species (option no.species)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "species.names"      ={ cat(sepLine,file=file,append=T)
             cat("# Species names, for information only. See file species_names.in \n# ",sp.names,"\n",file=file,append=T)
           },
           "first.age"         ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## first age all species (option first.age)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "rec.season"        ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## recruitment season (option rec.season). Use 1 for annual data\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "max.age.all"       ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## maximum age for any species(max.age.all)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           
           "species.info"      ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## various information by species\n",
                 "# 1. last age \n",
                 "# 2. first age where catch data are used (else F=0 assumed), First age if species is other predator\n", 
                 "# 3. plus group, 0=no plus group, 1=plus group\n",
                 "# 4. predator species, 0=no, 1=VPA predator, 2=Other predator \n",
                 "# 5. prey species, 0=no, 1=yes\n",
                 "# 6. Spawning season (not used yet, but set to 1)\n", 
                 "# 7. Additional data1\n",
                 "# 8. additional data2\n",
                 "# 9. additional data3\n",
                 "##\n",
                 file=file,append=T,sep="")
             wr.matrix.nice(slot(control,x),paste(1:length(sp.names),sp.names))
           } else wr.matrix(slot(control,x),x)
           },
           "SSB.R"       ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("# Stock Recruit relation\n", 
             "#      0_random walk 1=Ricker, 2=Beverton & Holt, 3=Geom mean,\n",
             "#      4= Hockey stick, 5=hockey stick with smoother,\n",
             "#      >100= hockey stick with known breakpoint (given as input)\n",file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),VPA.species) else wr.vector.nice(slot(control,x),VPA.species)
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "min.catch.CV"      ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## minimum CV of catch observation used in ML-estimation (option min.catch.CV)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "discard"           ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## Use proportion landed information in calculation of yield (option calc.discard)\n",
                 "#    0=all catches are included in yield\n",
                 "#    1=yield is calculated from proportion landed (file proportion_landed.in)\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),VPA.species) else wr.vector.nice(slot(control,x),VPA.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           
           "combined.catches"  ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## use seasonal or annual catches",inpOnly,  "(option combined.catches)\n",
                 "#    0=annual catches with annual time steps or seasonal catches with seasonal time steps\n",
                 "#    1=annual catches with seasonal time steps, assume seasonal catches=annual_catches/n_seasons\n",
                 "#    2=annual catches with seasonal time steps, input seasonal proportions at age from file xx (not implemented yet)\n",
                 file=file,append=T,sep="")   
             if (expand) wr.vector.expand(slot(control,x),VPA.species) else wr.vector.nice(slot(control,x),VPA.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           }, 
           "keyVarLogN" ={if (nice) {
             cat(sepLine,file=file,append=T)
             if (expand) wr.list.expand(slot(control,x),"# Process noise: number of age  groups with the same process N by species",
                                        "first age in sd group (option keyVarLogN)",VPA.species)
             if (!expand) wr.list.nice(slot(control,x),"Process noise: number of age  groups with the same process N by species",
                                       "first age in sd group (option keyVarLogN)",VPA.species)
           } else  wr.list(slot(control,x),"keyVarLogN",x)
           },

           "fModel"   ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## Model for F (default=1) (option fModel):\n",
                 "#   1=random walk (deafult SAM)\n",
                 "#   2=random walk used as year effect in separable F model\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),VPA.species) else wr.vector.nice(slot(control,x),VPA.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           
           "keyLogFsta" ={if (nice) {
             cat(sepLine,file=file,append=T)
             if (expand) wr.list.expand(slot(control,x),"# number of separate age groups for use for random walk F stages (option keyLogFsta)",
                                        "first ages in each group by species",VPA.species)
             if (!expand) wr.list.nice(slot(control,x),"number of separate age groups for use for random walk F stages (option keyLogFsta)",
                                       "first ages in each group by species",VPA.species)
           } else  wr.list(slot(control,x),"keyLogFsta",x)
           },
           
           "catch.s2.group" ={if (nice) {
             cat(sepLine,file=file,append=T)
             if (expand) wr.list.expand(slot(control,x),"catch observations: number of separate catch variance groups by species",
                                        "first age group in each catch variance group",VPA.species)
             if (!expand) wr.list.nice(slot(control,x),"catch observations: number of separate catch variance groups by species",
                                       "first age group in each catch variance group",VPA.species)
           } else  wr.list(slot(control,x),"n.catch.s2.group",x)
           },
           
           "firstAgeYearEffect" ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## First age in separable model (in case of fModel==2), (option firstAgeYearEffect))\n",
                 file=file,append=T,sep="")   
             if (expand) wr.vector.expand(slot(control,x),VPA.species) else wr.vector.nice(slot(control,x),VPA.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           }, 

           "catch.sep.age"  ={if (nice) {
             cat(sepLine,file=file,append=T)
             if (expand) wr.list.expand(slot(control,x),"number of age groups with the same age selection (used  in case of fModel==2)",
                                        "first age in each group (option catch.sep.age)",VPA.species)                                
             
             if (!expand) wr.list.nice(slot(control,x),"number of age groups with the same age selection (used  in case of fModel==2)",
                                       "first age in each group (option catch.sep.age)",VPA.species)
           } else  wr.list(slot(control,x),"catch.sep.age",x)
           },
           
           "catch.sep.year"  ={if (nice) {
             cat(sepLine,file=file,append=T)
             if (expand) wr.list.expand(slot(control,x),"number of year groups with the same age selection (used  in case of fModel==2)",
                                        "first year in each group (option catch.sep.year)",VPA.species)                                
             
             if (!expand) wr.list.nice(slot(control,x),"number of year groups with the same age selection (used  in case of fModel==2)",
                                       "first year in each group (option catch.sep.year)",VPA.species)
           } else  wr.list(slot(control,x),"catch.sep.year",x)
           },
           
           "use_rho"   ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## Correlation in catches (exploitation pattern) between age group (option useRho) (1=use, 0=no use)\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),VPA.species) else wr.vector.nice(slot(control,x),VPA.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           
           "avg.F.ages"        ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## first and last age in calculation of average F by species (option avg.F.ages)\n",
                 file=file,append=T)
             wr.matrix.nice(slot(control,x),VPA.species)
           } else wr.matrix(slot(control,x),x)
           },
          "catch.sep.year"    ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## \n",file=file,append=T)
             if (expand) wr.list.expand(slot(control,x),"catch observations: number of year groups with the same age and seasonal selection","first year in each group (please note #1 will always be changed to first model year)",VPA.species)
             if (!expand) wr.list.nice(slot(control,x),"catch observations: number of year groups with the same age and seasonal selection","first year in each group (please note #1 will always be changed to first model year)",VPA.species)
           } else  wr.list(slot(control,x),"catch.sep.year",x)
           },
            "zero.catch.year.season"={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## year season combinations with zero catch (F=0) (option zero.catch.year.season)\n",
                 "# 0=no, all year-seasons have catches,\n",
                 "# 1=yes there are year-season combinations with no catch.\n",
                 "#   Read from file zero_catch_seasons_ages.in\n",
                 "# default=0\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "zero.catch.season.age"={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## season age combinations with zero catch (F=0) (option zero.catch.season.ages)\n",
                 "# 0=no, all seasons have catches,\n",
                 "# 1=yes there are seasons with no catch. Read from file zero_catch_season_ages.in\n",
                 "# default=0\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
          "checkValueSingle" ={if (nice) {
            cat(sepLine,file=file,append=T)
            cat("## Check value, must be -999 \n",
                slot(control,x),"\n",file=file,append=T,sep="")
          } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
          },
          "checkValueMulti" ={if (nice) {
            cat(sepLine,file=file,append=T)
            cat("## Check value, must be -999 \n",
                slot(control,x),"\n",file=file,append=T,sep="")
          } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
          },         
          
           "incl.stom.all"     = {if (nice) {
             cat(sepLine,file=file,append=T)
             cat("#\n#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n#\n",
                 "# multispecies parameters\n#\n",
                 "#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n#\n",
                 file=file,append=TRUE,sep='')
             cat("# Exclude year, season and predator combinations where stomach data are not included,\n",
                 "#",inpOnly,"(option incl.stom.all)\n",
                 "#   0=no, all stomach data are used in likelihood\n",
                 "#   1=yes there are combinations for which data are not included in the likelihood.\n",
                 "#      Read from file: incl_stom.in\n",
                 "#   default(0)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },                                
           "use.Nbar"        = { if (nice)  {                              
             cat(sepLine,file=file,append=TRUE) 
             cat("##  N in the beginning of the period or N bar for calculation of M2 (option use.Nbar)\n",
                 "#  0=use N in the beginning of the time step (default)\n",
                 "#  1=use N bar\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },      
           "M2.iterations"   = {if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("## Maximum M2 iterations (option M2.iterations) in case of use.Nbar=1\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },      
           "max.M2.sum2"      = {if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("## convergence criteria (option max.M2.sum2) in case of use.Nbar=1\n",
                 "#  use max.M2.sum2=0.0 and M2.iterations=7 (or another high number) to make Hessian\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                            
           "stom.likelihood"   = {if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("## likelihood model for stomach content observations (option stom.likelihood)\n",
                 "#  1 =likelihood from prey weight proportions only (see option below)\n",
                 "#  2 =likelihood from prey weight proportions and from prey numbers to estimate size selection\n",
                 "#  3 =Gamma distribution for prey absolute weight and size selection from prey numbers\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },  
           "stomach.variance"   = {if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("# Variance used in likelihood model for stomach contents as prey weight proportion (option stomach.variance)\n",
                 "#  0 = not relevant, \n",
                 "#  1 = log normal distribution, \n",
                 "#  2 = normal distribution,\n",
                 "#  3 = Dirichlet distribution (with estimation of 'concentration' parameter)\n",
                 "#  4 = Dirichlet distribution (with input of 'concentration' parameter)\n",
                 "#  5 = Dirichlet distribution (with input of maximum 'concentration' parameter)\n",
                 
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                             
           
           "simple.ALK"         = {if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("## Usage of age-length-keys for calc of M2 (option simple.ALK))\n",
                 "#  0=Use only one size group per age (file lsea.in or west.in)\n",
                 "#  1=Use size distribution per age (file ALK_all.in)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },   
           "consum"            = {if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("## Usage of food-rations from input values or from size and regression parameters (option consum)\n",
                 "#",inpOnly,"\n",
                 "#  0=Use input values by age (file consum.in)\n",
                 "#  1=use weight at age (file west.in) and regression parameters (file consum_ab.in)\n",
                 "#  2=use length at age (file lsea.in), l-w relation and regression parameters (file consum_ab.in)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                        
           "size.select.model"   ={if (nice)  {
             cat(sepLine,file=file,append=TRUE)
             cat("## Size selection model based on (option size.select.model)\n",
                 "#  1=length:\n",
                 "#      M2 calculation:\n",
                 "#         Size preference:\n",
                 "#           Predator length at age from file: lsea.in\n",
                 "#           Prey     length at age from file: lsea.in\n",
                 "#         Prey mean weight is weight in the sea from file: west.in\n",
                 "#      Likelihood:\n",
                 "#         Size preference:\n",
                 "#           Predator mean length per length group (file: stom_pred_length_at_sizecl.in) \n",
                 "#           Prey mean length per ength group (file stomlen_at_length.in \n",
                 "#         Prey mean weight from mean weight per prey length group (file: stomweight_at_length.in \n",
                 "#  2=weight:\n",
                 "#      M2 calculation:\n",
                 "#         Size preference:\n",
                 "#           Predator weight at age from file: west.in\n",
                 "#           Prey     weight at age from file: west.in\n",
                 "#         Prey mean weight is weight in the sea from file: west.in\n",
                 "#      Likelihood:\n",
                 "#         Size preference\n",
                 "#           Predator mean weight is based on mean length per predator length group (file: stom_pred_length_at_sizecl.in)\n",
                 "#              and l-w relation (file: length_weight_relations.in), \n",
                 "#           Prey mean weight per prey length group (file: stomweight_at_length.in) \n",
                 "#         Prey mean weight from mean weight per prey length group (file: stomweight_at_length.in \n",file=file,append=T,sep="")
             cat("#  3=weight:\n",
                 "#       M2 calculation: Same as option 2\n",
                 "#       Likelihood:\n",
                 "#         Size preference:\n",
                 "#           Predator mean weight is based on mean length per predator length group (file: stom_pred_length_at_sizecl.in)\n",
                 "#              and l-w relation (file: length_weight_relations.in), \n",
                 "#           Prey mean weight per prey length group (file: stomlen_at_length.in) and l-w relation (file:length_weight_relations.in)\n",
                 "#         Prey mean weight from prey mean length per prey length group (file: stomlen_at_length.in) and l-w relation (file: length_weight_relations.in) \n",
                 "#  4=weight:\n",
                 "#       M2 calculation:\n",
                 "#         Size preference:\n",
                 "#           Predator mean weight from file lsea.in (length in the sea) and l-w relation (file: length_weight_relations.in) \n",
                 "#           Prey mean weight from file lsea.in (length in the sea) and l-w relation (file: length_weight_relations.in) \n",
                 "#       Likelihood:  Same as option 3\n",
                 "#  5=weight in combination with simple.ALK=1:\n",
                 "#       M2 calculation:\n",
                 "#         Size preference:\n",
                 "#           Predator weight based on length from file ALK_all.in (length distribution at age) and l-w relation (file: length_weight_relations.in) \n",
                 "#           Prey     weight based on length from file ALK_all.in (length distribution at age) and l-w relation (file: length_weight_relations.in) \n",
                 "#         Prey mean weight based on length from file ALK_all.in (length distribution at age) and l-w relation (file: length_weight_relations.in) \n",
                 "#       Likelihood: Same as for option 2\n",
                 "#  6=weight in combination with simple.ALK=1:\n",
                 "#       M2 calculation: Same as option 5\n",
                 "#       Likelihood: Same as option 3\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           }
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },                                            
         "size.selection"        ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## spread of size selection (option size.selection)\n",
                 "#   0=no size selection, predator/preys size range defined from observations\n",
                 "#   1=normal distribution size selection\n",
                 "#   3=Gamma distribution size distribution\n",
                 "#   4=no size selection, but range defined by input min and max regression parameters (file pred_prey_size_range_param.in)\n",
                 "#   5=Beta distributed size distribution, within observed size range\n",
                 "#   6=log-Beta size distributed, within observed size range\n",
                 "#\n",
                 "# by predator\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           }, 
           "sum.stom.like"     ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## sum stomach contents over prey size for use in likelihood for prey weight proportions (option sum.stom.like)\n",
                 "#   0=no, use observations as they are; 1=yes, sum observed and predicted stomach contents before used in likelihood for prey weight proportions\n",
                 "#\n",
                 "# by predator\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           }, 
           "stom.obs.var"     ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## # Use estimated scaling factor to link number of observation to variance for stomach observation likelihood (option stom_obs_var)\n",
                 "#    0=no, do not estimate factor (assumed=1);  1=yes, estimate the factor;  2=equal weight (1) for all samples\n",
                 "#\n",
                 "# by predator\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           }, 
           "stom.max.sumP"     ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## # Upper limit for Dirichlet sumP. A low value (e.g. 10) limits the risk of overfitting. A high value (e.g. 100) allows a full fit. (option stom_max_sumP)\n",
                 "# by predator\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           }, 
           "var.scale.stom"    ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## Scaling factor (to bring parameters close to one) for relation between no of stomachs sampling and variance\n",
                 "#  value=0: use default values i.e. 1.00 for no size selection and otherwise 0.1 (option var.scale.stom)\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           
           "size.other.food.suit"={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## other food suitability size dependency  (option size.other.food.suit)\n",
                 "#  0=no size dependency\n",
                 "#  1=yes, other food suitability is different for different size classes\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "min.stom.cont"     ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## Minimum observed relative stomach contents weight for inclusion in ML estimation (option min.stom.cont)\n",
                 "#",inpOnly,'\n',
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "max.stom.sampl"     ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## Upper limit for no of samples used for calculation of stomach observation variance (option max.stom.sampl)\n",
                 "#",inpOnly,'\n',
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },
           "prey.pred.size.fac" ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## Max prey size/ pred size factor for inclusion in M2 calc (option max.prey.pred.size.fac)\n",
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },                                 
           "stom.type.include"  ={if (nice) {
             cat(sepLine,file=file,append=T)
             cat("## inclusion of individual stomach contents observations in ML for weight proportions (option stom.type.include)\n",
                 "#",inpOnly,'\n',
                 "# 1=Observed data\n",
                 "# 2= + (not observed) data within the observed size range (=fill in)\n",
                 "# 3= + (not observed) data outside an observed size range. One obs below and one above (=tails)\n",
                 "# 4= + (not observed) data for the full size range of a prey species irrespective of predator size (=expansion)\n", 
                 file=file,append=T,sep="")
             if (expand) wr.vector.expand(slot(control,x),pred.species) else wr.vector.nice(slot(control,x),pred.species)
             
           } else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
           },            
           
           "use.overlap"       ={if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("## use overlap input values by year and season (use.overlap)\n",
                 "#   0: overlap assumed constant or estimated within the model \n",
                 "#   1: overlap index from file overlap.in (assessment only, use overlap from last year in forecast)\n",
                 "#   2: overlap index from file overlap.in (assessment and forecast)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                              
           "phase.vulnera"     ={if (nice)  {                               
             cat(sepLine,file=file,append=TRUE) 
             cat("## parameter estimation phases for predation parameters\n",
                 "#  the number gives the phase, -1 means no estimation\n#\n",
                 "#  vulnerability (default=2) (phase phase.vulnera)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                               
           "phase.other.suit.slope"  ={if (nice)  {                                
             cat("# other food suitability slope (default=-1) (option phase.other.suit.slope)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                          
           "phase.pref.size.ratio"  ={if (nice)  {                                
             cat("# preferred size ratio (default=2) (option phase.pref.size.ratio)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                          
           "phase.pref.size.ratio.correction"={if (nice)  {                                
             cat("# predator size ratio adjustment factor (default=-1) (option phase.pref.size.ratio.correction))\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                 
           "phase.prey.size.adjustment" ={if (nice)  {                                
             cat("# prey species size adjustment factor (default=-1) (option phase.prey.size.adjustment)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                      
           "phase.var.size.ratio" =  {if (nice)  {                                
             cat("# variance of preferred size ratio (default=2) (option phase.var.size.ratio)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                              
           "phase.season.overlap"   ={if (nice)  {                                
             cat("# season overlap (default=-1) (option phase.season.overlap)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                             
           "phase.stom.var"     ={if (nice)  {                                
             cat("# Stomach variance parameter (default=2) (option phase.Stom.var)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },  
           "phase.mesh.adjust"     ={if (nice)  {                                
             cat("# Mesh size selection of stomach age length key (default=-1) (option phase.mesh.adjust)\n",
                 slot(control,x),"\n",file=file,append=T,sep="")
             cat(sepLine,file=file,append=T)
           } 
             else  cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE) 
           },                                                                                                               
           #other wise                       
           # cat(slot(control,x),"\t#",x,"\n",file=file,append=TRUE)
    )  #end switch
  }
  
  if (writeSpNames) {
    sp.names<-slot(control,"species.names")
    sp.names<-substr(paste(sp.names,"___________",sep=''),1,11)
    sp.names<-gsub(' ','_',sp.names)
    write(sp.names,file=file.path(path,"species_names.in"))
    cat("12345678901\nPlease note. exactly 11 charaters for species names !!!!\n",file=file.path(path,"species_names.in"),append=TRUE) 
  }
  setwd(old.path)
}

read.RSMS.control<-function(dir='.',file="rsms_new.dat",test=FALSE) {       
  #print(file.path(dir,file))
  opt<-scan(file=file.path(dir,file), comment.char = "#",quiet=T) 
  
  n<-1
  control<-new("RSMS.control")
  n.<-slotNames(control)
  n.
  for (x in n.) {
    switch(x,
       "no.species"          = {slot(control,x)<-as.integer(opt[n]); nsp<-as.integer(opt[n]); n<-n+1;
       species.names<-readLines(file.path(dir,"species_names.in"), n=control@no.species)
       species.names<-gsub('_',' ',species.names)
       species.names<-sub('[[:space:]]+$', '', species.names)
       },
       "species.names"       = { slot(control,x)<-species.names; },
       
       "species.info"        = {    ncols<-9
           tmp<-matrix(opt[n:(n-1+nsp*ncols)],ncol=ncols,nrow=nsp,byrow=TRUE,
                       dimnames<-list(species.names,c("last-age","first-age F>0","+group","predator","prey","SpawningQ","add1","add2","add3")));
           slot(control,x)<-tmp;
           n<-n+nsp*ncols;
           n.prey<-0; n.pred<-0; n.oth.pred<-0
           for (s in (1:nsp)) { 
             if(tmp[s,'predator']>=1) n.pred<-n.pred+1;
             if(tmp[s,'predator']==2 | tmp[s,'predator']==3) n.oth.pred<-n.oth.pred+1;
           }
           first.VPA<-n.oth.pred+1
           n.VPA.sp<-nsp-first.VPA+1
       },       
       "SSB.R"                ={slot(control,x)<-as.vector(opt[n:(n-1+n.VPA.sp)]); n<-n+n.VPA.sp},
       
       "min.catch.CV"        = {slot(control,x)<-as.numeric(opt[n]); n<-n+1},
       "discard"             ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp},
       "combined.catches"    ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp},
       "keyVarLogN"      = {v<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp;
           group<-vector("list", length=n.VPA.sp); 
           for (j in 1:n.VPA.sp) {group[[j]]<-as.integer(opt[n:(n-1+v[j])]); n<-n+v[j]; } 
            slot(control,x)<-group},
       "fModel"             ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp},
       "keyLogFsta"      = {v<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp;
         group<-vector("list", length=n.VPA.sp); 
         for (j in 1:n.VPA.sp) {group[[j]]<-as.integer(opt[n:(n-1+v[j])]); n<-n+v[j]; } 
         slot(control,x)<-group},
       "firstAgeYearEffect"= {slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp},
       
      
       "catch.s2.group"      = {v<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp;
             group<-vector("list", length=n.VPA.sp); 
             for (j in 1:n.VPA.sp) {group[[j]]<-as.integer(opt[n:(n-1+v[j])]); n<-n+v[j]; } 
             slot(control,x)<-group},
       "seasonal.catch.s2"   = {slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp},
      "catch.sep.age"    = {v<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp;
         group<-vector("list", length=n.VPA.sp); 
         for (j in 1:n.VPA.sp) {group[[j]]<-as.integer(opt[n:(n-1+v[j])]); n<-n+v[j]; } 
         slot(control,x)<-group }, 
      "catch.sep.year"    = {v<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp;
      group<-vector("list", length=n.VPA.sp); 
      for (j in 1:n.VPA.sp) {group[[j]]<-as.integer(opt[n:(n-1+v[j])]); n<-n+v[j]; } 
      slot(control,x)<-group }, 
      
      
       "use_rho"             ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp},
       
       "avg.F.ages"          = {slot(control,x)<-matrix(opt[n:(n-1+n.VPA.sp*2)],ncol=2,nrow=n.VPA.sp,byrow=TRUE,
                                                        dimnames=list(species.names[first.VPA:nsp],c("first-age","last-age"))); n<-n+n.VPA.sp*2},
       "catch.sep.year"      = {v<-as.vector(as.integer(opt[n:(n-1+n.VPA.sp)])); n<-n+n.VPA.sp;
       group<-vector("list", length=n.VPA.sp); 
       for (j in 1:n.VPA.sp) {group[[j]]<-as.integer(opt[n:(n-1+v[j])]); n<-n+v[j]; } 
       slot(control,x)<-group },
      "checkValueSingle"   = {slot(control,x)<-as.numeric(opt[n]); n<-n+1
                       cat('Check value:',slot(control,x),'\n')
                       stopifnot('Check value should be -999. Something is wrong\n'=slot(control,x) ==  -999) },
      
      "checkValueMulti"   = {slot(control,x)<-as.numeric(opt[n]); n<-n+1
      cat('Check value:',slot(control,x),'\n')
      stopifnot('Check value should be -999. Something is wrong\n'=slot(control,x) ==  -999) },
      
      
      "est.calc.sigma"      = {slot(control,x)<-as.vector(opt[n:(n-1+3)]); n<-n+3},
      "size.selection"      ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.pred)])); n<-n+n.pred},
      "sum.stom.like"       ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.pred)])); n<-n+n.pred},
      "stom.obs.var"        ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.pred)])); n<-n+n.pred},
      "stom.max.sumP"        ={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.pred)])); n<-n+n.pred},
      "var.scale.stom"      ={slot(control,x)<-as.vector((opt[n:(n-1+n.pred)])); n<-n+n.pred},
      "size.other.food.suit"={slot(control,x)<-as.vector(as.integer(opt[n:(n-1+n.pred)])); n<-n+n.pred},
      "min.stom.cont"       ={slot(control,x)<-as.vector(opt[n:(n-1+n.pred)]); n<-n+n.pred},
      "max.stom.sampl"       ={slot(control,x)<-as.vector(opt[n:(n-1+n.pred)]); n<-n+n.pred},
      "prey.pred.size.fac"  ={slot(control,x)<-as.vector(opt[n:(n-1+n.pred)]); n<-n+n.pred},
      "stom.type.include"   ={slot(control,x)<-as.vector(opt[n:(n-1+n.pred)]); n<-n+n.pred},                          
      
      
       # otherwise                        
       {slot(control,x)<-opt[n];n<-n+1}                                 
    )
    if (test) {
      cat(x,'\n')
      #cat(class(slot(control,x)),'\n')
      if (class(slot(control,x))!='list') cat(slot(control,x),'\n') else print(slot(control,x))
    }
  }
  control
}



## show (a replacement for print of S3 classes)
setMethod("show", signature(object="RSMS.control"),
          function(object){
            n.<-slotNames(object)
            for (i in 1:length(n.)) {
              cat(n.[i]) 
              for (j in nchar(n.[i]):25) cat(" ")
              
              if (is.integer(slot(object,n.[i]))) cat(slot(object,n.[i]),"\n")
              
              else if (is.vector(slot(object,n.[i]))) {          
                if (is.list(slot(object,n.[i]))) { 
                  for (k in 1:length(slot(object,n.[i]))) cat("\n\t",slot(object,n.[i])[[k]])
                  cat("\n")
                }
                else cat(slot(object,n.[i]),"\n")
              }
              else if (is.matrix(slot(object,n.[i]))) {
                m<-slot(object,n.[i])
                cat("\n")
                print(m)
              }
            }
          }
)



if (FALSE) {
  RSMS<-new("RSMS.control")    
  RSMS<-read.RSMS.control(dir=data.path,file='rsms.dat')
  write.RSMS.control(RSMS,file=file.path(data.path,"trsms.dat"),write.multi=T,nice=T)
  
  write.FLSMS.control(SMS.dat,file=file.path(data.path,"trSMS.dat"),write.multi=T,nice=F)
  
  write.FLSMS.control(SMS.dat,file=file.path(data.path,"trSMS.dat"),write.multi=T,nice=F)
  
  print(validFLSMS.control(SMS.dat))
  a<-FLSMS.control(first.year=1976,last.year=2000,no.species=3,
                   no.other.predators=1,no.VPA.predators=1,species.names=c('Bird','Cod','Herring'))
  write.FLSMS.control(a,file=file.path(data.path,'tSMS.dat'))
}  
