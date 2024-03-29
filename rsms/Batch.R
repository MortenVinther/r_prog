# First you have to run: run ini.r in SMS dir, and  _init_rsms.R
#
collectData<-TRUE
testSp<-1L:12L  #species combinations
#testSp<-c(2L,6L)


load(file=file.path(rsms.root,"rsms_input_all.Rdata"),verbose=TRUE)


if (collectData) {
  inp<-pick_species(ps=testSp, inp=inp_all)
  allPars<-inp[['parameters']]
  allData<-inp[['data']]
}



Annual<- FALSE  # annual or quarterly data

### init
source(file.path(rsms.root.prog,"make_rsms_data_function.R"))
source(file.path(rsms.root.prog,"pick_species.R"))
source(file.path(rsms.root.prog,"rsms_function.R"))


### Extract data from SMS
SMSenv<-"rsms_input"


for (sp in testSp) {
  # select a combination of species from the (full) data set
  inp<-pick_species(ps=sp, inp=inp_all) 
  
  #  transform quarterly data into to annual data (testing)
  if (Annual) inp<-into_annual(inp)

  
  #### prepare to run
  
  data<-inp[['data']]
  cat(data$spNames)
  parameters<-inp[['parameters']]
  
  #data$stockRecruitmentModelCode[]<-1L
  data$Debug<-1L
  
  
  #### Run rsms
  
  
  #Going back to our objective function f, first step is to check that you can evaluate the function as a normal R function:
  # func(parameters)   # KALDET VIL PÅVIRKE KALDET TIL MakeAdFun !!!?
  #An error at this point is obviously not due to RTMB.
  
  # adjust if the are species/year combination with zero catches (or assumed F is very low and highly uncertain)
  if (data$zeroCatchYearExists==1) {
    UfMap<-matrix(1L:(dim(parameters$Uf)[[1]]*dim(parameters$Uf)[[2]]),nrow=sum(data$nlogF),byrow=TRUE)
    for (s in 1:data$nSpecies) if (length(data$zeroCatchYear[[s]]) >0 ) {
      zy<-data$zeroCatchYear[[s]]
      fromTo<-data$nlogFfromTo[s,]
      UfMap[fromTo[1]:fromTo[2],zy]<-NA
      parameters$Uf[fromTo[1]:fromTo[2],zy]<-log(0.001)
    }
    UfMap<-factor(UfMap)
  }
  
  if (data$zeroCatchYearExists==1) my.map<-list(Uf=UfMap) else my.map=list()
  random=c("Un","Uf")
  obj <- MakeADFun(func, parameters, random,silent=TRUE,map=my.map)
  
  source(file.path(rsms.root.prog,"lowerUpper.R"))
  
  opt <- nlminb(obj$par, obj$fn, obj$gr, lower=lower, upper=upper)
  save(obj,opt,data,file=file.path(rsms.root,paste0('B',sp,Annual,'.Rdata')))
  
  cat("\nspecies comb:",sp," objective:",opt$objective,"  convergence:",opt$convergence,'\n\n') 
  
  if (collectData) {
    sdrep <- sdreport(obj); cat('Hesssian:',sdrep$pdHess,'\n')
    x<-as.list(sdrep, what="Est")
    if (sp == testSp[1]) newPar<-x else {
      for (i in names(newPar)) {
        if (i %in% random) newPar[[i]]<-rbind(newPar[[i]],x[[i]]) else newPar[[i]]<-c(newPar[[i]],x[[i]])
      }
    }
  }
  
} 

if (collectData) {
  inp[['parameters']]<-newPar
  str(newPar)
  inp[['data']]<-allData
  
  # you can now re-run all species in the same run
}

convert_var<-function(x) {
  xx<-lapply(x,function(x) {
    as.data.frame(x)  %>% mutate(year=as.numeric(rownames(x))) %>% 
      pivot_longer(!year, names_to = "Age_idx", values_to = "N")
  })
  for (s in names(x)) xx[[s]]$species<-s
  xx<-do.call(rbind,xx)
  xx$Age<-as.integer(matrix(unlist(strsplit(xx$Age_idx,' ')),ncol=2,byrow=TRUE)[,2]) -data$off.age
  xx  
}


sms<-Read.summary.table(dir=file.path(root,SMSenv),read.init.function=TRUE) %>% select(Species,Year,Rec,SSB,mean.F)
sms$source='sms'

tab<-NULL
for (sp in testSp) {
  load(file.path(rsms.root,paste0('B',sp,Annual,'.Rdata')))
  sdrep <- sdreport(obj); 
  cat(sp,data$spNames,"objective:",opt$objective,"  convergence:",opt$convergence, " Hessian: ",sdrep$pdHess,'\n')
  rep<-obj$report()
  tab<-rbind(tab,cbind(rep$nlls,all=rowSums(rep$nlls)))
}  

tab

for (sp in testSp) {
  load(file.path(rsms.root,paste0('B',sp,Annual,'.Rdata')))
  
   rep<-obj$report()
  

  if (data$nSeasons==1) N<-lapply(rep$logNq,function(x) (exp(x[,1,])))
  if (data$nSeasons==4) N<-lapply(rep$logNq,function(x) (exp(x[,3,])))
  Recruit<-convert_var(N) %>% filter(Age==0) %>% rename(Species=species,Year=year)
  
  sdrep <- sdreport(obj)
  cat("SP :",sp,data$spNames,'\n')
  cat("objective:",opt$objective,"  convergence:",opt$convergence, " Hessian:",sdrep$pdHess ,"\n")
  
  
  x<-as.list(sdrep, "Est", report=TRUE)
  
  ssb<-x$ssb
  colnames(ssb)<-data$years
  rownames(ssb)<-data$spNames
  ssb<-array2DF(ssb); colnames(ssb)<-c('Species','Year','SSB')
  ssb<- ssb %>% mutate(Year=as.numeric(Year))
  
  x<-as.list(sdrep, "Est", report=FALSE)
  ff<-exp(x$Uf)
  
  avg_F<-NULL
  for (s in (1:data$nSpecies)) {
    i<-data$nlogFfromTo[s,]
    key<-data$keyLogFsta[s,][data$keyLogFsta[s,]>0]
    fff<-ff[i[1]:i[2],,drop=FALSE][key,]
    data$fbarRange[s,]
    ageF<-data$fbarRange[s,]+data$off.age-data$info[s,'faf']+1
    
    avg_F<-rbind(avg_F,data.frame(Year=data$year,Species=data$spNames[s],Species.n=s,mean.F=apply(fff[ageF[1]:ageF[2],],c(2),mean)))
  }

  names(ssb);names(Recruit);names(avg_F)
  rsms<-left_join(left_join(ssb,Recruit),avg_F) %>% select(Species, Year, SSB, N, mean.F) %>% rename(Rec=N) %>% mutate(source='rsms')

  rsp<-unique(rsms$Species)
 
  b<-rbind(rsms,filter(sms,Species %in% rsp )) %>% filter(Year %in% data$years) %>% mutate(Rec=Rec/1000)
  b<-pivot_longer(b,cols=c(Rec,SSB,mean.F),names_to='variable') %>% mutate_if(is_character, as.factor)
  
  
  for (s in (rsp)) {
     bb=filter(b,Species==s)
  
    p<-ggplot(data=bb, aes(x=Year, y=value, group=source)) +
      geom_line(aes(linetype=source,col=source))+
      geom_point(aes(shape=source,col=source))+
      facet_grid(variable ~ ., scales="free_y")+
      ggtitle(unlist(bb[1,'Species']))
    print(p)
    cat('press return to see the next plot:')
    readLines(n=1)
    cat('\n')
  }

  }


