# First you have to run: run ini.r in SMS dir, and  _init_rsms.R

combSp<-c("S16_S21","S16_S17_S21","S16","S17","S18","S19","S20","S21","s22","S23","S24","S25","S26","S27")


combSp<-c("S16_S17_S21")

Annual<- FALSE  # annual or quarterly data

# species combinations


Batch<-TRUE
my.comb<-"S16"

source(file.path(rsms.root.prog,"make_rsms_data_function.R"))

for (my.comb in combSp) {

  make_rsms_data(dir=my.comb,annual=FALSE)
  dat<-make_rsms_data(dir=my.comb,annual=F,outDir=rsms.root)
  # makes  save(data,parameters,file=file.path(rsms.root,"rsms_input.Rdata"))
 
   source(file.path(root.prog,"r_prog","rsms","rsms.R")) 
  save(obj,opt,data,file=file.path(rsms.root,paste0(my.comb,'.Rdata')))
}

for (my.comb in combSp) {
  load(file.path(rsms.root,paste0(my.comb,'.Rdata')))
  cat("Comb :",my.comb,'\n')
  cat("objective:",opt$objective,"  convergence:",opt$convergence, "  # 0 indicates successful convergence\n")
  
  #rep<-obj$report()
  #print(rep$nlls)
  
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


for (my.comb in combSp) {

  my.comb<-"S16_S21"
  load(file.path(rsms.root,paste0(my.comb,'.Rdata')),verbose=T)
  cat("Comb :",my.comb,'\n')
  cat(data$spNames,'\n')
  cat("objective:",opt$objective,"  convergence:",opt$convergence, "  # 0 indicates successful convergence\n")
  
  rep<-obj$report()
  

  if (data$nSeasons==1) N<-lapply(rep$logNq,function(x) (exp(x[,1,])))
  if (data$nSeasons==4) N<-lapply(rep$logNq,function(x) (exp(x[,3,])))
  Recruit<-convert_var(N) %>% filter(Age==0) %>% rename(Species=species,Year=year)
  
  sdrep <- sdreport(obj)
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
    fff<-ff[i[1]:i[2],][key,]
    data$fbarRange[s,]
    ageF<-data$fbarRange[s,]+data$off.age-data$info[s,'faf']
    
    avg_F<-rbind(avg_F,data.frame(Year=data$year,Species=data$spNames[s],Species.n=s,mean.F=apply(fff[ageF[1]:ageF[2],],c(2),mean)))
  }

  names(ssb);names(Recruit);names(avg_F)
  rsms<-left_join(left_join(ssb,Recruit),avg_F) %>% select(Species,Species.n, Year, SSB, N, mean.F) %>% rename(Rec=N) %>% mutate(source='rsms')

  
  
  sms<-Read.summary.table(dir=file.path(root,my.comb),read.init.function=TRUE) %>% select(Species.n,Year,Rec,SSB,mean.F)
  sms$source='sms'
  sms$Species=data$spNames[sms$Species.n]

  b<-rbind(rsms,sms) %>% filter(Year %in% data$years) %>% mutate(Rec=Rec/1000)
  b<-pivot_longer(b,cols=c(Rec,SSB,mean.F),names_to='variable') %>% mutate_if(is_character, as.factor)
  
  
  for (s in (1:data$nSpecies)) {
     bb=filter(b,Species.n==s)
  
    p<-ggplot(data=bb, aes(x=Year, y=value, group=source)) +
      geom_line(aes(linetype=source,col=source))+
      geom_point(aes(shape=source,col=source))+
      facet_grid(variable ~ ., scales="free_y")+
      ggtitle(unlist(bb[1,'Species']))
    print(p)
    cat('press return\n')
    readLines(n=1)
  }
}
