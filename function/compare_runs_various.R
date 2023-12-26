compare_runs_various<-function(
  paper=TRUE,   # graphs on file (paper=TRUE) or on screen (paper=FALSE)
  do.log=FALSE,  # show values as log values
  first.year.on.plot=1974,
  last.year.on.plot=2070,
  vari=c("M1","M2", "M", "F", "Z", "N", "c.obs","propmat", "west", "weca","ration","Z",'prop.in')[1],
  maxAge=5,     # age > maxage are put in a separate plot 
  nonFish=c('Fulmar','Gannet','GBB. Gull','Grey seal','Guillemot','H. porpoise','Her. Gull','Kittiwake','Puffin','Razorbill'),
  makeAllGraphs=FALSE, # make plots for HTML output
  compare.dir=data.path,
  dirs,
  labels
) {

for (dir in dirs) {
  if ( file.access(file.path(root,dir,"sms.dat"), mode = 0)!=0)  stop(paste('Directory',dir,'does not exist'))
} 

Init.function() # get SMS.control object  including sp.names
#SMS.control@combined.catches

if (SMS.control@zero.catch.year.season==1) {
  a<-scan(file.path(data.path,"zero_catch_year_season.in"),comment='#')
  a<-head(a,-1)
  yq<-data.frame(expand.grid(Quarter=1:SMS.control@last.season,Year=SMS.control@first.year:SMS.control@last.year,Species.n=first.VPA:nsp),yq=a)
  
  yq_all <- yq %>% group_by(Species.n,Quarter) %>% summarize(yq_all=!all(yq==0)) %>% ungroup()
  yq<-yq %>% mutate(yq= (yq==1))
}

a2<-NULL
for (dir in dirs) {
     Init.function(dir=file.path(root,dir)) # get SMS.contol object  including sp.names
    a<-Read.summary.data(dir=file.path(root,dir),read.init.function=F)
    a<-subset(a,(Year>=first.year.on.plot & Year<=last.year.on.plot ))
    a$label<-labels[ which(dirs==dir)]
    if (dir==dirs[1]) usedVars<-names(a) else a<- a %>% dplyr::select(all_of(usedVars))
    if (dir==dirs[1]) a2<-a else a2<-rbind(a2,a)
}

a2<-droplevels(subset(a2,!(Species %in% nonFish)))

isPredVPA<-data.frame(Species=names(SMS.control@species.info[,'predator']),
                   isPred=ifelse(SMS.control@species.info[,'predator']>=1,TRUE,FALSE),
                   isVPA=ifelse(SMS.control@species.info[,'predator']!=2,TRUE,FALSE))
a2<-merge(x=a2,y=isPredVPA,all.x=TRUE)
if (vari %in% c("M1","M2","M","F","c.obs","weca",'propmat','Z')) a2<-subset(a2,isVPA)

if (makeAllGraphs) dev<-'png'

if (makeAllGraphs) my.dir<-compare.dir

an_q<-data.frame(Species.n=first.VPA:nsp,annual=SMS.control@combined.catches )


do_plot_val<-function(vari='M1',my.dir) {

  ######################
  if (vari=="M1") a2$value<-a2$M1
  if (vari=="M2") a2$value<-a2$M2
  if (vari=="M") a2$value<-a2$M
  if (vari=="F") a2$value<-a2$F
  if (vari=="N") a2$value<-a2$N/1000
  if (vari=="c.obs") a2$value<-a2$C.obs/1000
  if (vari=="west") a2$value<-a2$west
  if (vari=="weca") a2$value<-a2$weca
  if (vari=="propmat") a2$value<-a2$propmat
  if (vari=="ration") a2$value<-a2$ration
  if (vari=="Z") a2$value<-a2$Z
  if (vari=="prop.in") a2$value<-a2$prop.in
  
  
  if (vari=="c.obs" | vari=="weca") {
   a2<-left_join(a2,an_q,by = join_by(Species.n)) %>% filter(isVPA)
   a2<-a2 %>% filter((annual==1 & Quarter==3) | annual==0)  
   a2[a2$annual==1,'Quarter']<-9
   
   if (SMS.control@zero.catch.year.season==1) {
     if (FALSE) {
       a2<-left_join(a2,yq_all,join_by(Quarter, Species.n))
       a2[is.na(a2$yq_all),'yq_all'] <-TRUE
       a2<-filter(a2,yq_all)
     }
     a2<-left_join(a2,yq,join_by(Year,Quarter, Species.n))
      a2[is.na(a2$yq),'yq'] <-TRUE
     a2<-filter(a2,yq)
     #filter(a2,Species %in% c('Cod','Herring') & Year==2020 & Age==3)
   }
  
   
  }
  if (vari=="prop.in"){
    
    a2<-a2 %>% filter(value>=0) %>% group_by(Species.n) %>% mutate(allIn=all(value==1)) %>% filter(!allIn)
  }  
  
  
  if (vari=="propmat"){
    a2<-filter(a2,Quarter==1)
  }  
  a3<-aggregate(value~label+Species+Year+Quarter+Age+isVPA+isPred,data=a2,sum) 
 
  if (vari %in% c("ration")) a3<-subset(a3,isPred)
  
  if (do.log) a3$value<-log(a3$value)
  
  pp<-unique(a3$Species)
  
  if (paper) cleanup()
  
  for (p in pp){
    a<-droplevels(subset(a3,Species==p))
    a<-tapply(a$value,list(a$label,a$Species,a$Year,a$Quarter,a$Age),sum) # to get the same number of quarter in all ages
    nq<-dim(a)[[4]]
    print(nq)
    a<-arr2dfny(a)
    colnames(a)<-c('label','Species','Year','Quarter','Age','value')
    a$Year<-as.numeric(as.character(a$Year))
    a$Age<-as.numeric(as.character(a$Age))
  
    a$yo<-ifelse(a$Age<=maxAge,'young','older')
    
    if (length(unique(a$label))>1) for (ages in c('young','older')) {
      b<-droplevels(subset(a,yo==ages))  
      if (dim(b)[[1]]>0 ) { 
        nn<-max(b$Age)-min(b$Age)+1 
        print(paste(p,ages))
        filename=file.path(paste0('Compare_',vari,'_',p,'_',ages,'.png'))
        if (makeAllGraphs) filename=file.path(my.dir,filename) 
      
        if (!paper) trellis.device(device = "windows", color = T, width=9, height=9,pointsize = 2,new = TRUE, retain = FALSE)
        if (paper)trellis.device(device='png',file=filename,width = 1000, height = 800/4*nq+50)
        trellis.par.set(my.trellis.settings)
        if (do.log) yylab<-paste0('log(',vari,')') else yylab=vari
        print(xyplot( value~Year |paste(Species,' Q',Quarter,' Age:',formatC(Age,width=2,digits=0,format='f'),sep=''),groups=label, data=b,
          type='b',lwd=2 , layout=c(nn,nq), ylab=yylab,
           strip = strip.custom( bg='white'),par.strip.text=list(cex=1.25, lines=1.7),
           auto.key = list(space = "bottom", points = T, lines = F,cex=1.5, columns = 2) ,
            scales = list(x = list( cex=1), y= list(cex=1),alternating = 1,relation='same')     # ,relation='free' 'same'
        ))
    }}
  }
  if (paper) cleanup()

  Init.function()
  }
do_plot_val(vari,my.dir)
}

