
stomDirs<-c(file.path(root,"NS_2023_02__Simple_0001_haul_as_observedCopy"),file.path(root,"NS_2023_02__Simple_0001_haul_as_observed"))  
stomLabels<-c("2020-key","2023")

#stomDirs<-c(file.path(root,"NS_2023_02__Boots_0500_haul_hb_as_observed"),file.path(root,"NS_2023_02__Boots_0500_sample_id_hb_as_observed"))  
#stomLabels<-c("haul","sample")

data.frame(stomDirs,stomLabels)

stom<-Read.stomach.data.start(dir=stomDirs[1])
filter(stom,is.na(size.ratio))

a<-select(stom,Predator.no,Predator,Prey.no,Prey,Quarter,Year, Predator.length, Predator.length.mean, ,Prey.length.class, Prey.length.mean,  Prey.size, N.haul,stomcon.input,size.ratio,alfa0) 
  
summary(a)


my.colors<-c('red','green','plum','blue','cyan','yellow','coral','skyblue','purple','magenta','limegreen','pink' )

#my.colors<-my.colors[1:length(new.names)]
palette(my.colors)


stom<-Read.stomach.data.start(dir=stomDirs[2])
filter(stom,is.na(size.ratio))
b<-select(stom,Predator.no,Predator,Prey.no,Prey,Quarter,Year, Predator.length, Predator.length.mean, Prey.length.class, Prey.length.mean,  Prey.size, N.haul,stomcon.input,size.ratio,alfa0)
summary(b)

# alpha0
aa<-a%>% select(Predator.no,Predator,Year,Quarter,Predator.length,alfa0) %>% unique() %>% rename(alfa0_a=alfa0)
bb<-b%>% select(Predator.no,Predator,Year,Quarter,Predator.length,alfa0) %>% unique() %>% rename(alfa0_b=alfa0)
ab<-full_join(aa,bb) 
summary(ab)
a0<-by(ab,list(paste(formatC(ab$Predator.no,w=2,flag='0'),ab$Predator)),function(x)summary(x$alfa0_a))
do.call(rbind,a0)
b0<-by(ab,list(paste(formatC(ab$Predator.no,w=2,flag='0'),ab$Predator)),function(x)summary(x$alfa0_b))
do.call(rbind,b0)


ab



#without prey size.class to avoid mess with that 
aa<-a%>% group_by(Predator.no,Predator,Year,Quarter,Predator.length,Prey) %>% summarize(stomcon_a=sum(stomcon.input)) %>% ungroup()
bb<-b%>% group_by(Predator.no,Predator,Year,Quarter,Predator.length,Prey) %>% summarize(stomcon_b=sum(stomcon.input)) %>% ungroup()
ab<-full_join(aa,bb) %>% filter(Predator=='Cod' & Year==1981 & Quarter=="Q1") 
print(ab,n=20)



# special for 2020 keyrun
a[a$Year==1981 & a$Prey.no>0,"Prey.length.class"] <-a[a$Year==1981 & a$Prey.no>0,"Prey.length.class"]-6
a[a$Year>=1991 & a$Prey.no>0,"Prey.length.class"] <-a[a$Year>=1991 & a$Prey.no>0,"Prey.length.class"]-8

#differences, predator information
aa<-a %>% select(Predator.no,Predator, Quarter, Year, Predator.length ) %>% mutate(source.a=TRUE) %>% unique()
bb<-b %>% select( Predator.no,Predator, Quarter, Year, Predator.length ) %>% mutate(source.b=TRUE) %>% unique()

anti_join(x=aa,y=bb) #all rows from x without a match in y.
anti_join(x=bb,y=aa) #all rows from x without a match in y.


in_both<-inner_join(aa,bb) %>% dplyr::select( Predator.no,Predator,Quarter,Year,Predator.length) %>% unique()
dim(in_both);dim(aa);dim(bb)

a<-right_join(a,in_both) %>% mutate(source.a=TRUE)
b<-right_join(b,in_both)%>% mutate(source.b=TRUE)



#differences prey species
aa<-a %>% select(Predator, Quarter, Year, Predator.length,Prey)  %>% unique()
bb<-b %>%select(Predator, Quarter, Year, Predator.length,Prey)  %>% unique()


anti_join(x=aa,y=bb) #all rows from x without a match in y.
anti_join(x=bb,y=aa) #all rows from x without a match in y.

ab<-full_join(aa,bb) %>% arrange(  Predator,Year, Quarter, Predator.length, Prey )
head(ab)

#differences prey information
aa<-select(a, Predator.no,Predator, Quarter, Year, Predator.length,Prey.no,Prey,Prey.length.class,Prey.length.mean,source.a,stomcon.input,Prey.size,size.ratio) %>% 
     rename(stom_a=stomcon.input,preyl_a=Prey.length.mean,Prey.size_a=Prey.size,size.ratio_a=size.ratio) %>% unique()
bb<-select(b, Predator.no,Predator, Quarter, Year, Predator.length,Prey.no,Prey,Prey.length.class,Prey.length.mean,source.b,stomcon.input,Prey.size,size.ratio) %>% 
  rename(stom_b=stomcon.input,preyl_b=Prey.length.mean,Prey.size_b=Prey.size,size.ratio_b=size.ratio) %>% unique()

# special key 2020
aa<-filter(aa,Year %in% c(1981,1991))
bb<-filter(bb,Year %in% c(1981,1991))


a1<-anti_join(x=aa,y=bb) #all rows from x without a match in y.
a1
a2<-anti_join(x=bb,y=aa) #all rows from x without a match in y.
a2

aa
ab<-full_join(aa,bb) %>%arrange( Predator.no,  Predator,Year, Quarter, Predator.length, Prey.no,Prey ,Prey.length.class) %>% as_tibble()

dim(ab)
dim(a1);dim(a2)

ab %>% mutate(prey_size_ab=Prey.size_a/Prey.size_b,source.a=NULL,source.b=NULL,size.ratio_a=NULL,size.ratio_b=NULL) %>% filter(preyl_a == preyl_b) %>% arrange(desc(prey_size_ab)) 
#ab %>% mutate(prey_size_ab=Prey.size_a/Prey.size_b,source.a=NULL,source.b=NULL) %>% filter(preyl_a != preyl_b) %>% arrange(desc(prey_size_ab)) 

ab2<-ab %>% mutate(prey_size_a=NULL,prey_size_b=NULL,source.a=NULL,source.b=NULL) %>% filter(Prey !='Other' & Predator.no>=16) 
summary(ab2$size.ratio_a)
filter(ab2,is.na(size.ratio_a))
filter(ab2,is.na(size.ratio_b))
ss<-ab2 %>% filter(!is.na(size.ratio_a) & !is.na(size.ratio_b)) %>% group_by(Predator.no, Predator,Prey) %>% 
  summarize(min_size_a=min(size.ratio_a),min_size_b=min(size.ratio_b),max_size_a=max(size.ratio_a),max_size_b=max(size.ratio_b), min_ab=min_size_a /min_size_b, max_ab=max_size_a /max_size_b)
print(ss,n=50)





head(ab)
dim(ab)
tail(ab)
filter(ab,Predator.no==20 & Year==1991 &Quarter=="Q4" & Predator.length==400)
filter(ab,Predator.no==16 & Year==1981 &Quarter=="Q2" & Predator.length==400)


# compare  relative stomach content by predator quarter and year from two sources
compare_stoms<-function(stom) {
  cleanup()
  
  stom<-ab
  dev<-"png"
  #dev<-"screen"
  nox<-2
  noy<-3
  
  i<-0
  b<- tapply(stom$stom_a,list(stom$Prey.no),sum)
  all.prey.col<-sort(as.numeric(dimnames(b)[[1]]),decreasing = TRUE)
  all.names<-rep('aaa',length(all.prey.col))
  for (s in (1:length(all.prey.col))) all.names[s]<-sp.other.names[all.prey.col[s]+1]
  
  oldY<-0
  if (length(all.names<5)) my.cex <- 1.5 else my.cex<-1
  
  if (TRUE) by(stom,list(stom$Quarter,stom$Year,stom$Predator.no),function(x) {
    newY<<-x[1,'Year']
    b<- tapply(x$stom_a,list(x$Prey.no,x$Predator.length),sum,na.rm=TRUE)
    b[is.na(b)]<-0
    
    c<- tapply(x$stom_b,list(x$Prey.no,x$Predator.length),sum,na.rm=TRUE)
    c[is.na(c)]<-0
    
    b<-rbind(b,c)
    prey.no<-as.numeric(dimnames(b)[[1]])
    prey.names<-rep('aaa',length(prey.no))
    for (s in (1:length(prey.names))) prey.names[s]<-sp.other.names[prey.no[s]+1]
    length.names<-dimnames(b)[[2]]
    if (oldY != newY) {
      if (dev=='png') cleanup()
      newplot(dev,nox,noy,Portrait=F,filename=paste('comp_stmch_',x[1,]$Predator,x[1,]$Year,sep='-'),dir=data.path);
      par(mar=c(3,5,3,2))
      if (dev=="wmf" ) par(mar=c(2,4,2,2))
      
      # text(x=0.0,y=0.07,"lower: observed",pos=4)
      plot.new();
      title(main=list(paste0("upper:",stomLabels[2],"\nlower:",stomLabels[1]),cex=my.cex*1.0))
      legend("center",all.names,fill=all.prey.col,cex=my.cex,ncol=2)
      oldY<<-newY   
    }
    barplot(b,names=length.names,col=prey.no)
    title(main=paste(x[1,]$Year,x[1,]$Quarter," Pred:",x[1,]$Predator))
    #title(main=paste(x[1,]$Year,x[1,]$Quarter))
    
    abline(h=1,lwd=2)
  })
  cleanup()
}

compare_stoms(stom=ab)

aa<-filter(ab, source.a & is.na(source.b));aa;dim(aa)
bb<-filter(ab, source.b & is.na(source.a));bb;dim(bb)
X11()

#First, it would be easier to convert sn to a factor.

df$sn <- factor(df$sn)
                
p<-select(ab,Prey.no,Prey) %>% unique() %>% arrange(Prey.no)
p
ab$Prey<-factor(ab$Prey,levels=p$Prey,ordered=TRUE)

                
ggplot(ab, aes(x=stom_a, y=stom_b,group=Prey,shape=Prey, color=Prey)) + geom_point()+
  scale_shape_manual(values=1:nlevels(ab$Prey)) +
  facet_grid(Prey ~ .)+ xlab('2019 keyrun')+ylab('2022 keyrun')

ggplot(ab, aes(x=stom_a, y=stom_b,group=Prey,shape=Prey, color=Prey)) + geom_point()+
  scale_shape_manual(values=1:nlevels(ab$Prey)) +
  facet_wrap(Prey.length.class ~ ., ncol=3)+geom_abline(slope=1,intercept=0,color='black')+
  xlab('2019 keyrun')+ylab('2022 keyrun')

ggplot(ab, aes(x=stom_a, y=stom_b,group=Prey,shape=Prey, color=Prey)) + geom_point()+
  facet_wrap(Prey.length.class ~ ., ncol=3)+geom_abline(slope=1,intercept=0,color='black')+
  xlab('2019 keyrun')+ylab('2022 keyrun')

ggplot(ab, aes(x=stom_a, y=stom_b,group=Prey,shape=Prey, color=Prey)) + geom_point()+
  scale_shape_manual(values=1:nlevels(ab$Prey)) +
  facet_grid(Prey.length.class ~ Prey)+geom_abline(slope=1,intercept=0,color='black')+
  xlab('2019 keyrun')+ylab('2022 keyrun')

ggplot(filter(ab,Prey !='Other'), aes(x=stom_a, y=stom_b,group=Prey,shape=Prey, color=Prey)) + 
  scale_shape_manual(values=1:nlevels(ab$Prey)) +geom_point()+
  facet_grid(Prey.length.class ~ Prey)+geom_abline(slope=1,intercept=0,color='black')+
  xlab('2019 keyrun')+ylab('2022 keyrun')

ggplot(filter(ab,Prey=='Herring'), aes(x=stom_a, y=stom_b)) + geom_point()+
  scale_shape_manual(values=1:nlevels(ab$Prey)) +
  facet_wrap(Prey.length.class ~ ., ncol=3)+geom_abline(slope=1,intercept=0,color='black')+ geom_smooth(method = "lm", se = FALSE)+
  xlab('2019 keyrun')+ylab('2022 keyrun')+ggtitle('Herring')
ggplot(filter(ab,Prey=='Sprat'), aes(x=stom_a, y=stom_b)) + geom_point()+
  scale_shape_manual(values=1:nlevels(ab$Prey)) +
  facet_wrap(Prey.length.class ~ ., ncol=3)+geom_abline(slope=1,intercept=0,color='black')+ geom_smooth(method = "lm", se = FALSE)+
  xlab('2019 keyrun')+ylab('2022 keyrun')+ggtitle('Sprat')

ggplot(filter(ab,Prey=='Other'), aes(x=stom_a, y=stom_b)) + geom_point()+
  scale_shape_manual(values=1:nlevels(ab$Prey)) +
  facet_wrap(Predator ~ ., ncol=3)+geom_abline(slope=1,intercept=0,color='black')+ geom_smooth(method = "lm", se = FALSE)+
  xlab('2019 keyrun')+ylab('2022 keyrun')+ggtitle('Sprat')





# sum within preys
ab<-ab %>% group_by(Predator,Year, Quarter, Predator.length, Prey) %>% summarize(stom_a=sum(stom_a,na.rm=TRUE),stom_b=sum(stom_b,na.rm=TRUE))
newplot(dev='png',nox=1,noy=1,Portrait=TRUE,filename='stom_diffrence',dir=data.path,w8=8,w11=10,pointsize=12,doGgplot=TRUE)
ggplot(ab, aes(x=stom_a, y=stom_b,group=Prey,shape=Prey, color=Prey)) + geom_point()+
  facet_grid(Prey ~ .)+geom_abline(slope=1,intercept=0,color='red')+
  xlab('2019 keyrun')+ylab('2022 keyrun')
cleanup()

X11()

ggplot(ab, aes(x=stom_a, y=stom_b,group=Prey,shape=Prey, color=Prey)) + geom_point()+
  facet_grid(Prey ~ .)+geom_abline(slope=1,intercept=0,color='red')+
  xlab('2019 keyrun')+ylab('2022 keyrun')


