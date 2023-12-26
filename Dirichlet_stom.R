stom<-Read.stomach.data()

stom
read_diri<-function(dir){
stom<-Read.stomach.data(dir=dir,read.init.function=T) %>% rename(Area=SMS.area)

a<-stom %>% select(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter,N.haul, N.samples, Prey,Diri.sum.p ) %>% unique() %>%
  group_by(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter,N.haul, N.samples,Diri.sum.p) %>% summarize(n_prey_sp=dplyr::n()) %>%
  arrange(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter) %>%
  mutate(pred=paste(formatC(Predator.no,flag='0',width=2),Predator)) %>% ungroup()
return(a)
}

dir1<-"NS_2023_04__Simple_0001_haul_as_observed"
dir2<-"NS_2023_04__Boots_0500_haul_mu"
my.labels=c('default','alpha prey')
my.types<-c('default','bootstrap')

# Diri.sum.p ################################
# long format
a<-read_diri(dir=file.path(root,dir1)) %>% select(pred,Predator.length,Year,Quarter,Diri.sum.p) %>% rename(alfa0=Diri.sum.p) %>% mutate(type='dir1') 
b<-read_diri(dir=file.path(root,dir2)) %>%  select(pred,Predator.length,Year,Quarter,Diri.sum.p) %>% rename(alfa0=Diri.sum.p) %>% mutate(type='dir2')
ab<-rbind(a,b)
head(ab)

ab %>% group_by(pred,type) %>% summarize(medi=median(alfa0))

abw<-pivot_wider(ab,names_from=type,values_from=alfa0)


on<-abw %>% group_by(pred) %>% summarize(medi_dir1=median(dir1),medi_dir2=median(dir2))

write.csv(on,file=file.path(data.path,'old_new.csv'))

sp<-sort(unique(ab$pred))
sp<-data.frame(sp,no=as.numeric(substr(sp,1,2)),species=substr(sp,4,30)) %>%
  mutate(species=if_else(no<=8,'Birds',species))

ab$pred<-factor(ab$pred,levels=sp$sp,labels=sp$species)

abw$pred<-factor(abw$pred,levels=sp$sp,labels=sp$species)

ab<- ab %>% mutate(group=cut(ab$alfa0,breaks=c(0,1,2,3,5,7,10,20,30,50,100,500)),n=1) %>% 
   group_by(pred,type,group) %>% summarize(n=sum(n)) 


ab<-ab %>% mutate(type=factor(type,levels=c('dir1','dir2'),labels=my.labels))

 ggplot(ab, aes(fill=type, y=n, x=group)) + 
   geom_bar(position=position_dodge2(preserve='single',padding=0.05), stat="identity")+
   facet_grid(vars(pred),scales = "free")+
   theme(strip.text.y = element_text(angle = 0))+
   #theme_minimal() +
   theme( panel.grid.major = element_line(linetype = "blank"),
          #panel.grid.minor = element_line(linetype = "blank"),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
   )+
   labs(x='alpha0 group',y='count')
 

head(abw)

ggplot(filter(abw,!(pred %in% c('Birds','Hake','H.porpoise',"Grey.seal"))))+
  geom_point(aes(x=old, y=new))+
  facet_wrap(vars(pred),scales = "free")+
   labs(x='Old alpha0',y='New alpha0')+
  geom_abline(intercept=0,slope=1,color='red',lwd=1 )



ab<-ab %>% mutate(type=factor(type,levels=c('old','new'),labels=my.labels))


x<- ggplot(filter(ab,pred %in% c('Cod','Whiting')), aes(fill=type, y=n, x=group)) + 
  geom_bar(position=position_dodge2(preserve='single',padding=0.05), stat="identity")+
  facet_wrap(vars(pred),scales = "free")+
  theme(strip.text.y = element_text(angle = 0))+
  #theme_minimal() +
  theme( panel.grid.major = element_line(linetype = "blank"),
         #panel.grid.minor = element_line(linetype = "blank"),
         axis.text.x = element_text(angle = 90, vjust = 0.5),
  )+
  labs(x='alpha0',y='count',fill='SMS config.')
x

png(filename=file.path(data.path,'cod_whi_alpha_prey.png'),width=700,height=500,pointsize=25)
print(x)
cleanup()


######################################################


read_diri.like<-function(dir){
  stom<-Read.stomach.data(dir=dir,read.init.function=T) %>% rename(Area=SMS.area)
  
  a<-stom %>% select(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter,N.haul, N.samples, Prey,Diri.like ) %>% unique() %>%
    group_by(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter,N.haul, N.samples,Diri.like) %>% summarize(n_prey_sp=dplyr::n()) %>%
    arrange(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter) %>%
    mutate(pred=paste(formatC(Predator.no,flag='0',width=2),Predator)) %>% ungroup()
  return(a)
}

# Diri.like ################################
# long format
a<-read_diri.like(dir=file.path(root,dir1)) %>% select(pred,Predator.length,Year,Quarter,Diri.like) %>% rename(loglike=Diri.like) %>% mutate(type='new') 
b<-read_diri.like(dir=file.path(root,dir2)) %>%  select(pred,Predator.length,Year,Quarter,Diri.like) %>% rename(loglike=Diri.like) %>% mutate(type="old")
ab<-rbind(a,b)
head(ab)

ab %>% group_by(pred,type) %>% summarize(medi=median(loglike))

abw<-pivot_wider(ab,names_from=type,values_from=loglike)

on<-abw %>% group_by(pred) %>% summarize(medi_old=median(old),medi_new=median(new))
write.csv(on,file=file.path(data.path,'old_new_like.csv'))

sp<-sort(unique(ab$pred))
sp<-data.frame(sp,no=as.numeric(substr(sp,1,2)),species=substr(sp,4,30)) %>%
  mutate(species=if_else(no<=8,'Birds',species))

ab$pred<-factor(ab$pred,levels=sp$sp,labels=sp$species)

abw$pred<-factor(abw$pred,levels=sp$sp,labels=sp$species)
summary(ab$loglike)

ab<- ab %>% filter(loglike !=0) %>% mutate(group=cut(loglike,breaks=c(-30,-25,-20,-15,-10,-5,0,5,10,15,20,25,30)),n=1) %>% 
  group_by(pred,type,group) %>% summarize(n=sum(n)) 

ggplot(ab, aes(fill=type, y=n, x=group)) + 
  geom_bar(position=position_dodge2(preserve='single',padding=0.05), stat="identity")+
  facet_grid(vars(pred),scales = "free")+
  theme(strip.text.y = element_text(angle = 0))+
  #theme_minimal() +
  theme( panel.grid.major = element_line(linetype = "blank"),
         #panel.grid.minor = element_line(linetype = "blank"),
         axis.text.x = element_text(angle = 90, vjust = 0.5),
  )+
  labs(x='alpha0 group',y='count')


head(abw)

ggplot(filter(abw,!(pred %in% c('Birds','Hake','H.porpoise',"Grey.seal","A.radiata","W.horse.mac","N.horse.mac"))))+
  geom_point(aes(x=old, y=new))+
  facet_wrap(vars(pred),scales = "free")+
  labs(x='default method',y='alpha prey method')+
  geom_abline(intercept=0,slope=1,color='red',lwd=1 )



ab<-ab %>% mutate(type=factor(type,levels=c('old','new'),labels=my.labels))


x<- ggplot(filter(ab,pred %in% c('Cod','Whiting')), aes(fill=type, y=n, x=group)) + 
  geom_bar(position=position_dodge2(preserve='single',padding=0.05), stat="identity")+
  facet_wrap(vars(pred),scales = "free")+
  theme(strip.text.y = element_text(angle = 0))+
  #theme_minimal() +
  theme( panel.grid.major = element_line(linetype = "blank"),
         #panel.grid.minor = element_line(linetype = "blank"),
         axis.text.x = element_text(angle = 90, vjust = 0.5),
  )+
  labs(x='negative log likelihood',y='count',fill='SMS config.')
x

png(filename=file.path(data.path,'cod_whi_like.png'),width=700,height=500,pointsize=25)
print(x)
cleanup()
