# *Total biomass 205 kt over the years 1977 -2013;
# %let fyear=1977;
# *%let lyear=1986;
# %let lyear=2013;
# %let totBio=205;
#


fy=1977;
#lyear=1986;
ly=2013;
totBio=205;


species2<-'GUR'
sp<-'G.gurnard'
fdir<-file.path( root,"Data_NorthSea","input_NS_2023",'ByStock',species2)


b<-read_csv(file.path(fdir,"CPUE per length per haul per hour_2023-09-15 13_19_11.csv")) %>%
  rename(l=LngtClass,round=Area,no=CPUE_number_per_hour, square=SubArea) %>%
  mutate(quarter=lubridate::quarter(DateTime),l_10='XXXX') %>%
  filter(round %in% (1:7) & Year>=1974 & Year<=2022 & Species=="Eutrigla gurnardus") %>%
  filter(quarter %in% c(1,3) | (quarter %in% c(1,2,3,4) & Year %in% (1991:1997))) %>%
  filter(!(quarter==3 & Year<1991)) %>%
  mutate(AphiaID=NULL, Sex=NULL,DateofCalculation=NULL, DateTime=NULL,Species=NULL)


b[b$l==0,'l_10']<-'Minus'
b[b$l>0 &b$l<100,'l_10']<-'L00a10'
b[b$l>=100 &b$l<200,'l_10']<-'L10a20'
b[b$l>=200 &b$l<300,'l_10']<-'L20a30'
b[b$l>=300,'l_10']<-'L30a99'



b<-aggregate(no~Year+ quarter+Ship+HaulNo+Gear+round+square+l_10,data=b,FUN=sum)

head(b)
bb<- as_tibble(b) %>%  tidyr::spread(l_10, no)
bb

bb[is.na(bb)]<-0
bb$Minus<-NULL
bb


bbb<-bb %>% mutate(all=L00a10+L10a20+L20a30+L30a99) %>%group_by(Year,quarter,Ship) %>% mutate(sumN=sum(all)) %>% ungroup() %>%
       mutate(status=if_else(sumN==0,'DEL','OK'))

crit<-(bbb$Ship=="ISI" & bbb$Year==1989) | bbb$Gear=="HOB"
bbb[crit,'status']<-'DEL'

bbb[substr(bbb$Gear,1,3)=='DHT','Gear'] <-'DHT'
bbb[substr(bbb$Gear,1,2)=='HT','Gear'] <-'HT'
d<-filter(bbb,status=='OK')


d<-d %>% select(-all,-sumN,-status)

d<-d %>% pivot_longer(cols=starts_with("L"))

info<-data.frame(name=c("L00a10","L10a20","L20a30","L30a99"),
           size=c("00-10","10-20","20-30","30-99"),
           w=   c(4.240,43.542,124.202, 347.154),
           age=c(1,2,3,4)
  )

d<-left_join(d,info,by = join_by(name))


library(mgcv)

# proc genmod data=model;
# by size;
# class year quarter  gear round;
# ods output ParameterEstimates=Parms;
# model cpue= quarter gear round year / d=p link=log type3 SCALE=DEVIANCE obstats;
# *ods output ObStats=obsstat;  *obsstat er kun brugt for at overbevise Anna om at man kan få det samme resultat ud fra kun parameter estimatet.
# *model  cpue= quarter gear round year / d=p link=log type3 SCALE=DEVIANCE ;
#
# run;
#

bb<-d %>% filter(Year>=1991 & Year<=1997 & !is.na(value)) %>% mutate(bio=value*w) %>%
  group_by(Year,quarter,size) %>% summarize(bio=mean(bio),value=mean(value)) %>%
  group_by(     quarter,size) %>% summarize(bio=mean(bio),value=mean(value))

bbb<-bb %>% group_by(size)%>% summarize(biodist=sum(bio)) %>% ungroup()
bbb

d$value<-round(d$value)

dp1<-expand.grid(Gear='GOV',round='6',Year=as.character(min(d$Year):max(d$Year)), quarter=as.character(1:4))  #,  L00a10=NA, L10a20=NA,L20a30=NA,L30a99=NA)
summary(dp1)
dp1

d<-d %>% mutate(Year=as.character(Year),quarter=as.character(quarter) ,round=as.character(round))

o<-by(d,d$name,function(x) {
   a20<-glm(value ~ quarter+Year+Gear+round,family = poisson(link = "log"),data=x)

   d2<-cbind(dp1,predicted=predict(a20, type = "terms",newdata=dp1),name=x[1,'name'])
   d2<-left_join(d2,info)

})
head(o[[1]])


all<-do.call(rbind,o) %>%as_tibble()
all
all<-  all  %>% mutate(y.effect=exp(predicted.Year),q.effect=exp(predicted.quarter), predicted.N=y.effect*q.effect,  bio=predicted.N*w) %>%
          mutate(Year=as.numeric(as.character(Year)))

ay<-all %>% select(Year,size,y.effect,bio) %>% unique() %>% arrange(size,Year)
X11()
ggplot(ay, aes( Year, y.effect,group=1)) +
  geom_point()+geom_line()+facet_wrap(vars(size))


bbb<-bbb %>% mutate(biodist=biodist/sum(biodist))
# find scaling factor
 ny<-ly-fy+1;nq<-4
a2<-all %>% filter(Year>=fy & Year<=ly) %>%
  mutate(Gear=NULL,round=NULL,predicted.Year=NULL,predicted.quarter=NULL,predicted.Gear=NULL,predicted.round=NULL) %>%
  group_by(size) %>%  summarize(sumbio=sum(bio)/ny/nq) %>% ungroup()
a2

target<-totBio*1000000 # kg
bbbb<-bbb %>% mutate(target_size=biodist*target)

factor<-left_join(a2,bbbb) %>% mutate(factor=target_size/sumbio)

factor
sum(factor$target_size)

sms<-left_join(all,factor,by = join_by(size))%>%
  transmute(year=Year,size,quarter,age,N=predicted.N*factor,WSEA=w/1000,bio=N*WSEA,sub_area='A',species='GUR') %>%
  arrange(year,size,quarter)
sms


tst<-filter(sms,year>=fy & year<=ly) %>% group_by(size) %>% summarize(bio=sum(bio)/ny/nq)
tst
sum(tst$bio)

old<-read.table(file=file.path(fdir,"gur2020.dat"),sep=' ',header=T) %>%as_tibble() %>% arrange(year,size,quarter)
tst<-filter(old,year>=fy & year<=ly)
tst<- tst %>% group_by(size) %>% summarize(bio=sum(bio)/ny/nq)
tst
sum(tst$bio)



old$type<-'old'
sms$type<-'new'
old;sms

both<-rbind(old,sms)
ggplot(filter(both,quarter==1), aes( year, bio,group=type,color=type)) +
  geom_point()+geom_line()+facet_wrap(vars(size),,scales='free_y')

ggplot(filter(both,quarter==3), aes( year, bio,group=type,color=type)) +
  geom_point()+geom_line()+facet_wrap(vars(size),scales='free_y')

sms$type<-NULL
sms$quarter=as.numeric(sms$quarter)

write.table(sms,file=file.path(fdir,'gur2023.dat'),row.names=FALSE)
#
# PROC EXPORT DATA= WORK.newdata
# OUTFILE= "&newOut\gur2020.dat"
# DBMS=DLM REPLACE;
# DELIMITER=' ';
# PUTNAMES=YES;
# RUN;
#
