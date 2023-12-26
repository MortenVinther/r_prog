
###################################################
# Other predators

a<-read.csv(file=file.path(root,exchangeDir,'other_sp_key_2020.csv'))


sort(unique(a$species))
a<-subset(a,!(species %in% c('N_M','W_M')),select=-species.n)
a<-subset(a,!(species %in% c('HAK','HKE')))

# just add the new years, re-using old data
b<-subset(a,year==lastYold)
for (y in ((lastYold+1):lastY)) {
  b$year<-y
  a<-rbind(a,b)
}

#replace Grey seal number with 2023 estimates
gse_old<-filter(a,species=='GSE' & N>0) %>% mutate(type="old")
names(gse_old)
summary(gse_old)
gse<-read.csv(file=file.path(root,exchangeDir,'GSE_stock_numbers_from_Vanessa.csv')) %>% filter(species=='GSE') %>% mutate(type="new",species.n=NULL)
names(gse)
gg<-rbind(gse,gse_old) %>% mutate(time=year+(quarter-1)*0.25+0.125)
# plot(gg$time,y=gg$N)

a<-subset(a,species !='GSE')
names(b)
gse$type<-NULL
a<-rbind(a,gse)

a$species.n<-a$sub_area<-NULL  # added later on;


# birds
birds<-c(  "FUL",           "GBG",       "GLT",              "GNT",          "HEG",         "KTW" ,                  "PUF",           "RAZ" )
species<-c("NorthernFulmar","GBB_Gull",  "CommonGuillemot", "NorthernGannet","HerringGull", "BlackLeggedKittiwake", "AtlanticPuffin", "Razorbill")
bb<-data.frame(birds,species)
species2<-'Birds'
sp<-'XXX'
fdir<-file.path( root,"Data_NorthSea","input_NS_2023",'ByStock',species2)
b<-read.csv(file=file.path(fdir,'cont_prop_sms.csv')) # file from Michael Spencer
b<-b %>% mutate(X=NULL) %>% rename(year=Year,quarter=Quarter) %>%
  pivot_longer(names_to="species", values_to="Nnew",cols=NorthernFulmar:Razorbill)
b<-left_join(b,bb,by = join_by(species)) %>% mutate(species=birds,birds=NULL,age=1)
head(filter(a,species=='FUL'))
head(filter(b,species=='FUL'))

a<-left_join(a,b,by = join_by(year, species, age, quarter)) %>% mutate(N=if_else(is.na(Nnew),N,Nnew)) %>% mutate(Nnew=NULL)


# Hake
species2<-'HKE'
sp<-'Hake'
fdir<-file.path( root,"Data_NorthSea","input_NS_2023",'ByStock',species2)
b<-read.csv(file=file.path(fdir,'HAKE_N_WSEA.csv'))

a<-subset(a,species!=species2)
a<-rbind(a,subset(b,select=names(a)))
sort(unique(a$species))



 # Rays
  sp<-'RAJ'
  spName<-'Amblyraja radiata'
  old<-subset(a,species==sp)
  old$biomass<-old$N*old$WSEA
  t1<-tapply(old$biomass,list(old$year,old$quarter),sum)
  aver<-apply(t1,1,mean)
  old<-round(cbind(t1,aver))
  old
  mean(old[as.character(1982:2013),5])



  fdir<-file.path( root,"Data_NorthSea","input_NS_2023",'ByStock',sp)

  b<-read.table(file=file.path(fdir,"raj2023.dat"),header=TRUE)
  head(b)
  b<-subset(b,year<=lastY)
  t1<-tapply(b$N*b$WSEA,list(b$year,b$quarter),sum)
  ftable(t1)

  aver<-apply(t1,1,mean)
  new<-round(cbind(t1,aver))
  new
  mean(new[as.character(1982:2013),5])  # 100000
  dim(old);dim(new)
  new[as.character(1974:2019),]/old[as.character(1974:2019),]


  a<-subset(a,species!=sp)
  a<-rbind(a,subset(b,select=names(a)))



  # grey gurnard
  sp<-'GUR'
  spName<-'grey gurnard'

  old<-subset(a,species==sp)
  old$biomass<-old$N*old$WSEA
  t1<-tapply(old$biomass,list(old$year,old$quarter),sum)
  aver<-apply(t1,1,mean)
  old<-round(cbind(t1,aver))
  old
  mean(old[as.character(1977:2013),5])

  fdir<-file.path( root,"Data_NorthSea","input_NS_2023",'ByStock',sp)


  b<-read.table(file=file.path(fdir,"gur2023.dat"),header=TRUE)
  head(b)
  b<-subset(b,year<=lastY)
  t1<-tapply(b$N*b$WSEA,list(b$year,b$quarter),sum)
  aver<-apply(t1,1,mean)
  new<-round(cbind(t1,aver))
  new
  mean(new[as.character(1977:2013),5])
  dim(old);dim(new)
  new/old


  a<-subset(a,species!=sp)
  a<-rbind(a,subset(b,select=names(a)))



# horse mackerel
sp<-'W_H'
spName<-'Western Horsemackerel'
old<-subset(a,species==sp)
old$biomass<-old$N*old$WSEA
old<-round(tapply(old$biomass,list(old$year,old$quarter),sum))
old


a<-subset(a,species!=sp)

fdir<-file.path( root,"Data_NorthSea","input_NS_2023",'ByStock',sp,"SS3")
b <- read.csv(file=file.path(fdir,"horse_mac_N_West_2023.csv"),header=TRUE)
b<-subset(b,species==sp)
head(b)
new<-round(tapply(b$N*b$WSEA,list(b$year,b$quarter),sum))
new
old
a<-rbind(a,subset(b,select=names(a)))


sp<-'N_H'
spName<-'North Sea  Horsemackerel'
a<-subset(a,species!=sp)

b <- read.csv(file=file.path(fdir,"horse_mac_N_West_2023.csv"),header=TRUE)
b<-subset(b,species==sp)
a<-rbind(a,subset(b,select=names(a)))


if (dim(a[a$N>0 & a$WSEA<=0,])[[1]]>0) stop('Other predator N>0 and Wsea <=0')

round(tapply(a$N*a$WSEA,list(a$year,a$species),sum)/4)
####

BIO.02<-a

subset(a,species=='W_H' & year==1986)

save(BIO.02,file=file.path(root,exchangeDir,'BIO_02.Rdata'))

########################
old<-read.csv(file=file.path(root,exchangeDir,'other_sp_key_2020.csv'))
load(file=file.path(root,exchangeDir,'BIO_02.Rdata'),verbose=TRUE)
new<-BIO.02

new$oldNew<-'new'
old$oldNew<-'old'
old$sub_area<-NULL
old$species.n<-NULL

head(old)
head(filter(old,species=='HAK'))
old[old$species=='HAK','species']<-'HKE'
head(filter(new,species=='HKE'))

if (makeGraphs) {
 sort(names(new));sort(names(old))
  new$len<-NULL
  a<-rbind(new,old)
  a$bio<-a$N*a$WSEA
  a<-subset(a,N>0)

  b<-aggregate(bio~species+year+quarter+oldNew,sum,na.rm=T,data=a)

  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Other_bio_',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(bio~as.numeric(year)|
                   paste(species,'Q:',quarter),groups=oldNew,data=x,type='a',lwd=2,
                 scales = list(x=list(relation="same"),y=list(relation='free')),xlab='year',ylab='biomass',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })


  b<-aggregate(N~species+year+quarter+age+oldNew,sum,na.rm=T,data=a)


  by(b,b$species,function(x) {
    trellis.device(device='png',file=file.path(OutputDir,paste0('Other_N_',x[1,'species'],'.png')),width = 1000, height = 800)
    print(xyplot(N~as.numeric(year)|
                   paste(species,"Q:",quarter,"age:",formatC(age, digits = 0, width = 2, format = "f", flag = "0")),groups=oldNew,data=x,type='a',lwd=2,
                 scales = "free",xlab='year',ylab='Number',
                 auto.key =   list(space = "right", points = FALSE, lines = TRUE)))
    cleanup()
  })
}
