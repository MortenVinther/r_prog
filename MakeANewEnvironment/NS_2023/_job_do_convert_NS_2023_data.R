# this script makes
#  1. write a number of of files on a spreadsheet like format with data extracted from the current (old)  environment
#  2. read the ICES single species assessment input/output
#  3. update the old data set with new data
#  4. convert files produced under 1) to files with data on SMS format

# between 1) and 2) you are supposed to update data !

# start by submitting init.R, (for the old keyrun or source) so file paths are established

do.1<-FALSE
do.2<-TRUE
do.stock.dist.basic<-FALSE
add_cod_Mmagic<-TRUE

newEnv<-'NS_2023'
exchangeDir<-file.path( "Data_NorthSea","input_NS_2023")   # directory with data  files on spread sheet format
RexchangeDir<-"MakeANewEnvironment"     # directory with R scripts to convert files etc.
finalExchangeDir<-file.path(root,"Data_NorthSea","final_input_NS_2023")

addfn<-'key2023'  # a label for files
old_key_label<-'key_2020'

firstY<-1974    # first year in key run
lastY<-2022     # last year in (new) key run
newYear<-firstY:lastY
lastYold<-2019  # last year in previous key run
OutputDir<-file.path(root,'Output')
makeGraphs<-TRUE

new.code.name<-c("FUL","GLT","HEG","KTW","GBG","GNT","PUF","RAZ","RAJ","GUR","W_H","N_H","GSE","HBP",'HKE','COD','WHG','HAD','POK','MAC','HER','NSA','SSA','NOP','SPR','PLE','SOL')
options(stringsAsFactors = FALSE)


#  first 1.

if (do.1) {  # write a number of of files on a spreadsheet like format with data extracted from the current (old) environment.
  # You have first to run init.R with the old environment where data come from
  source(file.path(prog.path,RexchangeDir,newEnv,'from_sms_format_to_list.r'))
  my.code.name<- new.code.name
  #my.code.name<-c("MAC")
  From_SMS_format_to_list(otherPredExist=T,catchMultiplier=1,code.name=my.code.name,exchangeDir=file.path(root,exchangeDir),addfn=old_key_label)
  old_ages<-data.frame(species=new.code.name,plusage=SMS.control@species.info[,'last-age'],use_plus=SMS.control@species.info[,'+group']==1)
  save(old_ages,file=file.path(root,exchangeDir,'old_ages.Rdata'))
}

load(file=file.path(root,exchangeDir,'old_ages.Rdata'),verbose=TRUE)

###   Change directory to new key-run
my.stock.dir<-'NS_2023_04'
data.path<-file.path(root,my.stock.dir)

#####################################################################################
# Then 2, read ICES assessment data and update SMS data
if (do.2) {
   source(file.path(prog.path,RexchangeDir,newEnv,"read_update_assessment_data_NS_2023.R"))
  #catch  round(tapply(CAT.01$CATCHN*CAT.01$WCATCH,list(CAT.01$year,CAT.01$species),sum))
}

######################################################################################

# stock distribution
if (do.stock.dist.basic) source(file.path(prog.path,RexchangeDir,newEnv,'stock_distribution_NS_2023sp.r')) # Stock distribution by stock, produces maps and cod Rdata

# Other predators
source(file.path(prog.path,RexchangeDir,newEnv,'Horse_mack.R'))
source(file.path(prog.path,RexchangeDir,newEnv,'hake.R'))

source(file.path(prog.path,RexchangeDir,newEnv,'other_predators.R'))

################
# mean length
source(file.path(prog.path,RexchangeDir,newEnv,'mean_length.R'))

source(file.path(prog.path,RexchangeDir,newEnv,'stock_distribution_NS_2023.r')) # Stock distribution by stock


################
# consum
source(file.path(prog.path,RexchangeDir,newEnv,'consum.R'))
#

#correct where presently used age range differs from previously used

load(file=file.path(root,exchangeDir,'old_ages.Rdata'),verbose=T)
new_ages<-data.frame(species=new.code.name,plusage=SMS.control@species.info[,'last-age'],use_plus=SMS.control@species.info[,'+group']==1)
old_ages==new_ages

if (!all(old_ages==new_ages)) {
  # start with a cleanup
  years=firstY:lastY
  load(file.path(root,exchangeDir,'BIO_01.Rdata'),verbose=T)
  BIO.01<-subset(BIO.01,year %in% years,select=c(species,year,age,quarter,sub_area,WSEA,PROPMAT,M,M1,PROP_M2))
  a<-left_join(BIO.01,old_ages,by = join_by(species)) %>%as_tibble() %>% mutate(clean=age>plusage,notBorn=age==0 & quarter<=2)
  crit<-a$clean | a$notBorn
  a[crit,'WSEA']<- NA
  a[crit,'PROPMAT']<- NA
  a[crit,'M']<- NA
  a[crit,'M1']<- NA
  a[crit,'PROP_M2']<- NA
  a<-a %>% mutate(plusage=NULL,clean=NULL,notBorn=NULL, use_plus=NULL)

  # new ages
  a<-left_join(a,new_ages,by = join_by(species))
  a<-a %>% mutate(w_fac=if_else(age>plusage,exp(0.5*(plusage-age)),1))
  filter(a,age>plusage & !is.na(WSEA) &age>6 &year==1977)

  b<- a %>% mutate(age=if_else(age>plusage,plusage,age)) %>%
          group_by(species,year,age, quarter,sub_area) %>%
          summarize(WSEA=weighted.mean(x=WSEA,y=w_fac,na.rm=TRUE),
              M=  weighted.mean(x=M,y=w_fac,na.rm=TRUE),
              PROPMAT=weighted.mean(x=PROPMAT,y=w_fac,na.rm=TRUE),
              M1=weighted.mean(x=M1,y=w_fac,na.rm=TRUE),
              PROP_M2=weighted.mean(x=PROP_M2,y=w_fac,na.rm=TRUE)) %>% ungroup()


  filter(b,species=='WHG' & year==1974 & quarter==1)
  BIO.01<-b %>% mutate(across(where(is.numeric), ~replace(., is.na(.), -1)))
  save( BIO.01,file=file.path(root,exchangeDir,'BIO_01.Rdata'))


  # CAT 01

  # start with a cleanup
  years=firstY:lastY
  load(file.path(root,exchangeDir,'CAT_01.Rdata'),verbose=T)
  head(CAT.01)
  CAT.01<-subset(CAT.01,year %in% years,select=c(species,year,age,quarter,CATCHN,WCATCH,PROP_CAT,oldNew))
  a<-left_join(CAT.01,old_ages,by = join_by(species)) %>%as_tibble() %>% mutate(clean=age>plusage,notBorn=age==0 & quarter<=2)
  crit<-a$clean | a$notBorn
  a[crit,'WCATCH']<- NA
  a[crit,'CATCHN']<- NA
  a[crit,'PROP_CAT']<- NA
  filter(a,species=='WHG' & year==1974 & quarter==1)
  a<-a %>% mutate(plusage=NULL,clean=NULL,notBorn=NULL, use_plus=NULL)

  # new ages
  a<-left_join(a,new_ages,by = join_by(species))
  tst<-filter(a,age>plusage & !is.na(WCATCH))
  tst
  sort(unique(tst$species))
  filter(a,age>plusage & !is.na(WCATCH) &age>6 &year==1974)

  b<- a %>% mutate(age=if_else(age>plusage,plusage,age)) %>%
    group_by(species,year,age, quarter,oldNew) %>%
    summarize(WCATCH=weighted.mean(x=WCATCH,y=CATCHN,na.rm=TRUE),
              PROP_CAT=  weighted.mean(x=PROP_CAT,y=CATCHN,na.rm=TRUE),
              CATCHN=sum(CATCHN,na.rm=TRUE)) %>% ungroup()

  if (FALSE) {
    dim(b)
    tst<- b

    dups<- duplicated(select(tst,species,year,age,quarter))
    aa<-tst[dups,]
    aa
    tst$dups<-dups
    tst[dups,]
    dim(tst); dim(unique(tst))
    dups<-tst[duplicated(tst),]
    xtabs(~ year+species+age,data=dups)
  }


  filter(b,species=='WHG' & year==1974 & quarter==1)
  CAT_01<-b %>% mutate(across(where(is.numeric), ~replace(., is.na(.), -1)))
  save( CAT_01,file=file.path(root,exchangeDir,'CAT_01.Rdata'))
}

#write.table(BIO.01,file=file.path(finalExchangeDir,paste0('VPA_Bi01.IN')),row.names = F,quote = T,sep=',')



writeLists<-function(years=firstY:lastY) {
  ages<- new_ages %>% select(species, plusage)

  # BIO 01
  load(file.path(root,exchangeDir,'BIO_01.Rdata'),verbose=T)
  sort(unique(BIO.01$species))
  BIO.01<-left_join(new_ages,BIO.01,by = join_by(species)) %>% filter(age<=plusage)
  BIO.01<-subset(BIO.01,year %in% years,select=c(species,year,age,quarter,sub_area,WSEA,PROPMAT,M,M1,PROP_M2))
  write.table(BIO.01,file=file.path(finalExchangeDir,paste0('VPA_Bi01.IN')),row.names = F,quote = T,sep=',')

  # BIO 02
  load(file.path(root,exchangeDir,'BIO_02.Rdata'),verbose=T)
  sort(unique(BIO.02$species))
  BIO.02<-left_join(new_ages,BIO.02,by = join_by(species)) %>% filter(age<=plusage)
  BIO.02<-subset(BIO.02,year %in% years)
  write.table(BIO.02,file=file.path(finalExchangeDir,paste0('VPA_Bi02.IN')),row.names = F,quote = T,sep=',')

  # CAT 01
  load(file.path(root,exchangeDir,'CAT_01.Rdata'),verbose=T)
  #round(tapply(CAT.01$CATCHN*CAT.01$WCATCH,list(CAT.01$year,CAT.01$species),sum))
  #ftable(round(tapply(CAT.01$CATCHN*CAT.01$WCATCH,list(CAT.01$year,CAT.01$species,CAT.01$quarter),sum)))
  CAT.01<-left_join(new_ages,CAT.01,by = join_by(species)) %>% filter(age<=plusage)
  CAT.01<-subset(CAT.01,year %in% years,select=c(species,year,age,quarter,CATCHN,WCATCH,PROP_CAT)) %>%
    mutate(PROP_CAT=if_else(CATCHN>0,PROP_CAT,0), WCATCH=if_else(CATCHN>0,WCATCH,0))
  filter( CAT.01,)
  write.table(CAT.01,file=file.path(finalExchangeDir,paste0('VPA_Ca01.IN')),row.names = F,quote = T,sep=',')
}
writeLists()


########################################################################################
if (FALSE) {
s<-subset(SS,species=='W_H')
head(s)
subset(s,year==1986)
tapply(s$WSEA,list(s$year,s$quarter,s$species,s$age),sum)
WW[as.character(1986:1990),,'W_H',as.character(0:2)]
setwd(data.path)
source(file.path(prog.path,RexchangeDir,newEnv,'write_surveys.R'))
}

