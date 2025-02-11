
###################################################
# consum (ration)

a<-read.csv(file=file.path(root,exchangeDir,'consum_key_2020.csv'))
head(a)
sort(unique(a$species))
#a[a$species=='HAK','species']<-'HKE'

# just add the new years, re-using old data
b<-subset(a,year==lastYold)
for (y in ((lastYold+1):lastY)) {
  b$year<-y
  a<-rbind(a,b)
}
head(a)

b<-filter(a,species=='W_H' & age<6 & year==1974)
xtabs(CONSUM~quarter+age,data=b)
head(a)




new_ages<-data.frame(species=new.code.name,plusage=SMS.control@species.info[,'last-age'],use_plus=SMS.control@species.info[,'+group']==1)
old_ages==new_ages

if (!all(old_ages==new_ages)) {
  # start with a cleanup
  years=firstY:lastY
  a<-subset(a,year %in% years)
  a<-left_join(a,old_ages,by = join_by(species)) %>%as_tibble() %>% mutate(clean=age>plusage,notBorn=age==0 & quarter<=2)
  crit<-a$clean | a$notBorn
  a[crit,'CONSUM']<- NA
  a<-a %>% mutate(plusage=NULL,clean=NULL,notBorn=NULL, use_plus=NULL)

  # new ages
  a<-left_join(a,new_ages,by = join_by(species))
  a<-a %>% mutate(w_fac=if_else(age>plusage,exp(0.5*(plusage-age)),1))
  filter(a,age>plusage & !is.na(CONSUM) &age>6 &year==1977)

  b<- a %>% mutate(age=if_else(age>plusage,plusage,age)) %>%
    group_by(SMS_area,species,species.n,year,age, quarter) %>%
    summarize(CONSUM=weighted.mean(x=CONSUM,y=w_fac,na.rm=TRUE)) %>% ungroup()


  filter(b,species=='WHG' & year==1974 & quarter==1)
  CONSUM<-b %>% mutate(across(where(is.numeric), ~replace(., is.na(.), -1)))
}
save(CONSUM,file=file.path(root,exchangeDir,'CONSUM.Rdata'))


write.table(CONSUM,file=file.path(finalExchangeDir,paste0('consum.dat')),row.names = F,quote = T,sep=',')

if (FALSE) { # do not update

    ### consum ab
    #REMEMBER TO CHECK IF IT INCLUDES THE RIGHT SPECIES

    #file.copy(from=file.path(data.path,'consum_ab.in'),to=file.path(data.path,'consum_ab_OLD.in'),overwrite=TRUE)

    a<-scan(file=file.path(root,exchangeDir,'consum_ab_key_2017.in'),comment.char = "#")
    # test a<-1:length(a)
    b<-array(a,dim=c(2,4,npr))
    dimnames(b)<-list(c('a','b'),paste0('quarter',1:4),new.code.name[1:npr])
    round(ftable(b),1)


    # read new values
    bb<-read.csv(file=file.path(root,exchangeDir,'consum_AB_2017.csv'))

    for (i in (1:dim(bb)[[1]])) {
      b[1,bb[i,'quarter'],bb[i,'Species']]<-bb[i,'FOODA']
      b[2,bb[i,'quarter'],bb[i,'Species']]<-bb[i,'FOODB']
    }


    # use W.horse mac values for the North Sea horse mac
    b[,,'N_H']<-b[,,'W_H']

    b[,"quarter2",'N_H']<-b[,"quarter3",'N_H']
    b['a',"quarter2",'N_H']<-b['a',"quarter3",'N_H']*0.7 # guessing


    # write new values
    out<-file.path(data.path,'consum_ab.in')
    cat("# paramter a and b for consumption=a*weight^b\n",file=out)
    for (s in (1:npr)) {
      cat("# ",sp.names[s],'\n',file=out,append=TRUE)
      for (q in (1:4)) {
        cat(round(b[1,q,s],5), round(b[2,q,s],5), '  # quarter ',q,'\n',file=out,append=TRUE)
      }
    }
    cat("-999 # check",'\n',file=out,append=TRUE)
}
