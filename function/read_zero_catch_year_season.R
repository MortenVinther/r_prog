a<-scan(file.path(data.path,"zero_catch_year_season.in"),comment='#')
a<-head(a,-1)
tail(a)
yq<-data.frame(expand.grid(Quarter=1:SMS.control@last.season,Year=SMS.control@last.year:SMS.control@first.year,Species.n=first.VPA:nsp),yq=a) 
summary(yq)           

yq<-yq %>% group_by(Species.n,Quarter) %>% summarize(yq=all(yq==0))

print(yq,n=40)

#filter(yq,Species.n==23)

if (SMS.control@zero.catch.year.season==1)
