
library(tidyverse)
a<-expand.grid(y=1993:1998,q=1:4,pred=1:3,predAge=1:5,prey=1:2,preyAge=1:4) %>% 
  filter(predAge>preyAge) %>%
  mutate(predW=predAge*1.2,preyW=1.5) %>% as_tibble()

dim(a)
a
a1<-a %>% nest(data=c(pred,predAge,predW,prey,preyAge,preyW))
a1
str(a1,1)

a1[1,'data'][[1]][[1]]

rm(a2)
a2<-a1 %>% rowwise() %>% mutate(data=list(data %>% nest(data=c(prey,preyAge,preyW))))
a2

a2[,'y']
a2[1,'data'][[1]][[1]]
a2[1,'data'][[1]][[1]][['data']]
a2[1,'data'][[1]][[1]][['data']][[2]]

# the same but condensed

r<-a %>% nest(data=c(pred,predAge,predW,prey,preyAge,preyW)) %>% 
  rowwise() %>% mutate(data=list(data %>% nest(data=c(prey,preyAge,preyW))))


allRec<-0
nyq<-dim(r)[[1]]
for (yy in (1:nyq)) {
  y<-unlist(r[yy,'y']) 
  q<-unlist(r[yy,'q'])
  cat('y;',y,' q:',q,'\n')
  npp<-dim(r[yy,'data'][[1]][[1]])[[1]]
  for (pp in (1:npp)) {
    pred<-unlist(r[yy,'data'][[1]][[1]][pp,'pred'])
    predAge<-unlist(r[yy,'data'][[1]][[1]][pp,'predAge'])
    predW<-unlist(r[yy,'data'][[1]][[1]][pp,'predW'])
    cat('     pred:',pred,' predAge:',predAge,'\n')
    diet<-r[yy,'data'][[1]][[1]][pp,'data'][[1]][[1]]
    #print(diet)
    allRec<-allRec+dim(diet)[[1]]
  }
}
cat(allRec,dim(a)[[1]])


