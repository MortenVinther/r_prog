a<-Read.summary.data(dir=file.path(root,"NS_2023_04")) %>% 
   select(Species,Year,Quarter,Species.n,Age,C.obs,weca  ) %>% rename(C.obs.a=C.obs, weca.a=weca)

b<-Read.summary.data(dir=file.path(root,"NS_2023_key_run")) %>% 
  select(Species,Year,Quarter,Species.n,Age,C.obs,weca  ) %>% rename(C.obs.b=C.obs, weca.b=weca)

ab<-left_join(a,b,by = join_by(Species, Year, Quarter, Species.n, Age))

filter(ab, C.obs.a != C.obs.b)

tst<-filter(ab, weca.b== -9) 
filter(tst,weca.a>0)


# SOP
a<-Read.summary.table(dir=file.path(root,"NS_2023_04")) %>% 
  select(Species,Year,Species.n,SOP ) %>% rename(SOP.a=SOP)

b<-Read.summary.table(dir=file.path(root,"NS_2023_key_run")) %>% 
  select(Species,Year,Species.n,SOP ) %>% rename(SOP.b=SOP)


ab<-left_join(a,b,by = join_by(Species, Year, Species.n))

filter(ab, SOP.a != SOP.b) %>% mutate(difPercent=round((SOP.a-SOP.b)/SOP.b*100,2))


