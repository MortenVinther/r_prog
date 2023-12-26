
a<-read.csv(file=file.path("C:","_C_drev","SMS-git","NS_2023_04__Boots_0500_haul_mu","Input_Output","Output","WhoEatsWhom","who_eats_whom_level3.csv")) %>% mutate(eatenW=eatenW/1000) %>% rename(mu_w=eatenW)


b<-read.csv(file=file.path("C:","_C_drev","SMS-git","NS_2023_04__Simple_0001_haul_as_observed","Input_Output","Output","WhoEatsWhom","who_eats_whom_level3.csv")) %>% mutate(eatenW=eatenW/1000) %>% rename(prey_w=eatenW)


ab<-left_join(a,b,by = join_by(Year, Predator, Prey, Prey.no)) %>%as_tibble()

birds<-c("Fulmar","Gannet", "GBB.Gull", "Guillemot", "Her.Gull",  "Kittiwake", "Puffin","Razorbill")  
mammals<-c( "Grey.seal","H.porpoise" )
ab<-ab %>% mutate(Predator=if_else(Predator %in% birds,'Birds',Predator)) %>%mutate(Predator=if_else(Predator %in% mammals,'Mammals',Predator)) %>%
  filter( Prey %in% c('Cod','Herring') )
ab<- ab %>% group_by(Year,Predator,Prey) %>% summarize(mu_w=sum(mu_w), prey_w=sum(prey_w))

ab<-ab  %>% mutate(ratio=mu_w/prey_w)



x<- ggplot(filter(ab,Prey=='Cod'), aes(x=prey_w, y=mu_w,group=Predator)) +
  geom_point()+
  geom_abline(col='red')+facet_wrap(~Predator, scales='free')+
  labs(title='Cod',x = "Eaten biomass (kt) non bootstrap",y="Eaten biomass (tonnes) from bootstrap")
png(filename=file.path(data.path,'cod.png'),width=700,height=700,pointsize=25)
print(x);cleanup()


x<- ggplot(filter(ab,Prey=='Herring'), aes(x=prey_w, y=mu_w,group=Predator)) +
  geom_point()+
  geom_abline(col='red')+facet_wrap(~Predator, scales='free')+
  labs(title='Herring',x = "Eaten biomass (kt) non bootstrap",y="Eaten biomass (tonnes) from bootstrap")
png(filename=file.path(data.path,'herring.png'),width=700,height=700,pointsize=25)
print(x);cleanup()


x<- ggplot(filter(ab,Prey=='Cod'), aes(x=Year, y=mu_w/prey_w,group=Predator)) +
  geom_point()+
  geom_abline(col='red')+facet_wrap(~Predator, scales='free')+
  labs(x = "non bootstrap weight proportion",y="proportion from alpha")
png(filename=file.path(output_dir,'bias_11b.png'),width=700,height=700,pointsize=25)
print(x);cleanup()

