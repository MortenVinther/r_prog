plot_alpha0_Dirichlet<-function(dir=data.path,   output.dir=data.path,inclSp=NULL)  {
  stom<-Read.stomach.data(dir=dir,read.init.function=T) %>% rename(Area=SMS.area)
  stom_st<-Read.stomach.data.start(dir=dir,read.init.function=T) 
  
  a<-stom %>% select(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter,N.haul, N.samples,Diri.sum.p ) %>% unique()
  b<-stom_st %>% select(Predator.no, Predator, Predator.length.class, Predator.length,Year,Quarter,alfa0 ) %>% unique()
  a<-left_join(a,b,by = join_by(Predator.no, Predator, Predator.length.class, Predator.length, Year, Quarter))
 

  if (!is.null(inclSp)) a<-filter(a,Predator %in% inclSp)
  
  allalfa<-a %>% group_by(Predator) %>% summarize(inputalfa=all(alfa0>0))
  a<-a %>% group_by(Predator.no) %>% mutate(inputalfa=all(alfa0>0))
  
  b<-read.table(file.path(data.path,'pred_format.dat'),header=TRUE) %>% rename(Predator=old)
  b<-left_join(b,allalfa,by = join_by(Predator)) %>% mutate(inputalfa=if_else(is.na(inputalfa),FALSE,inputalfa),new2=if_else(inputalfa,paste0(new,'*'),new))
  
  bf<-b %>% select(new,new2,newno) %>% unique %>% filter(newno<99)%>% arrange(newno)
  
  a<-left_join(a,b) %>% mutate(Predator=new2) %>%
    group_by( Predator,newno,inputalfa)  %>% mutate(meanAlpha0=mean(Diri.sum.p),medianAlpha0=median(Diri.sum.p),alfa0=NULL) %>% 
    ungroup() %>%
    arrange(newno,Predator,Predator.length,Year,Quarter) %>%rename(alpha0=Diri.sum.p) %>% mutate(Length=(Predator.length))

  a$Predator<-factor(a$Predator,levels=bf$new2)
  l<-data.frame(Length=c(100,  120,  150,  200,  250,  300,  350,  400,  500,  600,  700,  800, 1000, 1200),
                lNew=  c(10,  10,  10,  20,  20,  30,  30,  40,  50,  50,  70,  70, 100, 100))
  
  a<-left_join(a,l,by = join_by(Length)) %>% mutate(Length=factor(lNew)) 
  

  png(filename=file.path(output.dir,'alpha0.png'),width=600,height=800)
  pp<-  ggplot() + geom_histogram(data = a, aes(x = alpha0,  color=Length,fill=Length)) +
    #geom_vline(data =a, mapping = aes(xintercept = meanAlpha0) ,color="blue", linetype="dashed", linewidth=1)+
    geom_vline(data =a, mapping = aes(xintercept = medianAlpha0) ,color="blue", linetype="dashed", linewidth=1)+
    facet_grid(vars(Predator),scales = "free")+
    theme(strip.text.y = element_text(angle = 0))
  print(pp)
  cleanup()
}

#plot_alpha0_Dirichlet(dir=data.path,   output.dir=data.path,inclSp=c("Cod", "Whiting", "Haddock","Saithe","Mackerel","N.horse.mac","A.radiata","G.gurnards","W.horse.mac"))  

#plot_alpha0_Dirichlet(dir=data.path,   output.dir=data.path)  
