plot_alk<-function(dir=data.path,   output.dir=data.path, my.year=2010){
  

a<-Read.ALK.stom(dir=dir) %>% filter(year==my.year) 
b<-Read.stom.size.classes(dir =dir )  %>% filter(year==my.year) %>% 
  rename(Species.no=species.no,LengthGroup=size.class) %>% mutate(species=NULL)
a<-left_join(a,b,by = join_by(year, quarter, Species.no, LengthGroup)) %>% 
  mutate(quarter=paste0("Q:",quarter),size=substr(size,1,3))


sort(unique(a$Species))
head(a)


make_plots<-T
if (make_plots) {
  age<-sort(unique(a$Age))
  age<-data.frame(age,Age2=paste('Age',age))
  a$age= factor(a$Age,levels=age$age,labels=age$Age)
  xx<-by(a,list(a$Species.no),function(x) {
    prey_name<-x[1,'Species']
    p<-ggplot(data=x, aes(x=size, y=proportion)) +
      geom_bar(stat="identity", color='red',width=0.5)+
      facet_grid(rows=vars(age),cols=vars(quarter),scales = "free_y")+
      labs(x='Length class (cm)',y='Proportion',title=prey_name) +
      #theme_minimal() +
      theme( #panel.grid.major = element_line(linetype = "blank"),
        #panel.grid.minor = element_line(linetype = "blank"),
        axis.text.x = element_text(angle = 90, vjust = 0.5))
    
    
    png(filename=file.path(output.dir,paste0('ALK_',prey_name,'.png')),width=800,height=1000,pointsize=35)
    print(p)
    cleanup()
    return(NULL)
  })
}
}


#plot_alk(dir=data.path,   output.dir=data.path, my.year=2010)
  
  
