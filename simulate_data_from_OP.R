who_eats_whom(
    first.year=1974,                #first year on plot, negative value means value defined by data
    last.year=2050,                  #last year on plot
    OperatingModel=TRUE,   # include data from forecast (the OP )
    op.dir=data.path,
    makeAllGraphs=FALSE,
    output.dir=data.path, 
    my.colors=c('red','green','plum','blue','cyan','yellow','coral','skyblue','purple','magenta','limegreen','pink' )
)

d<-read_csv(file=file.path(data.path,"who_eats_whom_level1.csv"))  

size<-Read.OP.size(dir=data.path,infile='op_size.in')
head(size)
pred_l<-filter(size,type=="length") %>%mutate(type=NULL) %>% rename(pred_l=size,Predator=Species)
pred_l
str(pred_l)

d<-left_join(d,pred_l)
