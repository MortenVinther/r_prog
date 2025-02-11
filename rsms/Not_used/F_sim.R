Nini<-1000
CN<-c(9,126,34,0)
Z=-19

cn<-data.frame(q=1L:4L,CN)
scen<-expand.grid(M=seq(0.1,0.4,by=0.01),q=1L:4L)
scen<-left_join(scen,cn,by = join_by(q))%>% arrange(M,q)
head(scen)

Ffunc<-function(x,N,CN,M) abs(CN- (N*(1-exp(-(x+M))) /(x+M)*x))

aa<-lapply(1:dim(scen)[[1]],function(s) {
  q<-scen[s,'q']; 
  if (q==1) N<<-Nini else N<<-N*exp(-Z)
  cat(N,'\n')
  M<-scen[s,'M']; 
  CN<-scen[s,'CN']
   
  ff<-optimize(Ffunc,N=N,M=M,CN=CN,interval=c(0,4))
  FF<-ff$minimum
  Z<<-M+FF
  out<-data.frame(M,q,FF,CN,Z,N)
  out
})
a<-do.call(rbind,aa)

b<-a %>% mutate(Nbar=N*(1-exp(-Z))/Z,DF=Nbar*FF,DM=Nbar*M,Nnext=N*exp(-Z)) %>% as_tibble() %>%
  group_by(M) %>% mutate(pCN=CN/sum(CN),pFF=FF/sum(FF),rattoCNf=pCN/pFF,q=factor(q),cumM=cumsum(M),cumCN=cumsum(pCN) )
b

ggplot(data=filter(b,CN>2), aes(x=M, y=rattoCNf,col=q,shape=q,group=q)) +
    geom_line()+geom_point()
