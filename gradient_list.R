g<-Read.gradient.dat(gfile="run_ms2.grd")
g$parNo<-1:dim(g)[[1]]
g<-g[order(abs(g$Gradient)),]
print(tail(g,10))

if (FALSE) {
  file<-file.path(data.path,"par_exp.out")
  s<-read.table(file,header=TRUE)
  head(s)
  
  a<-merge(x=g,y=s,all.x=T)
  head(a)
  a$Sp<-sp.names[a$species]
  a<-a[order(abs(a$Gradient)),]
  print(tail(a,10))
  
}
