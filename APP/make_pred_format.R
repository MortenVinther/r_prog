s<-Read.species.names()
s<-tail(s,-1)
s

s<-data.frame(no=(1:length(s)),old=s,new=s)
s

if (my.area=='North Sea') {
  i<-s$no %in% ((1:8))
  s[i,"new"]<-'Sea birds'
  
  
  i<-s$no %in% ((10))
  s[i,"new"]<-'Grey gurnard'
  
  i<-s$no %in% ((11:12))
  s[i,"new"]<-'Horse mackerel'
  
  i<-s$no %in% ((14))
  s[i,"new"]<-'Habour porpoise'
  
  #i<-s$no %in% ((22:23))
  #s[i,"new"]<-'Sandeel'
  
  
  #i<-s$no %in% ((24))
  #s[i,"new"]<-'Norway pout'
  
  a<-!duplicated(s$new)
  s$new_no<-unlist(lapply(1:length(a),function(x)sum(a[1:x])))
  s$no<-NULL
s
}

s$group<-"Other predators"
if (my.area=='North Sea') s[first.VPA:npr,'group']<-'VPA.pred'
s[(npr+1):nsp,'group']<-'VPA.prey'
if (my.area=='North Sea') s[(nsp+1-sum(SMS@species.info[,'prey']==0 & SMS@species.info[,'predator']==0)):nsp,'group']<-'flat'
s

adds<-data.frame(old=rep('NONE',2),new=c('Humans','Residual mortality'), new_no=c(0,-1), group=rep('NONE',2))
if (my.area=='Baltic Sea') {s$new_no<-1:dim(s)[[1]]; s$no<-NULL}

s<-rbind(adds,s)

write.csv(s,file=file.path(out_op,'pred_format.csv'),row.names=FALSE)
write.csv(s,file=file.path(data.path,'pred_format.csv'),row.names=FALSE)       

