fname<-'zero_catch_year_season.in'
z<-head(scan(file.path(dir,fname),comment.char='#',quiet = TRUE),-1)
bz<-expand.grid(s=(first.VPA:nsp) +off.species,y=1:nYears,q=1:nSeasons)
bz<-bz[order(bz$s,bz$y,bz$q),]
bz$z<-z 

yield<-inp_all[['data']]$catchNumber # structure
for (i in seq_along(yield)) yield[[i]]<-inp_all[['data']]$catchMeanWeight[[i]]*inp_all[['data']]$catchNumber[[i]]

lapply(yield,function(x) apply(x,c(1,2),sum))

       
