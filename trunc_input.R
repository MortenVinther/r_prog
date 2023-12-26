trunc_input<-function(outDir=file.path(root,exchangeDir),inclVPA=c(16,22)) {
  la<-SMS.control@max.age.all
  fa<-SMS.control@first.age
  years<-c(1,1)
  years[1]<-SMS.control@first.year
  years[2]<-SMS.control@last.year
  ny<-years[2]-years[1]+1
  npr<-sum(SMS.control@species.info[,'predator']>=1)
  nproth<-sum(SMS.control@species.info[,'predator']>=2)
  nsp<-SMS.control@no.species
  nq<-SMS.control@last.season
  noAreas<-SMS.control@no.areas
  
  #############  catch data
  
  canum<-head(scan(file.path(data.path,'canum.in'),comment.char='#'),-1)
  weca<-head(scan(file.path(data.path,'weca.in'),comment.char='#'),-1)
  Prop.landed<-head(scan(file.path(data.path,'proportion_landed.in'),comment.char='#'),-1)
  Prop.landed<-Prop.landed[1:length(weca)]
  b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:years[2],quarter=1:nq,age=fa:la)
 
  b$species<-sp.names[b$species.n]
  b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
  b<-data.frame(b,CATCHN=canum,WCATCH=weca,PROP_CAT=Prop.landed)
  bb<-subset(b,species.n %in% inclVPA,select=c(year,species.n,species,quarter,sub_area,age,WCATCH,CATCHN,PROP_CAT))
  
  ############## bio data
  
  wsea<-head(scan(file.path(data.path,'west.in'),comment.char='#'),-1)
  if (nproth>0)  wsea<-wsea[((nproth)*noAreas*ny*(la-fa+1)*nq+1):length(wsea)]
  propmat<-head(scan(file.path(data.path,'propmat.in'),comment.char='#'),-1)
  length(wsea)
  m<-head(scan(file.path(data.path,'natmor.in'),comment.char='#'),-1)
  m1<-head(scan(file.path(data.path,'natmor1.in'),comment.char='#'),-1)
  prop_m2<-head(scan(file.path(data.path,'n_proportion_m2.in'),comment.char='#'),-1)
  
  b<-expand.grid(sub_area=1:noAreas,species.n=first.VPA:nsp,year=years[1]:(years[2]),quarter=1:nq,age=fa:la)
  b$species<-sp.names[b$species.n]
  b<-b[order(b$sub_area,b$species.n,b$year,b$quarter,b$age),]
  
  b<-data.frame(b,WSEA=wsea, PROPMAT=propmat,M=m,M1=m1,PROP_M2=prop_m2)
  b<-subset(b,select=c(year,species,species.n,quarter,age,sub_area,WSEA,PROPMAT,M,M1,PROP_M2))
  aa<-subset(b,species.n %in% inclVPA,select=c(year,species.n,species,quarter,sub_area,age,WSEA,PROPMAT,M,M1,PROP_M2))
  
  bb<-left_join(aa,bb,by = join_by(year, species.n, species, quarter, sub_area, age))
  print(head(bb))
  checksum<-function(file='a'){
    cat("-999 # Checksum",file=file,append=TRUE)
  }
  
  attach(bb)
  transf<-function(item,file.name,dig) {

    CC<-tapply(item,list(year,quarter,species.n,age),sum)
    out<-file.path(outDir,file.name)
    unlink(out)
    
    for (sp in as.character(inclVPA)) {
      cat(sp,'\n')
      out1<<-CC[,,sp,]
      for (y in years[1]:years[2]){
        out2<-out1[as.character(y),,]
        out2[is.na(out2)]<--1
        cat(paste("#",sp.names[as.numeric(sp)] ,"year:",y,"\n"),file=out,append=TRUE)
        write.table(format(round(out2,dig),width=11),file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)    
      }
    }
    checksum(file=out)
  }
  
 transf(CATCHN,"canum.in",2)
 transf(WCATCH,"weca.in",5)
 
 transf(WSEA,"west.in",6)
 transf(M,"natmor.in",4)
 transf(M1,"natmor1.in",4)
 transf(PROPMAT,"propmat.in",4)
 transf(PROP_M2,"n_proportion_m2.in",2)
 
 detach(bb)
 
 fname<-'zero_catch_year_season.in'
 z<-head(scan(file.path(data.path,fname),comment.char='#'),-1)
 b<-expand.grid(species.n=first.VPA:nsp,year=years[1]:(years[2]),quarter=1:nq)
 b$species<-sp.names[b$species.n]
 b<-b[order(b$species.n,b$year,b$quarter),]
 b$z<-z 
 aa<-subset(b,species.n %in% inclVPA,select=c(year,species.n,species,quarter,z))

 CC<-tapply(aa$z,list(aa$year,aa$quarter,aa$species.n),sum)
 out<-file.path(outDir,fname)
 unlink(out)
 
 for (sp in as.character(inclVPA)) {
   cat(sp,'\n')
   out1<<-CC[,,sp]
   for (y in years[1]:years[2]){
     out2<-out1[as.character(y),]
     cat(paste("#",sp.names[as.numeric(sp)] ,"year:",y,"\n"),file=out,append=TRUE)
     write.table(format(round(t(out2),0),width=4),file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)    
   }
 }
 checksum(file=out)
 

 
 fname<-'zero_catch_season_ages.in'
 z<-head(scan(file.path(data.path,fname),comment.char='#'),-1)
 b<-expand.grid(species.n=first.VPA:nsp,quarter=1:nq,age=fa:la)
 b$species<-sp.names[b$species.n]
 b<-b[order(b$species.n,b$quarter,b$age),]
 b$z<-z 
 aa<-subset(b,species.n %in% inclVPA,select=c(species.n,species,quarter,age,z))
  CC<-tapply(aa$z,list(aa$quarter,aa$age,aa$species.n),sum)
 out<-file.path(outDir,fname)
 unlink(out)
 for (sp in as.character(inclVPA)) {
   cat(sp,'\n')
   out1<-CC[,,sp]
   cat(paste("#",sp.names[as.numeric(sp)] ,"\n"),file=out,append=TRUE)
   write.table(format(round(out1,0),width=4),file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)    
 }
 checksum(file=out)

 fname<-'recruitment_years.in'
 z<-head(scan(file.path(data.path,fname),comment.char='#'),-1)
 b<-expand.grid(species.n=first.VPA:nsp,year=years[1]:(years[2]))
 b$species<-sp.names[b$species.n]
 b<-b[order(b$species.n,b$year),]
 b$z<-z 
 aa<-subset(b,species.n %in% inclVPA,select=c(species.n,year,z))
 CC<-tapply(aa$z,list(aa$species.n,aa$year),sum)
 out<-file.path(outDir,fname)
 unlink(out)
 write.table(format(round(CC,0),width=3),file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)    
 checksum(file=out)
 
 
 
 fname<-'known_recruitment.in'
 z<-head(scan(file.path(data.path,fname),comment.char='#'),-1)
 b<-expand.grid(species.n=first.VPA:nsp,year=years[1]:(years[2]))
 b$species<-sp.names[b$species.n]
 b<-b[order(b$species.n,b$year),]
 b$z<-z 
 aa<-subset(b,species.n %in% inclVPA,select=c(species.n,year,z))
 CC<-tapply(aa$z,list(aa$species.n,aa$year),sum)
 out<-file.path(outDir,fname)
 unlink(out)
 write.table(format(round(CC,0),width=12),file=out,quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE)    
 checksum(file=out)
}




