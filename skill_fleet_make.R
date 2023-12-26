skill_fleet_make<-function(control,path=NULL,fleet.inf="fleet_info.dat",fleet.index="fleet_catch.in",
                        fleet.name="fleet_names.in",make_excl=TRUE,write_data=TRUE,exclFile='skill_exclude.in',adj) {
  
  old.wd<-getwd()
  if (!is.null(path)) setwd(path)

  if (missing(adj)) change<-FALSE else change<-TRUE

  firstYear<-slot(control,"first.year")
  lastYear<-slot(control,"last.year")
  years<- firstYear:(lastYear+1)
  nsp<-slot(control,"no.species")
  nq<-slot(control,"last.season")
  
  #count number of other predators    
  info<-slot(control,"species.info")[,"predator"]
  no.oth<-sum(info==2) 
  nsp<-nsp-no.oth  
  
  s<-readLines(fleet.name, n=1000)
  s<-gsub('_',' ',s)
  fl.names<-sub('[[:space:]]+$', '', s)
  
  info<-scan(fleet.inf,comment.char = "#",quiet=TRUE) 
  minCV<-info[1]
  i<-2
  n.fleet<-as.vector(info[i:(i-1+nsp)])
  i<-i+nsp
  sum.fleet<-sum(n.fleet)
  fl.info<-matrix(info[i:(i-1+sum.fleet*10)],ncol=10,nrow=sum.fleet,byrow=TRUE)

  i<-1
  sp.fl<-0
  fileName<-'skill_cpue.in'
  if (write_data){
    cat('# file for skill assessment with selection of fleet data. \n',file=fileName)
    cat('# 0 is no inclsion of obs in liklihood, 1  is inclusion, and -9 is otside range of fleet data\n',file=fileName,append=TRUE)
  }
  skill<-years
  names(skill)<-as.character(years)
  
  excl<-NULL
  
  for (sp in 1:nsp) {
    if (write_data) {
      cat('# ',sp.names[sp+ no.oth],'\n',file=fileName,append=TRUE)
      cat('# ',formatC(years,width=4),'\n',file=fileName,append=TRUE)
    }
    for (fl in 1:n.fleet[sp]) {
      skill[]<- -9
      
      sp.fl<-sp.fl+1
      fy<-fl.info[sp.fl,1]
      ly<-fl.info[sp.fl,2]
      skill[as.character(fy:ly)]<- 1
      
      if (change) {
        pot<-filter(adj,sp2==sp & fl2==fl & excl==0)
        if (dim(pot)[[1]]>0) {
          skill[as.character(pot$year2)]<- 0
        }
      }
      excl<-rbind(excl,data.frame(sp2=sp,fl2=fl,excl=1,year2=fy:ly))
    
      if (write_data) cat('  ',formatC(skill,width=4),' #', fl.names[sp.fl],'\n',file=fileName,append=TRUE)
    } 
  } 
  if (make_excl) write_csv(excl,file=exclFile)
  if (write_data) cat(-999,' # check value\n',file=fileName,append=TRUE)
  setwd(old.wd)

}
