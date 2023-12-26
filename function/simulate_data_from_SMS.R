make_survey_data<-function(smult=1,seed=5) {

  set.seed(seed)
  
  ofile<-file.path(data.path,paste0('survey_mult_',smult,'.in'))
  
  if (smult < 0) {
    cat("multiplier is negative, fleet_catch.in is just copied\n")
    file.copy(from=file.path(data.path,'fleet_catch.in'),to=ofile,overwrite = TRUE)
  } else {
  
    indx<-SMS2FLIndices(control=SMS.control)
    effort<-lapply(indx,function(x)as.vector(x@effort))
    fnames<-lapply(indx,function(x)as.vector(x@desc))
    
    
    file<-file.path(data.path,'catch_survey_residuals.out')
    res<-read.table(file,comment.char = "#",header=T) %>%
      filter(data=='survey') %>% mutate(data=NULL) # %>% mutate(s=if_else(s2>0,sqrt(s2),0))
    res$s<-0
    res[res$s2>0,'s']<- sqrt(res[res$s2>0,'s2'])
    # dplyr::select(res,Species.n,Quarter,Age,s) %>% unique()  %>% arrange(Species.n,Quarter,Age)
    
    
    noise<-function(x,s2){
      xx1<- x>0
      x[xx1]<-log(x[xx1])+rnorm(length(x[xx1]), mean=0, sd=sqrt(s2))
      x[xx1]<-exp(x[xx1])
      return(x)
    }
    cpue<-res %>% group_by(Species.n,Quarter,fleet,s)%>% mutate(obs.noise=noise(x=model,s2=(s*smult)^2)) %>% ungroup() %>%
      mutate(ID=paste(formatC(Species.n,width=2,flag='0'),formatC(fleet,,width=2,flag='0'),sep='-')) 
    
  
  
    ids<-sort(unique(cpue$ID))
  
    b<-lapply(ids,function(x)  {
      a<-subset(cpue,ID==x)
      a<-tapply(a$obs.noise,list(a$Year,a$Age),sum)
      return(a)
    })
  
#lapply(b,dimnames)


    cat("# Cpue with noise factor",smult,"\n",file=ofile)
    for (i in (1:length(ids))) {
     # out<-cbind(effort[[i]], b[[i]])
      out<-cbind(1, b[[i]])
      out[is.na(out)]<- -99
      cat(" #",fnames[[i]],"\n",file=ofile,append=TRUE)
      write.table(round(out,2),file=ofile,col.names=FALSE,row.names=FALSE,append=TRUE)
    }
    cat("-999 # Checksum\n",file=ofile,append=TRUE)
  }
}
################  catch

make_catch_data<-function(smult=1,seed=5) {

  set.seed(seed)
  
  ofile<-file.path(data.path,paste0('catch_mult_',smult,'.in'))
  
  if (smult < 0) {
    cat("multiplier is negative, canum.in is just copied\n")
    file.copy(from=file.path(data.path,'canum.in'),to=ofile,overwrite = TRUE)
  } else {
 
  # info<-data.frame(Species.n=first.VPA:nsp,combine_catches=SMS.control@combined.catches,combine_s2=SMS.control@seasonal.catch.s2)
  #  b<-left_join(b,info,by = "Species.n")%>% arrange(parNo)

    file<-file.path(data.path,'catch_survey_residuals.out')
    res<-read.table(file,comment.char = "#",header=T) %>%
       filter(data=='catch' & observed >0) %>% mutate(fleet=NULL,data=NULL,s=sqrt(s2))
    
    # dplyr::select(res,Species.n,Quarter,Age,s) %>% unique()  %>% arrange(Species.n,Quarter,Age)

    noise<-function(x,s2){
      xx1<- x>0
      x[xx1]<-log(x[xx1])+rnorm(length(x[xx1]), mean=0, sd=sqrt(s2))
      x[xx1]<-exp(x[xx1])
      return(x)
    }
    
    #summary(catch)
    #filter(catch,Age==0)
    catch<-dplyr::select(res,Species.n,Year,Quarter,Age,model,s) %>% group_by(Species.n,Quarter,s)%>% mutate(obs.noise=noise(x=model,s2=(s*smult)^2)) %>% ungroup()
    catch[catch$Quarter==9,"Quarter"]<-1L
    
 
    c2<-tapply(catch$obs.noise,list(catch$Species.n,catch$Year,catch$Quarter,catch$Age),sum)
    ages<-SMS.control@first.age:SMS.control@max.age.all
    firstF<-SMS.control@species.info[,"first-age F>0"]
      
    dm<-dimnames(c2)
    
    faca<-min(as.numeric(dm[[4]]  ))
    laca<-max(as.numeric(dm[[4]]  ))
  
   # c2["16","1974",,]
    
    cat("# Catch with noise factor",smult,"\n",file=ofile)
    for (s in dm[[1]]) {
      ff<-firstF[sp.names[as.numeric(s)]]
      for (y in dm[[2]]) {
         out<-c2[s,y,,] 
         cat("# ",sp.names[as.numeric(s)], " Year:",y,"\n",file=ofile,append=TRUE)
  
         for (a in (SMS.control@first.age:SMS.control@max.age.all)) {
           if (a<faca) out<-cbind(0,out)
           if (a>laca)out<-cbind(out,0)
         }
         out[is.na(out)]<- 0
         write.table(round(out,2),file=ofile,col.names=FALSE,row.names=FALSE,append=TRUE)
      }
    }
    cat("-999 # Checksum\n",file=ofile,append=TRUE)
  }
}

