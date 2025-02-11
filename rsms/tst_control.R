RSMS.control <- function(  
    RSMS=NULL,
    first.year=1900,
    first.year.model=1900,
    last.year=first.year+1, 
    last.season=1,  
    no.species=1,
    no.VPA.predators=0,
    no.other.predators=0,
    species.names=c("sp1"),
    first.age=0, 
    last.age=rep(9,no.species+no.other.predators),
    max.age.all=10)
{
  if (TRUE) {
    RSMS=NULL
    first.year=1976
    last.year=2000
    first.year.model=1900
    no.species=3
    last.season=1
    no.other.predators=1
    no.VPA.predators=1
    species.names=c('Bird','Cod','Herring')
    last.age=c(1,11,9)
  }
  if (is.null(RSMS)){
    tmplt<-new("RSMS.control")
    if (no.species == 1) {
      no.VPA.sp<-no.species
      no.predators<-0
      no.other.predators<-0
    }
    if (no.species>1) {
      if (no.other.predators>=no.species) stop("no.other.predators cannot be larger than no.species")
      no.VPA.sp<-no.species-no.other.predators
      no.predators<-no.other.predators+no.VPA.predators
    }
    if (species.names[1] != c("sp1") & length(species.names)!=no.species) 
      stop("no.species is diffrent from number of species names")
    
    species.info<-matrix(0,ncol=ncol(tmplt@species.info),nrow=no.species,dimnames=list(species.names,colnames(tmplt@species.info)))
    species.info[,1]<-last.age
    species.info[,2]<-first.age
    species.info[,3]<-0
    species.info[,7]<-c(rep(0,no.predators),rep(0,no.VPA.sp-no.VPA.predators))
    species.info[,8]<-c(rep(0,no.other.predators),rep(1,no.VPA.sp))
    species.info[,9]<-c(rep(0,no.other.predators),rep(3,no.VPA.sp))
    species.info[,10]<-c(rep(0,no.other.predators),rep(1,no.VPA.sp))
    species.info[,11]<-c(rep(0,no.other.predators),rep(0,no.VPA.sp))
    catch.s2.group<-vector("list", length=no.VPA.sp);
    for (j in 1:no.VPA.sp) catch.s2.group[[j]]<-as.integer(first.age:(first.age+2)) 
    catch.season.age<-vector("list", length=no.VPA.sp)
    for (j in 1:no.VPA.sp) {
      if (last.season>1) catch.season.age[[j]]<-as.integer(first.age:(first.age+1)) else 
        catch.season.age[[j]]<-as.integer(first.age)
    }
    catch.sep.year<-vector("list", length=no.VPA.sp)
    for (j in 1:no.VPA.sp) catch.sep.year[[j]]<-as.integer(first.year.model); 
    catch.spline.year<-vector("list", length=no.VPA.sp)
    for (j in 1:no.VPA.sp) catch.spline.year[[j]]<-as.integer(first.year.model);
    zero.catch.year.season<-as.integer(0)
    zero.catch.season.age<-as.integer(0)
    
    res <- new("RSMS.control",
               test.output=as.integer(0) ,,
               VPA.mode=as.integer(0),
               no.areas=as.integer(1),
               first.year=as.integer(first.year), 
               first.year.model=as.integer(first.year.model), 
               last.year=as.integer(last.year),
               last.year.model=as.integer(last.year),
               last.season=as.integer(last.season),
               no.species=as.integer(no.species),
               first.age=as.integer(first.age),
               max.age.all=as.integer(max.age.all),
               species.names=species.names,
               species.info=species.info,
               beta.cor        =rep(1E6,no.VPA.sp),
               SSB.R.year.first=rep(first.year.model,no.VPA.sp),
               SSB.R.year.last=rep(last.year,no.VPA.sp)  ,
                 seasonal.catch.s2=rep(0,no.VPA.sp),                                                 
               catch.s2.group=catch.s2.group,
               catch.season.age=catch.season.age,
               avg.F.ages      =matrix(0,ncol=2,nrow=no.VPA.sp,
                                       dimnames=list(species.names[1:no.VPA.sp],c("first-age","last-age"))),
               min.catch=rep(-5,no.VPA.sp),
               catch.sep.year=catch.sep.year,
               catch.spline.year=catch.spline.year,
               zero.catch.year.season=zero.catch.year.season,
               zero.catch.season.age=zero.catch.season.age,
               size.selection      =rep(0,no.predators),
               var.scale.stom      =rep(0,no.predators),
               size.other.food.suit=rep(0,no.predators),                                
               min.stom.cont       =rep(1E-4,no.predators),
               max.stom.sampl      =rep(1E4,no.predators),
               prey.pred.size.fac  =rep(0.5,no.predators),
               stom.type.include   =rep(1,no.predators)                               
               
    )
  }
  else {  # We re-use an RSMS.control object 
    if (!inherits(RSMS, "RSMS"))
      stop("RSMS must be an 'RSMS' object!")
    
    res <- RSMS@control
    
    # ... and possibly redefine some of its parameters
    if (!missing(test.output))
      res@test.output <- test.output
    if (!missing(VPA.mode))
      res@VPA.mode <- as.integer(VPA.mode)
    if (!missing(no.areas))
      res@no.areas <- as.integer(no.areas)
    if (!missing(first.year))
      res@first.year <- first.year
    if (!missing(first.year.model))
      res@first.year.model <- first.year.model
    if (!missing(last.year))
      res@last.year <- last.year
    if (!missing(last.year.model))
      res@last.year.model <- as.integer(last.year.model)
    if (!missing(last.season))
      res@last.season <- as.integer(last.season)
    if (!missing(last.season.last.year))
      res@last.season.last.year <- as.integer(last.season.last.year)
    if (!missing(no.species))
      res@no.species <- as.integer(no.species)
    if (!missing(first.age))
      res@first.age <- as.integer(first.age)
    if (!missing(rec.season))
      res@rec.season <- as.integer(rec.season)
    if (!missing(max.age.all))
      res@max.age.all <- as.integer(max.age.all)
    if (!missing(species.info))
      res@species.info <- as.matrix(species.info)
    
    # Verify that this object is valid
    test <- validObject(res)
    
    if (!test) stop("Invalid object:", test)
  }
  
  return(res)
}


a<-RSMS.control(first.year=1976,last.year=2000,no.species=3,
                 no.other.predators=1,no.VPA.predators=1,species.names=c('Bird','Cod','Herring'))
