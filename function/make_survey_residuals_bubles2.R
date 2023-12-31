
#Function to read and plot survey residuals
# parameter start.year: first year on X-axis, default=0 (defined from data)
# parameter end.year: end year on X-axis, default=0 (defined from data)
# 
# use over.all.max to set the maximum size of the reference buble. A value of 0 scales bubbles individually  


over.all.max<-6  # use over.all.max different from 0 to set bubble size  
over.all.max<-2.5

use.ref.dot<-TRUE


plot.survey.residuals2<-function(dev,nox=1,noy=1,Portrait=T,start.year=0,end.year=0,reverse.colors=F,use.ref.dot=TRUE,add.title=TRUE,over.all.max=1.5,my.species=NA,standardize=F,oDir=data.path,my.palette=palette("default"),...) {
  palette(my.palette) 
  
fleet.names<-Read.fleet.names()

file<-file.path(data.path,'fleet_info.dat') 
finfo<-scan(file,comment.char = "#") 

file<-file.path(data.path,'catch_survey_residuals.out')
res<-read.table(file,comment.char = "#",header=T)
res<-subset(res,data=='survey')
if (standardize) res$residual<- res$stand.residual
res[res$residual==-99.9 ,'residual']<-NA
if (reverse.colors) res$residual<- -res$residual

max.buble<-max(abs(res$residual),na.rm=TRUE)

Init.function() # get SMS.control object  including sp.names

nsp<-nsp-first.VPA+1

nox.noy<-nox*noy
plot.no<-0

years<-rep(0,2)
ages<-rep(0,2)
if (is.na(my.species)) my.species<-1:nsp

for (sp in 1:nsp) if (sp %in% my.species) {
  aa<-subset(res,Species.n==first.VPA+sp-1)
  nf<-finfo[sp+1] #no. of fleets
  sp.name<-sp.names[sp+first.VPA-1]

  print(paste(sp.name,sp,"Number of fleets:",nf))

  for (f in 1:nf) {

    print(paste("sp:",sp,"  fleet:",f))

    nyr<-years[2]-years[1]+1
    nag<-ages[2]-ages[1]+1
    plot.no<-plot.no+1

    if (plot.no%%nox.noy==0 || f==1){
     newplot(dev,nox,noy,Portrait=Portrait,filename=paste('Survey',sp.name,f,plot.no),dir=oDir,...)
      par(mar=c(3,4,3,2))
      if (dev=="wmf") par(mar=c(2,4,2,2))
      plot.no<-0
    }

    bb<-subset(aa,fleet==f)
    tmp<-tapply(bb$residual,list(age=bb$Age,year=bb$Year),sum,na.rm=T)
    tmp[tmp==-99.99]<-0

    
    xpos <- as.numeric(dimnames(tmp)[[2]]) # years
    ypos <- as.numeric(dimnames(tmp)[[1]]) #ages
    title<- paste(sp.name," fleet:",f,sep="")
    title<- fleet.names[sp,f]
    #title<- paste(sp.name,": ",fleet.names[f],sep="")

    if (length(ypos)==1) {
      
        tmp2<-tmp
        tmp2[]<-0
        if (ypos[1]>=1) {
          tmp<-rbind(tmp2,tmp,tmp2)
          ypos <- (ypos-1):(ypos+1)
        }
        else if (ypos[1]==0) {
          tmp<-rbind(tmp,tmp2)
          ypos <- ypos:(ypos+1)
        }      
    }
     #print(tmp)
    if (over.all.max>0) residplot(tmp,xpos,ypos,main=title,refdot=use.ref.dot,start.year=start.year,end.year=end.year,maxn=over.all.max)
    else residplot(tmp,xpos,ypos,main=title,refdot=use.ref.dot,start.year=start.year,end.year=end.year)

  }
}

cat("Max buble size=",max.buble,'\n')
if (reverse.colors) cat("log(Survey observed CPUE) - log(expected CPUE). 'Red' negative, 'White' positive\n")
if (! reverse.colors) cat("log(Survey observed CPUE) - log(expected CPUE). 'Red' positive, 'White' negative\n")

if (dev !='screen') cleanup()
}
