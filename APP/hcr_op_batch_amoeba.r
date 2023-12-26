
# make_data<-TRUE
# my.last.year<-2059

if (make_data) {
############################


batch<-T # with batch==T several runs are made with condensed output, while batch==F makes one run with detailed output

##################################################################################

# Stochastic recruitment.
# 0=No,  S-R parameters from OP_config.dat
# 1=Yes  S-R parameters and variance from OP_config.dat
# 2=Yes, S-R parameters from OP_config.dat and variance-covariance matrix from file covariance_Rec.in
# 3=Yes, S-R parameters from OP_config.dat and residuals from OP_SSB_rec_residuals.in

stochastic.recruitment<-0

HCR<-read.csv(file.path(data.path,'HCR_ini.csv'))
fsq<-subset(HCR,selectt=c(Species,Species.n,AMOEBA_F))

a<-tapply(fsq$AMOEBA_F,list(fsq$Species.n),sum)
species.R<-gsub(" ","",fsq$Species)
names(a)<-species.R

a<-data.frame(t(a))
a<-data.frame(a,change.name='all',change=1.00,change.n=0,stringsAsFactors = FALSE)

nsp<-dim(a)[[2]]-3


b<-a
a

do_change<-function(change=1,baseLine) {
  for (i in (1:nsp)) {
     b2<-baseLine
     b2[1,i]<-b2[1,i]*change; b2$change.name<-sp.names[first.VPA+i-1]; b2$change<-change; b2$change.n=first.VPA+i-1
     if (i==1)  b<-b2 else b<-rbind(b,b2)  
  }
  return(b)
}

if (batch) {
  change=1.1
  bb<-rbind(a,do_change(change=change,baseLine=a))
  b<-a; b[1,1:nsp]<-b[1,1:nsp]*change; b$change<-change
  bb<-rbind(bb,b)

  change=1.25
  bb<-rbind(bb,do_change(change=change,baseLine=a))
  b<-a; b[1,1:nsp]<-b[1,1:nsp]*change; b$change<-change
  bb<-rbind(bb,b)
  
  
  change=0.9
  bb<-rbind(bb,do_change(change=change,baseLine=a))
  b<-a; b[1,1:nsp]<-b[1,1:nsp]*change; b$change<-change
  bb<-rbind(bb,b)
  
  change=0.75
  bb<-rbind(bb,do_change(change=change,baseLine=a))
  b<-a; b[1,1:nsp]<-b[1,1:nsp]*change; b$change<-change
  bb<-rbind(bb,b)
  
    
  if (FALSE) {
    change=1.5
    bb<-rbind(bb,do_change(change=change,baseLine=a))
    b<-a; b[1,1:nsp]<-b[1,1:nsp]*change; b$change<-change
    bb<-rbind(bb,b)
    
    change=0.5
    bb<-rbind(bb,do_change(change=change,baseLine=a))
    b<-a; b[1,1:nsp]<-b[1,1:nsp]*change; b$change<-change
    bb<-rbind(bb,b)
    
  }
  
  bb$run<-1:dim(bb)[[1]]
  #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #bout<-subset(bb,select=c(COD,WHG,HAD,POK,HER,NSA,SSA,NOP,SPR,PLE,SOL))
  bout<-subset(bb,select=species.R)
} else bout<-subset(b,select=species.R)

cat(dim(bout)[[1]],'\n',file='op_multargetf.in')
write.table(bout,file='op_multargetf.in',col.names=F,row.names=F,append=T)


if (OS=='windows') op.command<-'op.exe' else if (OS=='unix') op.command<-'op'
if (OS=='windows') HPC<-F

system(paste(file.path(data.path,op.command), " -nox -maxfn 0 -nohess"))

if (batch) {
   condensed<-read.table('op_average_val.out',header=F)
  allnames<-c('value','yield', 'CWsum', 'Fcomb','Fbar', 'SSB', 'TSB', 'recruit','belowBlim','belowBpa', 'Species.n', 'run','iteration',species.R)
  dimnames(condensed)[[2]]<-allnames
  # intersect(names(bb), names(condensed))
  
  a<-merge(subset(bb,select=c(change.name, change.n,change, run)),condensed)
  #a<-merge(bb,condensed)
  dim(bb);dim(condensed);dim(a)
   subset(a,run %in% c(2,3))
  
  a<-a[order(a$run,a$Species.n),] 
  a$Species<-sp.names[a$Species.n]
  a$value<-a$belowBlim<-a$belowBpa<-a$iteration<-a$Fcomb<-a$TSB<-NULL
  head(a,20)
  
  print(xyplot(yield~Fbar| Species,data=a,scales = "free"))
  print(xyplot(SSB~Fbar| Species,data=a,scales = "free"))
  
  aa<-rbind(tapply(a$SSB,list(paste(a$Species.n,a$Species)),min),
            tapply(a$SSB,list(paste(a$Species.n,a$Species)),max),
            tapply(a$SSB,list(paste(a$Species.n,a$Species)),mean))
  rownames(aa)<-c('min','max','mean')
  round(t(aa))
  a<-subset(a,select=c("run",	"change.name","change.n",	"change",	"yield",	"CWsum",	"Fbar",	"SSB",	"recruit",	"Species.n",		"Species"))
  write.csv(a,file=file.path(data.path,paste('Amoebae_',my.last.year,'.csv',sep='')),row.names=F)
}
} #make_data

### data are now made and stored in csv file
# start using them
a<-read.csv(file=file.path(data.path,paste('Amoebae_',my.last.year,'.csv',sep='')))
a$Species<-paste(formatC(a$Species.n,width=2),a$Species,sep='_')
a$run<-paste(formatC(a$run,width=2,flag="0"),formatC(a$change*100,width=3),a$change.name,sep='_')

n.species<-length(unique(a$Species))
n.fleet<-n.species



# basis change factor used for approximations and predictions
my.change<-1.1

do_amoeba_data<-function (x,my.change=1.1,scenario_names){
  b<-subset(x,(change==1.0 | change==my.change) & !(change != 1.0 & change.n==0))
  
  SSB<-tapply(b$SSB,list(b$run,b$Species),sum)
  rownames(SSB)<-scenario_names
  
  recruit<-tapply(b$recruit,list(b$run,b$Species),sum)
  rownames(recruit)<-scenario_names
 
  scenario_no<-length(scenario_names)
  # make change matrix 
  emat<-matrix(nrow=scenario_no,ncol=scenario_no)
  emat[,]<-1
  diag(emat)<-my.change
  emat[1,1]<-1
 
  # einv is its inverse calculated as
  xinv <- solve(emat)

  ssbparm <- xinv%*%SSB
  rownames(ssbparm)<-scenario_names
  
  recparm <- xinv%*%recruit
  rownames(recparm)<-scenario_names
  
  yield<-tapply(b$yield,list(b$run,b$Species),sum)
  rownames(yield)<-scenario_names
  
  
 
  rr<-sort(unique(b$run))
  ss<-sort(unique(b$Species))
  ff<-paste(ss,'Fleet',sep='_')
  ff<-expand.grid(run=rr,Species=ss,Fleet=ff)
  ff<-merge(x=b,y=ff,all.y=TRUE)
  
  ff$del<-sapply(1:dim(ff)[[1]],function(x) grepl(ff[x,"Species"],ff[x,"Fleet"]))
  ff[ !ff$del,'yield']<-0
  yy<-tapply(ff$yield,list(ff$Fleet,ff$run,ff$Species),sum)

 # yy[,1,]  # status quo
#  yy[,2,]  # 10% increase in cod fleet  
 # yy[,2,] -yy[,1,] 

  #estimates yield parameters for each fleet, 
  yieldparm<-lapply(1:dim(yy)[[1]],function(fltn){
    cpue<-yy[fltn,,]
    cpue[1 + fltn, ] <- cpue[1 + fltn,]/my.change
    params <- xinv %*% cpue
    rownames(params) <-rownames(cpue)
    return(params)
  })
  
  list(SSB=SSB,yield=yield,recruit=recruit,xmat=emat,xinv=xinv,ssbparm=ssbparm,recparm=recparm,yieldparm=yieldparm,yieldYY=yy,yieldff=ff)  
}




parms<-do_amoeba_data(x=a,my.change=my.change,scenario_names=c('Fsq,',paste('F',HCR$Species)))

#output data for AMOEBA
ssb<-parms$SSB
colnames(ssb)<-HCR$Species
write.table(parms$SSB,file="ssb110_2020.csv",sep=',',row.names = FALSE)

fn<-rownames(parms$SSB)
fl<-data.frame(order=0:(length(fn)-1),fleets=fn,fleetNames=fn)
write.table(fl,file="fleetNames_2020.csv",sep=',',row.names = FALSE)

yy<-parms$yieldYY
ysq<-yy[,1,]  # status quo
colnames(ysq)<-HCR$Species
write.table(ysq,file="statusq_2020.csv",sep=',',row.names = FALSE)

ff<-parms$yieldff
ff<-data.frame(id=ff$run,fleetName=ff$Fleet,species=ff$Species,	Yield=ff$yield)
ff<-ff[order(ff$id,ff$fleetName,ff$species),]
write.table(ff,file="yield110.csv",sep=',',row.names = FALSE)


## end output


xchange.no<-dim(parms[["SSB"]])[[2]]

#xchange is vector of x  multipliers
predict_ssb<-function(parms,xchange) {
   xnew<-c(1,xchange) # first value is status quo
   return(t(parms[["ssbparm"]]) %*% xnew)
}

#xchange is vector of x  multipliers
predict_yield<-function(parms,xchange) {
  xnew<-c(1,xchange) # first value is status quo
  yield_params<-parms[["yieldparm"]]
  #catch by all fleets   
  yield <- matrix(nrow = n.species, ncol = n.fleet)

  for (j in 1:n.fleet) {
    cpuehat <- t(parms[["yieldparm"]][[j]]) %*% xnew
    yhat <- cpuehat * xnew[j+1]
    yield[, j ] <- yhat
  }
  round(cbind(yield,all=rowSums(yield,na.rm=TRUE)))
}



# check, no change in x 

xchange<-rep(1,xchange.no) 
ssb<-predict_ssb(parms,xchange)
#should be the same as SSB status quo 
ssb-parms[['SSB']][1,]


yield<-predict_yield(parms,xchange)
#should be the same as yield status quo 
tail(t(yield),1) -parms[['yield']][1,]

############ make plots of "real" and approximated values, using different change factors
##  SSB

b<-subset(a,change.name=='all')
b$type<-"model"
b<-subset(b,select=c(change,Species,SSB,type))
babs<-subset(b,select= -type); babs$absSSB<-babs$SSB; babs$SSB<-NULL

# SMS estimates
ggplot(b, aes(x=change, y=SSB,group=type)) +
  facet_wrap(~ Species, ncol=3,scales='free_y')+
  geom_point(size=2)

all_changes<-sort(unique(a$change))

##


for (ch in all_changes ) if (ch !=1) {
  cat(ch,'\n')
  parms<-do_amoeba_data(x=a,my.change=ch,scenario_names=c('Fsq,',paste('F',HCR$Species)))
  for (x in all_changes) {
    xchange<-rep(x,xchange.no) 
    SSB<-predict_ssb(parms,xchange)
    bb<-data.frame(SSB)
    bb$Species<- row.names(ssb)
    bb$change<-x
    bb$type<-paste('est',ch)
    b<<-rbind(b,bb)
  }
}

b$SSB<-b$SSB/1000
# SMS estimates
ggplot(b, aes(x=change, y=SSB,group=type,shape=type, color=type, size=type)) +
  facet_wrap(~ Species, ncol=3,scales='free_y')+
  geom_point(size=2)+
labs(x="factor used for model approximation", y = "SSB (000' tonnes)")

b<-merge(b,babs)
b$SSBrel<-b$SSB/b$absSSB*1000
X11()
ggplot(b, aes(x=change, y=SSBrel,group=type,shape=type, color=type, size=type)) +
  facet_wrap(~ Species, ncol=3,scales='free_y')+
  geom_point(size=2)+
  geom_hline(yintercept = 1)+
  labs(x="Relative F", y = "Approximated SSB relative to SMS SSB")
ggsave(filename = "Compar_ssb.png")



#### and now yield

b<-subset(a,change.name=='all')
b$type<-"model"
b<-subset(b,select=c(change,Species,yield,type))
babs<-subset(b,select= -type); babs$absyield<-babs$yield; babs$yield<-NULL

# SMS estimates
ggplot(b, aes(x=change, y=yield/1000,group=type)) +
  facet_wrap(~ Species, ncol=3,scales='free_y')+
  geom_point(size=2)

all_changes<-sort(unique(a$change))


for (ch in all_changes ) if (ch !=1) {
  cat(ch,'\n')
  parms<-do_amoeba_data(x=a,my.change=ch,scenario_names=c('Fsq,',paste('F',HCR$Species)))
  for (x in all_changes) {
    xchange<-rep(x,xchange.no)  
    yield<-predict_yield(parms,xchange)
    yield<-t(t(yield[,dim(yield)[[2]]]))
    bb<-data.frame(yield)
    bb$Species<-colnames(parms[['yield']])
    bb$change<-x
    bb$type<-paste('est',ch)
    b<<-rbind(b,bb)
  }
}

b$yield<-b$yield/1000
# SMS estimates
ggplot(b, aes(x=change, y=yield,group=type,shape=type, color=type, size=type)) +
  facet_wrap(~ Species, ncol=3,scales='free_y')+
  geom_point(size=2)+
  labs(x="factor used for model approximation", y = "Yield (000' tonnes)")

b<-merge(b,babs)
b$yieldrel<-b$yield/b$absyield*1000

ggplot(b, aes(x=change, y=yieldrel,group=type,shape=type, color=type, size=type)) +
  facet_wrap(~ Species, ncol=3,scales='free_y')+
  geom_point(size=2)+
  geom_hline(yintercept = 1)+
  labs(x="Relative F", y = "Approximated yield relative to SMS yield")

ggsave(filename = paste("Compar_yield.png"))


if (FALSE) {
  eff_chg<-1
  effnew <- c(1, rep(eff_chg,n.fleet)) # effnew is vector of effort multipliers, one per fleet
  ssb <- t(ssbparm) %*% effnew
  ssb- ssbtab[1,] #Fsq should be the same as first row of ssbparam
  
  effnew[n.fleet+1]<-1.1
  ssb <- t(ssbparm) %*% effnew
  ssb- ssbtab[n.fleet+1,] # should be the same as last row of ssbparam
}

  