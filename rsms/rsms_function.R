
func <- function(parameters) {
  square <- function(x){x*x}
  ### Parameter transform
  # f <- function(x) {2/(1 + exp(-2 * x)) - 1}
 
  getAll(data, parameters,warn=TRUE)
  
  ## Optional (enables extra RTMB features) 
  #obs<-OBS(obs) #statement tells RTMB that obs is the response. This is needed to enable automatic simulation and residual calculations from the model object.
  logCatchObs<-OBS(logCatchObs)
  logSurveyObs<-OBS(logSurveyObs)
 
  # first we have all the N states
  logN<-list()
  for (s in 1:nSpecies) logN<-c(logN, list(Un[(nlogNfromTo[s,1]:nlogNfromTo[s,2]),,drop=FALSE]))

  # and then the F states
  # with seasonal data, logF is redefined as the sum of seasonal F values (which is not the same as the annual F)
  
  logF<-list()
  for (s in 1:nSpecies) logF<-c(logF, list(Uf[(nlogFfromTo[s,1]:nlogFfromTo[s,2]),,drop=FALSE]))
  
  timeSteps <- nYears
  stateDimF <- nlogF
  stateDimN <- nlogN
  sdLogFsta <- exp(logSdLogFsta)
  varLogN <- exp(logSdLogN*2)
  varLogObsCatch <- exp(logSdLogObsCatch*2)
  maxAgePlusGroup <-info[,'+group']
  varLogObsSurvey = exp(logSdLogObsSurvey*2)
  
  makeVar2<-function() {
    d<-lapply(1:nSpecies,function(x) {
      matrix(0,nrow=info[x,'la'],ncol=timeSteps,dimnames=list(a=paste0("a",1:info[x,'la']),y=years))
    })
    names(d)<-spNames
    d
  }
  
  makeVar3<-function() {
    d<-lapply(1:nSpecies,function(x) {
      array(0,dim=c(timeSteps,nSeasons,info[x,'la']),dimnames=list(y=years,q=1:nSeasons,a=paste("a",1:info[x,'la'])))
    })
    names(d)<-spNames
    d
  }
  
  logNq<-makeVar3()
  logNbarq<-makeVar3()
  Zq<-makeVar3()
  MM<-makeVar3()  # natural mortality
  for (i in (1:length(MM))) MM[[i]]<-natMor[[i]]
  M2<-makeVar3()  # predation mortality (M2)
  
  
  #Cq<-makeVar3()
  Chat<-makeVar2()
  predN<-makeVar2()
  #outF<-makeVar2()
  ssb<-matrix(0,nrow=nSpecies,ncol=nYears)
  fq<-1L;     # first season (e.g. quarter)
  lq<-nSeasons #  use 1 for annual data
  nlls<-matrix(0,ncol=4,nrow=nSpecies,dimnames=list(species=spNames,nll=c("catch","F","N","survey")))
  
  
  # survey initialization
  predSurveyObs<-numeric(length(logSurveyObs))
  surveyType<-keySurvey.overview[,'type']
  mina<-keySurvey.overview[,'mina']
  maxa<-keySurvey.overview[,'maxa']
  techCreepIdx<-keySurvey.overview[,'techCreep']
  
  
  ## Initialize joint negative log likelihood
  ans <- 0
  
  # needed for simulation  and consistencyCheck ? does not work
  # i<-1;
  # for (s in 1:nSpecies) {faf=info[s,'faf']; for (a in nlogFfromTo[s,1]:nlogFfromTo[s,2]) ans<- ans -dnorm(logF[[s]][a,i], sdLogFsta[keyLogFstaSd[s,a+faf-1L]], log=TRUE)}
  # for (s in 1:nSpecies) {for (a in nlogNfromTo[s,1]:nlogNfromTo[s,2]) ans<- ans -dnorm(logN[[s]][a,i]-predN[[s]][a,i],  log=TRUE)}
  

  ##########################################################################################
  
  SSB_R<-function(s,y,a=1) {
    if (stockRecruitmentModelCode[s]==0 | !recruitYears[s,y]){    ## straight RW
      rec = logN[[s]][a, y-1]
    } else {
      if (stockRecruitmentModelCode[s]==1){ ## Ricker
        rec<-rec_loga[s]+log(ssb[s,y-recAge])-exp(rec_logb[s])*ssb[s,y-recAge]
      } else {
        if(stockRecruitmentModelCode[s]==2){  ## B&H
          rec<-rec_loga[s]+log(ssb[s,y-recAge])-log(1+exp(rec_logb[s])*ssb[s,y-recAge])
        } else {
          if(stockRecruitmentModelCode[s]==3){  ## GM
            rec<-rec_loga[s]
          } else {
            stop(paste0("SR model code ",stockRecruitmentModelCode[s]," not recognized"))
          }
        }
      }
    }
  }

 
  #function for Fishing mortality
  fiMo<-function(s,y,a,expIt=TRUE) {
    ff<-logF[[s]][keyLogFsta[s,a],y]  
    if (fModel==2) {
     # cat(s,y,a,fSepar,'\n')
      if (a >=fSepar)  ff<-ff+logSeparF[keyLogSeparF[[s]][y,a]] 
    }
    #outF[[s]][a,y]<<-ff
    if (expIt) return(exp(ff)) else return(ff)
  }

  
  # food suitability
  suitability <-function(pred,prey,predSize,preysize,observed) { 
   #   FUNCTION dvariable suit(int y, int q, int d, int pred,int prey,double pred_size,double prey_size, int observed)
  #  vul=season_overlap(d,pred,q,prey)*vulnera(pred_prey_comb(d,pred,prey));
 
  vul<-vulnera[vulneraIdx[prey,pred]]  
  #if (use_overlap==1 && y<=lyModel) vul=vul*overlap(d,pred,y,q,prey);
  #else if (use_overlap==2 && y>lyModel) vul=vul*overlap_forecast(d,pred,y,q,prey); 
  
  opt<-info[pred,'sizeSlct']
   
  if (opt==0) {         # (uniform) no size selection
    if (observed) return(vul) else {
      ratio<-pred_size/prey_size;  
      if (ratio>=predPreySize['lower_intc',prey,pred] & ratio<=predPreySize['upper_intc',prey,pred]) return(vul) else return(0.0)
    }
  }
  else if (opt==4) {        # (confined uniform) no size selection, but within limits
    if (observed==1) return(vul) else { 
      ratio<-pred_size-prey_size; # sizes are on log scale  
      if (ratio>=(predPreySize['lower_intc',prey,pred] + predPreySize['lower_slope',prey,pred]*pred_size) & 
          ratio<=(predPreySize['upper_intc',prey,pred] + predPreySize['upper_slope',prey,pred]*pred_size)) return(vul) else return(0.0)
    }
  }
  } #end function
   
  # else {           //  size selection
  #   if (ratio >= min_pred_prey_size_ratio(pred,prey) &&  ratio <= max_pred_prey_size_ratio(pred,prey)){    
  #     
  #     if (size_selection(pred)==1  || size_selection(pred)==2 ) { //normal distribution or asymmetric normal distribution
  #       tmp=log(ratio)-(pref_size_ratio(pred)*prey_size_adjustment(prey)+pref_size_ratio_correction(pred)*log(pred_size));     
  #       return vul*exp(-square(tmp)/(2.0* var_size_ratio(pred)));
  #     }
  #     else if (size_selection(pred)==3) {      //Gamma     1.0/(b^a*gamma(a))*x^(a-1)*exp(-x/b)
  #       // use gammln function: exp(-(log(b)*a+log(gamma(a)))+log(x)*(a-1)-x/b)
  #       tmp=log(ratio);
  #       // cout<<"gamma var:"<<var_size_ratio(pred)<<" pref: "<<pref_size_ratio(pred)<<" tmp:"<<tmp<<endl;                                                   
  #       size_sel= exp(-(log(var_size_ratio(pred))*pref_size_ratio(pred)+gammln(var_size_ratio(pred)))+log(tmp)*(pref_size_ratio(pred)-1.0)-tmp/pref_size_ratio(pred));
  #       return vul*size_sel;
  #     }
  #     else if (size_selection(pred)==5 || size_selection(pred)==6) {      //beta     gamma(a+b)/gamma(a)/gamma(b)*x^(a-1)*(1-x)^(b-1)
  #       // use gammln function: exp(lgamma(a+b)-lgamma(a)-lgamma(b)+log(x)*(a-1)+log(1-x)*(b-1))                                             
  #       //rescale ratio to [0;1]                                             
  #       ratio=(log(ratio)-all_min_pred_prey_size_ratio(pred))/all_range_pred_prey_size_ratio(pred); 
  #       size_sel=exp(gammln(pref_size_ratio(pred)+var_size_ratio(pred))-gammln(pref_size_ratio(pred))-gammln(var_size_ratio(pred))+log(ratio)*(pref_size_ratio(pred)-1.0)+log(1.0-ratio)*(var_size_ratio(pred)-1));
  #       //cout<<"pred:"<<pred<<" prey:"<<prey<<setprecision(3)<<" ratio:"<<ratio<<" size_sel:"<<size_sel<<" suit:"<<vul*size_sel<<endl;
  #       return vul*size_sel;
  #     }
  #     
  #     else return -1000.0;  // error
  #   }
  #   else return 0.0;
  # }  
  # 
  
 # Other food suitability
 otherSuit<-function(pred, pred_size,y, q) {

  # tmp=exp(stl_other_suit_slope(pred)*log(pred_size/AV_other_food_size(pred)))*season_overlap(d,pred,q,0);
  # 
  # if (use_overlap>1 && y<=lyModel) tmp*=overlap(d,pred,y,q,0); 
  # else if (use_overlap==2 && y>lyModel) tmp*=overlap_forecast(d,pred,y,q,0); 
  # 
  # return tmp;
  # 
  
  return(1.0)
 }
    
    
  ###################  now we begin 
  for (s in 1:nSpecies) {
  
    cat('Species:',s,spNames[s],'\n')
    fModel <- as.integer(info[s,'fModel']) 
    fSepar <- as.integer(info[s,'fSepar'])
    
    ## First take care of F
    fsd <- sdLogFsta[keyLogFstaSd[s,keyLogFstaSd[s,]>0]]
    if (useRho[s]) {
      fcor <- outer(1:stateDimF[s],
                    1:stateDimF[s],
                    function(i,j)(i!=j)*rho[s] + (i==j))
      fvar <- outer(1:stateDimF[s],
                    1:stateDimF[s],
                    function(i,j)fcor[cbind(i,j)]*fsd[i]*fsd[j])
    } else {
      fvar <- diag(fsd[1:stateDimF[s]]*fsd[1:stateDimF[s]],ncol=stateDimF[s], nrow=stateDimF[s])
    }
    
    if (zeroCatchYearExistsSp[s]==1)  ans <- ans - sum(dmvnorm( diff( t(logF[[s]][,-zeroCatchYear[[s]]])),mu=0, Sigma=fvar, log=TRUE)) else {
      ans <- ans - sum(dmvnorm( diff( t(logF[[s]])),mu=0, Sigma=fvar, log=TRUE))
    }
    
    if (Debug) nlls[s,"F"]<-  -sum(dmvnorm( diff( t(logF[[s]])) , 0, Sigma=fvar, log=TRUE))
    
   # Spawning Stock Biomass

  # # SSB 1. January !!   
  # for (s in 1:nSpecies) { 
  #   for (y in (1:nYears)) {
  #     ssb[s,y]<-sum(exp(logN[[s]][,y]) *propMat[[s]][y,q,]*stockMeanWeight[[s]][y,q,])  #HUSK AT SÆTTE PropMat for ikke rekruterede aldre til 0
  #   }
  # }  
  
  # SSB 1. January !!
  #  Class 'simref' does not allow operation 'sum' (argument must be fully simulated)
    for (y in 1:nYears) {
      for (a in 1:stateDimN[s]) {
        ssb[s,y]<- ssb[s,y] + exp(logN[[s]][a,y])*propMat[[s]][y,spawnSeason,a]*stockMeanWeight[[s]][y,spawnSeason,a]
      }
    }

    
    ## Now take care of N
    nvar <- outer(1:stateDimN[s], 1:stateDimN[s],
                  function(i,j) (i==j)*varLogN[ keyVarLogN[s,i]])
   
    #recruits first year given recruitment at age 0 (later in the same year)
    if (stockRecruitmentModelCode[s] >=1 & recAge==0 & recruitYears[s,1]){
      predN[[s]][1,1]<-SSB_R(s,y=1);
      ans <- ans - dmvnorm(logN[[s]][1,1], predN[[s]][1,1], nvar[1,1], log=TRUE) ## N-Process likelihood
    }
    
   for(i in 2:timeSteps) { 
      predN[[s]][1,i]<-SSB_R(s,y=i);
  
      for(j in 2:stateDimN[s]) {
          if (keyLogFsta[s,j-1]>0) {  # to take account for ages with no F
             predN[[s]][j,i]=logN[[s]][j-1,i-1]-fiMo(s,i-1,j-1)-sum(MM[[s]][i-1,,j-1]) 
          } else { 
            predN[[s]][j,i]=logN[[s]][j-1,i-1]-sum(MM[[s]][i-1,,j-1]) 
          }
       }
       if(maxAgePlusGroup[s]==1){
          predN[[s]][stateDimN[s],i] = log( exp(logN[[s]][stateDimN[s]-1,i-1]-fiMo(s,i-1,stateDimN[s]-1)-sum(MM[[s]][i-1,,stateDimN[s]-1])) +
                                       exp(logN[[s]][stateDimN[s],i-1]  - fiMo(s,i-1,stateDimN[s])  -sum(MM[[s]][i-1,,stateDimN[s]]))  )
       }
      ans <- ans - dmvnorm(logN[[s]][,i], predN[[s]][,i], nvar, log=TRUE) ## N-Process likelihood
      if (Debug) nlls[s,"N"]<- nlls[s,"N"]  - dmvnorm(logN[[s]][,i], predN[[s]][,i], nvar, log=TRUE) 
    }

  #We have now calculated logN (1 quarter)
 
    for(i in 1:timeSteps) {
      for(j in 1:stateDimN[s]) {
        
        #first season (q=fq)
        if (j >1 | recSeason==1 ) {  # first age might not have enter the model yet
          logNq[[s]][i,fq,j]<-logN[[s]][j,i]
          if((keyLogFsta[s,j])>0) {
            Zq[[s]][i,fq,j] <- fiMo(s,i,j)*seasFprop[[s]][i,fq,j] + MM[[s]][i,fq,j]
          } else {
            Zq[[s]][i,fq,j] <- MM[[s]][i,fq,j]
          }
          logNbarq[[s]][i,fq,j] <- logNq[[s]][i,fq,j]-log(Zq[[s]][i,fq,j]) + log(1.0 -exp(-Zq[[s]][i,fq,j]))
          if (keyLogFsta[s,j]>0) Chat[[s]][j,i] <- exp(logNbarq[[s]][i,fq,j]+fiMo(s,i,j,expIt=FALSE)+log(seasFprop[[s]][i,fq,j]))
        } else Chat[[s]][j,i] = 0.1; #SNYD for at undgå log(0), værdien bruges ikke, så OK
        
        # recruits within the year, need simpler code
        if (j ==1 & recSeason>1) {
           logNq[[s]][i,recSeason,j]<-logN[[s]][j,i] 
           if (keyLogFsta[s,j]>0) {
             Zq[[s]][i,recSeason,j] <- fiMo(s,i,j)*seasFprop[[s]][i,recSeason,j] + MM[[s]][i,recSeason,j]
           } else {
             Zq[[s]][i,recSeason,j] <- MM[[s]][i,recSeason,j]
           }
           logNbarq[[s]][i,recSeason,j] <- logNq[[s]][i,recSeason,j]-log(Zq[[s]][i,recSeason,j]) +log(1.0 -exp(-Zq[[s]][i,recSeason,j]))
           if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <-  exp(logNbarq[[s]][i,recSeason,j]+fiMo(s,i,j,expIt=FALSE)+log(seasFprop[[s]][i,recSeason,j]))
           } else { Chat[[s]][j,i]<-0.1  }
        }
        
        #remaining seasons
        if (nSeasons>1) for (q in 2:lq) if (j >1 | q>recSeason ) {
          logNq[[s]][i,q,j]<- logNq[[s]][i,q-1,j]-Zq[[s]][i,q-1,j]
          if (keyLogFsta[s,j]>0) {
            Zq[[s]][i,q,j] <- fiMo(s,i,j)*seasFprop[[s]][i,q,j] + MM[[s]][i,q,j]
          } else {
            Zq[[s]][i,q,j] <- MM[[s]][i,q,j]
          }
          logNbarq[[s]][i,q,j] <- logNq[[s]][i,q,j]-log(Zq[[s]][i,q,j]) +log(1.0 -exp(-Zq[[s]][i,q,j]))
          if (keyLogFsta[s,j]>0)  {Chat[[s]][j,i] <- Chat[[s]][j,i]+ exp(logNbarq[[s]][i,q,j]+fiMo(s,i,j,expIt=FALSE)+log(seasFprop[[s]][i,q,j]))}  else { Chat[[s]][j,i]<-0.1  }
        }
      }
    }

  
  Chat[[s]]<-log(Chat[[s]])
  
  
  ##  match to observations
  
  # first catches
    #cat('Species: ',s,'\n')
    idx<-keyCatch[,"s"]==s
    key.v<-keyCatch[idx,"keyVarObsCatch"]
    yy<-keyCatch[idx,"y"]
    aa<-keyCatch[idx,'a']
    ay<-cbind(aa,yy)
    obs<-logCatchObs[idx]
    predObs<-Chat[[s]][ay]
    for (i in 1:length(predObs)) {
      var <- varLogObsCatch[key.v[i]]
      ans <- ans - dnorm(obs[i],predObs[i],sqrt(var),log=TRUE)
      
      if (Debug==1)  nlls[s,'catch']<- nlls[s,'catch']   - dnorm(obs[i],predObs[i],sqrt(var),log=TRUE)
    }
    # if (Debug==1)  {  # the same and more efficient than above, but problems with convergence ??
    #    var <- varLogObsCatch[key.v]
    #    nlls[s,'catch']<- nlls[s,'catch']   - sum(dnorm(obs,predObs,sqrt(var),log=TRUE))
    #  }

  # and now surveys
  fleets<-keySurvey.overview[keySurvey.overview[,'s']==s,'f']
  
  for (fl in  fleets)  {
    cat("survey: ",fl,'\n')
    keys<-keySurvey[keySurvey[,"f"]==fl ,]
    q<-keys[1,"q"]
    
    if (surveyType[fl]==1) {
      for (a in mina[fl]:maxa[fl]) {
         keysA<-keys[keys[,"a"]==a  ,]
         flYears<-keysA[,'y']
         keyPowerQ<-keysA[1,"keyPowerQ"]
         keyCatchability<-keysA[1,"keyCatchability"]
         keyVarObsSurvey<-keysA[1,"keyVarObsSurvey"]
         obs.no<-keysA[,'obs.no']
         predSurveyObs[obs.no]<- logNq[[s]][flYears,q,a] - Zq[[s]][flYears,q,a]*sampleTimeWithinSurvey[fl]
         if(keyPowerQ>0) predSurveyObs[obs.no] <- predSurveyObs[obs.no]*exp(logQpow[keyPowerQ])
         predSurveyObs[obs.no]<- predSurveyObs[obs.no]+logCatchability[keyCatchability]
         var <- varLogObsSurvey[keyVarObsSurvey]
         ans <- ans - sum(dnorm(logSurveyObs[obs.no], predSurveyObs[obs.no],sqrt(var),log=TRUE))
         if (Debug==1)  {  
          nlls[s,'survey']<- nlls[s,'survey']   - sum(dnorm(logSurveyObs[obs.no], predSurveyObs[obs.no],sqrt(var),log=TRUE))
        }
       } 
    } else if (surveyType[fl]==2) {  # exploitable biomass (assumed all age with F>0)
      flYears<-keys[,'y']
      faf<-info[s,'faf']; laf<-info[s,'last-age']
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      # does not work     predObs<-sapply(flYears,function(y) log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] )
      for (f in 1:length(obs.no))  {
        y<-flYears[f]
        predSurveyObs[obs.no[f]]<- log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] 
      }
      var <- varLogObsSurvey[keyVarObsSurvey]
      ans <- ans - sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      if (Debug==1)  {  
        nlls[s,'survey']<- nlls[s,'survey']   - sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      }
    } else if (surveyType[fl]==3) {  # SSB index
      flYears<-keys[,'y']
      faf<-1L; laf<-info[s,'la']
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      for (f in 1:length(obs.no))  {
        y<-flYears[f]
        predSurveyObs[obs.no[f]]<- log(sum(exp(logNq[[s]][y,q,faf:laf] - Zq[[s]][y,q,faf:laf]*sampleTimeWithinSurvey[fl])*stockMeanWeight[[s]][y,q,faf:laf]*propMat[[s]][y,q,faf:laf]))+ logCatchability[keyCatchability] 
      }
      var <- varLogObsSurvey[keyVarObsSurvey]
      ans <- ans - sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      if (Debug==1)  {  
        nlls[s,'survey']<- nlls[s,'survey']   - sum(dnorm(logSurveyObs[obs.no],predSurveyObs,sqrt(var),log=TRUE))
      } 
    } else if (surveyType[fl]==4) {  # effort index, one "catchability" by age group
      for (a in mina[fl]:maxa[fl]) {
        keysA<-keys[keys[,"a"]==a, ]
        flYears<-keysA[,'y']
        keyPowerQ<-keysA[1,"keyPowerQ"]
        keyCatchability<-keysA[1,"keyCatchability"]
        keyVarObsSurvey<-keysA[1,"keyVarObsSurvey"]
        techCreepNo<-techCreepIdx[fl]
        obs.no<-keysA[,'obs.no']
        predSurveyObs[obs.no]<- fiMo(s,flYearss,a,expIt=FALSE)+log(seasFprop[[s]][flYears,q,a]) + logCatchability[keyCatchability] 
        if (techCreepNo>0) for (n in 1:length(obs.no))  {
            predSurveyObs[obs.no[n]]<- predSurveyObs[obs.no[n]]+ n*logTechCreep[techCreepNo] 
        }
        var <- varLogObsSurvey[keyVarObsSurvey]
        ans <- ans - sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
        if (Debug==1) nlls[s,'survey']<- nlls[s,'survey']  - sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      }
    } else if (surveyType[fl]==5) {  # effort used as index for Fbar
      flYears<-keys[,'y']
      faf<-fbarRange[s,1]; laf<-fbarRange[s,2]; naf<-laf-faf+1
      obs.no<-keys[,'obs.no']
      keyCatchability<-keys[1,"keyCatchability"]
      keyVarObsSurvey<-keys[1,"keyVarObsSurvey"]
      techCreepNo<-techCreepIdx[fl]
      for (n in 1:length(obs.no))  {
        y<-flYears[n]
        predSurveyObs[obs.no[n]]<-0.0;
        for(j in faf:laf) {
          # cat("n:",n," y:",y," q:",q," j:",j," F",logF[[s]][keyLogFsta[s,j],y]); cat( "seasFprop[[s]][y,q,j]",seasFprop[[s]][y,q,j],'\n')
          predSurveyObs[obs.no[n]]<- predSurveyObs[obs.no[n]]+ fiMo(s,y,j) * seasFprop[[s]][y,q,j] #sum F within Fbar range
        }
        predSurveyObs[obs.no[n]]<- log(predSurveyObs[obs.no[n]]/naf) + logCatchability[keyCatchability] 
        if (techCreepNo>0) predSurveyObs[obs.no[n]]<- predSurveyObs[obs.no[n]]+ n*logTechCreep[techCreepNo] 
      }
      var <- varLogObsSurvey[keyVarObsSurvey]
      ans <- ans - sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      if (Debug==1)  {  
        nlls[s,'survey']<- nlls[s,'survey']  - sum(dnorm(logSurveyObs[obs.no],predSurveyObs[obs.no],sqrt(var),log=TRUE))
      } 
    }  else stop(paste("s:",s,"  fleet:",fl,'  fleet type:',surveyType[fl],' is not known\n'))
    
  } # end fleet loop
  } #end species loop

  if (sms.mode>0) {
    MM<-makeVar3()  # natural mortality
    for (i in (1:length(MM))) MM[[i]]<-natMor[[i]]
    
    
  }
  ADREPORT(ssb)
  REPORT(logNq)
  #REPORT(logF)
  #REPORT(outF)
  REPORT(predN)
  REPORT(Zq)
  REPORT(Chat)
  REPORT(predSurveyObs)
  REPORT(nlls)
  
  ans
}

