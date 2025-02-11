plotSurveyResiduals<-function(surv,inclFleets=1:100, inclSp=1:100,fileLabel='survey_residuals',
                              outFormat=c('screen','pdf','png')[1],longSpNames=TRUE,standardized=FALSE) {

  survResid<-surv$rep$resSurv  %>% 
    transmute(s,species,year,Age=factor(age),quarter=q,f,fleet=fleetName,predict=logPred,obserded=logObs,Residual=resid,variance)
    if (standardized) survResid$Residual<- survResid$Residual/sqrt(survResid$variance)
    survResid<-survResid %>% mutate(cols=if_else(Residual<0,'negativ','positiv')) %>% filter(f %in% inclFleets & s %in% inclSp)
  
  plt<-by(survResid,survResid$s,function(x){
    if (longSpNames) sp<-surv$data$allSpNamesLong[unlist(x[1,'s'])] else sp<<-x[1,'species']
    pl<-x %>% ggplot(aes(year, Age,color= cols,fill=cols,size = sqrt(abs(Residual)))) +
      geom_point(shape = 21,alpha=0.75) +
      scale_color_manual(values = c("blue", "red")) +
      labs(x = "", y = "Age",title=paste('Survey residuals, ',sp)) +
      # facet_grid(rows =vars(fleet), scales="free_y")
      facet_wrap(vars(fleet),ncol=1)+
      my_theme() + 
      theme(strip.text = element_text(size = 10,margin = margin(0.01,0,0.01,0, "cm")))
    if (outFormat=='screen'){
      print(pl)
    } else {
      ggexport(pl, filename = paste0(fileLabel,'_',sp,'.',outFormat),width = 900,height = 600, pointsize = 16)
      cleanup()
    }
  })
  invisible(NULL)
}

plotSurveyResiduals(surv=sms,inclFleets=1:300, inclSp=1:100, outFormat=c('screen','pdf','png')[1],standardized=TRUE) 

