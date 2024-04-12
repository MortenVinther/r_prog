
#fleet with commercial effort (assumed to cover the full international fishery)
efl<-data$keySurvey.overview[data$keySurvey.overview[,'type']==4,]

for (i in (1:dim(efl)[[1]])) {
  x<-filter(data.frame(data$keySurvey),f==efl[i,'f'] & s==efl[i,'s']) %>% select(obs.no,f,s,y,q) a q keyVarObsSurvey keyCatchability keyPowerQ
  print(x)
}
data$keySurvey
