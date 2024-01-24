nl<-names(obj$par)

lower <- obj$par*0-Inf
upper <- obj$par*0+Inf
lower[nl=='rho'] <- 0.01
upper[nl=='rho'] <- 0.99

grep('rho',nl)

lower[nl=="logSdLogObsSurvey"]<-rep(log(0.1),length(parameters$logSdLogObsSurvey))
upper[nl=="logSdLogObsSurvey"]<-rep(log(2.0),length(parameters$logSdLogObsSurvey))

lower[nl=="logSdLogObsCatch"]<-rep(log(0.15),length(parameters$logSdLogObsCatch))
upper[nl=="logSdLogObsCatch"]<-rep(log(2.0),length(parameters$logSdLogObsCatch))


# N.sandeel
#upper[nl=="logSdLogN"]<-log(c(20,0.1))


# mac
#upper[nl=="logSdLogN"]<-log(c(0.10))
