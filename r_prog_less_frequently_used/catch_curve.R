library(stockassessment)
setwd(file.path("c:","_c_drev","sms-git","XSA",'input'))


catch.no<-read.ices('canum.txt')


# catch curves
a<-catch.no
a<-data.frame(v=as.vector(a),year=as.numeric(dimnames(a)[[1]]),age=rep(as.numeric(dimnames(a)[[2]]),each=dim(a)[[1]]))
subset(a,year==2015)
a$cohort<-a$year-a$age  

library(lattice) 
trellis.device(device='png',filename=file.path('catch_curves.png'),width=1000,height=600)


ttl <- list(label="", cex=1.5)
yttl <- list(label="log(catch numbers)", cex=1.0)
xttl <- list(cex=1.0)
stripttl <- list(cex=1.0)
ax <- list(cex=0.7)
b<--0.6    # slope

print(xyplot(log(v)~age| factor(cohort), data=a, subset=(cohort>1979 & cohort<2015),
             main=ttl, ylab=yttl, xlab=xttl, scales=ax, par.strip.text=stripttl, layout = c(7, 5),
             panel = function(x, y) {
               panel.xyplot(x, y,col="blue", type="l",lwd=2)
               panel.abline(a =10, b=b, lwd=0.7, lty=3, col='grey')  # der må være en nemmere måde
               panel.abline(a =12, b=b, lwd=0.7, lty=3, col='grey')
               panel.abline(a =14, b=b, lwd=0.7, lty=3, col='grey')
               panel.abline(a =16, b=b, lwd=0.7, lty=3, col='grey')
               panel.abline(a =18, b=b, lwd=0.7, lty=3, col='grey')
               panel.abline(a =20, b=b, lwd=0.7, lty=3, col='grey')
               panel.abline(a =22, b=b, lwd=0.7, lty=3, col='grey')
               panel.abline(a =24, b=b, lwd=0.7, lty=3, col='grey')
             }
))

 cleanup()




