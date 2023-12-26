library(DirichletReg)


x<-c(5,35,60)
preys<-c('Cod','Sprat','Other')
names(x)<-preys

x<-x/sum(x)

n<-100
pSum<-20
alfa<-x*pSum
xx<-rdirichlet(n,alfa)
xx
xx<-as.data.frame(xx)
names(xx)<-preys

AL <- DR_data(xx)
AL
cbind(AL,toSimplex(AL))


a1<-DR_data(t(matrix(x,dimnames=list(preys,NULL))))
a1
a11<-toSimplex(a1)
cbind(a1,a11)

plot(AL, cex = 0.5, a2d = list(colored = F, c.grid = T))
points(x=a11[,1],y=a11[,2],pch = 16,col='red',cex=2)
points(a11,pch = 16,col='blue',cex=2)

text(x=a11[,1],y=a11[,2],label='ddd',col='red',cex=2)
text(a11,label='dds',col='blue',cex=2)
text(x=1:10,y=1:10,labels=paste(1:10,1:10))


text(x=1:10,y=1:10,labels=paste(1:10,1:10))
text(toSimplex(t(matrix(x))),col='red',label='XXXX',cex=8)
text(x=0.7,y=0.06,col='red',label='XXXXXX',cex=8)

cbind(AL,toSimplex(AL))


toSimplex()
round(x,3)
round(xx,3)
apply(xx,1,sum)

ddirichlet(xx, x*pSum, log = T, sum.up = FALSE)
ddirichlet(xx, x*pSum, log = T, sum.up = T)

one<-xx[1,]
sum(one)

logSumP<-lgamma(sum(alfa))
loggamma<-sum(lgamma(one))
pPowStom<-sum((alfa-1)*log(one))

logLike<-logSumP-loggamma+pPowStom

logLike

-logLike
