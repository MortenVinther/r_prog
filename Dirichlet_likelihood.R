
multibeta <- function(alphas, log = TRUE) {
  lmb <- sum(lgamma(alphas)) - lgamma(sum(alphas))
  if (log) lmb else exp(lmb)
}

obs<-c(0.01,0.09,0.10,0.15,0.65)
sum(obs)

obs_bias<-c(0.1,0.20,0.10,0.15,0.65)
obs_bias<-obs_bias/sum(obs_bias)
sum(obs_bias)

obs;round(obs_bias,2)

a0<- c(2,5,10,100)

multibeta(alphas=2*obs)
multibeta(alphas=5*obs)
multibeta(alphas=10*obs)

mul<-function(a0= c(2,5,10), obs) {
  beta<-data.frame(a0,beta=0)
  for (i in 1:length(a0)) {
    beta[i,'beta']<-multibeta(a0[i]*obs,log=TRUE)
  }
  return(beta)
}
mul(a0,obs)


productx_a<-function(a0,obs,dolog=TRUE) {
  alpha<-a0*obs -1
  y<-sum(log(obs)*alpha)
  if (!dolog) y=exp(y) 
  y
}
productx_a(a0=2,dolog=T,obs)
productx_a(a0=2,dolog=T,obs)
productx_a(a0=20,dolog=T,obs)
productx_a(a0=100,dolog=T,obs)


library(DirichletReg)

doDir<-function(a0,obs,n=100) {
  X1 <- DirichletReg::rdirichlet(n=n, alpha=a0*obs)
  head(X1)
  DirichletReg::ddirichlet(x=(X1), alpha=a0*obs,log=TRUE, sum.up = T)
}

doDir(a0=2,obs)
doDir(a0=10,obs)
doDir(a0=20,obs)
doDir(a0=50,obs)
doDir(a0=100,obs)

alpha0<-10
# + log_likelihood should all be the same
-multibeta(alphas=alpha0*obs)+productx_a(a0=alpha0,dolog=T,obs)
-mul(alpha0,obs)['beta']             +productx_a(a0=alpha0,dolog=T,obs)
ddirichlet(x=t(matrix(obs)), alpha=alpha0*obs,log=TRUE, sum.up = T)


#negative log likelihood becomes better with increasing factor on alpha0 (and the same observations)
-ddirichlet(x=t(matrix(obs)), alpha=alpha0*obs*0.1,log=TRUE, sum.up = T)
-ddirichlet(x=t(matrix(obs)), alpha=alpha0*obs*0.2,log=TRUE, sum.up = T)
-ddirichlet(x=t(matrix(obs)), alpha=alpha0*obs*0.5,log=TRUE, sum.up = T)
-ddirichlet(x=t(matrix(obs)), alpha=alpha0*obs*1.0,log=TRUE, sum.up = T)
-ddirichlet(x=t(matrix(obs)), alpha=alpha0*obs*10,log=TRUE, sum.up = T)


#### simulation

obs

n=10
bias<-1
a0<-c(3,10,100)

mul<-function(a0, obs,bias=1,n=10) {
  for (i in 1:length(a0)) {
   if (i==1) {
   X<-matrix(rep(obs,n),byrow=T,nrow=n) 
   X<-data.frame(X)
   colnames(X)<-paste0('prey',1:length(obs))
   X$alpha<-"The Truth"
   X$alpha0<- -9
   X$sample<-factor(1:n)
   X$nsamp<-n
   }  
   X1 <- DirichletReg::rdirichlet(n=n, alpha=a0[i]*obs*bias)
   X1<-data.frame(X1)
   colnames(X1)<-paste0('prey',1:length(obs))
  
   X1$alpha<-paste("alpha:",a0[i])
   X1$alpha0<- a0[i]
   X1$sample<-factor(1:n)
   X1$nsamp<-n
   
   X<-rbind(X,X1)
  }
  return(X)
}
xx<-mul(a0,obs,n=n)
head(xx)

x2<-pivot_longer(xx,cols=prey1:prey5) %>% rename(prey=name)


if (n<=10) {
ggplot(x2, aes(fill=prey, y=value, x=sample)) + 
  geom_bar(position="stack", stat="identity")+  
  facet_wrap(vars(alpha),scales = "free")+
  theme_minimal() 
}


x3<-filter(xx,alpha0>0)
by(x3,list(x3$alpha0),function(x){
  alpha0<-x[1,'alpha0']
  nsamp<-x[1,'nsamp']
 cat(alpha0,nsamp,'\n')
  x$alpha<-x$alpha0 <-  x$sample <-  x$nsamp <-NULL
  xx<-as.matrix(x)
  #print(round(xx,3))
  list(no_process_error=ddirichlet(x=xx,alpha=alpha0*obs,log=TRUE, sum.up = T)/nsamp,
       with_process_error=ddirichlet(x=xx,alpha=alpha0*obs_bias,log=TRUE, sum.up = T)/nsamp)
})


X1 <- rdirichlet(100, c(5, 5, 10))
ddirichlet(X1, 1*c(5, 5, 10) ,sum.up = TRUE,log=T)
ddirichlet(X1, 2*c(5, 5, 10) ,sum.up = TRUE,log=T)
ddirichlet(X1, 10*c(5, 5, 10) ,sum.up = TRUE,log=T)

# new alpha deviate from alpha used to create data (process error)
ddirichlet(X1, c(5, 5, 10) ,sum.up = TRUE,log=T)
ddirichlet(X1, c(4, 5, 10) ,sum.up = TRUE,log=T)
ddirichlet(X1, c(1, 5, 10) ,sum.up = TRUE,log=T)

