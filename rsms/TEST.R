abc<-function(a,b,c) {
 
 # if (c>3) c<-78
  res<-a+b+c
  return(res)
}
n=10
x<-data.frame(no=1:10,a=1:n,b=rep(8,n),c=1:n+100)
abc(a=x$a,b=x$b,c=x$c)


for (i in 1:n) print(abc(a=x[i,'a'],b=x[i,'b'],c=x[i,'c']))
lapply(1:n,function(i) abc(a=x[i,'a'],b=x[i,'b'],c=x[i,'c']))

x<-x %>% tibble()  %>% nest_by(no)
x

y<-x %>% mutate(model=list(abc(a=a,b=b,c=c)))

y
str(y)
