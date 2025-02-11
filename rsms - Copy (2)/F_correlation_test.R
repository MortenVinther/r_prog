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

n=5
s=1
rho=0.6
fsd<-c(0.3,0.4,0.5,0.1,0.7)
fsd^2

fcor <- outer(1:n,
              1:n,
              function(i,j)(i!=j)*rho[s] + (i==j))
fcor
fvar <- outer(1:n,
              1:n,
              function(i,j)fcor[cbind(i,j)]*fsd[i]*fsd[j])
fvar


ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 1:n - 1)
  rho^exponent
}

fcor
ar1_cor(n,rho)

rho^abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 1:n - 1)


fcor <- outer(1:n,
              1:n,
              function(i,j)(i!=j)*rho[s] + (i==j))


 outer(1:n,
              1:n,
              function(i,j)(i!=j) + abs(i-j))
 
 
