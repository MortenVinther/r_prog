install.packages('DirichletReg')
library(DirichletReg)

ddirichlet_R<-function (x, alpha, log = FALSE, sum.up = FALSE) 
{
  if (is.null(dim(x))) 
    stop("x must be a matrix")
  if (is.vector(alpha)) {
    if (ncol(x) != length(alpha)) 
      stop("alpha must be a vector/matrix fitting to the data in x")
    alpha <- matrix(rep(alpha, nrow(x)), nrow(x), byrow = TRUE)
  }
  if (any(dim(alpha) != dim(x))) 
    stop("check if x and alpha are correctly specified")
  if (any(alpha <= 0)) {
    warning("all values in alpha must be > 0")
    if (sum.up) 
      return(NaN)
    else return(rep(NaN, nrow(x)))
  }
  res <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) + 
    rowSums((alpha - 1) * log(x))
  if (sum.up) {
    if (log) 
      return(sum(res))
    else return(exp(sum(res)))
  }
  else {
    if (log) 
      return(res)
    else return(exp(res))
  }
}

##################
library(MCMCpack)
MCMCpack::ddirichlet

ddirichlet_R
function (x, alpha, log = FALSE, sum.up = FALSE) 
{
  if (is.null(dim(x))) 
    stop("x must be a matrix")
  if (is.vector(alpha)) {
    if (ncol(x) != length(alpha)) 
      stop("alpha must be a vector/matrix fitting to the data in x")
    alpha <- matrix(rep(alpha, nrow(x)), nrow(x), byrow = TRUE)
  }
  if (any(dim(alpha) != dim(x))) 
    stop("check if x and alpha are correctly specified")
  if (any(alpha <= 0)) {
    warning("all values in alpha must be > 0")
    if (sum.up) 
      return(NaN)
    else return(rep(NaN, nrow(x)))
  }
  res <- lgamma(rowSums(alpha)) - rowSums(lgamma(alpha)) + 
    rowSums((alpha - 1) * log(x))
  if (sum.up) {
    if (log) 
      return(sum(res))
    else return(exp(sum(res)))
  }
  else {
    if (log) 
      return(res)
    else return(exp(res))
  }
}


##################
library(MCMCpack)

MCMCpack::ddirichlet<-function (x, alpha) 
{
  dirichlet1 <- function(x, alpha) {
    logD <- sum(lgamma(alpha)) - lgamma(sum(alpha))
    s <- sum((alpha - 1) * log(x))
    exp(sum(s) - logD)
  }
  if (!is.matrix(x)) 
    if (is.data.frame(x)) 
      x <- as.matrix(x)
  else x <- t(x)
  if (!is.matrix(alpha)) 
    alpha <- matrix(alpha, ncol = length(alpha), nrow = nrow(x), 
                    byrow = TRUE)
  if (any(dim(x) != dim(alpha))) 
    stop("Mismatch between dimensions of x and alpha in ddirichlet().\n")
  pd <- vector(length = nrow(x))
  for (i in 1:nrow(x)) pd[i] <- dirichlet1(x[i, ], alpha[i, 
  ])
  pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
  pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
  return(pd)
}
