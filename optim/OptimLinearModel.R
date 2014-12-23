
Y <- read.table("http://www.biostat.umn.edu/~sudiptob/pubh5485/Y.txt", header=F)
Y <- Y[,1]

X <- as.matrix(read.table("http://www.biostat.umn.edu/~sudiptob/pubh5485/X.txt", header=F))

##Only to compare with standard linear model:

simple.lm <- lm(Y ~ X[,2]+X[,3]+X[,4])
simple.lm.summary <- summary(simple.lm)

no.parameters = ncol(X)+1

## Different ways of writing the negative OLS likelihood:

ols.lf1 <- function(theta, y, X) {
  beta <- theta[-1]
  sigma2 <- theta[1]
  if (sigma2 <= 0) return(NA)
  n <- nrow(X)
  e <- y - X%*%beta                                  # t() = matrix transpose
  logl <- ((-n/2)*log(2*pi)) - ((n/2)*log(sigma2)) - ((t(e)%*%e)/(2*sigma2))
  return(-logl) # since optim() does minimisation by default.
}

## Pretty code, by Douglas Bates --
ols.lf3 <- function(theta, y, X) {
  if (theta[1] <= 0) return(NA)
  -sum(dnorm(y, mean = X %*% theta[-1], sd = sqrt(theta[1]), log = TRUE))
}         ## dnorm(..., log=TRUE) is generally superior to log(dnorm())


## This is the gradient of ols.lf1
## It is written in the (mathematical) notation of ols.lf1()
ols.gradient <- function(theta, y, X) {
  beta <- theta[-1]
  sigma2 <- theta[1]
  e <- y - X%*%beta
  n <- nrow(X)

  g <- numeric(length(theta))
  g[1] <- (-n/(2*sigma2)) + (t(e)%*%e)/(2*sigma2*sigma2) # d logl / d sigma
  g[-1] <- (t(X) %*% e)/sigma2                           # d logl / d beta

  return(-g)
}


## Let us use some in-built optimization functions in R to see how well they compute the mle's and the var-cov matrices

##We will first illustrate the optim function

cat("\nGradient-free (constrained optimisation) --\n")
p1 <- optim(c(0.5, rep(1, times=ncol(X))), method="L-BFGS-B", fn=ols.lf1, lower=c(1e-6,rep(-Inf, times=ncol(X))), upper=rep(Inf,no.parameters), y=Y, X=X, hessian=TRUE)

inverted1 <- solve(p1$hessian)
results1 <- cbind(p1$par[-1], sqrt(diag(inverted1[-1,-1])), p1$par[-1]/sqrt(diag(inverted1[-1,-1])))
colnames(results1) <- c("Coefficient", "Std. Err.", "t")
rownames(results1) <- names(simple.lm$coefficients)
cat("MLE results for beta's --\n")
print(results1)
cat("MLE results for sigma^2 --\n")
results1.sigmasq <- cbind(p1$par[1], sqrt(diag(inverted1)[1]))
colnames(results1.sigmasq) <- c("Estimate", "Std. Err.")
row.names(results1.sigmasq) <- c("sigma^2")
print(results1.sigmasq)
cat("Compare with the OLS results --\n")
print(simple.lm.summary)


cat("\nUsing the gradient (constrained optimisation) --\n")
p2 <- optim(c(0.5, rep(1, times=ncol(X))), method="L-BFGS-B", fn=ols.lf1, gr=ols.gradient, lower=c(1e-6,rep(-Inf, times=ncol(X))), upper=rep(Inf,no.parameters), y=Y, X=X, hessian=TRUE)

inverted2 <- solve(p2$hessian)
results2 <- cbind(p2$par[-1], sqrt(diag(inverted2[-1,-1])), p2$par[-1]/sqrt(diag(inverted2[-1,-1])))
colnames(results2) <- c("Coefficient", "Std. Err.", "t")
rownames(results2) <- names(simple.lm$coefficients)
cat("MLE results --\n")
print(results2)
cat("MLE results for sigma^2 --\n")
results2.sigmasq <- cbind(p2$par[1], sqrt(diag(inverted2)[1]))
colnames(results2.sigmasq) <- c("Estimate", "Std. Err.")
row.names(results2.sigmasq) <- c("sigma^2")
print(results2.sigmasq)
cat("Compare with the OLS results --\n")
print(simple.lm.summary)
cat("Compare with the OLS results --\n")
print(simple.lm.summary)


## Compare the variance-covariance matrices from the different methods:
cat("\n Covariance matrix from the gradient-free method --\n")
print(inverted1[-1,-1])
cat("\n Covariance matrix from the method using gradients --\n")
print(inverted2[-1,-1])
cat("\n Covariance matrix from the lm function --\n")
print(vcov(simple.lm))


## We now find the posterior mode using optim on (beta, log(sigma^2))

lposterior <- function (theta, a, b)
{
  beta <- theta[1:ncol(X)]
  sigmasq <- exp(theta[length(theta)])
  
  out <- -(a+1+length(Y)/2.0) * log(sigmasq) - (1/sigmasq)*(b + 0.5*t(Y-X%*%beta)%*%(Y-X%*%beta))

  return(-out)
}

p.mode <- optim(c(rep(0,times=ncol(X)), 0.0), method="BFGS", fn=lposterior, a=0.001, b=0.001, hessian=TRUE)
posterior.mode <- p.mode$par
posterior.vcov <- chol2inv(chol(p.mode$hessian))
cat("\n Posterior modes for (beta, log(sigma^2)) --\n")
print(posterior.mode)
cat("\n Posterior variance-covariances --\n")
print(posterior.vcov)

## Sampling from the approximate normal posterior: Bayesian CLT

NITER <- 500
library(mvtnorm)
library(coda)

posterior.samples <- rmvnorm(NITER, posterior.mode, posterior.vcov)
posterior.samples[,no.parameters] <- exp(posterior.samples[,no.parameters])
posterior.samples <- as.mcmc(posterior.samples)
summary(posterior.samples)

plot(posterior.samples, trace=FALSE)


## A digression: Getting the posterior modes for (beta,sigma^2)
## First we use a delta-method approach. We already have var-cov of (beta,z),
## where z= log(sigma^2). Therefore, sigma^2 = g(z) = e^z, so g'(z)=e^z.
## Thus, variance var(g(z)) = (g'(z))Var(z)(g'(z)).
## MLE is obtained by replacing z with the mle of z.

posterior.mode.with.sigmasq <- c(p.mode$par[1:ncol(X)], exp(p.mode$par[-(1:ncol(X))]))
g.prime <-  diag(c(rep(1, times=ncol(X)), exp(p.mode$par[-(1:ncol(X))])))
posterior.vcov2 <- g.prime%*%posterior.vcov%*% t(g.prime)
cat("\n Posterior modes for (beta, sigma^2) --\n")
print(posterior.mode.with.sigmasq)
cat("\n Posterior variance-covariances for (beta, sigma^2) using delta-method --\n")
print(posterior.vcov2)

## Alternative: can we find the posterior var-cov of (beta,sigma^2) directly?

lposterior <- function (theta, a, b)
{
  beta <- theta[1:ncol(X)]
  sigmasq <- theta[length(theta)]
  
  out <- -(a+1+length(Y)/2.0) * log(sigmasq) - (1/sigmasq)*(b + 0.5*t(Y-X%*%beta)%*%(Y-X%*%beta))

  return(-out)
}

## Now we need to perform constrained maximization:
## We need to define bounds

inits <- c(rep(0, times=ncol(X)), 1.0)
lower.bd <- c(rep(-Inf, times=ncol(X)), 0.00001)
upper.bd <- rep(Inf, times=no.parameters)

p.mode <- optim(inits, method="L-BFGS-B", fn=lposterior, lower=lower.bd, upper=upper.bd, a=0.001, b=0.001, hessian=TRUE)
posterior.mode <- p.mode$par
posterior.vcov <- chol2inv(chol(p.mode$hessian))
cat("\n Posterior modes for (beta, sigma^2) using constrained optimization --\n")
print(posterior.mode)
cat("\n Posterior variance-covariances for (beta, sigma^2) using constrained optimization --\n")
print(posterior.vcov)
