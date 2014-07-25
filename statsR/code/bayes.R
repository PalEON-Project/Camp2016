
## @knitr flatprior
y <- 35; n <- 100
thetaSeq <- seq(0, 1, length = 200)
lik <- dbinom(y, n, thetaSeq)
prior <- rep(1, 200)
alphaPost <- y+1; betaPost <- n-y+1
post <- dbeta(thetaSeq, alphaPost, betaPost)
plot(thetaSeq, post, type = 'l')
lines(thetaSeq, prior, col = 'red')
lines(thetaSeq, lik*80, col = 'blue') # scaling of likelihood is arbitrary -- it's not a density for theta
legend(0, 8, legend = c('prior', 'lik', 'post'), col = c('red','blue','black'), lty = rep(1,3), bty = 'n')


## @knitr infPrior
y <- 35; n <- 100
thetaSeq <- seq(0, 1, length = 200)
lik <- dbinom(y, n, thetaSeq)
alphaPrior <- 3; betaPrior <- 12
prior <- dbeta(thetaSeq, alphaPrior, betaPrior) 
alphaPost <- y+alphaPrior; betaPost <- n-y+betaPrior
post <- dbeta(thetaSeq, alphaPost, betaPost)
plot(thetaSeq, post, type = 'l')
lines(thetaSeq, prior, col = 'red')
lines(thetaSeq, lik*80, col = 'blue')  # scaling arbitrary
legend(0, 8, legend = c('prior', 'lik', 'post'), col = c('red','blue','black'), lty = rep(1,3), bty = 'n')


## @knitr posteriorInference
y <- 35; n <- 100
alphaPost <- y+1; betaPost <- n-y+1
qbeta(c(.025, .5, .975), alphaPost, betaPost)
postMean <- alphaPost/(alphaPost + betaPost)
postSD <- sqrt(alphaPost*betaPost / ((alphaPost + betaPost)^2 * (alphaPost + betaPost + 1)))
postMean
postSD


## @knitr posteriorInferenceWithSample
####### PAY NO ATTENTION TO THE MAN BEHIND THE CURTAIN ######
fakeMCMCsample <- rbeta(2000, alphaPost, betaPost) 
#############################################################

# now use 'fakeMCMCsample' as if it came from an actual MCMC, 

# ignoring the fact that theta has a beta distribution

quantile(fakeMCMCsample, c(.025, .5, .975))
mean(fakeMCMCsample)
sd(fakeMCMCsample)


## @knitr bivarNorm
nIts <- 1000
nPlot <- 40
mu1 <- 0; mu2 <- 1; sig1 <- 1; sig2 <- 1; rho <- .7 
store <- matrix(NA, nrow = nIts, ncol = 2)
luckyStartingValues <- c(.1, 1.2)

store[1, ] <- luckyStartingValues
for(it in 2:nIts){
	curTheta2 <- store[it-1, 2]
	curTheta1 <- rnorm(1, mu1 + sig1*rho*(curTheta2 - mu2)/sig2, sig1*sqrt(1-rho^2))
	curTheta2 <- rnorm(1, mu2 + sig2*rho*(curTheta1 - mu1)/sig1, sig2*sqrt(1-rho^2))
	store[it, ] <- c(curTheta1, curTheta2)
}
# we know the 'true' 'posterior', so let's plot its contours to see how well we are doing
library(mvtnorm)
gridLength <- 40
grid1d <- seq(-8, 3, length = gridLength)
grid2d <- expand.grid(grid1d, grid1d)
densValues <- dmvnorm(grid2d, c(mu1, mu2), matrix(c(sig1^2, rho*sig1*sig2, rho*sig1*sig2, sig2^2), 2, 2))
contour(grid1d, grid1d, matrix(densValues, gridLength))
lines(store[1:nPlot, 1], store[1:nPlot, 2], col = 'blue')

store2 <- matrix(NA, nrow = nIts, ncol = 2)

unluckyStartingValues <- c(-6, -8)
store2[1, ] <- unluckyStartingValues
for(it in 2:nIts){
	curTheta2 <- store2[it-1, 2]
	curTheta1 <- rnorm(1, mu1 + sig1*rho*(curTheta2 - mu2)/sig2, sig1*sqrt(1-rho^2))
	curTheta2 <- rnorm(1, mu2 + sig2*rho*(curTheta1 - mu1)/sig1, sig2*sqrt(1-rho^2))
	store2[it, ] <- c(curTheta1, curTheta2)
}
lines(store2[1:nPlot, 1], store2[1:nPlot, 2], col = 'red')


## @knitr convergence
par(mfrow = c(2,2))
subset <- 1:100
tsplot <- function(vec, ...){
	plot(1:length(vec), vec, xlab = 'index', ylab = '', ...)
}
tsplot(store[subset, 1], main = 'theta1, good start', type = 'l')
tsplot(store[subset, 2], main = 'theta2, good start', type = 'l')
tsplot(store2[subset, 1], main = 'theta1, bad start', type = 'l')
tsplot(store2[subset, 2], main = 'theta2, bad start', type = 'l')


## @knitr usingPosterior
library(coda) # a package for manipulating MCMC output
nonBurn <- 100:nIts  # conservative in this case
use <- store[nonBurn, ]
ess1 <- effectiveSize(use[, 1])
ess2 <- effectiveSize(use[, 2])
postMean <- colMeans(use)
postSD <- apply(use, 2, sd)
uncIntervals <- apply(use, 2, quantile, c(.025, .975))
postCorr <- cor(use[, 1], use[ , 2])
postMean; postSD; uncIntervals; postCorr


## @knitr posteriorFunctionals
derivedQuant <- exp(use[ , 1])
postMeanDer <- mean(derivedQuant)
postCIsDer <- quantile(derivedQuant, c(.025, .975))
postMeanDer; postCIsDer


## @knitr mixedModelData
theta <- 2
sig <- 2
nMoms <- 10
#### IGNORE THE MAN BEHIND THE CURTAIN ####
mus <- rgamma(nMoms, shape = 2, scale = 1)
nPups <- rbinom(nMoms, 15, prob = .4) 
momIds <- rep(1:nMoms, times = nPups)
y <- rnorm(sum(nPups), mus[momIds], sig)
n <- sum(nPups)
###########################################

### model function (note this looks like R code but is really JAGS code for defining a model
mixedModel <- function(){
  # (hyper)priors
  sig ~ dunif(0, 100)
  tau ~ dunif(0, 100)
  theta ~ dnorm(0,.00001)
  # some transformations to get precisions from sd's
  tau2Inv <- 1/(tau^2)
  sig2Inv <- 1/(sig^2)
  # latent value layer
  for(i in 1:nMoms){
    mus[i] ~ dnorm(theta, tau2Inv)
  }
  # likelihood
  for(j in 1:n)
    {
      y[j] ~ dnorm(mus[momIds[j]], sig2Inv)
    }
}


## @knitr mixedModelRun
library(R2jags) # or BRugs or R2WinBUGS or R2OpenBUGS
out <- jags(data = list(y = y, n = n, nMoms = nMoms, momIds = momIds),
  parameters.to.save = c('theta','sig','tau', 'mus'), n.chains = 1,
  n.iter = 2000, n.burnin = 1000, model.file = mixedModel, DIC = FALSE)


## @knitr mixedModelPostProcess
out.mcmc <- as.mcmc(out)
plot(out)
crosscorr.plot(out.mcmc)
mu1 <- c(out.mcmc[ , 2])
sig <- c(out.mcmc[ , 12])
tau <- c(out.mcmc[ , 13])
theta <- c(out.mcmc[ , 14])
plot(mu1, theta)
plot(sig, tau)


## @knitr fakeDLMdata
nT <- 100
n <- 100
### IGNORE THE MAN BEHIND THE CURTAIN #######
theta0 <- 0.2
tau <- 0.2
tau2Inv <- 1/(tau^2)
logitTheta <- p <- y <- rep(NA, nT)
logitTheta0 <- log(theta0/(1-theta0))
logitTheta[1] <- rnorm(1, logitTheta0, tau)
theta[1] <- exp(logitTheta[1])/(1+exp(logitTheta[1]))
y[1] <- rbinom(1, n, theta[1])
for(i in 2:nT){
	logitTheta[i] <- rnorm(1, logitTheta[i-1], tau)
	theta[i] <- exp(logitTheta[i])/(1+exp(logitTheta[i]))
	y[i] <- rbinom(1, n, theta[i])
}
##############################################

# JAGS model definition
dlmFake <- function(){
  # (hyper)priors
  tau ~ dunif(0, 100)
  logitTheta0 ~ dnorm(0, .000001)	
  # deterministic transformation
  tau2Inv <- 1/(tau^2)
    
  # latent process evolution and likelihood
  logitTheta[1] ~ dnorm(logitTheta0, tau2Inv)
  theta[1] <- exp(logitTheta[1])/(1+exp(logitTheta[1]))
  y[1] ~ dbin(theta[1], n)
  for(i in 2:nT){
    logitTheta[i] ~ dnorm(logitTheta[i-1], tau2Inv)
    theta[i] <- exp(logitTheta[i])/(1+exp(logitTheta[i]))
    y[i] ~ dbin(theta[i], n)
  }
}



## @knitr fakeDLMfit
out <- jags(data = list(nT = nT, n = n, y = y), parameters.to.save = c('tau', 'theta'), 
	n.chains = 1, n.iter = 10000, n.burnin = 2000, model.file = dlmFake, DIC = FALSE)
out.mcmc <- as.mcmc(out)
plot(1:nT, theta, type = 'l')
lines(1:nT, colMeans(out.mcmc)[1:nT], col = 'blue')
points(1:nT, y/n, col = 'red')
plot(1:nT, theta, type = 'l')
quants <- apply(out.mcmc[ , 1:nT], 2, quantile, c(.025, .975))
# lines(1:nT, quants[1, ], col = 'blue', lty = 2)
# lines(1:nT, quants[2, ], col = 'blue', lty = 2)
polygon(cbind(c(1:nT, nT:1, 1), c(quants[1, ], rev(quants[2, ]), quants[1, 1])), border = NA, col = 'lightblue')
lines(1:nT, theta)
lines(1:nT, colMeans(out.mcmc)[1:nT], col = 'blue')
points(1:nT, y/n, col = 'red')
plot(1:nrow(out.mcmc), out.mcmc[ , nT+1]) # tau


## @knitr fakeDLMstartingValues
inits <- function(){
	list('tau' = runif(1, .01, 5), 'logitTheta' = rnorm(nT, 0, 6)) 
}
# one can also specify inits as a list of fixed starting values
# see the documentation for jags()
out <- jags(data = list(nT = nT, n = n, y = y), parameters.to.save = c('tau', 'theta'), 
	n.chains = 3, n.iter = 5000, n.burnin = 1000, model.file = dlm, DIC = FALSE)
traceplot(out)
out.mcmc <- as.mcmc(out)
par(mfrow = c(3, 1))
for(i in 1:3)
	plot(1:nrow(out.mcmc[[1]]), out.mcmc[[i]][ , nT+1], type = 'l') #tau


## @knitr readData
data <- read.csv('data/newEngl/PollenTimeSeries.csv')
names(data)[4] <- "calAge"
sub <- subset(data, sitename == 'aino')
sub$calAge <- -sub$calAge
sub <- sub[order(sub$calAge), ]
sub$n <- rowSums(sub[ , 5:14])
nT <- nrow(sub)
plot(sub$calAge, sub$oak/sub$n)
y <- sub$oak
n <- round(sub$n)

# JAGS model definition
dlm <- function(){
  logitTheta0 ~ dnorm(0, .000001)	
  tau ~ dunif(0, 100)
  tau2Inv <- 1/(tau^2)
  
  logitTheta[1] ~ dnorm(logitTheta0, tau2Inv)
  theta[1] <- exp(logitTheta[1])/(1+exp(logitTheta[1]))
  y[1] ~ dbin(theta[1], n[1])
  for(i in 2:nT){
    logitTheta[i] ~ dnorm(logitTheta[i-1], tau2Inv)
    theta[i] <- exp(logitTheta[i])/(1+exp(logitTheta[i]))
    y[i] ~ dbin(theta[i], n[i])
  }
}

## @knitr fitModel
out <- jags(data = list(nT = nT, n = n, y = y), parameters.to.save = c('tau', 'theta'), 
	n.chains = 1, n.iter = 10000, n.burnin = 2000, model.file = dlm, DIC = FALSE)
out.mcmc <- as.mcmc(out)
thetaHat <- sub$oak/sub$n
plot(sub$calAge, thetaHat, col = 'red')
lines(sub$calAge, colMeans(out.mcmc)[1:nT], col = 'blue')
quants <- apply(out.mcmc[ , 1:nT], 2, quantile, c(.025, .975))
polygon(cbind(c(sub$calAge, rev(sub$calAge), sub$calAge[1]), c(quants[1, ], 
	rev(quants[2, ]), quants[1, 1])), border = NA, col = 'lightblue')
lines(sub$calAge, thetaHat, col = 'red')
lines(sub$calAge, colMeans(out.mcmc)[1:nT], col = 'blue')
title('oak at Aino Pond, MA')
abline(v = 1650 - 1950, col = 'grey')


