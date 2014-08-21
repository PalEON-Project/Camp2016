### rIntro

## Challenge: boxplot split by taxa

boxplot(full[ , taxaCols], xlab = 'taxon', ylab = 'proportion')
boxplot(full$oak, full$hickory)

## Challenge: extract even rows

even <- seq(1, nrow(full), by = 2)
full[even, ]

## Challenge: find median across all taxa by pond

lapply(full[ , taxaCols], median)

## Creating a function
## example: plotting proportions on a map

## basic code for a single case

require(maps)
x <- full$lon
y <- full$lat
z <- full$oak
propBreaks = c(0, 0.01, 0.05, 0.10, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1)
col <- rev(terrain.colors(length(propBreaks)))
cats <- cut(z, breaks = propBreaks)
cats <- as.integer(cats)
plot(x, y, xlim = c(-73.5, -70), ylim = c(41.3, 43.5), col = col[cats], pch = 16)
map("state", add = TRUE)

### rPractice

## 1)

setwd('data/upperMidwest')
poll <- read.csv('settlementPollen.csv')
pls <- read.csv('pls.csv')

## 2)
library(maps)

poll[1, ]
pls[1, ]

par(mfrow = c(1, 2))
plot(poll$long, poll$lat, pch = 16)
map('state', add = TRUE)

plot(pls$x, pls$y)
points(poll$x, poll$y, col = 'red')

plsSub <- pls[pls$x < max(poll$x) & pls$x > min(poll$x) &
pls$y < max(poll$y) & pls$y > min(poll$y), ]

plsTaxaCols <- 3:ncol(plsSub)
pollTaxaCols <- 12:ncol(poll)

identical(sort(names(plsSub)[plsTaxaCols]),
  sort(names(poll)[pollTaxaCols]))


taxa <- names(pls)[plsTaxaCols]

plsTotal <- colSums(pls[ , taxa], na.rm = TRUE)
pollTotal <- colSums(poll[ , taxa], na.rm = TRUE)

plsProp <- plsTotal / sum(plsTotal)
pollProp <- pollTotal / sum(pollTotal)

plot(plsProp, pollProp)

which(plsProp > .19)
which(pollProp > .3)

## 3)

library(fields)
dists <- rdist(pls[ , c('x', 'y')], poll[ , c('x', 'y')])

dim(dists)

matches <- apply(dists, 2, which.min)
plsMatch <- pls[matches, ]

## 4)

plsMatch[ , taxa] <- plsMatch[ , taxa] / rowSums(plsMatch[ , taxa], na.rm = TRUE)
poll[ , taxa] <- poll[ , taxa] / rowSums(poll[ , taxa], na.rm = TRUE)

nTaxa <- length(taxa)

par(mfrow = c(4, 5), mai = c(.4, .4, .3, .1))
for(taxon in taxa) {
plot(plsMatch[ , taxon], poll[ , taxon], main = taxon,  xlab = 'pls', 
     ylab = 'pollen')
abline(0, 1)
}



## 5)

pollTime <- read.csv('pollenTimeSeries.csv')

ponds <- unique(pollTime[ , c('sitename', 'lat', 'long', 'dataset', 'x', 'y')]) 
# notice the repetition - multiple cores per lake?

## 6)

glimmer <- pollTime[pollTime$sitename == "Glimmerglass Lake", ]
head(glimmer)
glimmerTaxaCols <- 12:ncol(glimmer)
taxa <- names(glimmer)[glimmerTaxaCols]
glimmer[ , taxa] <- glimmer[ , taxa] / rowSums(glimmer[ , taxa])

par(mfrow = c(2, 3))
for(taxon in c('birch', 'hemlock', 'oak', 'pine', 'fir', 'spruce')) {
  plot(-glimmer$age, glimmer[[taxon]], 
       main = taxon, type = "l")
  points(-glimmer$age, glimmer[[taxon]])
}

## 7)

dists <- rdist(ponds[ponds$sitename == "Glimmerglass Lake", c('x', 'y')], ponds[ , c('x', 'y')])

ord <- order(dists)
nearby <- ponds$sitename[ord[1:10]]


pollTime[ , taxa ] <- pollTime[ , taxa ] / rowSums(pollTime[ , taxa])

site <- nearby[1]
sub <- pollTime[pollTime$sitename == site, ]
plot(-sub$age, sub$hemlock, xlim = range(-sub$age), 
     type = 'l')

cols <- c('black', 'darkgrey', 'lightgrey', 'red1', 'orange', 'red2', 'green', 'purple', 'lightblue', 'darkblue')
cnt <- 2
for(site in nearby[2:length(nearby)]) {
  sub <- pollTime[pollTime$sitename == site, ]
  lines(-sub$age, sub$hemlock, col = cols[cnt])
  cnt <- cnt + 1
}

## 8)
library(mgcv)

# get glimmer again, as need counts not proportions
pollTime <- read.csv('pollenTimeSeries.csv')
glimmer <- pollTime[pollTime$sitename == "Glimmerglass Lake", ]

glimmer$total <- rowSums(glimmer[ , glimmerTaxaCols])
y <- cbind(glimmer[[taxon]], glimmer$total - glimmer[[taxon]])
x <- -glimmer$age
mod <- gam(y ~ s(x), family = 'binomial')

newAges <- seq(min(x), max(x), length = 200)
preds <- predict(mod, newdata = list(x = newAges), 
                 type = 'response', se.fit = TRUE)

par(mfrow = c(1, 1))
taxon <- 'hemlock'
plot(-glimmer$age, glimmer[[taxon]] / glimmer$total, 
     main = taxon, type = "l")
points(-glimmer$age, glimmer[[taxon]] / glimmer$total)
lines(newAges, preds$fit, col = 'red')
lines(newAges, preds$fit + 2*preds$se.fit, col = 'blue')
lines(newAges, preds$fit - 2*preds$se.fit, col = 'blue')

### forwardModel stats practice problems

setwd('data')
draws <- read.csv('drawsOak.csv')
dim(draws)

## 1)

par(mfrow = c(2, 3))
for(i in seq_len(ncol(draws)))
  hist(draws[ , i])
for(i in seq_len(ncol(draws)))
  plot(seq_len(nrow(draws)), draws[ , i], type = 'l')

summary(draws)

## 2)

plot(draws)

## 3)

diff <- draws$loc1.1750 - draws$loc2.1750
c(mean(diff), sd(diff))
ci_diff <- quantile(diff, c(.025, .975))

ci_1 <- quantile(draws$loc1.1750, c(.025, .975))
ci_2 <- quantile(draws$loc2.1750, c(.025, .975))


### forwardModel simulate 

## 1)
logit <- function(x) log(x/(1-x))
expit <- function(x) exp(x)/(1+exp(x))

simRW <- function(nTimes = 20, n = 40, tau = 1, theta0 = 0.5) {
  if(length(n) == 1) n <- rep(n, nTimes)
  y <- logitTheta <- rep(0, nTimes)
  logitTheta[1] <- rnorm(1, logit(theta0), tau)
  y[1] <- rbinom(1, n[1], expit(logitTheta[1]))
  for(i in 2:nTimes) {
    logitTheta[i] <- rnorm(1, logitTheta[i-1], tau)
    y[i] <- rbinom(1, n[i], expit(logitTheta[i]))    
  }
  return(list(n = n, y = y, theta = expit(logitTheta)))
}

nT <- 40
set.seed(0)
vals <- simRW(nT)

plot(seq_len(nT), vals$y)
lines(seq_len(nT), vals$y)

plot(seq_len(nT), vals$theta, ylim = c(0,1), type = 'l')

## generate data series for students to assimilate

set.seed(999)
vals <- simRW(nTimes = 20, n = 50, tau = 0.2, theta0 = 0.2)
write.csv(data.frame(n = vals$n, y = vals$y), file = 'tsData.csv', row.names = FALSE)
write(vals$theta, file = 'tsHidden.csv', sep = ',', ncol = 20)

## 2)

set.seed(0)
nReps <- 30
thetas <- matrix(0, nrow = nT, ncol = nReps)
for(j in seq_len(nReps)) {
  vals <- simRW(nT, tau = 0.2, theta0 = 0.2)
  thetas[ , j] <- vals$theta
}

plot(seq_len(nT), thetas[ , 1], type = 'l', ylim = c(0,1))
for(j in 2:nReps)
  lines(seq_len(nT), thetas[ , j], col = j)

## 3)

data <- read.csv('tsData.csv')
nT <- nrow(data)

logLikRW <- function(theta, n, y) {
  logLik <- dbinom(y, n, theta, log = TRUE)
  return(sum(logLik))
}


set.seed(0)
nReps <- 10000
thetas <- matrix(0, nrow = nT, ncol = nReps)
for(j in seq_len(nReps)) {
  vals <- simRW(nT, tau = 0.2, theta0 = 0.2)
  thetas[ , j] <- vals$theta
}

logLiks <- apply(thetas, 2, logLikRW, data$n, data$y)
ord <- order(logLiks, decreasing = TRUE)

plot(seq_len(nT), data$y/data$n, col = 'orange')
for(j in seq_len(5))
  lines(seq_len(nT), thetas[ , ord[j]], col = 'blue')
hidden <- scan('tsHidden.csv', sep = ',')
  lines(seq_len(nT), hidden, col = 'red')

plot(seq_len(nT), data$y/data$n, ylim = c(0,1))
for(j in seq_len(nReps))
  lines(seq_len(nT), thetas[ , j])
for(j in seq_len(5))
  lines(seq_len(nT), thetas[ , ord[j]], col = 'blue')
points(seq_len(nT), data$y/data$n, ylim = c(0,1), col = 'orange')

### bayesExercise

data <- read.csv(file.path('data','rings.csv'))
data <- data[data$SPP == "QURU", ]
ringCols <- grep("X[0-9]{4}", names(data))

ringModel <- function() {

  for(i in 1:N) {
    D[i, 1] ~ dunif(0,1000)
    # note that ~ dlnorm(0, .001) is very peaked near zero with 
    # long tail, so hard to get a 'non-informative' prior this way
    beta_i[i] ~ dnorm(beta0, tau_prec)
    
    for(j in 1:nT) {
      logDobs[i, j] ~ dnorm(log(D[i, j]), w_prec)
      logXobs[i, j] ~ dnorm(log(X[i, j]), v_prec)
      X[i, j] ~ dlnorm(beta_i[i] + beta_t[j], sigma_prec)
    }
    for(j in 2:nT) {
      D[i, j] <- D[i, j-1] + X[i, j-1]
    }
  }
  for(j in 1:nT) {
    beta_t[j] ~ dnorm(0, phi_prec)
  }
  
  beta0 ~ dnorm(0, .00001)
  w_prec <- 1 / (w_sd^2)
  v_prec <- 1 / (v_sd^2)
  sigma_prec <- 1 / (sigma_sd^2)
  phi_prec <- 1 / (phi_sd^2)
  tau_prec <- 1 / (tau_sd^2)
  w_sd ~ dunif(0, 100)
  v_sd ~ dunif(0, 100)
  sigma_sd ~ dunif(0, 100)
  phi_sd ~ dunif(0, 100)
  tau_sd ~ dunif(0, 100)
}


N <- nrow(data)
nT <- length(ringCols)
logDobs <- log(cbind(matrix(NA, N, nT-2), data$DBH11, data$DBH12))
logXobs <- log(as.matrix(data[ , ringCols]))

# run MCMC

require(R2jags, quietly = TRUE) 
out <- jags(data = list(logDobs = logDobs, logXobs = logXobs, N = N, nT = nT),
  parameters.to.save = c('D[1,1]', 'D[1,30]','w_sd', 'v_sd', 'tau_sd', 'phi_sd', 'sigma_sd','beta0', 'beta_i', 'beta_t'), inits = list(list(w_sd=1,v_sd=1, tau_sd=1, phi_sd=1,sigma_sd=1)), 
  n.chains = 1,
  n.iter = 5000, n.burnin = 1000, model.file = ringModel, DIC = FALSE)
out_noInit <- jags(data = list(logDobs = logDobs, logXobs = logXobs, N = N, nT = nT),
            parameters.to.save = c('D[1,1]', 'D[1,30]','w_sd', 'v_sd', 'tau_sd', 'phi_sd', 'sigma_sd','beta0', 'beta_i', 'beta_t'), 
            n.chains = 1,
            n.iter = 5000, n.burnin = 1000, model.file = ringModel, DIC = FALSE)

# MCMC diagnostics
tsplot <- function(x) plot(seq_along(x), x, type = 'l')

out.mcmc <- as.mcmc(out)[[1]]

tsplot(out.mcmc[ , 'beta0'])

beta <- out.mcmc[ , paste('beta[', 1:N, ']', sep = '')]

par(mfrow = c(4, 4))
for(i in 1:N)
  tsplot(beta[ , i])

beta_t <- out.mcmc[ , paste('beta_t[', 1:nT, ']', sep = '')]

par(mfrow = c(5, 6))
for(i in 1:nT)
  tsplot(beta_t[ , i])

par(mfrow = c(2,3))
tsplot(out.mcmc[ , 'w_sd'])
tsplot(out.mcmc[ , 'v_sd'])
tsplot(out.mcmc[ , 'tau_sd'])
tsplot(out.mcmc[ , 'phi_sd'])
tsplot(out.mcmc[ , 'sigma_sd'])

id <- 1
Xvals <- out.mcmc[ , paste('X[', id, ',', 1:nT, ']', sep = '')]
plot(1:nT, exp(logXobs[id, ]), type = 'l')
lines(1:nT, colMeans(Xvals), col = 'red')
lines(1:nT, apply(Xvals, 2, quantile, .025), col = 'blue')
lines(1:nT, apply(Xvals, 2, quantile, .975), col = 'blue')



