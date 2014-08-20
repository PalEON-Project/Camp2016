require(R2jags, quietly = TRUE) 
library(ggplot2)

#  We're only going to use the last 6000 years:
allData <- allData[allData$age < 6000, ]

#  Now we need to bin the data into 200 year bins:
library(reshape2)
melted.pollen <- melt(allData, value.name='count', id=colnames(allData)[1:9])
melted.pollen$agebin <- round(melted.pollen$age, -2)

good.dates <- melted.pollen$date.type %in% c('Calibrated radiocarbon years BP', 'Varve years BP')

melted.pollen <- melted.pollen[good.dates,]

up.pollen <- dcast(melted.pollen, agebin ~ variable, value.var = 'count', fun.aggregate = sum, na.rm=TRUE)

hem.model <- data.frame(age = -up.pollen$agebin,
                        n = rowSums(up.pollen[ , 2:ncol(up.pollen)], na.rm=TRUE),
                        hemlock = up.pollen$Tsuga)

hem.model$pct <- hem.model$hemlock / hem.model$n

hem.model <- hem.model[order(hem.model$age), ]
hem.model <- na.omit(hem.model)

ggplot(hem.model, aes(x = age, y = pct)) +
  geom_point() + geom_smooth()

hmm <- function(){
  mu ~ dunif(-5, 5)
  rho ~ dunif(-1, 1)
  tau ~ dunif(0, 100)
  logitTheta0 ~ dnorm(0, .000001)	
  tau2Inv <- 1/(tau^2)
  
  # latent process evolution and likelihood
  logitTheta[1] ~ dnorm(mu + rho*(logitTheta0 - mu), tau2Inv)
  theta[1] <- exp(logitTheta[1])/(1+exp(logitTheta[1]))  
  y[1] ~ dbin(theta[1], n[1])
  for(i in 2:nT){
    logitTheta[i] ~ dnorm(mu + rho*(logitTheta[i-1] - mu), tau2Inv)
    theta[i] <- exp(logitTheta[i])/(1+exp(logitTheta[i]))
    y[i] ~ dbin(theta[i], n[i])
  }
}

#  And let's run it:
nT <- nrow(hem.model)
y <- hem.model$hemlock
n <- round(hem.model$n)

out <- jags(data = list(nT = nT, n = n, y = y), parameters.to.save = c('rho', 'tau', 'theta'), 
            n.chains = 1, n.iter = 10000, n.burnin = 2000, model.file = hmm, DIC = FALSE)

out.mcmc <- as.mcmc(out)[[1]]
thetaHat <- y / n

colNames <- dimnames(out.mcmc)[[2]]
whichTheta <- grep('theta', colNames)
thetaPost <- out.mcmc[201:1000 , whichTheta]
thetaNames <- dimnames(thetaPost)[[2]]
index <- gsub("theta\\[", "", thetaNames)
index <- as.numeric(gsub("\\]", "", index))
thetaPost <- thetaPost[ , order(index)]

par(mfrow = c(1,1))
plot(hem.model$age, thetaHat, col = 'red')
quants <- apply(thetaPost, 2, quantile, c(.025, .975))
polygon(cbind(c(hem.model$age, rev(hem.model$age), hem.model$age[1]), 
              c(quants[1, ], rev(quants[2, ]), quants[1, 1])), 
        border = NA, col = 'lightblue')
lines(hem.model$age, thetaHat, col = 'red')
points(hem.model$age, thetaHat, col = 'red')
lines(hem.model$age, colMeans(thetaPost), col = 'blue')
title('hemlock in the Upper Peninsula')
abline(v = 1650 - 1950, col = 'grey')
