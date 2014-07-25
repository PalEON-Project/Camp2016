require(R2jags)

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

mixedModel <- function(){
  for(j in 1:nMoms){
    tau2Inv[j,j] <- 1/(tau^2)
    thetas[j] <- theta
  }
  
  for(j in 2:n){
    for(k in 1:(j-1)){
        tau2Inv[j,k] <- 0
        tau2Inv[k,j] <- 0
    }
  }
  mus[1:nMoms] ~ dmnorm(thetas[1:nMoms], tau2Inv[1:nMoms,1:nMoms])

  for(j in 1:n){
    y[j] ~ dnorm(mus[momIds[j]], sig2Inv)
  }
  sig2Inv <- 1/(sig^2)
  theta ~ dnorm(0,.00001)
  sig ~ dunif(0, 100)
  tau ~ dunif(0, 100)
}
  
out <- jags(data = list(y = y, n = n, nMoms = nMoms, momIds = momIds),
  parameters.to.save = c('theta','sig','tau', 'mus'), n.chains = 1,
  n.iter = 2000, n.burnin = 1000, model.file = mixedModel, DIC = FALSE)



# x ~ dmulti(probs, n)
# x should be vector counts
