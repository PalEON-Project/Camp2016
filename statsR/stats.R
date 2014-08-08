
## @knitr distns
y <- seq(-3, 3, length = 300)
par(mfrow = c(1,1))
plot(y, dnorm(y), type = 'l')
pnorm(-2)
pnorm(-1.96)
qnorm(0.975)
n <- 1000000
y <- rnorm(n)
# empirical analogs of the distributional quantities above
par(mfrow = c(1,2))
hist(y); hist(y, probability = TRUE)
mean(y)
mean(y <= -2)  # sum(y <= 2)/n
sort(y)[round(.975*n)] # empirical 97.5th percentile


## @knitr branching
x = 4; y = 7
if(x > 3 && y < 5){
	print("success")
	print("doing some stuff")
} else{
	print("failure")
	print("doing something different")
}


## @knitr forloop
nSteps <- 100
track <- matrix(NA, nr = nSteps, nc = 2)
track[1, ] <- c(0,0)
bound <- 20
plot(track[1, 1], track[1, 2], xlim = c(-bound, bound), ylim = c(-bound, bound), 
  	pch = 16, xlab = 'x', ylab = 'y')
for(it in 2:nSteps){
	track[it, ] <- track[it-1, ] + sample(c(-1,1), 2, replace = TRUE)
	arrows(track[it-1, 1], track[it-1, 2], track[it, 1], track[it, 2],
		length = .05, angle = 20)
}


## @knitr whileloop
maxSteps <- 1000
track <- matrix(NA, nr = maxSteps, nc = 2)
track[1, ] <- c(0,0)
bound <- 20
plot(track[1, 1], track[1, 2], xlim = c(-bound, bound), ylim = c(-bound, bound), 
  	pch = 16, xlab = 'x', ylab = 'y')
it <- 1
while(max(abs(track[it, ])) < bound && it < maxSteps){
	it <- it + 1
	track[it, ] <- track[it-1, ] + sample(c(-1,1), 2, replace = TRUE)
	arrows(track[it-1, 1], track[it-1, 2], track[it, 1], track[it, 2],
		length = .05, angle = 20)
}
title(main = "Computational modern art in R")
cat("Number of steps taken: ", it, ".\n", sep = '')


