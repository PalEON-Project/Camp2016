
## @knitr reading files
getwd()
setwd('~/Desktop/camp2012')
pol <- read.csv('data/newEngl/modPol.csv', header = TRUE)
ponds <- read.csv('data/newEngl/pondInfo.csv') 
help(read.csv) 
help(read.table) 
names(pol) 
head(pol) 
class(pol)


## @knitr objects
class(pol) 
x <- c(7.3, 5, 11) 
class(x) 
is.data.frame(x) 
class(pol$sitename) 
pol$sitename 
as.character(pol$sitename) 
class(pol$SPRUCE) 
object.size(x) 
x <- rnorm(10000000) 
object.size(x) 
print(object.size(x), units = "Mb")  
# that code is a bit opaque and I'm lazy, so let's write our own function: 
size <- function(x){   
  print(object.size(x), units = "Mb") 
} 
size(x) 
ls() 
rm(x) 
ls() 


## @knitr dataframeInfo
names(pol) 
class(names(pol)) 
names(pol) <- tolower(names(pol)) 
names(ponds) <- tolower(names(ponds)) 
cor(pol$oak, pol$spruce) 
dim(pol) 
taxaCols <- 3:12
pol$total <- rowSums(pol[ , taxaCols]) 
pol2 <- pol 
pol2[ , taxaCols] <- pol2[ , taxaCols] / pol2$total # careful here: I know that R will divide each numerator column the denominator column
head(pol2) 
full <- merge(ponds, pol2, by.x = "site", by.y = "sitename", all.x = FALSE, all.y = TRUE) 
head(full) 
summary(full) 
taxaCols <- 9:18
cor(full[ , taxaCols]) 
require(fields, quietly = TRUE)
image.plot(1:10, 1:10, cor(full[ , taxaCols])) # I should change the z-scale and the color scheme


## @knitr initialPlots
require(maps)
plot(full$lon, full$lat, xlim = c(-73.5, -70), ylim = c(41.3,43.5), col = 'red', pch = 16) 
map('state', add = TRUE) 
par(mfrow=c(3,4)) 
for(i in taxaCols){   
  hist(full[ , i], main = names(full)[i], xlab = "proportion") 
} 


## @knitr subsetting
sub <- full[1:5, ] 
full$lon < (-72.5) 
west <- full[full$lon < -72.5, ] 
full[full$lon < -72.5, c("oak", "pine")]
indices <- which(full$lon < -72.5)
indices
full[indices, c(15,16)]
west2 <- subset(full, lon < -72.5, select = c("oak", "pine"))


## @knitr sorting
full <- full[order(full$lat, full$lon), ]
head(full)[ , 1:8]


## @knitr vectorized operations
log(full$oak)
full$coldSpp <- full$beech + full$hemlock + full$spruce 
full$west <- full$lon < -72.5 
full$north <- full$lat > 42.5 
full$north <- FALSE # example of recyling 
full$north[full$lat > 42.5] <- TRUE  
full$nw <- full$west & full$north 
round(full$lat, digits = 1) 


## @knitr with
full <- within(full, warmSpp <- oak + hemlock)
full$warmSpp[1:8]


## @knitr recycling
mat <- matrix(1:2, nrow = 4, ncol = 4) # what do you think will happen? 
mat <- matrix(1:3, nrow = 4, ncol = 4) 


## @knitr lists
myList <- list(a = 7, b = c(8, 9, 11), d = "wtfit?", e = list(first = 7, second = 8)) 
myList$b 
myList[[2]] 
myList$e$first
is.list(full)
myModel = lm(spruce ~ lat, data = full) 
is.list(myModel) 
names(myModel) 
myModel$coefficients 


## @knitr apply
myList = list(1:3, 4:8, 101, 111:150) 
lapply(myList, max) # vs. 
out = rep(NA, length(myList)) 
for(i in 1:length(myList))   
   out[i] <- max(myList[[i]]) 


## @knitr pie
help(pie)
pie(full[1, taxaCols]) 
pie(full[1, taxaCols] + 1e-9) 
pie(full[1, taxaCols] + .1) 
class(full[1, taxaCols]) 
class(c(full[1, taxaCols])) 
class(unlist(full[1, taxaCols])) 
vals <- unlist(full[1, taxaCols]) 
pie(vals) 
pie(vals, col = rainbow(length(vals))) 
help(polygon) 
pie(vals, col = rainbow(length(vals)), border = "NA", main = paste(full$sitename[1], 'Pond'), cex.main = 3) 
pie 


## @knitr pieMap
pieMap <- function(proportions, centers, xlim = c(629000-6000, 809000+6000), ylim = c(4621000-6000, 4801000+6000), radius = NULL, scale = 1, xlab = 'x', ylab = 'y', ...){
# plots multiple pie composition charts as a map
	centers <- as.matrix(centers)
	proportions <- as.matrix(proportions)
	if(is.null(xlim)){
		rg <- range(centers[,1])
		df <- (scale-1) * diff(range(centers[,1]))
		xlim <- c(rg[1]-df, rg[2]+df)
	}
	if(is.null(ylim)){
		rg <- range(centers[ ,2])
		df <- (scale-1) * diff(range(centers[,2]))
		ylim <- c(rg[1]-df, rg[2]+df)
	}
	plot(centers, type = 'n',xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab,...)
	n <- length(centers[,1])
	if(is.null(radius))
		radius <- .025 * diff(range(centers[,1]))
	minVal <- min(proportions)
	if(minVal < 0)
		stop("Some proportions are less than zero.")
	if(length(radius) == 1)
		radius = rep(radius,n)
	for(i in 1:n){
		if(sum(proportions[i,]) > 0){
			minVal <- min(proportions[i,])
			if(minVal == 0){
				warning("Some proportions are zero; proceeding with jittering.")
				eps <- 1e-10 * max(proportions[i,])
				proportions[i,] <- proportions[i,] + eps
			}
			pieAdd(as.vector(proportions[i,]), as.vector(centers[i,]), radius = radius[i] ,...)
		} else{
			points(centers[i,], pch = 'X')
		}
	}
	map('state', add = TRUE)
}

pieAdd <- function (x, center, labels = names(x), edges = 200, radius = 0.8, density = NULL, angle = 45, col = NULL, border = NULL, lty = NULL, ...) # modified from the pie() function in R
{
	if (!is.numeric(x) || any(is.na(x) | x <= 0)) 
		stop("pie: `x' values must be positive.")
	if (is.null(labels)) 
		labels <- rep("",length(x))
	x <- c(0, cumsum(x)/sum(x))
	dx <- diff(x)
	pin <- par("pin")
	nx <- length(dx)
	if (is.null(col)) 
		col <- if (is.null(density)) 
			c("white", "black","lightblue", "red","darkblue","yellow",
				"purple","orange","lightgreen","darkgreen")
		else par("fg")
	col <- rep(col, length.out = nx)
	border <- rep(border, length.out = nx)
	lty <- rep(lty, length.out = nx)
	angle <- rep(angle, length.out = nx)
	density <- rep(density, length.out = nx)
	for (i in 1:nx) {
		n <- max(2, floor(edges * dx[i]))
		t2p <- 2 * pi * seq(x[i], x[i + 1], length = n)
		xc <- c(cos(t2p), 0) * radius + center[1]
		yc <- c(sin(t2p), 0) * radius + center[2]
		polygon(xc, yc, density = density[i], angle = angle[i], 
		border = border[i], col = col[i], lty = lty[i],...)
		t2p <- 2 * pi * mean(x[i + 0:1])
		xc <- cos(t2p) * radius + center[1]
		yc <- sin(t2p) * radius + center[2]
		if (!is.na(lab <- labels[i]) && lab != "") {
			lines(c(1, 1.05) * xc, c(1, 1.05) * yc)
			text(1.1 * xc, 1.1 * yc, lab, xpd = TRUE, 
				adj = ifelse(xc < 0, 1, 0), ...)
		}
	}
	invisible(NULL)
}


## @knitr advanced plots
pdf('treeComp.pdf', height = 4, width = 7)
par(mfrow = c(2,5), mai = c(.25,.25,.3,.1), omi = c(.3, .3, .4, 0))
for(i in taxaCols)
	hist(full[ , i], xlim = c(0, 1), xlab = '', ylab = '', main = names(full)[i])
title(ylab = 'frequency', outer = TRUE, line = 1, cex.lab = 1.5)
title(xlab = 'proportion of pollen', outer = TRUE, line = 1, cex.lab = 1.5)
title(main = paste('Pollen composition of ', nrow(full), ' ponds in New England', sep = ''),
	outer = TRUE, line = 1, cex.main = 1.5)
dev.off()
help(par)
par()$mfrow


## @knitr packages
# install.packages('foreign') 
library(foreign) 
library(help = foreign) 


## @knitr source
source('silly.R') 


## @knitr save
save(full, pol, ponds, file = 'pollenAnalData.RData') 
save.image('pollenFullAnal.RData')

rm(list = ls()) 
ls() 
load('pollenAnalData.RData') 
ls() 


