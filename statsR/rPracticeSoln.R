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
