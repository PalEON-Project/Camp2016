# code to format dbh and ring data from Mike for bayesExercise

trees <- read.csv('../ForestPlots/treecores2014.csv', stringsAsFactors = FALSE)
trees$DBH14 <- as.numeric(trees$DBH14) * 10  # now in mm dbh

source('matchInventoryRings.R')
source('read_tucson.R')
source('extract.stringCode.R')
rings <- Read_Tucson("../ForestPlots/Tucson")

# ok through here

 ## Match observations & format for JAGS
combined <- matchInventoryRings(trees,rings,nyears=30,coredOnly=TRUE)

ringCols <- 10:ncol(combined)
# transform to mm of diameter not radius
combined[ , ringCols] <- 2 * combined[ , ringCols]

write.csv(combined, file = 'data/rings.csv', row.names = FALSE, quote = FALSE)
