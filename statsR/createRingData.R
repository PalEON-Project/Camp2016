# code to format dbh and ring data from Mike for bayesExercise

trees <- read.csv('../PecanInputs/H 2012 Adult Field Data.csv', stringsAsFactors = FALSE)
trees$DBH11[trees$DBH11 == "missing data"] <- NA
trees$DBH12[trees$DBH12 == "missing data"] <- NA

trees$DBH11 <- as.numeric(trees$DBH11) / 10  # now in mm dbh
trees$DBH12 <- as.numeric(trees$DBH12) / 10


source('matchInventoryRings.R')
source('Read_Tuscon.R')
source('extract.stringCode.R')
rings <- Read_Tuscon("../PecanInputs/Revised 2/")
      
 ## Match observations & format for JAGS
combined <- matchInventoryRings(trees,rings,nyears=30,coredOnly=TRUE)

# restrict to (0, 15) x (0, 15) as that seems to be where sampling is concentrated
combined <- combined[combined$X > 0 & combined$X < 15 & combined$Y < 15 & combined$Y > 0, ]

write.csv(combined, file = 'rings.csv', row.names = FALSE, quote = FALSE)
