#### Code for preparing Upper Midwest PLS and pollen datasets for Camp Peon 2016
#### Chris Paciorek and Simon Goring
#### July 2016

setwd(file.path('data', 'upperMidwest'))

### PLS data

## preparatory step: download western_comp_v0.6-2.csv from Wiki

dat <- read.csv('western_comp_v0.6-2.csv')

coreTaxa <- c('Ash', 'Basswood', 'Beech', 'Birch', 'Cedar.juniper', 'Elm', 'Fir', 'Hemlock', 'Hickory', 'Ironwood', 'Maple', 'Oak', 'Pine', 'Poplar.tulip.poplar', 'Spruce', 'Tamarack' , 'Walnut', 'Other.hardwood')

other <- c('Alder', 'Black.gum.sweet.gum', 'Cherry', 'Dogwood', 'Buckeye', 'Hackberry', 'Locust', 'Mulberry', 'Unknown.tree', 'Sycamore', 'Willow')

## omit some columns (bald cypress because we have no "other conifer" category
omit <- c("Bald.cypress", "No.tree")

taxa <- names(dat)
taxa <- taxa[!taxa %in% c('x', 'y')]
taxa <- taxa[!taxa %in% omit]

dat <- dat[ , c('x', 'y', taxa)]

## combine other hardwood taxa
other <- taxa[!taxa %in% coreTaxa]
taxa <- taxa[!taxa %in% other]
dat[ , "Other.hardwood"] <- dat[ , "Other.hardwood"] + rowSums(dat[ , other])
dat <- dat[ , c('x', 'y', taxa)]

write.table(dat, file = 'pls.csv', quote = FALSE, sep = ",", row.names = FALSE)

### Pollen data from Neotoma

#  Simon Goring - July 25, 2014, updated by Chris Paciorek July 21, 2016
#
#  1.  Get all sites (and assoc. metadata) within a bounding box around Tower, bounding box should be big 
#      enough to get ~30 sites.
#  2.  Download data and convert to 'paleon' taxonomy
#  3.  Generate pollen diagrams for each site going back through time (highlighting the set era 
#      assemblage?)
#  4.  Create data.frame with one settlement era pollen assemblage per site (as rows), including 
#      a. site name
#      b. age
#      c. lat/long and albers x/y 
#      d. taxa as counts (or proportions?)

library(neotoma) 
library(rgdal)

towerLoc <- c(-(86 + 2/60 + 13/3600), 46 + 32/60 + 31/3600)

##  These are the sites around UNDERC and Tower Hill Lake.
##  This gives us 46 records (on 25 July, 2014) and 72 records (on 21 July, 2016)
boundingBox <- c(-92, 44, -85, 48)

boxSites <- get_dataset(loc = boundingBox, datasettype = 'pollen')

allSites <- get_download(boxSites)

allCompiled <- compile_taxa(allSites, list.name = 'P25')

allData <- compile_downloads(allCompiled)

coordinates(allData) <- c("long", "lat")
proj4string(allData) <- CRS("+proj=longlat")
tmp <- spTransform(allData, CRSobj=CRS('+init=epsg:3175'))

pollen <- as.data.frame(allData)
pollen$x <- coordinates(tmp)[ , 1]
pollen$y <- coordinates(tmp)[ , 2]

transTable <- read.csv(file.path('..', '..', 'pollenTranslation.csv'), stringsAsFactors = FALSE)

omit <- transTable$latin[transTable$common == 'omit']

colnames <- names(pollen)
pollen <- pollen[ , colnames[!colnames %in% omit]]

others <- transTable$latin[transTable$common == "Other.hardwood"]

pollen$Other.hardwood <- rowSums(pollen[ , others], na.rm = TRUE)

pollen <- pollen[ , !names(pollen) %in% others] 

toConvert <- names(pollen)[names(pollen) %in% transTable$latin]

for(nm in toConvert) {
  common <- transTable$common[transTable$latin == nm]
  pollen[[nm]][is.na(pollen[[nm]])] <- 0
  if(!common %in% names(pollen)) {
    pollen[[common]] <- pollen[[nm]]
  } else {
    pollen[[common]] <- pollen[[common]] + pollen[[nm]]
  }
  pollen[nm] <- NULL
}

## downloaded July 21 2016 11:45 am PDT
write.table(pollen, file = 'pollenTimeSeries.csv', quote = FALSE, sep = ",", row.names = FALSE)

### Pollen at settlement time

pollenBySite <- split(pollen, pollen$site.name)

getSettlement <- function(x) {
  x <- x[x$age > 1950 - 1830 & x$age < 1950 - 1700, ]
  if(nrow(x)) {
    wh <- which.min(abs(x$age - (1950 - 1800)))
    return(x[wh, ])
  } else return(NULL)
}


settlementPollen <- sapply(pollenBySite, getSettlement)
settlementPollen <- do.call(rbind, settlementPollen)

write.table(settlementPollen, file = 'settlementPollen.csv', quote = FALSE, sep = ",", row.names = FALSE)


