#### Code for preparing Upper Midwest PLS and pollen datasets for Camp Peon 2014
#### Chris Paciorek and Simon Goring
#### July 2014

setwd(file.path('data', 'upperMidwest'))

### PLS data

## preparatory step: download westerncompv0.3.csv from Wiki

dat <- read.csv('westerncompv0.3.csv')

coreTaxa <- c('ash', 'basswood', 'beech', 'birch', 'cedar.juniper', 'elm', 'fir', 'hemlock', 'hickory', 'ironwood', 'maple', 'oak', 'pine', 'poplar.tulip.poplar', 'spruce', 'tamarack' , 'walnut', 'other.hardwood')

other <- c('alder', 'black.gum.sweet.gum', 'cherry', 'dogwood', 'buckeye', 'hackberry', 'locust', 'mulberry', 'unknown.tree', 'sycamore', 'willow')

## omit some columns (bald cypress because we have no "other conifer" category
omit <- c("bald.cypress", "no.tree", "unknown")

taxa <- names(dat)
taxa <- taxa[!taxa %in% c('x', 'y')]
taxa <- taxa[!taxa %in% omit]

dat <- dat[ , c('x', 'y', taxa)]

## combine other hardwood taxa
other <- taxa[!taxa %in% coreTaxa]
taxa <- taxa[!taxa %in% other]
dat[ , "other.hardwood"] <- dat[ , "other.hardwood"] + rowSums(dat[ , other])
dat <- dat[ , c('x', 'y', taxa)]

write.table(dat, file = 'pls.csv', quote = FALSE, sep = ",", row.names = FALSE)

stop()
### Pollen data from Neotoma

#  Simon Goring - July 25, 2014
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

towerLoc <- c(-(86 + 2/60 + 13/3600), 46 + 32/60 + 31/3600)


library(neotoma)
library(rgdal)
library(neotoma)

##  These are the sites around UNDERC and Tower Hill Lake.
##  This gives us 46 records (on 25 July, 2014)
boundingBox <- c(-92, 44, -85, 48)

boxSites <- get_dataset(loc = boundingBox, datasettype = 'pollen')

allSites <- get_download(sapply(boxSites, function(x)x$DatasetID))

allCompiled <- lapply(allSites, compile_taxa, list.name = 'P25')

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

others <- transTable$latin[transTable$common == "other.hardwood"]

pollen$other.hardwood <- rowSums(pollen[ , others], na.rm = TRUE)

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

## downloaded July 31 10:30 am PDT
write.table(pollen, file = 'pollenTimeSeries.csv', quote = FALSE, sep = ",", row.names = FALSE)

### Pollen at settlement time

pollenBySite <- split(pollen, pollen$sitename)

getSettlement <- function(x) {
  x <- x[x$age > 1950 - 1830 & x$age < 1950 - 1700, ]
  if(nrow(x)) {
    wh <- which.min(abs(x$age - (1950 - 1800)))
    return(x[wh, ])
  } else return(NULL)
}

getSettlement_oldVersion <- function(x) {
  settlementDate <- 1950 - 1820
  wh <- which.min(abs(x$age - settlementDate))
  if(abs(x$age[wh] - settlementDate) < 50)
    return(x[wh, ]) else return(NULL)
}

settlementPollen <- sapply(pollenBySite, getSettlement)
settlementPollen <- do.call(rbind, settlementPollen)

write.table(settlementPollen, file = 'settlementPollen.csv', quote = FALSE, sep = ",", row.names = FALSE)


## older code from Simon

#  We need to get the chroncontrol tables so we can assign pre-settlement:
if(FALSE) {
  all.chronIDs <- sapply(all.sites, function(x)x$sample.meta$ChronologyID[1])

  chron.tables <- lapply(all.chronIDs, function(x)try(get_chroncontrol(x)))
  
  find_european <- function(x){
    settlement <- which(x$chron.control$ControlType == 'European')
  }
}

