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

install.packages("devtools") 
require(devtools) 
install_github("neotoma", "SimonGoring") 
require(neotoma) 
library(rgdal)
library(ggplot2)

#  This is the explicit location of Tower Hill Lake, that we're not using, but whatever.
towerLoc <- c(-(86 + 2/60 + 13/3600), 46 + 32/60 + 31/3600)

##  These are the sites around UNDERC and Tower Hill Lake.
##  This gives us 46 records (on 25 July, 2014)
boundingBox <- c(-92, 44, -85, 48)

boxSites <- get_dataset(loc = boundingBox, datasettype = 'pollen')

allSites <- get_download(boxSites)

allCompiled <- compile_taxa(allSites, list.name = 'P25')

allData <- compile_downloads(allCompiled)
