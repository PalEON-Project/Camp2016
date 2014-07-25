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

library(neotoma)

#  These are the sites around UNDERC and Tower Hill Lake.
#  This gives us 46 records (on 25 July, 2014)
bounding.box <- c(-92, 44, -85, 48)

box.sites <- get_dataset(loc = bounding.box, datasettype = 'pollen')

all.sites <- get_download(sapply(box.sites, function(x)x$DatasetID))

all.compiled <- lapply(all.sites, compile_taxa, list.name = 'P25')

all.data.frame <- compile_downloads(all.compiled)

#  write.csv(all.data.frame, 'Neotoma_Pollen/data/output/all_site_df.csv')

#  We need to get the chroncontrol tables so we can assign pre-settlement:
all.chronIDs <- sapply(all.sites, function(x)x$sample.meta$ChronologyID[1])

chron.tables <- lapply(all.chronIDs, function(x)try(get_chroncontrol(x)))

find_european <- function(x){
  settlement <- which(x$chron.control$ControlType == 'European')
}
