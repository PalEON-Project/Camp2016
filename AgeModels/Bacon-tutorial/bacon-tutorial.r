# update neotoma
install.packages("devtools")
library(devtools)
install_github("ropensci/neotoma")
library(neotoma)

# load packages
library(neotoma)
library(ggplot2)
library(reshape2)
#library(mapdata)

# set your working directory!

# load Bacon
source('Bacon.R')

# function to compute posterior age estimates for specified depths
bacon_age_posts <- function(d, b.depths, out, thick)
{ 
  its=out[,1:(ncol(out)-1)]
  
  elbows <- cbind(its[,1])
  accs <- its[,2:(ncol(its)-1)]
  for(i in 2:ncol(accs))
    elbows <- cbind(elbows, elbows[,ncol(elbows)] + (thick * accs[,i-1]))
  
  if (d %in% b.depths)
    ages <- elbows[,which(b.depths == d)] 
  else
  {
    maxd <- max(which(b.depths < d))
    ages <- elbows[,maxd] + ((d-b.depths[maxd]) * accs[,maxd])
  }
  ages
}

# we are going to look for irwin smith
# we can get the dataset id from neotoma explorer, or using get_site

# get the data
pol  <- get_download(13047)

site <- get_site(pol)

# if you have the package mapdata installed
# A crude way of making the oceans blue.
# par(mfrow=c(1,1))
# plot(1, type = 'n',
#      xlim=range(site$long)+c(-30, 30),
#      ylim=range(site$lat)+c(-30, 30),
#      xlab='Longitude', ylab = 'Latitude')
# rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "lightblue")
# map('world',
#     interior=TRUE,
#     fill=TRUE,
#     col='gray',
#     xlim=range(site$long)+c(-30, 30),
#     ylim=range(site$lat)+c(-30, 30),
#     add=TRUE)
# points(site$long, site$lat, pch=19, cex=2, col='red')

# get the geochron
geochron <- get_geochron(get_site(pol))

# there are two chronologies
pol[[1]]$chronologies
pol[[1]]$chronologies[[1]]$chronology.id[1]
pol[[1]]$chronologies[[2]]$chronology.id[1]

# get the chronological controls
cc1 = get_chroncontrol(pol[[1]]$chronologies[[1]]$chronology.id[1])
# cc2 = get_chroncontrol(pol[[1]]$chronologies[[2]]$chronology.id[1])

# let's plot the controls
plot_chroncontrols <- function(x) {
  
  dat = data.frame(x$chron.control)
  dat$age.young
  
  dat$vals = dat$control.type %in% c('Radiocarbon')
  
  ggplot(dat) + geom_point(aes(x=depth, y=age, colour=vals)) + 
    geom_errorbar(data=dat, aes(x=depth, ymax=age.old, ymin=age.young, colour=vals)) + 
    xlab('Age') + ylab('Depth') + scale_color_manual("Control type", values=c('#F8766D', '#00BF7D'))
 
}

plot_chroncontrols(cc1)
# plot_chroncontrols(cc2)

write_agefile(pol[[1]], chronology = 1, path = ".",
              corename = "IRWIN1", cal.prog = 'Bacon') 
# write_agefile(pol[[1]], chronology = 2, path = ".",
#               corename = "IRWIN2", cal.prog = 'Bacon') 


# fix plotting limits?
thick = 50
core_bac <- Bacon('IRWIN1', acc.mean = 10, thick = thick, #plot.pdf = FALSE, 
                  depths.file = TRUE, suggest = FALSE, ask = FALSE)

# Add the Bacon model back into the object:
pol[[1]] <- read_bacon(x = "IRWIN1", sections = 5, add = TRUE, download = pol[[1]], chron_name = "Bacon")
pol[[1]] <- compile_taxa(pol[[1]], "P25")

trees_shrubs <- pol[[1]]$taxon.list$compressed[which(pol[[1]]$taxon.list$ecological.group == "TRSH")]

irwin_pol <- counts(pol[[1]]) / rowSums(counts(pol[[1]]))
irwin_pol <- irwin_pol[,colnames(irwin_pol) %in% c("Pinus", "Betula", "Tsuga", "Fagus")]

output <- data.frame(age = c(pol[[1]]$chronologies$Bacon$age,
                             pol[[1]]$chronologies$`Booth et al. 2012`$age),
                    chron = rep(c("Bacon", "Booth"), each = nrow(irwin_pol)),
                    rbind(irwin_pol, irwin_pol))
          
long_table <- melt(output, id.vars = c("age", "chron"))

ggplot(long_table, aes(y = value, x = age)) + 
  geom_line(aes(color = chron)) + 
  facet_wrap(~variable, nrow = 1) + 
  coord_flip() +
  theme_bw() +
  xlab('Age (YBP)') + 
  ylab('Proportion') + scale_x_reverse()


# sometimes we might need a break in the core - this is a hiatus
# now let's add a hiatus
acc.mean.mod  = 3.02
acc.mean.old  = 15
acc.shape.mod = 0.53
acc.shape.old = 0.9
mem.strength  = 4
mem.mean      = 0.7
hiatus.depth  = 16.5
thick = 5

# hiatus gamma prior
hiatus.shape = 0.1
hiatus.mean  = 1/100

par(mfrow=c(1,1))
x = seq(0.00001,20,length=30)
y = dgamma(x, shape=hiatus.shape, rate=hiatus.shape/hiatus.mean)
plot(x,y, xlim=c(0,20), ylim=c(0,8), type='l')

core_bac_hiatus <- Bacon('IRWIN', 
                         acc.mean      = c(acc.mean.mod, acc.mean.old), 
                         acc.shape     = c(acc.shape.mod, acc.shape.old),
                         #                  acc.shape     = acc.shape.val,
                         mem.strength  = mem.strength,
                         mem.mean      = mem.mean,
                         thick         = thick,
                         ask           = FALSE,
                         suggest       = FALSE,
                         # depths.file   = FALSE, # we could pass one
                         hiatus.shape  = hiatus.shape,
                         hiatus.mean   = hiatus.mean,
                         hiatus.depths = hiatus.depth)




# # Add the Bacon model back into the object:
# pol[[1]] <- read_bacon(x = "IRWIN", sections = 45, add = TRUE, download = pol[[1]], chron_name = "Bacon2")
# pol[[1]] <- compile_taxa(pol[[1]], "P25")
# 
# trees_shrubs <- pol[[1]]$taxon.list$compressed[which(pol[[1]]$taxon.list$ecological.group == "TRSH")]
# 
# irwin_pol <- counts(pol[[1]]) / rowSums(counts(pol[[1]]))
# irwin_pol <- irwin_pol[,colnames(irwin_pol) %in% c("Pinus", "Betula", "Tsuga", "Fagus")]
# 
# output <- data.frame(age = c(pol[[1]]$chronologies$Bacon$age,
#                              pol[[1]]$chronologies$Bacon2$age,
#                              pol[[1]]$chronologies$`Booth et al. 2012`$age),
#                      chron = rep(c("Bacon", "Bacon-hiatus", "Booth"), each = nrow(irwin_pol)),
#                      rbind(irwin_pol, irwin_pol, irwin_pol))
# 
# long_table <- melt(output, id.vars = c("age", "chron"))
# 
# ggplot(long_table, aes(y = value, x = age)) + 
#   geom_line(aes(color = chron)) + 
#   facet_wrap(~variable, nrow = 1) + 
#   coord_flip() +
#   theme_bw() +
#   xlab('Age (YBP)') + 
#   ylab('Proportion') + scale_x_reverse()
