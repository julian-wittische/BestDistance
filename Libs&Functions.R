################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

# This is largely based on work by William Peterman and Kristopher Winiarski

# Loading necessarey packages
library(corMLPE)
library(ResistanceGA)
library(radish)
library(RandomFields)
library(raster)
library(gdistance)
library(adegenet)
library(readr)
library(spatialEco)
library(sp)
library(EcoGenetics)
library(RClone)

source("cdpop_from_R_function_JW.R")
source("ResGA_CDPOP_empir_simFunction_JW.R")

geosites <- SpatialPoints(read.table("Data/geo.txt", header=TRUE)[,2:3],
                          CRS(SRS_string = "EPSG:3044"))

catraster <- raster("Data/all5.asc", crs = "EPSG:3044")
catraster <- reclassify(catraster, cbind(c(3,4),c(4,3))) 

# PURPOSE: Landscape is a bit too big and it would make sense to have more than
# 1 individual per cell
# Let us try to keep linear features intact but have a coarser resolution

# Create garbage raster with 120m res
garbage_120m <- aggregate(catraster, fact=3, fun=max, na.rm=TRUE) #max
# Resample original using ngb
set.seed(1) #ngb is random so final result may differ between runs
catraster_resamp <- resample(catraster, garbage_120m, method = "ngb") 

catbrick <- layerize(catraster, falseNA=FALSE)
catbrick <- catbrick*c(1,2,3,4,5) 
catbrick_agg <- aggregate(catbrick, fact=3, fun=max, na.rm=TRUE)

catbrick_agg$X0[catbrick_agg$X0==0] <- NA
catbrick_agg$X1[catbrick_agg$X1==0] <- NA
catbrick_agg$X2[catbrick_agg$X2==0] <- NA
catbrick_agg$X3[catbrick_agg$X3==0] <- NA
catbrick_agg$X4[catbrick_agg$X4==0] <- NA
catbrick_agg <- catbrick_agg-1
plot(merge(catraster_resamp, catbrick_agg$X3))
lol <- mosaic(catbrick_agg$X3, catbrick_agg$X4, fun=max)

catraster_coarser <- merge(lol,catraster_resamp)

################################################################################
####### Checking step: Is the DYI downsampling  meaningful for analysis? #######
# plot(catraster)
# plot(catraster_coarser)
################################################################################
     
catraster_coarser_cropped <- crop(catraster_coarser,
                                  extent(geosites)+c(-360,360,-360,360))
catraster_coarser_cropped <- reclassify(catraster_coarser_cropped,
                                        cbind(c(3,4),c(4,3)))
plot(catraster_coarser_cropped)
plot(geosites, add=TRUE, pch = 19)
plot(catraster)
plot(geosites, add=TRUE, pch = 19)

# Find the empirical distribution of lizard abundance at our new res
lizpercell <- as.matrix(table(extract(catraster_coarser_cropped, geosites,
                      cellnumbers = TRUE,
                      fun=length)[,1]))

# Use a narrow smoother on it and get a reasonable number of individuals on it
lizgrid <- catraster_coarser_cropped
lizgrid[] <- 0
lizcellrowcol <- rowColFromCell(lizgrid, as.numeric(names(lizpercell[,1])))
lizgrid[lizcellrowcol] <- lizpercell
lizgridsmooth <- raster.gaussian.smooth(lizgrid, sigma=1, n=5)*200

################################################################################
########### Checking step: Is the smoothing  meaningful for analysis? ##########
# plot(lizgridsmooth)
# sum(na.omit(getValues(lizgridsmooth)))
################################################################################

# BELOW: Same but smaller study area
# CDPOP does not allow for a K limit PER cell so I have to change my strategy
catraster_SA <- raster("Data/study_area5.asc", crs = "EPSG:3044")
catraster_SA <- reclassify(catraster_SA, cbind(c(3,4),c(4,3))) 
plot(catraster_SA)
plot(geosites, add=TRUE, pch = 19)
plot(catraster_SA, colNA="blue")

garbage_120m <- aggregate(catraster_SA, fact=3, fun=max, na.rm=TRUE) #max
# Resample original using ngb
set.seed(1) #ngb is random so final result may differ between runs
catraster_SA_resamp <- resample(catraster_SA, garbage_120m, method = "ngb") 

catbrick_SA <- layerize(catraster_SA, falseNA=FALSE)
catbrick_SA <- catbrick_SA*c(1,2,3,4,5)
catbrick_SA_agg <- aggregate(catbrick_SA, fact=3, fun=max, na.rm=TRUE)

catbrick_SA_agg$X0[catbrick_SA_agg$X0==0] <- NA
catbrick_SA_agg$X1[catbrick_SA_agg$X1==0] <- NA
catbrick_SA_agg$X2[catbrick_SA_agg$X2==0] <- NA
catbrick_SA_agg$X3[catbrick_SA_agg$X3==0] <- NA
catbrick_SA_agg$X4[catbrick_SA_agg$X4==0] <- NA
catbrick_SA_agg <- catbrick_SA_agg-1
plot(merge(catraster_SA_resamp, catbrick_SA_agg$X3))
lol_SA <- mosaic(catbrick_SA_agg$X3, catbrick_SA_agg$X4, fun=max)

catraster_SA_coarser <- merge(lol_SA,catraster_SA_resamp)

catraster_SA_coarser_cropped <- crop(catraster_SA_coarser,
                                  extent(geosites)+c(-360,360,-360,360))
catraster_SA_coarser_cropped <- reclassify(catraster_SA_coarser_cropped,
                                        cbind(c(3,4),c(4,3)))
plot(catraster_SA_coarser_cropped, colNA="red")
plot(geosites, add=TRUE, pch = 19)
plot(catraster)
plot(geosites, add=TRUE, pch = 19)

sum(getValues(!is.na(catraster_SA_coarser_cropped))) #1887
