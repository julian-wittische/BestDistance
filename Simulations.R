################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

### Acknowledgements:
# This is based on previous work by William Peterman and Kristopher Winiarski
# Thanks to Dr Erin Landguth who welcomed me for a short visit in her lab

source("Libraries.R")
################################################################################
### Sourcing two modified functions to get CDPOP to run from R
source("cdpop_from_R_function_JW.R")
source("ResGA_CDPOP_empir_simFunction_JW.R")

# Geographical coordinates
geosites <- SpatialPoints(read.table("Data/geo.txt", header=TRUE)[,2:3],
                          CRS(SRS_string = "EPSG:3044"))

# CDPOP does not allow for a K limit PER cell so I have to change my strategy
catraster_SA <- raster("Data/study_area5.asc", crs = "EPSG:3044")
catraster_SA <- reclassify(catraster_SA, cbind(c(3,4),c(4,3))) 
# plot(catraster_SA)
# plot(geosites, add=TRUE, pch = 19)
# plot(catraster_SA, colNA="blue")

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
lol_SA <- mosaic(catbrick_SA_agg$X3, catbrick_SA_agg$X4, fun=max)

catraster_SA_coarser <- merge(lol_SA,catraster_SA_resamp)

catraster_SA_coarser_cropped <- crop(catraster_SA_coarser,
                                     extent(geosites)+c(-360,360,-360,360))
catraster_SA_coarser_cropped <- reclassify(catraster_SA_coarser_cropped,
                                           cbind(c(3,4),c(4,3)))

empir.sim(catraster = catraster_SA_coarser_cropped,
          geosites = geosites,
          parallel=3,
          iters=1,
          loci=16,
          alleles=12,
          n_ind = 3500,
          habitat=0.5)

