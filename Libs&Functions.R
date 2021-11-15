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

source("cdpop_from_R_function.R") # another version is available on WP's GitHub

source("ResistanceGA_CDPOP_simFunction_JW.R")

source("ResGA_CDPOP_empir_simFunction_JW.R")

memory.limit()

geosites <- SpatialPoints(read.table("Data/geo.txt", header=TRUE)[,2:3],
                          CRS(SRS_string = "EPSG:3044"))

catraster <- raster("Data/all5.asc", crs = "EPSG:3044")
catraster <- crop(catraster, extent(geosites)+120)
plot(catraster)

empir.sim(catraster = catraster,
          geosites = geosites,
          parallel=3,
          iters=1,
          loci=17,
          alleles=20)
1
# Line 332 to 341 of this function will not dynamically function with a
# different number of rasters and have values different from defaults