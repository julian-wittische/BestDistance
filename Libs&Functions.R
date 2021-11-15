################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

# This is largely based on work by William Peterman and Kristopher Winiarski

# Loading necessarey packages
library(corMLPE)
library(ResistanceGA)
library(radish)
library(RandomFields)

source("cdpop_from_R_function.R") # another version is available on WP's GitHub

source("ResistanceGA_CDPOP_simFunction_JW.R")

ResistanceGAMS.sim(parallel=3, iters=1, loci=20, alleles=20)
1
# Line 332 to 341 of this function will not dynamically function with a
# different number of rasters and have values different from defaults