################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

### Acknowledgements:
# This is based on previous work by William Peterman and Kristopher Winiarski
# Thanks to Dr Erin Landguth who welcomed me for a short visit in her lab

source("Loading_Scripts&Data.R")

empir.sim(catraster = catraster_SA_coarser_cropped,
          geosites = geosites,
          parallel = 3,
          iters = 1,
          loci = 16,
          alleles = 12,
          n_ind = 4000,
          habitat = 0.5,
          matemoveno = 5,
          matemoveparA = 1,
          matemoveparB = 1.05,
          gridformat = "cdpop",
          Res_rem = 1,
          Res_for = 1,
          Res_urb = 4,
          Res_riv = 8,
          Res_inf = 1)
