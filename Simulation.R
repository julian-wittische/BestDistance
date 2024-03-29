################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

### Acknowledgements:
# This is based on previous work by William Peterman and Kristopher Winiarski
# Thanks to Dr Erin Landguth who welcomed me for a short visit in her lab

source("Loading_Scripts&Data.R")

empir.sim(catraster = catraster_SA_coarser_cropped,
          geosites = geosites,
          parallel = 6,
          iters = 1,
          loci = 16,
          alleles = 12,
          n_ind = 3500,
          habitat = 0.5,
          matemoveno = 5,
          matemoveparA = 1,
<<<<<<< HEAD
          matemoveparB = 3,
          gridformat = "cdpop",
          Res_rem = 3,
          Res_for = 2,
          Res_urb = 20,
          Res_riv = 30,
          Res_inf = 1,
          # JULIA_HOME = "C:/Users/Utilisateur/AppData/Local/Programs/Julia-1.7.1/bin/",
          # CDPOP.py = 'C:/Users/Utilisateur/Desktop/Projects/BestDistance/CDPOP-master/src/CDPOP.py',
          # sim_dir = "C:/Users/Utilisateur/Desktop/Projects/BestDistance/cdpop_sim_TEST/",
          seed = 1)
=======
          matemoveparB = 0.25,
          gridformat = "cdpop",
          Res_rem = 3,
          Res_for = 2,
          Res_urb = 15,
          Res_riv = 30,
          Res_inf = 1,
          JULIA_HOME = "C:/Users/Utilisateur/AppData/Local/Programs/Julia-1.7.1/bin/",
          CDPOP.py = 'C:/Users/Utilisateur/Desktop/Projects/BestDistance/CDPOP-master/src/CDPOP.py',
          sim_dir = "C:/Users/Utilisateur/Desktop/Projects/BestDistance/cdpop_sim_TEST/",
          seed = 1111)
>>>>>>> d8711714f3819e4cc44250ff2cc8b366b75c37fc
