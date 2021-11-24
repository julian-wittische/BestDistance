# Line 332 to 341 of this function will not dynamically function with a
# different number of rasters and have values different from defaults

source("Libs&Functions.R")

empir.sim(catraster = catraster_SA_coarser_cropped,
          geosites = geosites,
          parallel=3,
          iters=1,
          loci=16,
          alleles=10,
          n_ind = 5000)
1

cdpop_sim1 <- readRDS("cdpop_sim_TEST/Results/iter__1/cdpop_sim.rds")

cdpop_sim1$grid_list$gen_101@other$xy

a_tab <- adegenet::tab(cdpop_sim1$grid_list$gen_101)
pc <- prcomp(a_tab)
pc_dist <- as.matrix(dist(pc$x[,1:10]))
pc_dist
geo_dist <- as.matrix(dist(cdpop_sim1$grid_list$gen_101@other$xy))
plot(log(geo_dist), pc_dist)
mantel.randtest(dist(pc$x[,1:10]), dist(cdpop_sim1$grid_list$gen_101@other$xy))

# To calculate Loiselle's kinship I need to use the function eco.kin.Loiselle() from the EcoGenetics package
# First I need to convert genind in ecogen using ecogen2genind()
devtools::install_github("leandroroser/EcoGenetics-devel")
