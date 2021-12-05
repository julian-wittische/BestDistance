################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

source("Libraries.R")

################################################################################
### Loading the simulation data
cdpop_sim1 <- readRDS("cdpop_sim_TEST/Results/iter__1/cdpop_sim.rds")

################################################################################
### Loading the Podarcis muralis data

# Geographical coordinates
geosites <- SpatialPoints(read.table("Data/geo.txt", header=TRUE)[,2:3],
                          CRS(SRS_string = "EPSG:3044"))
empir_geo_dist <- as.matrix(dist(read.table("Data/geo.txt", header=TRUE)[,2:3]))


# Genetic data
lizgen <- read.csv("Data/TR_NA_header.csv", row.names="ID", na.strings="NA")
lizgen[lizgen==95] <- "095"
lizgen[lizgen==97] <- "097"
lizgen[lizgen==99] <- "099"
lizgen.df <- lizgen[,1:(ncol(lizgen)/2)]


for (i in seq(1,ncol(lizgen),2)){
  lizgen.df[ ,(i+1)/2] <- apply( lizgen[ , i:(i+1) ] , 1, paste , collapse = "/" )
}

colnames(lizgen.df) <- paste0("L", 1:ncol(lizgen.df))
lizgen.df[lizgen.df=="NA/NA"] <- NA
lizgen.genind <- df2genind(lizgen.df, NA.char=NA,
                           ploidy=2, type="codom", sep = "/", check.ploidy=TRUE)

empirLoiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(lizgen.genind))

mantel.randtest(as.dist(empir_geo_dist), as.dist(1-empirLoiselle_EcoGenetics))

empir_geo_dist2 <- empir_geo_dist
empir_geo_dist2[empir_geo_dist2==0] <- NA
IBD <- lm(c(as.dist(empirLoiselle_EcoGenetics))~log(c(as.dist(empir_geo_dist2))))
summary(IBD)
plot(log(empir_geo_dist2), empirLoiselle_EcoGenetics)
abline(IBD, col="red")

################################################################################
################################################################################


################################################################################
sim_geo_dist <- as.matrix(dist(cdpop_sim1$grid_list$gen_101@other$xy))
sim_geo_dist[sim_geo_dist==0] <- NA

simLoiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(cdpop_sim1$grid_list$gen_101))

mantel.randtest(as.dist(sim_geo_dist), as.dist(1-simLoiselle_EcoGenetics))
IBDsim <- lm(c(as.dist(simLoiselle_EcoGenetics))~log(c(as.dist(sim_geo_dist))))
summary(IBDsim)
plot(log(sim_geo_dist), simLoiselle_EcoGenetics)
abline(IBDsim, col="red")

################################################################################
# Find the empirical distribution of lizard abundance at our new res

lizpercell <- as.matrix(table(extract(catraster_SA_coarser_cropped, geosites,
                                      cellnumbers = TRUE,
                                      fun=length)[,1]))

lizgrid <- catraster_SA_coarser_cropped
lizgrid[] <- 0
lizcellrowcol <- rowColFromCell(lizgrid, as.numeric(names(lizpercell[,1])))
lizgrid[lizcellrowcol] <- lizpercell

plot(lizgrid)
################################################################################
# Let us try to subsample to get a similar sampling as the empirical dataset
sim_geosites <- SpatialPoints(cdpop_sim1$grid_list$gen_101@other$xy,
                          CRS(SRS_string = "EPSG:3044"))

lizgridno0 <- lizgrid 
lizgridno0[lizgridno0==0] <- NA

#subs <- extract(lizgridno0, sim_geosites, cellnumber-TRUE, sp=TRUE)

library(spatialEco)
lizgridno0poly <- rasterToPolygons(lizgridno0)
plot(lizgridno0poly)
subs <- erase.point(sim_geosites, lizgridno0poly, inside=FALSE)



"%IN%" <- function(x, y) interaction(x) %in% interaction(y)

index <- cdpop_sim1$grid_list$gen_101@other$xy  %IN% as.data.frame(subs)
sim_subs_genind <- cdpop_sim1$grid_list$gen_101[index]
sim_subs_genind

sim_subs_geo_dist <- as.matrix(dist(sim_subs_genind@other$xy))
sim_subs_geo_dist[sim_geo_dist==0] <- NA

sim_subs_Loiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(sim_subs_genind))

mantel.randtest(as.dist(sim_subs_geo_dist), as.dist(1-sim_subs_Loiselle_EcoGenetics))
IBDsim_subs <- lm(c(as.dist(sim_subs_Loiselle_EcoGenetics))~log(c(as.dist(sim_subs_geo_dist))))
summary(IBDsim_subs)
plot(log(sim_subs_geo_dist), sim_subs_Loiselle_EcoGenetics)
abline(IBDsim_subs, col="red")

par(mfrow=c(1,2))

plot(log(sim_subs_geo_dist), sim_subs_Loiselle_EcoGenetics, xlim=c(0,10), ylim=c(-0.3,0.45))
abline(IBDsim_subs, col="red")

plot(log(empir_geo_dist2), empirLoiselle_EcoGenetics, xlim=c(0,10), ylim=c(-0.3,0.45))
abline(IBD, col="red")

saveRDS(sim_subs_genind,"cdpop_sim_TEST/Results/iter__1/sim_subs_genind.rds")
lol <- readRDS("cdpop_sim_TEST/Results/iter__1/sim_subs_genind.rds")
lol
plot(lol@other$xy)