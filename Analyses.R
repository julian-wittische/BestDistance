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
# plot(log(empir_geo_dist2), empirLoiselle_EcoGenetics)
# abline(IBD, col="red")

################################################################################
################################################################################


################################################################################
# sim_geo_dist <- as.matrix(dist(cdpop_sim1$grid_list$gen_101@other$xy))
# sim_geo_dist[sim_geo_dist==0] <- NA
# 
# simLoiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(cdpop_sim1$grid_list$gen_101))
# 
# mantel.randtest(as.dist(sim_geo_dist), as.dist(1-simLoiselle_EcoGenetics))
# IBDsim <- lm(c(as.dist(simLoiselle_EcoGenetics))~log(c(as.dist(sim_geo_dist))))
# summary(IBDsim)
# plot(log(sim_geo_dist), simLoiselle_EcoGenetics)
# abline(IBDsim, col="red")

################################################################################
# Find the empirical distribution of lizard abundance at our new res

lizpercell <- as.matrix(table(extract(catraster_SA_coarser_cropped, geosites,
                                      cellnumbers = TRUE,
                                      fun=length)[,1]))

lizgrid <- catraster_SA_coarser_cropped
lizgrid[] <- 0
lizcellrowcol <- rowColFromCell(lizgrid, as.numeric(names(lizpercell[,1])))
lizgrid[lizcellrowcol] <- lizpercell

# plot(lizgrid)
################################################################################
# Let us try to subsample to get a similar sampling as the empirical dataset
sim_geosites <- SpatialPoints(cdpop_sim1$grid_list$gen_101@other$xy,
                          CRS(SRS_string = "EPSG:3044"))

lizgridno0 <- lizgrid 
lizgridno0[lizgridno0==0] <- NA

#subs <- extract(lizgridno0, sim_geosites, cellnumber=TRUE, sp=TRUE)

library(spatialEco)
lizgridno0poly <- rasterToPolygons(lizgridno0)
# plot(lizgridno0poly)
subs <- erase.point(sim_geosites, lizgridno0poly, inside=FALSE)

"%IN%" <- function(x, y) interaction(x) %in% interaction(y)

index <- cdpop_sim1$grid_list$gen_101@other$xy  %IN% as.data.frame(subs)
sim_subs_genind <- cdpop_sim1$grid_list$gen_101[index]
# This is to have exactly the same number of individuals as in the empirical dataset
sim_subs_genind <- sim_subs_genind[sample(1:nrow(sim_subs_genind@tab),
                                          nrow(lizgen.genind@tab),
                                          replace = FALSE),]

#sim_subs_genind <- readRDS("sim_subs_genind.rds")

sim_subs_geo_dist <- as.matrix(dist(sim_subs_genind@other$xy))
sim_subs_geo_dist[sim_subs_geo_dist==0] <- NA

sim_subs_Loiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(sim_subs_genind))

mantel.randtest(as.dist(sim_subs_geo_dist), as.dist(1-sim_subs_Loiselle_EcoGenetics))
IBDsim_subs <- lm(c(as.dist(sim_subs_Loiselle_EcoGenetics))~log(c(as.dist(sim_subs_geo_dist))))
summary(IBDsim_subs)
# plot(log(sim_subs_geo_dist), sim_subs_Loiselle_EcoGenetics)
# abline(IBDsim_subs, col="red")

plot(lizgrid)
points(sim_subs_genind@other$xy)

par(mfrow=c(1,2))

plot(log(sim_subs_geo_dist), sim_subs_Loiselle_EcoGenetics, xlim=c(0,10), ylim=c(-0.3,0.45))
abline(IBDsim_subs, col="red")

plot(log(empir_geo_dist2), empirLoiselle_EcoGenetics, xlim=c(0,10), ylim=c(-0.3,0.45))
abline(IBD, col="red")

sim_subs_genind

# saveRDS(sim_subs_genind,"sim_subs_genind_3.2.15.30.1RES_LOWIBD3.rds")
# 
get_slope <- function(sim_genind_object){
  sim_genind_object_geo_dist <- as.matrix(dist(sim_genind_object@other$xy))
  sim_genind_object_geo_dist[sim_genind_object_geo_dist==0] <- NA
  sim_genind_object_Loiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(sim_genind_object))
  IBDsim <- lm(c(as.dist(sim_genind_object_Loiselle_EcoGenetics))~log(c(as.dist(sim_genind_object_geo_dist))))
  return(summary(IBDsim))
}


AL1 <- readRDS("sim_subs_genind_3.2.15.30.1RES_LOWIBD1.rds")
AL2 <- readRDS("sim_subs_genind_3.2.15.30.1RES_LOWIBD2.rds")
AL3 <- readRDS("sim_subs_genind_3.2.15.30.1RES_LOWIBD3.rds")
AL1
AL2
AL3
get_slope(AL1)
get_slope(AL2)
get_slope(AL3)

AM1 <- readRDS("sim_subs_genind_3.2.15.30.1RES_MODIBD1.rds")
AM2 <- readRDS("sim_subs_genind_3.2.15.30.1RES_MODIBD2.rds")
AM3 <- readRDS("sim_subs_genind_3.2.15.30.1RES_MODIBD3.rds")
AM1
AM2
AM3
get_slope(AM1)
get_slope(AM2)
get_slope(AM3)

AH1 <- readRDS("sim_subs_genind_3.2.15.30.1RES_HIGHIBD1.rds")
AH2 <- readRDS("sim_subs_genind_3.2.15.30.1RES_HIGHIBD2.rds")
AH3 <- readRDS("sim_subs_genind_3.2.15.30.1RES_HIGHIBD3.rds")
AH1
AH2
AH3
get_slope(AH1)
get_slope(AH2)
get_slope(AH3)

TL1 <- readRDS("sim_subs_genind_1.1.15.30.1RES_LOWIBD1.rds")
TL2 <- readRDS("sim_subs_genind_1.1.15.30.1RES_LOWIBD2.rds")
TL3 <- readRDS("sim_subs_genind_1.1.15.30.1RES_LOWIBD3.rds")
TL1
TL2
TL3
get_slope(TL1)
get_slope(TL2)
get_slope(TL3)

TM1 <- readRDS("sim_subs_genind_1.1.15.30.1RES_MODIBD1.rds")
TM2 <- readRDS("sim_subs_genind_1.1.15.30.1RES_MODIBD2.rds")
TM3 <- readRDS("sim_subs_genind_1.1.15.30.1RES_MODIBD3.rds")
TM1
TM2
TM3
get_slope(TM1)
get_slope(TM2)
get_slope(TM3)

TH1 <- readRDS("sim_subs_genind_1.1.15.30.1RES_HIGHIBD1.rds")
TH2 <- readRDS("sim_subs_genind_1.1.15.30.1RES_HIGHIBD2.rds")
TH3 <- readRDS("sim_subs_genind_1.1.15.30.1RES_HIGHIBD3.rds")
TH1
TH2
TH3
get_slope(TH1)
get_slope(TH2)
get_slope(TH3)