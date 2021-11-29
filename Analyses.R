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
# Find the empirical distribution of lizard abundance at our new res

lizpercell <- as.matrix(table(extract(catraster_SA_coarser_cropped, geosites,
                                      cellnumbers = TRUE,
                                      fun=length)[,1]))

lizgrid <- catraster_SA_coarser_cropped
lizgrid[] <- 0
lizcellrowcol <- rowColFromCell(lizgrid, as.numeric(names(lizpercell[,1])))
lizgrid[lizcellrowcol] <- lizpercell

plot(lizgrid)

# ################################################################################
sim_geo_dist <- as.matrix(dist(cdpop_sim1$grid_list$gen_101@other$xy))
sim_geo_dist[sim_geo_dist==0] <- NA


# a_tab <- adegenet::tab(cdpop_sim1$grid_list$gen_101)
# pc <- prcomp(a_tab)
# pc_dist <- as.matrix(dist(pc$x[,1:10]))
# pc_dist

# RClonedf <- genind2df(cdpop_sim1$grid_list$gen_101)
# RCloneobj <- convert_GC(RClonedf, 2, "")
# simLoiselle_RClone <- kinship_Loiselle_core(RCloneobj)
simLoiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(cdpop_sim1$grid_list$gen_101))

mantel.randtest(dist(pc$x[,1:10]), dist(cdpop_sim1$grid_list$gen_101@other$xy))

plot(log(geo_dist), simLoiselle_EcoGenetics)
