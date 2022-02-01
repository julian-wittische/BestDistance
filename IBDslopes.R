get_slope <- function(sim_genind_object){
  sim_genind_object_geo_dist <- as.matrix(dist(sim_genind_object@other$xy))
  sim_genind_object_geo_dist[sim_genind_object_geo_dist==0] <- NA
  sim_genind_object_Loiselle_EcoGenetics <- eco.kin.loiselle(genind2ecogen(sim_genind_object))
  IBDsim <- lm(c(as.dist(sim_genind_object_Loiselle_EcoGenetics))~log(c(as.dist(sim_genind_object_geo_dist))))
  return(summary(IBDsim))
}

highlow <- readRDS("sim_subs_genind_HIGHRES_LOWIBD.rds")
highlow2 <- readRDS("sim_subs_genind_HIGHRES_LOWIBD2.rds")
highlow3 <- readRDS("sim_subs_genind_HIGHRES_LOWIBD3.rds")
get_slope(highlow)
get_slope(highlow2)
get_slope(highlow3)

highhigh <- readRDS("sim_subs_genind_HIGHRES_HIGHIBD.rds")
highhigh2 <- readRDS("sim_subs_genind_HIGHRES_HIGHIBD2.rds")
highhigh3 <- readRDS("sim_subs_genind_HIGHRES_HIGHIBD3.rds")
get_slope(highhigh)
get_slope(highhigh2)
get_slope(highhigh3)

lowlow <- readRDS("sim_subs_genind_LOWRES_LOWIBD.rds")
lowlow2 <- readRDS("sim_subs_genind_LOWRES_LOWIBD2.rds")
lowlow3 <- readRDS("sim_subs_genind_LOWRES_LOWIBD3.rds")
get_slope(lowlow)
get_slope(lowlow2)
get_slope(lowlow3)

lowhigh <- readRDS("sim_subs_genind_LOWRES_HIGHIBD.rds")
lowhigh2 <- readRDS("sim_subs_genind_LOWRES_HIGHIBD2.rds")
lowhigh3 <- readRDS("sim_subs_genind_LOWRES_HIGHIBD3.rds")
get_slope(lowhigh)
get_slope(lowhigh2)
get_slope(lowhigh3)