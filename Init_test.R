library(corMLPE)
library(ResistanceGA)
library(radish)
library(RandomFields)

r.dim = 100
RMexpvar_1 = 1
RMexpscale_1 = 1
RMexpvar_2 = 1
RMexpscale_2 = 2
transformation_1 = 3
r.tran_shape_1 = 3
r.tran_max_1 = 100
transformation_2 = 8
r.tran_shape_2 = 3
r.tran_max_2 = 100
RMexpvar_r = 1
RMexpscale_r = 15
RMexpvar_r2 = 1
RMexpscale_r2 = 10
habitat = 0.5
n_ind = 500
n_samplepoints = 100
start = 10
seed = 1
iters = 1
parallel = 3
method = 'standard'
maxiter = 100                          
JULIA_HOME = "C:/Users/jwittische/AppData/Local/Programs/Julia-1.6.3/bin/"
CDPOP.py = 'C:/Users/jwittische/Desktop/Projects/BestDistance/CDPOP-master/src/CDPOP.py'
sim_name = 'output_'
sim_dir = "C:/Users/jwittische/Desktop/Projects/BestDistance/cdpop_sim_TEST"
looptime = 101
output_years = 100
gridformat = 'cdpop'
loci = 30
alleles = 30
matemoveno = 2 ## 1 = Linear, 2 = Inv sq; 9 = custom prob matrix
matemovethresh = 0.025
MeanFecundity = 4
n_axes = 64

resist_mat = NULL
agefilename = NULL
mcruns = 1
looptime = 400
output_years = 50
gridformat = 'cdpop'
cdclimgentime = 0
matemoveno = 2
matemoveparA = 0
matemoveparB = 0
matemoveparC = 0
matemovethresh = 'max'
output_matedistance = 'N'
sexans = 'Y'
Freplace = 'N'
Mreplace = 'N'
philopatry = 'N'
multiple_paternity = 'N'
selfans = 'N'
Fdispmoveno = NULL
FdispmoveparA = 0
FdispmoveparB = 0
FdispmoveparC = 0
Fdispmovethresh = NULL
Mdispmoveno = NULL
MdispmoveparA = 0
MdispmoveparB = 0
MdispmoveparC = 0
Mdispmovethresh = NULL
offno = 2
MeanFecundity = 5
Femalepercent = 50
EqualsexratioBirth = 'N'
TwinningPercent = 0
popModel = 'exp'
r = 1

subpopmortperc = 0
muterate = 0.0005
mutationtype = 'random'
loci = 1000
intgenesans = 'random'
allefreqfilename = 'N'
alleles = 2
mtdna = 'N'
startGenes = 0
cdevolveans = 'N'
startSelection = 0
betaFile_selection = 'N'
epistasis = 'N'
epigeneans = 'N'
startEpigene = 0
betaFile_epigene = 'N'
cdinfect = 'N'
transmissionprob = 0

surface_corrs <- function(A, orig, rep) {
  rep <- rep * A
  rep <- orig + rep
}

# >> Make Random Surfaces -------------------------------------------------

# Specify random model

#  * Surface 1 ------------------------------------------------------------

# Create original and correlated surfaces for surface 1.

z <- 1 #####

set.seed(z)
corr = c(0.9, 0.5, 0.1)
A <- sqrt((1/(corr^2)) - 1)  # calcualate A
coo <- 1:r.dim
model <- RMexp(var=RMexpvar_1, scale=RMexpscale_1)
rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
surface_orig <- as.matrix(rf.sim)

rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
surface_rep <- as.matrix(rf.sim)

corr_surfaces <- lapply(A, surface_corrs, orig = surface_orig, rep = surface_rep)

out_rasts <- lapply(c(list(surface_orig), corr_surfaces), function(r) raster(scale(r)))
names(out_rasts) <- c("orig", paste0("cor_", corr))

bb <- extent(0.5, r.dim + 0.5, 0.5, r.dim + 0.5)

extent(out_rasts$orig) <- bb
out_rasts$orig <- setExtent(out_rasts$orig, bb, keepres=TRUE)
orig <- out_rasts$orig

extent(out_rasts$cor_0.9) <- bb
out_rasts$cor_0.9 <- setExtent(out_rasts$cor_0.9, bb, keepres=TRUE)
cor_0.9 <- out_rasts$cor_0.9

extent(out_rasts$cor_0.5) <- bb
out_rasts$cor_0.5 <- setExtent(out_rasts$cor_0.5, bb, keepres=TRUE)
cor_0.5 <- out_rasts$cor_0.5

extent(out_rasts$cor_0.1) <- bb
out_rasts$cor_0.1 <- setExtent(out_rasts$cor_0.1, bb, keepres=TRUE)
cor_0.1 <- out_rasts$cor_0.1

orig_1 <- orig
cor_0.9_1 <- cor_0.9
cor_0.5_1 <- cor_0.5
cor_0.1_1 <- cor_0.1

#  * Surface 2 ------------------------------------------------------------

# Create original and correlated surfaces for surface 2.
corr = c(0.9, 0.5, 0.1)
A <- sqrt((1/(corr^2)) - 1)  # calcualate A
coo <- 1:r.dim
model <- RMexp(var=RMexpvar_2, scale=RMexpscale_2)
rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
surface_orig <- as.matrix(rf.sim)

rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
surface_rep <- as.matrix(rf.sim)


corr_surfaces <- lapply(A, surface_corrs, orig = surface_orig, rep = surface_rep)

out_rasts <- lapply(c(list(surface_orig), corr_surfaces), function(r) raster(scale(r)))
names(out_rasts) <- c("orig", paste0("cor_", corr))

extent(out_rasts$orig) <- bb
out_rasts$orig <- setExtent(out_rasts$orig, bb, keepres=TRUE)
orig <- out_rasts$orig

extent(out_rasts$cor_0.9) <- bb
out_rasts$cor_0.9 <- setExtent(out_rasts$cor_0.9, bb, keepres=TRUE)
cor_0.9 <- out_rasts$cor_0.9

extent(out_rasts$cor_0.5) <- bb
out_rasts$cor_0.5 <- setExtent(out_rasts$cor_0.5, bb, keepres=TRUE)
cor_0.5 <- out_rasts$cor_0.5

extent(out_rasts$cor_0.1) <- bb
out_rasts$cor_0.1 <- setExtent(out_rasts$cor_0.1, bb, keepres=TRUE)
cor_0.1 <- out_rasts$cor_0.1

#Create original and correlated surfaces for surface 2.
orig_2 <- orig
cor_0.9_2 <- cor_0.9
cor_0.5_2 <- cor_0.5
cor_0.1_2 <- cor_0.1

#  * Random surface ------------------------------------------------------------
corr = c(0.9, 0.5, 0.1)
A <- sqrt((1/(corr^2)) - 1)  # calcualate A
coo <- 1:r.dim
model <- RMexp(var=RMexpvar_r, scale=RMexpscale_r)
rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
surface_orig <- as.matrix(rf.sim)

rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
surface_rep <- as.matrix(rf.sim)

corr_surfaces <- lapply(A, surface_corrs, orig = surface_orig, rep = surface_rep)

out_rasts <- lapply(c(list(surface_orig), corr_surfaces), function(r) raster(scale(r)))
names(out_rasts) <- c("orig", paste0("cor_", corr))

extent(out_rasts$orig) <- bb
out_rasts$orig <- setExtent(out_rasts$orig, bb, keepres=TRUE)
orig <- out_rasts$orig

extent(out_rasts$cor_0.9) <- bb
out_rasts$cor_0.9 <- setExtent(out_rasts$cor_0.9, bb, keepres=TRUE)
cor_0.9 <- out_rasts$cor_0.9

extent(out_rasts$cor_0.5) <- bb
out_rasts$cor_0.5 <- setExtent(out_rasts$cor_0.5, bb, keepres=TRUE)
cor_0.5 <- out_rasts$cor_0.5

extent(out_rasts$cor_0.1) <- bb
out_rasts$cor_0.1 <- setExtent(out_rasts$cor_0.1, bb, keepres=TRUE)
cor_0.1 <- out_rasts$cor_0.1

#Create original and correlated surfaces for surface 2.
random_1 <- orig
cor_0.9_3 <- cor_0.9
cor_0.5_3 <- cor_0.5
cor_0.1_3 <- cor_0.1


# >> Create Directory -----------------------------------------------------
# SD <- ifelse(i == 0.5, '0_5', '1_25')
dir.create(paste0(sim_dir, "Results/",
                  'iter__', z),
           recursive = TRUE) 

out <- paste0(sim_dir, "Results/",
              'iter__', z, "/")
tmp <- paste0(tempdir(), '\\')


# >> Create Truth ---------------------------------------------------------
r.stack <- stack(orig_1,
                 orig_2,
                 random_1)
names(r.stack) <- c('cont1', 'cont2', 'rand')


# #Run 'GA.prep'
GA.inputs <- GA.prep(ASCII.dir = r.stack[[-3]],
                     select.trans = list('A','A'),
                     max.cat = 500,
                     max.cont = 500,
                     Results.dir = tmp,
                     quiet = TRUE)


#  * Julia Prep -----------------------------------------------------------
x.coord <- SpatialPoints(cbind(runif(n_samplepoints, 10, 40),
                               runif(n_samplepoints, 10, 40)))
jl.inputs <- jl.prep(n.Pops = n_samplepoints,
                     CS_Point.File = x.coord,
                     JULIA_HOME = JULIA_HOME)

PARM <- c(transformation_1, #Transformation for continuous surface 1
          r.tran_shape_1, #shape parameter for continuous surface 1
          r.tran_max_1, #max parameter for continuous surface 1
          transformation_2, #Transformation for continuous surface 1
          r.tran_shape_2, #shape parameter for continuous surface 1
          r.tran_max_2) #max parameter for continuous surface 2


# Combine resistance surfaces
Resist <- Combine_Surfaces(PARM = PARM,
                           jl.inputs = jl.inputs,
                           GA.inputs = GA.inputs,
                           out = NULL,
                           rescale = TRUE)#,
                           #p.contribution = TRUE)

print(Combine_Surfaces(PARM = PARM,
                       jl.inputs = jl.inputs,
                       GA.inputs = GA.inputs,
                       out = NULL,
                       rescale = TRUE,
                       p.contribution = TRUE))
plot(stack(r.stack,Resist))

# Generate Points ---------------------------------------------------------
pts <- unique(floor(cbind(runif(5000, 0.15 * r.dim, 0.85 * r.dim), 
                          runif(5000, 0.15 * r.dim, 0.85 * r.dim))))

sample.thresh <- as.numeric(quantile(Resist, habitat))

sample.extract <- extract(Resist, pts)
sample.suit <- pts[sample.extract <= sample.thresh,]
pts <- SpatialPoints(sample.suit[sample(nrow(sample.suit), n_ind, replace = F),])
K_env = length(pts)
### 


# >> Calculate cost distance -------------------------------------------------

jl.inputs <- jl.prep(n.Pops = length(pts),
                     CS_Point.File = pts,
                     JULIA_HOME = JULIA_HOME)

##  Pairwise resistance distances
r.dist <- Run_CS.jl(jl.inputs = jl.inputs,
                    r = Resist,
                    full.mat = T
)

plot(Resist)
plot(pts, add = T, pch = 19)

m_thresh <- quantile(lower(r.dist), matemovethresh)

################################################################################
if(!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
suppressWarnings(
  dir.create(paste0(sim_dir,"/data/"), recursive = TRUE)
)
data_dir <- paste0(sim_dir,"/data/")

if(matemoveno == 9){
  if(class(resist_rast) == "RasterLayer"){
    stop('Specify a probability matrix instead of a raster layer!')
  }
  write.table(resist_rast,
              paste0(data_dir, "move_prob.csv"), 
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  
  cdmat <- 'move_prob'
  # write.table(matemoveno, paste0(data_dir, 'DispProb.csv'),
  #             sep = ",",
  #             row.names = F)
  # matemoveno <- 'DispProb'  
}

if(is.null(Fdispmoveno)){
  Fdispmoveno <- matemoveno
}

if(is.null(Mdispmoveno)){
  Mdispmoveno <- matemoveno
}

if(is.null(Fdispmovethresh)){
  Fdispmovethresh <- matemovethresh
}

if(is.null(Mdispmovethresh)){
  Mdispmovethresh <- matemovethresh
}

if(is.null(agefilename)){
  age_df <- data.frame(`Age class` = c(0,1),
                       Distribution = c(0,1),
                       `Male Mortality` = c(0,100),
                       `Female Mortality` = c(0,100),
                       `Mean Fecundity` = c(0,MeanFecundity),
                       `Std Fecundity` = c(0,0),
                       `Male Maturation` = c(0,1),
                       `Female Maturation` = c(0,1),
                       check.names = F)
  write.table(age_df, paste0(data_dir, 'AgeVars.csv'),
              sep = ",",
              row.names = F)
  # age_file <- paste0(data_dir, 'AgeVars.csv')
  age_file <- 'AgeVars.csv'
  
}

xyFile_df <- data.frame(Subpopulation = rep(1, length(pts)),
                        XCOORD = pts@coords[,1],
                        YCOORD = pts@coords[,2],
                        ID = paste0('initial',1:length(pts) - 1),
                        sex = sample(c(0,1), 
                                     replace = T,
                                     size = length(pts)),
                        Fitness_AA = rep(0, length(pts)),
                        Fitness_Aa = rep(0, length(pts)),
                        Fitness_aa = rep(0, length(pts)),
                        Fitness_AABB = rep(0, length(pts)),
                        Fitness_AaBB = rep(0, length(pts)),
                        Fitness_aaBB = rep(0, length(pts)),
                        Fitness_AABb = rep(0, length(pts)),
                        Fitness_AaBb = rep(0, length(pts)),
                        Fitness_aaBb = rep(0, length(pts)),
                        Fitness_AAbb = rep(0, length(pts)),
                        Fitness_Aabb = rep(0, length(pts)),
                        Fitness_aabb = rep(0, length(pts))  
)

write.table(xyFile_df, paste0(data_dir, 'xyFile.csv'),
            sep = ',',
            row.names = F)
# xyFile <- paste0(data_dir, 'xyFile')
xyFile <- 'xyFile'

resist_rast = Resist
resist_mat = r.dist
sim_dir = out

if(matemoveno == 9){
  write.table(resist_mat,
              paste0(data_dir, "resist_mat.csv"), 
              sep = ",",
              row.names = FALSE,
              col.names = FALSE)
  
  cdmat <- 'resist_mat'
}

if(matemoveno != 9){
  if(!is.null(resist_mat)){
    write.table(resist_mat,
                paste0(data_dir, "resist_mat.csv"), 
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)
    
    cdmat <- 'resist_mat'
    
  } else {
    print("Calculating resistance distance with `gdistance`...")
    
    trans <- transition(x = resist_rast,
                        transitionFunction = function(x)  1 / mean(x),
                        directions = 8)
    
    trR <- geoCorrection(trans, "r", scl = T)
    resist_mat <- as.matrix(commuteDistance(trR, pts) / 1000)  
    
    ## Check file format, row/col names?
    write.table(resist_mat,
                paste0(data_dir, "resist_mat.csv"), 
                sep = ",",
                row.names = FALSE,
                col.names = FALSE)
    
    cdmat <- 'resist_mat'
  }  
}

cdpop_df <- data.frame(xyfilename = xyFile,
                       agefilename = age_file,
                       mcruns = mcruns,
                       looptime = looptime,
                       output_years = output_years,
                       gridformat = gridformat,
                       cdclimgentime = cdclimgentime,
                       matecdmat = cdmat,
                       dispcdmat = cdmat,
                       matemoveno = matemoveno,
                       matemoveparA = matemoveparA,
                       matemoveparB = matemoveparB,
                       matemoveparC = matemoveparC,
                       matemovethresh = matemovethresh,
                       output_matedistance = output_matedistance,
                       sexans = sexans,
                       Freplace = Freplace,
                       Mreplace = Mreplace,
                       philopatry = philopatry,
                       multiple_paternity = multiple_paternity,
                       selfans = selfans,
                       Fdispmoveno = Fdispmoveno,
                       FdispmoveparA = FdispmoveparA,
                       FdispmoveparB = FdispmoveparB,
                       FdispmoveparC = FdispmoveparC,
                       Fdispmovethresh = Fdispmovethresh,
                       Mdispmoveno = Mdispmoveno,
                       MdispmoveparA = MdispmoveparA,
                       MdispmoveparB = MdispmoveparB,
                       MdispmoveparC = MdispmoveparC,
                       Mdispmovethresh = Mdispmovethresh,
                       offno = offno,
                       Femalepercent = Femalepercent,
                       EqualsexratioBirth = EqualsexratioBirth,
                       TwinningPercent = TwinningPercent,
                       popModel = popModel,
                       r = r,
                       K_env = K_env,
                       subpopmortperc = subpopmortperc,
                       muterate = muterate,
                       mutationtype = mutationtype,
                       loci = loci,
                       intgenesans = intgenesans,
                       allefreqfilename = allefreqfilename,
                       alleles = alleles,
                       mtdna = mtdna,
                       startGenes = startGenes,
                       cdevolveans = cdevolveans,
                       startSelection = startSelection,
                       betaFile_selection = betaFile_selection,
                       epistasis = epistasis,
                       epigeneans = epigeneans,
                       startEpigene = startEpigene,
                       betaFile_epigene = betaFile_epigene,
                       cdinfect = cdinfect,
                       transmissionprob = transmissionprob,
                       check.names = F)

write.table(cdpop_df,
            paste0(data_dir, "CDPOP_inputs.csv"), 
            sep = ",",
            row.names = FALSE,
            col.names = TRUE,
            quote = F)

print("Running CDPOP...")

system(paste("python", CDPOP.py, data_dir, "CDPOP_inputs.csv", sim_name))

if(!is.null(python)){
  system(paste(python, CDPOP.py, data_dir, "CDPOP_inputs.csv", sim_name))
} else {
  system(paste("python", CDPOP.py, data_dir, "CDPOP_inputs.csv", sim_name))
}


################################################################################
cdpop_grid <- cdpop_sim$grid_list[[length(cdpop_sim$pop_list)]]
pops <- cdpop_sim$pop_list[[length(cdpop_sim$pop_list)]]

# Subsample ---------------------------------------------------------------

ind_samp <- sort(sample(1:nInd(cdpop_grid), n_samplepoints, replace = F))
pops_ <- pops[ind_samp]

s_pts <- pts[pops_]

s_pops <- data.frame(pop = pops_,
                     x = pts@coords[pops_,1],
                     y = pts@coords[pops_,2])

write.table(s_pops,
            paste0(out, "sampled_pops.csv"),
            sep = ",",
            row.names = FALSE,
            col.names = TRUE)

# PCA dist ---------------------------------------------------------------
Dps <- 1-propShared(cdpop_grid[ind_samp])
pca <- pca_dist(cdpop_grid[ind_samp], n_axes = n_axes)

## Quick plot select
par(mfrow = c(2,2))
plot(lower(Dps) ~ lower(r.dist[pops_,pops_]), main = "Dps ~ Resist")
plot(lower(pca) ~ lower(r.dist[pops_,pops_]), main = "PCA ~ Resist")
plot(lower(Dps) ~ c(dist(pts@coords[pops_,])), main = "Dps ~ Euc dist")
plot(lower(pca) ~ c(dist(pts@coords[pops_,])), main = "PCA ~ Euc dist")
par(mfrow = c(1,1))

graphics.off()


# Optimize: All Comb ----------------------------------------------------------------
if(method == 'island') {
  GA.inputs <- GA.prep(ASCII.dir = r.stack,
                       Results.dir = 'all_comb',
                       max.cont = 500,
                       select.trans = list('A', 'A', 4),
                       maxiter = 52,
                       run = 8,
                       island.pop = island.pop,
                       numIslands = 4,
                       migrationInterval = 4,
                       gaisl = T
  )
} else {
  GA.inputs <- GA.prep(ASCII.dir = r.stack,
                       Results.dir = 'all_comb',
                       max.cont = 500,
                       select.trans = list('A', 'A', 'A'),
                       maxiter = maxiter, #65
                       run = 15,
                       # pop.size = 5,
                       parallel = parallel
  )
}


jl.inputs <- jl.prep(n.Pops = n_samplepoints,
                     CS_Point.File = s_pts,
                     response = lower(pca),
                     JULIA_HOME = JULIA_HOME)

dir.create(paste0(out,'rga/'), recursive = T)
Rga_out <- all_comb(jl.inputs = jl.inputs,
                    GA.inputs = GA.inputs,
                    # iters = 10,
                    results.dir = paste0(out,'rga/')
)

write.csv(pca,
          file = paste0(out, "pca_dist.csv"))

writeRaster(Resist,
            paste0(out, "true_resist.asc"))

write.csv(r.dist[pops_,pops_],
          file = paste0(out, "true_ResistDist.csv"))