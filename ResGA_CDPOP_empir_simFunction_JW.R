################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

# This is based on previous work by William Peterman and Kristopher Winiarski
# Simulations based on empirical data

empir.sim <- function(catraster = catraster,
                      geosites = geosites,
                      habitat = 0.5,
                       RMexpvar_r = 1,
                       RMexpscale_r = 15,
                       n_ind = 10000,
                       n_samplepoints = 250,
                       start = 1,
                       seed = 1,
                       iters = 1,
                       parallel =3,
                       method = 'standard',
                       maxiter = 100,                               
                       JULIA_HOME = "C:/Users/jwittische/AppData/Local/Programs/Julia-1.6.3/bin/",
                       CDPOP.py = 'C:/Users/jwittische/Desktop/Projects/BestDistance/CDPOP-master/src/CDPOP.py',
                       sim_name = 'output_',
                       sim_dir = "C:/Users/jwittische/Desktop/Projects/BestDistance/cdpop_sim_TEST/",
                       looptime = 101,
                       output_years = 100,
                       gridformat = 'cdpop',
                       loci = 17,
                       alleles = 20,
                       matemoveno = 1, ## 1 = Linear, 2 = Inv sq; 9 = custom prob matrix
                       matemovethresh = 0.1,
                       MeanFecundity = 4,
                       n_axes = 64)
  {

  if(is.null(start)) {
    stop("Specify integer `start` value!!!")
  }

  # Main Function -----------------------------------------------------------
  for(z in start:iters){

    # >> Make Random Surfaces -------------------------------------------------
    
    #  * Empirical surface ------------------------------------------------------------

    # Load ordinal variable template
    orig <- catraster #not the same as Copernicus! (3035)
    
    #  * Random surface --------------------------------------------------------
    # coo <- dim(orig)
    # bb <- extent(orig)
    # model <- RMexp(var=RMexpvar_r, scale=RMexpscale_r)
    # rf.sim <- RFsimulate(model = model, x=1:nrow(orig), y=1:ncol(orig), grid=TRUE)
    # rand <- raster(scale(as.matrix(rf.sim)))
    # rand <- setExtent(rand, bb, snap= TRUE)
    # random_1 <- rand

    # >> Create Directory -----------------------------------------------------

    dir.create(paste0(sim_dir, "Results/",
                      'iter__', z),
               recursive = TRUE) 
    
    out <- paste0(sim_dir, "Results/",
                  'iter__', z, "/")

    # >> Create Truth ---------------------------------------------------------
    Resist <- catraster 
    Resist[Resist==0] <- 40
    Resist[Resist==2] <- 400
    Resist[Resist==3] <- 400
    Resist[Resist==4] <- 40
    # Load sampling sites ------------------------------------------------------
    pts <- unique(floor(cbind(runif(10000, extent(catraster)[1], extent(catraster)[2]), 
                              runif(10000, extent(catraster)[3], extent(catraster)[4]))))
    
    sample.thresh <- as.numeric(quantile(Resist, habitat))
    
    sample.extract <- extract(Resist, pts)
    sample.suit <- pts[sample.extract <= sample.thresh,]
    sample.suit <- na.omit(sample.suit)
    pts <- SpatialPoints(sample.suit[sample(nrow(sample.suit), n_ind, replace = F),])
    
    ### 
    
    # >> Calculate cost distance -----------------------------------------------
    
    jl.inputs <- jl.prep(n.Pops = length(pts),
                         CS_Point.File = pts,
                         JULIA_HOME = JULIA_HOME)
    
    ##  Pairwise resistance distances
    r.dist <- Run_CS.jl(jl.inputs = jl.inputs,
                        r = Resist,
                        full.mat = T)
    
    plot(Resist)
    plot(pts, add = T, pch = 19)
    
    m_thresh <- quantile(lower(r.dist), matemovethresh)
    
    cdpop_sim <- cdpopJW(CDPOP.py = "C:/Users/jwittische/Desktop/Projects/BestDistance/CDPOP-master/src/CDPOP.py",
                       sim_name = sim_name,
                       pts = pts,
                       resist_rast = Resist,
                       resist_mat = r.dist,
                       sim_dir = out,
                       looptime = looptime,
                       output_years = output_years,
                       gridformat = gridformat,
                       loci = loci,
                       alleles = alleles,
                       K_env = 10000,
                       matemoveno = matemoveno, ## 1 = Linear, 5 = Neg exp; 9 = custom prob matrix
                       matemovethresh = m_thresh,
                       MeanFecundity = MeanFecundity)
    
    saveRDS(cdpop_sim, paste0(out,"cdpop_sim.rds"))
    
    writeRaster(Resist,
                paste0(out, "true_resist.asc"), overwrite=TRUE)
  } #  end iteration loop (z)
} # end function
