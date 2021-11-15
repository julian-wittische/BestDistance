################################################################################
########## Julian Wittische - November 2021 - Simulating connectivity ##########
################################################################################

# This is based on previous work by William Peterman and Kristopher Winiarski
# Simulations based on empirical data



empir.sim <- function(catraster,
                      geosites,
                       RMexpvar_r = 1,
                       RMexpscale_r = 15,
                       habitat = 0.5,
                       n_ind = 15000,
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
                       matemovethresh = 1,
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
    coo <- dim(orig)
    bb <- extent(orig)
    model <- RMexp(var=RMexpvar_r, scale=RMexpscale_r)
    rf.sim <- RFsimulate(model = model, x=1:115, y=1:128, grid=TRUE)
    rand <- raster(scale(as.matrix(rf.sim)))
    rand <- setExtent(rand, bb, snap= TRUE)
    random_1 <- rand

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
    pts <- na.omit(unique(floor(cbind(runif(50000, 328302.5,333422.5), 
                              runif(50000, 5512494 , 5517094 )))))
    
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
    
    cdpop_sim <- cdpop(CDPOP.py = "C:/Users/jwittische/Desktop/Projects/BestDistance/CDPOP-master/src/CDPOP.py",
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
                       K_env = 30000,
                       matemoveno = matemoveno, ## 1 = Linear, 5 = Neg exp; 9 = custom prob matrix
                       matemovethresh = m_thresh,
                       MeanFecundity = MeanFecundity)
    
    
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
    r.stack <- stack(orig,
                     random_1)
    names(r.stack) <- c('cat', 'rand')
    
    GA.inputs <- GA.prep(ASCII.dir = r.stack,
                           Results.dir = 'all_comb',
                           max.cont = 2500,
                           select.trans = list('A', 'A'),
                           maxiter = maxiter, 
                           run = 25,
                            parallel = parallel)
    
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
    
  } #  end iteration loop (z)
} # end function
