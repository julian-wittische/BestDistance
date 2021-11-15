#ResistanceGA_sim function

ResistanceGAMS.sim <- function(r.dim = 100,
                               RMexpvar_1 = 1,
                               RMexpscale_1 = 1,
                               RMexpvar_2 = 1,
                               RMexpscale_2 = 2,
                               transformation_1 = 3,
                               r.tran_shape_1 = 3,
                               r.tran_max_1 = 100,
                               transformation_2 = 8,
                               r.tran_shape_2 = 3,
                               r.tran_max_2 = 100,
                               RMexpvar_r = 1,
                               RMexpscale_r = 15,
                               RMexpvar_r2 = 1,
                               RMexpscale_r2 = 10,
                               habitat = 0.5,
                               n_ind = 500,
                               n_samplepoints = 100,
                               start = NULL,
                               seed = NULL,
                               iters = 25,
                               parallel = 24,
                               method = 'standard',
                               maxiter = 100,                               
                               JULIA_HOME = "C:/Users/peterman.73/AppData/Local/Julia-1.3.1/bin/",
                               CDPOP.py = 'C:/Users/peterman.73/Box/R/CDPOP/CDPOP-master/src/CDPOP.py',
                               sim_name = 'output_',
                               sim_dir = "C:/cdpop_sim_TEST/",
                               looptime = 101,
                               output_years = 100,
                               gridformat = 'cdpop',
                               loci = 30,
                               alleles = 30,
                               matemoveno = 2, ## 1 = Linear, 2 = Inv sq; 9 = custom prob matrix
                               matemovethresh = 0.025,
                               MeanFecundity = 4,
                               n_axes = 64
){
  library(ResistanceGA)
  library(RandomFields)
  # source("U:/CDPOP/cdpop_function.R")
  # library(PopGenReport)
  
  if(is.null(start)) {
    stop("Specify integer `start` value!!!")
  }
  # Helper Function ---------------------------------------------------------
  
  surface_corrs <- function(A, orig, rep) {
    rep <- rep * A
    rep <- orig + rep
  }
  
  # Main Function -----------------------------------------------------------
  # for(i in sd) {
  for(z in start:iters){
    # seed <- seed + 1
    
    
    # >> Make Random Surfaces -------------------------------------------------
    
    # Specify random model
    
    #  * Surface 1 ------------------------------------------------------------
    
    # Create original and correlated surfaces for surface 1.
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
                               rescale = TRUE)
    plot(stack(r.stack,Resist))
    
    # Generate Points ---------------------------------------------------------
    pts <- unique(floor(cbind(runif(5000, 0.15 * r.dim, 0.85 * r.dim), 
                              runif(5000, 0.15 * r.dim, 0.85 * r.dim))))
    
    sample.thresh <- as.numeric(quantile(Resist, habitat))
    
    sample.extract <- extract(Resist, pts)
    sample.suit <- pts[sample.extract <= sample.thresh,]
    pts <- SpatialPoints(sample.suit[sample(nrow(sample.suit), n_ind, replace = F),])
    
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
    
    cdpop_sim <- cdpop(CDPOP.py = 'C:/CDPOP-master/src/CDPOP.py',
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
    
  } #  end iteration loop (z)
} # end function
