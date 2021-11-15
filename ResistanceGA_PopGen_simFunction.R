#ResistanceGA_sim function

ResistanceGAMS.sim <- function(r.dim = 50,
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
                               habitat = 0.4,
                               n_samplepoints = 75,
                               start = NULL,
                               seed = NULL,
                               iters = 50,
                               method = 'standard',
                               island.pop = 50, ## 40
                               JULIA_HOME = "C:/Users/peterman.73/AppData/Local/Julia-1.1.0/bin/",
                               out.dir = "G:/BoxSync/Research/Manuscripts/Manuscripts_in_Progress/Collaborations/Winiarski/resistance_ms/multisurface/",
                               n.pops = 100,  # total of plots to simulate (100)
                               n.ind = 30, #50  # Number of individuals per pop (50)
                               sex.ratio = 0.5,
                               n.loci = 15, # 15
                               n.allels = 15, # 15
                               steps = 200,  # Generations (250) - 100 generations generally gives same results, just quicker to run.
                               n.offspring = 2, # 10 [2]
                               mig.rate = 0.2, #0.05 # 0.10 [0.075] 0.25
                               disp.max = 0.125, # Percent of max (0.15) 0.125; 
                               disp.rate = 0.25, #0.25# 0.35 [0.25]
                               mut.rate = 0.0005 # 0.0005  
){
  
  library(ResistanceGA)
  library(RandomFields)
  library(PopGenReport)
  
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
    
    bb <- extent(0.5, 50.5, 0.5, 50.5)
    
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
    
    bb <- extent(0.5, 50.5, 0.5, 50.5)
    
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
    
    bb <- extent(0.5, 50.5, 0.5, 50.5)
    
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
    dir.create(paste0(out.dir, "Results/",
                      'PopGen/',
                      'iter__', z),
               recursive = TRUE) 
    
    out <- paste0(out.dir, "Results/",
                  'PopGen/',
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
    
    sample.thresh <- as.numeric(quantile(Resist, habitat))
    # Select cells away from the edge that are suitable habitat
    rast.df <- as.data.frame(Resist, xy = T) # convert raster to data frame for sampling
    mn.cell <- floor(min(rast.df$x) + (.20*max(rast.df$x)))
    mx.cell <- max(rast.df$x) - mn.cell
    
    habitat.rast <- which(rast.df$x > mn.cell &
                            rast.df$x < mx.cell & 
                            rast.df$y > mn.cell & 
                            rast.df$y < mx.cell & 
                            rast.df$layer <= sample.thresh)
    
    sample.sites <- sample(habitat.rast, size = n.pops, replace = F)
    sample.xy <- rast.df[sample.sites, 1:2]
    sp.coord <- SpatialPoints(sample.xy)
    
    
    # >> Calculate cost distance -------------------------------------------------
    jl.inputs <- jl.prep(n.Pops = n.pops,
                         CS_Point.File = sp.coord,
                         JULIA_HOME = JULIA_HOME)
    
    ##  Pairwise resistance distances
    r.dist <- Run_CS.jl(jl.inputs = jl.inputs,
                        r = Resist,
                        full.mat = T
    )
    
    plot(Resist)
    plot(sp.coord, add = T)
    
    
    # Conduct simulation ------------------------------------------------------
    pop.sim_init <- init.popgensim(n.pops = n.pops,
                                   n.ind = n.ind, 
                                   sex.ratio = sex.ratio,
                                   n.loci = n.loci,
                                   n.allels = n.allels, 
                                   n.cov = 3)
    
    sim.out <-  run.popgensim(simpops = pop.sim_init,
                              steps = steps,
                              cost.mat = r.dist,
                              n.offspring = n.offspring,
                              n.ind = n.ind,
                              mig.rate = mig.rate,
                              disp.max = disp.max * max(r.dist),
                              disp.rate = disp.rate,
                              n.allels = n.allels,
                              mut.rate = mut.rate)
    
    pops.gi <- pops2genind(sim.out)
    
    fst.out <- pairwise.fstb(pops.gi)
    
    ## Convert `genind` object to `genpop` object
    pop.gp <- genind2genpop(pops.gi)
    
    ## Calculate Dc from `genpop` object
    dc.out <- dist.genpop(pop.gp, method = 2, diag = T, upper = T) %>% as.matrix()
    
    ## Need to subsample here from dc.out
    select.pops <- sort(sample(length(sp.coord), 
                               size = n_samplepoints,
                               replace = F))
    
    dc.select <- dc.out[select.pops, select.pops]
    fst.select <- fst.out[select.pops, select.pops]
    Sample.coord <- sp.coord[select.pops]
    
    r.select <- r.dist[select.pops,select.pops]
    
    
    ## Quick plot select
    par(mfrow = c(2,2))
    plot(lower(dc.select) ~ lower(r.select))
    plot(lower(fst.select) ~ lower(r.select))
    plot(lower(dc.select) ~ c(dist(Sample.coord@coords)))
    plot(lower(fst.select) ~ c(dist(Sample.coord@coords)))
    par(mfrow = c(1,1))
    
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
                           maxiter = 65, #65
                           run = 15,
                           # pop.size = 5,
                           parallel = 4
      )
    }
    
    
    jl.inputs <- jl.prep(n.Pops = n_samplepoints,
                         CS_Point.File = Sample.coord,
                         response = lower(dc.select),
                         JULIA_HOME = JULIA_HOME)
    
    Rga_out <- all_comb(jl.inputs = jl.inputs,
                        GA.inputs = GA.inputs,
                        # iters = 10,
                        results.dir = out
    )
    
    write.csv(dc.select, 
              file = paste0(out, "dc_select.csv"))
    
    write.csv(dc.out, 
              file = paste0(out, "dc_all.csv"))
    
    
    writeRaster(Resist,
                paste0(out, "true_resist.asc"))
    
    sample.pts <- cbind((1:n_samplepoints),Sample.coord@coords) 
    colnames(sample.pts) <- c("pop","x","y")
    write.table(sample.pts, paste0(out,"sample_pts.txt"),
                sep = '\t',
                row.names = F,
                col.names = T)
    
    sample.pts_all <- cbind((1:n_samplepoints),sp.coord@coords) 
    colnames(sample.pts_all) <- c("pop","x","y")
    write.table(sample.pts_all, paste0(out,"sample_pts_all.txt"),
                sep = '\t',
                row.names = F,
                col.names = T)
    
  } #  end iteration loop (z)
} # end function
