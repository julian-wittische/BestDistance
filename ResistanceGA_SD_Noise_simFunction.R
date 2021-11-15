#ResistanceGA_sim multisurface sim function

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
                               sd = c(0.5, 1.25), ## Specify level(s) of noise to add to resitance distance
                               start = NULL,
                               seed = NULL,
                               iters = 50,
                               method = 'standard',
                               JULIA_HOME = "C:/Users/peterman.73/AppData/Local/Julia-1.1.0/bin/",
                               out.dir = "G:/BoxSync/Research/Manuscripts/Manuscripts_in_Progress/Collaborations/Winiarski/resistance_ms/multisurface/"){
  library(ResistanceGA)
  library(RandomFields)
  
  if(is.null(start)) {
    stop("Specify integer `start` value!!!")
  }
  
  # Helper Function ---------------------------------------------------------
  surface_corrs <- function(A, orig, rep) {
    rep <- rep * A
    rep <- orig + rep
  }
  
  # Main Function -----------------------------------------------------------
  for(i in sd) {
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
      SD <- ifelse(i == 0.5, '0_5', '1_25')
      dir.create(paste0(out.dir, "Results/",
                        'SD__', SD,"/",
                        'iter__', z),
                 recursive = TRUE) 
      
      out <- paste0(out.dir, "Results/",
                    'SD__', SD,"/",
                    'iter__', z, "/")
      tmp <- paste0(tempdir(), '\\')
      
      
      # >> Create Truth ---------------------------------------------------------
      r.stack <- stack(orig_1,
                       orig_2,
                       random_1)
      names(r.stack) <- c('cont1', 'cont2', 'rand')
      
      # #Run 'GA.prep'
      GA.inputs <- GA.prep(ASCII.dir = r.stack[[-3]],
                           # maxiter = 1,
                           select.trans = list('A','A'),
                           # select.trans = list(c(1,2,3,5,7,8,9),c(1,2,3,5,7,8,9)),
                           max.cat = 500,
                           max.cont = 500,
                           Results.dir = tmp,
                           quiet = TRUE)
      
      #habitat <- 0.4 # Could add this to the function call for easy manipulation
      sample.thresh <- as.numeric(quantile(orig_1, habitat))
      # Select cells away from the edge that are suitable habitat
      rast.df <- as.data.frame(orig_1, xy = T) # convert raster to data frame for sampling
      mn.cell <- floor(min(rast.df$x) + (.20*max(rast.df$x)))
      mx.cell <- max(rast.df$x) - mn.cell
      
      habitat.rast <- which(rast.df$x > mn.cell &
                              rast.df$x < mx.cell & 
                              rast.df$y > mn.cell & 
                              rast.df$y < mx.cell & 
                              rast.df$layer <= sample.thresh)
      sample.sites <- sample(habitat.rast, size = n_samplepoints, replace = F)
      sample.xy <- rast.df[sample.sites, 1:2]
      Sample.coord <- SpatialPoints(sample.xy)
      
      # #  * Julia Prep -----------------------------------------------------------
      # 
      # jl.inputs <- jl.prep(n.Pops = n_samplepoints,
      #                      CS_Point.File = Sample.coord,
      #                      JULIA_HOME = JULIA_HOME)
      # 
      # PARM <- c(transformation_1, #Transformation for continuous surface 1
      #           r.tran_shape_1, #shape parameter for continuous surface 1
      #           r.tran_max_1, #max parameter for continuous surface 1
      #           transformation_2, #Transformation for continuous surface 1
      #           r.tran_shape_2, #shape parameter for continuous surface 1
      #           r.tran_max_2) #max parameter for continuous surface 2
      # 
      # 
      # # Combine resistance surfaces
      # Resist <- Combine_Surfaces(PARM = PARM,
      #                            jl.inputs = jl.inputs,
      #                            GA.inputs = GA.inputs,
      #                            out = NULL,
      #                            rescale = TRUE)
      # 
      # rast.df <- as.data.frame(Resist, xy = T) # convert raster to data frame for sampling
      # mn.cell <- floor(min(rast.df$x) + (.20*max(rast.df$x)))
      # mx.cell <- max(rast.df$x) - mn.cell
      # 
      # #Cells buffered from edge
      # habitat.rast <- which(rast.df$x > mn.cell &
      #                         rast.df$x < mx.cell & 
      #                         rast.df$y > mn.cell & 
      #                         rast.df$y < mx.cell) 
      # 
      # sample.sites <- sample(habitat.rast)
      # sample.xy <- rast.df[sample.sites, 1:3]
      # sample.xy <- as.data.frame(sample.xy, xy = T)
      # 
      # rast_min <- min(sample.xy$layer)
      # rast_max <- max(sample.xy$layer)
      # 
      # for(j in 1:nrow(sample.xy)){
      #   sample.xy$random[j] <- runif(1,rast_min, rast_max)
      #   col   }
      # 
      # sample.xy$keep <- ifelse(sample.xy$random >= sample.xy$layer, 1,0)
      # 
      # sample.xy <- subset(sample.xy, keep == 1)
      # 
      # sample.cells <- dplyr::sample_n(sample.xy, n_samplepoints)
      # sample.cells <- sample.cells[,1:2]
      # #sample.random <- rast.df[sample.cells]
      # Sample.coord <- SpatialPoints(sample.cells)
      # 
      # 
      # #plot(Resist, col=resistcols)
      # #plot(Sample.coord, pch=15, col="black", add= TRUE)
      # 
      # 
      # # Re-run Julia
      # jl.inputs <- jl.prep(n.Pops = n_samplepoints,
      #                      CS_Point.File = Sample.coord,
      #                      JULIA_HOME = JULIA_HOME)
      # 
      # #Pairwise resistance distacnes
      # jl.response <- Run_CS.jl(jl.inputs = jl.inputs,
      #                       r = Resist,
      #                       full.mat = T
      #                       )
      # 
      # #Add noise
      # ## The idea here is to have it be a relative error
      # noise <- rnorm(length(lower(jl.response)), mean = lower(jl.response), sd = i * lower(jl.response))
      # noise0 <- rnorm(length(lower(jl.response)), mean = 0, sd = i * mean(lower(jl.response)))
      # 
      # jl.response_noise <- lower(jl.response) + noise
      # jl.response_noise0 <- lower(jl.response) + noise0
      # 
      # resp1 <- ResistanceGA:::SCALE.vector(jl.response_noise, 0, 1)
      # resp2 <- ResistanceGA:::SCALE.vector(jl.response_noise0, 0, 1)
      # 
      # par(mfrow = c(2,1))
      # plot(jl.response_noise ~ lower(jl.response))
      # plot(jl.response_noise0 ~ lower(jl.response))
      # par(mfrow = c(1,1))
      # 
      # par(mfrow = c(2,1))
      # plot(resp1 ~ lower(jl.response))
      # plot(resp2 ~ lower(jl.response))
      # par(mfrow = c(1,1))
      # 
      # 
      # n.obs <- jl.inputs$n.Pops
      # genetic.mat <- matrix(0, n.obs, n.obs)
      # genetic.mat[lower.tri(genetic.mat)] <- resp1
      # 
      # graphics.off()
      # 
      # 
      # 
      # #  * gdist Prep -----------------------------------------------------------
      # 
      # gdist.inputs <- gdist.prep(n.Pops = n_samplepoints,
      #                            samples = Sample.coord)
      # 
      # PARM <- c(transformation_1, #Transformation for continuous surface 1
      #           r.tran_shape_1, #shape parameter for continuous surface 1
      #           r.tran_max_1, #max parameter for continuous surface 1
      #           transformation_2, #Transformation for continuous surface 1
      #           r.tran_shape_2, #shape parameter for continuous surface 1
      #           r.tran_max_2) #max parameter for continuous surface 2
      # 
      # 
      # # Combine resistance surfaces
      # Resist <- Combine_Surfaces(PARM = PARM,
      #                            gdist.inputs = gdist.inputs,
      #                            GA.inputs = GA.inputs,
      #                            out = NULL,
      #                            rescale = TRUE)
      # 
      # rast.df <- as.data.frame(Resist, xy = T) # convert raster to data frame for sampling
      # mn.cell <- floor(min(rast.df$x) + (.20*max(rast.df$x)))
      # mx.cell <- max(rast.df$x) - mn.cell
      # 
      # #Cells buffered from edge
      # habitat.rast <- which(rast.df$x > mn.cell &
      #                         rast.df$x < mx.cell & 
      #                         rast.df$y > mn.cell & 
      #                         rast.df$y < mx.cell) 
      # 
      # sample.sites <- sample(habitat.rast)
      # sample.xy <- rast.df[sample.sites, 1:3]
      # sample.xy <- as.data.frame(sample.xy, xy = T)
      # 
      # rast_min <- min(sample.xy$layer)
      # rast_max <- max(sample.xy$layer)
      # 
      # for(j in 1:nrow(sample.xy)){
      #   sample.xy$random[j] <- runif(1,rast_min, rast_max)
      #   col   }
      # 
      # sample.xy$keep <- ifelse(sample.xy$random >= sample.xy$layer, 1,0)
      # 
      # sample.xy <- subset(sample.xy, keep == 1)
      # 
      # sample.cells <- dplyr::sample_n(sample.xy, n_samplepoints)
      # sample.cells <- sample.cells[,1:2]
      # #sample.random <- rast.df[sample.cells]
      # Sample.coord <- SpatialPoints(sample.cells)
      # 
      # 
      # #plot(Resist, col=resistcols)
      # #plot(Sample.coord, pch=15, col="black", add= TRUE)
      # 
      # 
      # # Re-run gdistance
      # gdist.inputs <- gdist.prep(n.Pops = n_samplepoints,
      #                            samples = Sample.coord)
      # 
      # #Pairwise resistance distacnes
      # gdist.response <- Run_gdistance(gdist.inputs = gdist.inputs,
      #                              r = Resist
      # )
      # 
      # gdist.response <- as.matrix(gdist.response)
      # 
      # #Add noise
      # ## The idea here is to have it be a relative error
      # noise <- rnorm(length(lower(gdist.response)), mean = lower(gdist.response), sd = i * lower(gdist.response))
      # noise0 <- rnorm(length(lower(gdist.response)), mean = 0, sd = i * mean(lower(gdist.response)))
      # 
      # gd.response_noise <- lower(gdist.response) + noise
      # gd.response_noise0 <- lower(gdist.response) + noise0
      # 
      # resp1 <- ResistanceGA:::SCALE.vector(gd.response_noise, 0, 1)
      # resp2 <- ResistanceGA:::SCALE.vector(gd.response_noise0, 0, 1)
      # 
      # par(mfrow = c(2,1))
      # plot(gd.response_noise ~ lower(gdist.response))
      # plot(gd.response_noise0 ~ lower(gdist.response))
      # par(mfrow = c(1,1))
      # 
      # par(mfrow = c(2,1))
      # plot(resp1 ~ lower(gdist.response))
      # plot(resp2 ~ lower(cs.response))
      # par(mfrow = c(1,1))
      # 
      # 
      # n.obs <- cs.inputs$n.Pops
      # genetic.mat <- matrix(0, n.obs, n.obs)
      # genetic.mat[lower.tri(genetic.mat)] <- resp1
      # 
      # graphics.off()
      # 
      # 
      # 
      
      #  * CS Prep -----------------------------------------------------------
      
      cs.inputs <- CS.prep(n.Pops = n_samplepoints,
                           CS_Point.File = Sample.coord)
      
      PARM <- c(transformation_1, #Transformation for continuous surface 1
                r.tran_shape_1, #shape parameter for continuous surface 1
                r.tran_max_1, #max parameter for continuous surface 1
                transformation_2, #Transformation for continuous surface 1
                r.tran_shape_2, #shape parameter for continuous surface 1
                r.tran_max_2) #max parameter for continuous surface 2
      
      
      # Combine resistance surfaces
      Resist <- Combine_Surfaces(PARM = PARM,
                                 CS.inputs = cs.inputs,
                                 GA.inputs = GA.inputs,
                                 out = NULL,
                                 rescale = TRUE)
      
      rast.df <- as.data.frame(Resist, xy = T) # convert raster to data frame for sampling
      mn.cell <- floor(min(rast.df$x) + (.20*max(rast.df$x)))
      mx.cell <- max(rast.df$x) - mn.cell
      
      #Cells buffered from edge
      habitat.rast <- which(rast.df$x > mn.cell &
                              rast.df$x < mx.cell & 
                              rast.df$y > mn.cell & 
                              rast.df$y < mx.cell) 
      
      sample.sites <- sample(habitat.rast)
      sample.xy <- rast.df[sample.sites, 1:3]
      sample.xy <- as.data.frame(sample.xy, xy = T)
      
      rast_min <- min(sample.xy$layer)
      rast_max <- max(sample.xy$layer)
      
      for(j in 1:nrow(sample.xy)){
        sample.xy$random[j] <- runif(1,rast_min, rast_max)
        col   }
      
      sample.xy$keep <- ifelse(sample.xy$random >= sample.xy$layer, 1,0)
      
      sample.xy <- subset(sample.xy, keep == 1)
      
      sample.cells <- dplyr::sample_n(sample.xy, n_samplepoints)
      sample.cells <- sample.cells[,1:2]
      #sample.random <- rast.df[sample.cells]
      Sample.coord <- SpatialPoints(sample.cells)
      
      # Re-run CS
      cs.inputs <- CS.prep(n.Pops = n_samplepoints,
                           CS_Point.File = Sample.coord)
      
      #Pairwise resistance distacnes
      cs.response <- Run_CS(CS.inputs = cs.inputs,
                            r = Resist,
                            full.mat = T
      )
      
      #Add noise
      ## The idea here is to have it be a relative error
      noise <- rnorm(length(lower(cs.response)), mean = lower(cs.response), sd = i * lower(cs.response))
      noise0 <- rnorm(length(lower(cs.response)), mean = 0, sd = i * mean(lower(cs.response)))
      
      cs.response_noise <- lower(cs.response) + noise
      cs.response_noise0 <- lower(cs.response) + noise0
      
      resp1 <- ResistanceGA:::SCALE.vector(cs.response_noise, 0, 1)
      resp2 <- ResistanceGA:::SCALE.vector(cs.response_noise0, 0, 1)
      
      par(mfrow = c(2,1))
      plot(cs.response_noise ~ lower(cs.response))
      plot(cs.response_noise0 ~ lower(cs.response))
      par(mfrow = c(1,1))
      
      par(mfrow = c(2,1))
      plot(resp1 ~ lower(cs.response))
      plot(resp2 ~ lower(cs.response))
      par(mfrow = c(1,1))
      
      n.obs <- cs.inputs$n.Pops
      genetic.mat <- matrix(0, n.obs, n.obs)
      genetic.mat[lower.tri(genetic.mat)] <- resp1
      
      graphics.off()
      
      # Optimize: All Comb ----------------------------------------------------------------
      
      if(method == 'island') {
        GA.inputs <- GA.prep(ASCII.dir = r.stack,
                             Results.dir = 'all_comb',
                             max.cont = 500,
                             select.trans = list('A', 'A', 'A'),
                             maxiter = 40,
                             run = 8,
                             island.pop = 40,
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
      
      # jl.inputs <- jl.prep(n.Pops = n_samplepoints,
      #                      CS_Point.File = Sample.coord,
      #                      response = resp1,
      #                      JULIA_HOME = JULIA_HOME)
      
      # cs.inputs <- CS.prep(n.Pops = n_samplepoints,
      #                      CS_Point.File = Sample.coord,
      #                      response = resp1)
      
      gdist.inputs <- gdist.prep(n.Pops = n_samplepoints,
                                 samples = Sample.coord,
                                 response = resp1)
      
      Rga_out <- all_comb(gdist.inputs = gdist.inputs,
                          GA.inputs = GA.inputs,
                          # iters = 10,
                          results.dir = out
      )
      
      write.csv(genetic.mat, 
                file = paste0(out, "jl.response_noise.csv"))
      
      writeRaster(Resist,
                  paste0(out, "true_resist.asc"))
      
      sample.pts <- cbind((1:n_samplepoints),Sample.coord@coords) 
      colnames(sample.pts) <- c("pop","x","y")
      write.table(sample.pts, paste0(out,"sample_pts.txt"),
                  sep = '\t',
                  row.names = F,
                  col.names = T)
      
    } #  end iteration loop (z)
    
  } ## End SD loop
  
} # end function
