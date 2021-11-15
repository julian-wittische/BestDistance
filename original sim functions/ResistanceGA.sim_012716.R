#ResistanceGA_sim function

  wd <- "C:/ResistanceGA_Sim/"  
  
  id <- 1
  r.dim <- 50
  n_samplepoints <- 25
  RMexpvar <- 1
  RMexpscale <- 5
  transformation <- "Monomolecular"
  r.tran_shape <- 3
  r.tran_max <- 50
  habitat <- 0.4
  sd <- .20
  Mid_smooth <- 4
  Low_smooth <- 1
  High_smooth <- 7
  parameters <- 15
  
  max.attempts <- 1
  
  out.df <- as.data.frame(matrix(0,length(max.attempts),parameters))
  colnames(out.df) <- c("id", 
                        "r.dim",
                        "n_samplepoints", 
                        "RMexpvar",
                        "RMexpscale",
                        "trans",
                        "transformation",
                        "r.tran_shape",
                        "r.tran_max",
                        "habitat",
                        "sd", 
                        "e_transformation", 
                        "e_shape", 
                        "e_max",
                        "pcorr",
                        "ks_D",
                        "ks_p")
  
  library(RandomFields)
  library(ResistanceGA)
  library(gridio)
  library(raster)
  
  CS.program <- paste('"C:/Program Files/Circuitscape/cs_run.exe"')
  
  #for(z in 1:max.attempts){
    # Specify random model
    ####################################################
    ####################################################
    # Create original and correlated surface
  corr = c(0.9, 0.7, 0.5, 0.3, 0.1)
  A <- sqrt((1/(corr^2)) - 1)  # calcualate A
  coo <- 1:50
  model <- RMexp(var=RMexpvar, scale=RMexpscale)
  rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
  surface_orig <- as.matrix(rf.sim)
  
  rf.sim <- RFsimulate(model = model, x=coo, y=coo, grid=TRUE)
  surface_rep <- as.matrix(rf.sim)
  
  surface_corrs <- function(A, orig, rep) {
    rep <- rep * A
    rep <- orig + rep
  }
  
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
  
  extent(out_rasts$cor_0.7) <- bb
  out_rasts$cor_0.7 <- setExtent(out_rasts$cor_0.7, bb, keepres=TRUE)
  cor_0.7 <- out_rasts$cor_0.7
  
  extent(out_rasts$cor_0.5) <- bb
  out_rasts$cor_0.5 <- setExtent(out_rasts$cor_0.5, bb, keepres=TRUE)
  cor_0.5 <- out_rasts$cor_0.5
  
  extent(out_rasts$cor_0.3) <- bb
  out_rasts$cor_0.3 <- setExtent(out_rasts$cor_0.3, bb, keepres=TRUE)
  cor_0.3 <- out_rasts$cor_0.3
  
  extent(out_rasts$cor_0.1) <- bb
  out_rasts$cor_0.1 <- setExtent(out_rasts$cor_0.1, bb, keepres=TRUE)
  cor_0.1 <- out_rasts$cor_0.1
  
  r.stack <- stack(orig, cor_0.9, cor_0.7, cor_0.5, cor_0.3, cor_0.1)
  names(r.stack) <- c("orig","cor_0.9", "cor_0.7", "cor_0.5", "cor_0.3", "cor_0.1")
  
  #Correlation between the two surfaces
  pairs(r.stack)
  plot(r.stack)
  
  #Get exact correlation between correlated surfaces and original
  pcorr_orig_0.9 <- layerStats(r.stack, stat = 'pearson')[[1]][1,2]
  pcorr_orig_0.7 <- layerStats(r.stack, stat = 'pearson')[[1]][1,3]
  pcorr_orig_0.5 <- layerStats(r.stack, stat = 'pearson')[[1]][1,4]
  pcorr_orig_0.3 <- layerStats(r.stack, stat = 'pearson')[[1]][1,5]
  pcorr_orig_0.1 <- layerStats(r.stack, stat = 'pearson')[[1]][1,6]
  
  #Save the surface in the wd
    writeRaster(orig,filename=paste0(wd,"cont.asc"),overwrite=TRUE)
    
    
    #Prep for resistanceGA
    GA.inputs <- GA.prep(ASCII.dir=wd,
                         maxiter = 5,
                         #select.trans = trans,
                         max.cat=500,
                         max.cont=500,
                         quiet=TRUE)
    
    #Transform continuous surface using the Resistance.tran function 
    r.tran <- Resistance.tran(transformation= transformation, shape=r.tran_shape, max=r.tran_max, r=orig) 
    plot(r.tran)
    ##################################################################################
    ##################################################################################
    
    ## To constrain points to fall only within similar habitat, use the code below
    ## Note that this will generate random points, not points evenly spaced on a grid
    
    #habitat <- 0.4 # Could add this to the function call for easy manipulation
    sample.thresh <- as.numeric(quantile(r.tran, habitat))
    
    # Select cells away from the edge that are suitable habitat
    rast.df <- as.data.frame(r.tran, xy = T) # convert raster to data frame for sampling
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
    samples <- cbind((1:n_samplepoints),Sample.coord@coords) 
    colnames(samples)<-c("pop","x","y")
    
    ## View random points within habitat
    plot(Sample.coord, pch=16, col="blue", add=TRUE)
    
    plot(r.tran <= sample.thresh)
    plot(Sample.coord, pch = 18, col = 'red', add = T)
    
    #Save samples.txt in wd
    write.table(samples,file=paste0(wd,"samples.txt"),
    sep="\t",col.names=F,row.names=F)
    
    ##################################################################################
    ##################################################################################
    
    CS.inputs <- CS.prep(n.POPS=n_samplepoints,
                         CS_Point.File=paste0(wd,"samples.txt"),
                         CS.program=CS.program)
    
    #Run the transformed resistance surface through CIRCUITSCAPE to get effective resistance between each pair of points. Run.CS
    #returns the lower half of the pairwise resistance matrix for use with the optimization prep functions.  This will be the response that is optimized on.
    
    #Create the true resistance/response surface
    CS.response <- Run_CS(CS.inputs=CS.inputs, 
                          GA.inputs=GA.inputs, 
                          r=r.tran)
    
    #Add noise
    ## The idea here is to have it be a relative error
    noise <- rnorm(length(CS.response), 0, sd) 
    
    CS.response_noise <- CS.response + noise
    
    #Rerun CS.prep including the newly created CS.response
    CS.inputs <- CS.prep(n.POPS=n_samplepoints,
                         response=CS.response_noise,
                         CS_Point.File=paste0(wd,"samples.txt"),
                         CS.program=CS.program)
    
    #Run the single surface optimization
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    ##Results
    #Grab the estimated transformation, shape and max for out.df
    e_transformation <- as.character(SS_RESULTS$ContinuousResults$Equation[1])
    
    #e_transformation <- as.character(e_transformation)
    e_shape <-   SS_RESULTS$ContinuousResults$shape
    e_max <- SS_RESULTS$ContinuousResults$max
    e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
    #Stack truth and optim
    ## Have to create the surface using the optimized parameters
    optim <- Resistance.tran(transformation = e_transformation, 
                             shape = e_shape,
                             max = e_max,
                             r=orig) 
    
    #Compare results with truth
    ## Keep in mind that you have "TRUTH" as well as "TRUTH w/ NOISE", which is what is being optimized
    ## As such, after optimization you can compare:
    ## 1) The correlation in the true and optimized resistance surfaces 
    r.stack <- stack(orig, optim)
    names(r.stack) <- c("Truth", "Optimized")
    #Correlation between the two surfaces
    pairs(r.stack)
    # This gets the correlation of the true raster and optimized raster
    pcorr_surfaces <- layerStats(r.stack, stat = 'pearson')[[1]][1,2]
    
    #KS test
    extract.r.tran <- as.data.frame(r.tran)
    extract.optim <- as.data.frame(optim)
    #plot(ecdf(extract.optim$layer),verticals=TRUE)
    #plot(ecdf(extract.r.tran$layer),verticals=TRUE,lty=3)
    ks <- ks.test(extract.optim$layer,extract.r.tran$layer)
    ks_D <- ks$statistic[1]
    ks_p <- ks$p.value[1]
    
    ## 2) The correlation in the noisy CS resistance (the response in the simulation) and the optimized CS resistance
    #Pairwise resistance distance for optim?
    ## You already calculated this above, no need to do it again
    # r.tran <- Resistance.tran(transformation= e_transformation, shape=e_shape, max=e_max, r=continuous)
    CS.response_optim <- Run_CS(CS.inputs=CS.inputs, GA.inputs=GA.inputs, r=optim)
    
    pcorr_rdistances <- cor(CS.response_optim, CS.response_noise)
    
    #Now want lower bandwith surface to fit. Need to bring in original surface
    writeRaster(orig,filename=paste0(wd,"cont.asc"),overwrite=TRUE)
    in.grid <- readasciigrid(path=paste(wd,"cont.asc", sep=""))
    gridinit()
    setwindow(in.grid)
    SD <- Mid_smooth
    smooth_midband <- gaussiansmooth(in.grid,sd=SD,max.r=3*SD,kernel.dim=Inf)
    plot(smooth_midband)
    
    #Save the kernel surface in the wd
    writeasciigrid(smooth_midband,paste0(wd,"cont.asc"))
    
    #Measure correlation between lowband and true surface
    smooth_midband <- raster(x=smooth_midband$m, xmx=xmax(orig), ymx=ymax(orig), xmn=xmin(orig), ymn=ymin(orig))
    r.stack <- stack(true_smooth, smooth_midband)
    names(r.stack) <- c("true_smooth", "smooth_midband")
    #Correlation between the two surfaces
    pairs(r.stack)
    # This gets the correlation of the true raster and optimized raster
    corr_true_midband <- layerStats(r.stack, stat = 'pearson')[[1]][1,2]
    
    #Rerun the optimization on the lowband kernel surface
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    KSmid_e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
    ## Smooth the higher bandwith to fit.  Need to bring in original continuous surface again
    writeRaster(orig,filename=paste0(wd,"cont.asc"),overwrite=TRUE)
    in.grid <- readasciigrid(path=paste(wd,"cont.asc", sep=""))
    gridinit()
    setwindow(in.grid)
  
    SD <- High_smooth
    smooth_highband <- gaussiansmooth(in.grid,sd=SD,max.r=3*SD,kernel.dim=Inf)
    plot(smooth_highband)
    
    #Save the kernel surface in the wd
    writeasciigrid(smooth_highband,paste0(wd,"cont.asc"))
    
    #Measure correlation between lowband and true surface
    smooth_highband <- raster(x=smooth_highband$m, xmx=xmax(orig), ymx=ymax(orig), xmn=xmin(orig), ymn=ymin(orig))
    r.stack <- stack(orig, smooth_highband)
    names(r.stack) <- c("true_smooth", "smooth_highband")
    #Correlation between the two surfaces
    pairs(r.stack)
    # This gets the correlation of the true raster and optimized raster
    corr_true_highband <- layerStats(r.stack, stat = 'pearson')[[1]][1,2]
    
    #Rerun the optimization on the kernel surface
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    KShigh_e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
  #Now fit correlated surfaces beginning with 0.90
    #Save the correlated surface in the wd
    writeRaster(cor_0.9,filename=paste0(wd,"cont.asc"),overwrite=TRUE)
    
    #Rerun the optimization on the kernel surface
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    Cor_0.9_e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
    #Now fit correlated surfaces beginning with 0.70
    #Save the kernel surface in the wd
    writeasciigrid(correlated_0.7,paste0(wd,"cont.asc"))
    
    #Rerun the optimization on the kernel surface
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    Cor_0.7_e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
    #Now fit correlated surfaces beginning with 0.50
    #Save the correlated surface in the wd
    writeRaster(cor_0.5,filename=paste0(wd,"cont.asc"),overwrite=TRUE)
    #Save the kernel surface in the wd
    writeasciigrid(correlated_0.5,paste0(wd,"cont.asc"))
    
    #Rerun the optimization on the kernel surface
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    Cor_0.5_e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
    #Now fit correlated surfaces beginning with 0.30
    #Save the correlated surface in the wd
    writeRaster(cor_0.3,filename=paste0(wd,"cont.asc"),overwrite=TRUE)
    #Save the kernel surface in the wd
    writeasciigrid(correlated_0.3,paste0(wd,"cont.asc"))
    
    #Rerun the optimization on the kernel surface
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    Cor_0.3_e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
    #Now fit correlated surfaces beginning with 0.10
    #Save the correlated surface in the wd
    writeRaster(cor_0.1,filename=paste0(wd,"cont.asc"),overwrite=TRUE)
    #Save the kernel surface in the wd
    writeasciigrid(correlated_0.1,paste0(wd,"cont.asc"))
    
    #Rerun the optimization on the kernel surface
    SS_RESULTS <- SS_optim(CS.inputs=CS.inputs,
                           GA.inputs=GA.inputs)
    
    Cor_0.1_e_AIC <- SS_RESULTS$ContinuousResults$AICc
    
    
    out <- data.frame(e_transformation = e_transformation,
                      e_shape = e_shape,
                      e_max = e_max,
                      pcorr_surfaces = pcorr_surfaces,
                      pcorr_rdistance = pcorr_rdistance,
                      ks_D = ks_D,
                      ks_p = ks_p, 
                      id=id,
                      r.dim=r.dim,
                      n_samplepoints=n_samplepoints,
                      RMexpvar=RMexpvar,
                      RMexpscale=RMexpscale,
                      habitat=habitat,
                      transformation=transformation,
                      r.tran_shape=r.tran_shape,
                      r.tran_max=r.tran_max,
                      habitat=habitat,
                      sd=sd)
    
    out.dir <- "C:/ResistanceGA_Sim/Results/"
    out.file <- paste(out.dir, "ResistanceGA_sim.csv", sep = "")
    if(!file.exists(out.file)) {
      write.table(out, out.file, sep = ",", row.names = F, col.names = T)
    } else {
      write.table(out, out.file, sep = ",", append = T, row.names = FALSE, col.names = F)
    }
    
    
                    