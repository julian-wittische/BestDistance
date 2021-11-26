#' @author Bill Peterman
#' @description Function to run CDPOP from R
#' 
#' @param CDPOP.py
#' @param sim_name Name for simulation results. Defaults to 'output'
#' @param pts Spatial points object
#' @param sim_dir Directory where simulation results will be written
#' @param resist_rast Resistance surface
#' @param agefilename Path to age file. Default will create and use a non-overlapping generations file.
#' @param mcruns Default = 1
#' @param looptime Default = 401; Number generations to conduct simulation
#' @param output_years Default = 50; Interval to write out simulation results
#' @param gridformat Default = 'genepop'; c('genepop', 'genalex', 'structure', 'cdpop')
#' @param cdclimgentime Default = 0. To initiate the CDClimate module, this is the generation/year that the next effective distance matrix will be read in at. You can specify multiple generations by separating each generation to read in the next cost distance matrix by ‘|’.  Then in the following surface columns, a separate file can be given for each generation.
#' @param matemoveno Default = 2; Uses Inverse Square (1 / (Cost Distance^2)). This function gets rescaled to min and threshold of the inverse square cost distance.
#' @param matemoveparA Not used with inverse square movement
#' @param matemoveparB Not used with inverse square movement
#' @param matemoveparC Not used with inverse square movement
#' @param matemovethresh Default = 'max'; The maximum movement is the maximum resistance distance
#' @param output_matedistance
#' @param sexans Default = 'N'; No selfing
#' @param Freplace Default = 'N'; Females mate without replacement
#' @param Mreplace Default = 'N'; Males mate without replacement
#' @param philopatry Default = 'N';
#' @param multiple_paternity Default = 'N'; No philopatry of Males or Females
#' @param selfans Default = 'N'; No selfing
#' @param Fdispmoveno Default = NULL. Will be set equal to matemoveno 
#' @param FdispmoveparA Not used with inverse square movement
#' @param FdispmoveparB Not used with inverse square movement
#' @param FdispmoveparC Not used with inverse square movement
#' @param Fdispmovethresh Default = NULL. Will be set to matemovethresh
#' @param Mdispmoveno Default = NULL. Will be set equal to matemoveno 
#' @param MdispmoveparA Not used with inverse square movement
#' @param MdispmoveparB Not used with inverse square movement
#' @param MdispmoveparC Not used with inverse square movement
#' @param Mdispmovethresh Default = NULL. Will be set to matemovethresh
#' @param offno Default = 2; Poisson draw around ‘mean fecundity’
#' @param MeanFecundity Default = 5; Specifies mean fecundity in age variable file.
#' @param Femalepercent Default = 50
#' @param EqualsexratioBirth Default = 'N'
#' @param TwinningPercent Default = 0
#' @param popModel Default = 'exp'
#' @param r Population growth rate. No applicable when using exponential growth rate
#' @param K_env Equal to the number of individuals simulated
#' @param subpopmortperc Default = 0|0|0|0; Not using subpopulation features
#' @param muterate Default = 0.0005
#' @param mutationtype Default = 'forward'
#' @param loci Default = 1000; For simulating SNP-like markers
#' @param intgenesans Default = 'random'; Random initiation of alleles
#' @param allefreqfilename Default = 'N'
#' @param alleles Default = 2; For simulating SNP-like markers
#' @param mtdna Default = 'N'
#' @param startGenes Default = 0; 
#' @param cdevolveans Default = 'N'; No loci are under selection
#' @param startSelection Default = 0; No selection
#' @param betaFile_selection Default = 'N'; No selection
#' @param epistasis Default = 'N'; No epigenetics
#' @param epigeneans Default = 'N'; No epigenetics
#' @param startEpigene Default = 0; No epigenetics
#' @param betaFile_epigene Default = 'N'; No epigenetics
#' @param cdinfect Default = 'N'; No epigenetics
#' @param transmissionprob Default = 0; No epigenetics
#' 
#' 
#' 
cdpopJW <- function(CDPOP.py,
                  sim_name = 'output_',
                  pts,
                  sim_dir,
                  resist_rast,
                  resist_mat = NULL,
                  agefilename = NULL,
                  mcruns = 1,
                  looptime = 400,
                  output_years = 50,
                  gridformat = 'genepop',
                  cdclimgentime = 0,
                  matemoveno = 2,
                  matemoveparA = 0,
                  matemoveparB = 0,
                  matemoveparC = 0,
                  matemovethresh = 'max',
                  output_matedistance = 'N',
                  sexans = 'Y',
                  Freplace = 'N',
                  Mreplace = 'N',
                  philopatry = 'N',
                  multiple_paternity = 'N',
                  selfans = 'N',
                  Fdispmoveno = NULL,
                  FdispmoveparA = 0,
                  FdispmoveparB = 0,
                  FdispmoveparC = 0,
                  Fdispmovethresh = NULL,
                  Mdispmoveno = NULL,
                  MdispmoveparA = 0,
                  MdispmoveparB = 0,
                  MdispmoveparC = 0,
                  Mdispmovethresh = NULL,
                  offno = 2,
                  MeanFecundity = 5,
                  Femalepercent = 50,
                  EqualsexratioBirth = 'N',
                  TwinningPercent = 0,
                  popModel = 'exp',
                  r = 1,
                  K_env = length(pts),
                  subpopmortperc = 0,
                  muterate = 0.0005,
                  mutationtype = 'random',
                  loci = 1000,
                  intgenesans = 'random',
                  allefreqfilename = 'N',
                  alleles = 2,
                  mtdna = 'N',
                  startGenes = 0,
                  cdevolveans = 'N',
                  startSelection = 0,
                  betaFile_selection = 'N',
                  epistasis = 'N',
                  epigeneans = 'N',
                  startEpigene = 0,
                  betaFile_epigene = 'N',
                  cdinfect = 'N',
                  transmissionprob = 0){
  
  
  # Install / Load Libraries ------------------------------------------------
  
  list.of.packages <- c("gdistance",
                        "adegenet",
                        "readr",
                        "raster")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  library(raster)
  library(gdistance)
  library(adegenet)
  library(readr)
  
  
  # Create directories ------------------------------------------------------
  
  if(!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
  suppressWarnings(
    dir.create(paste0(sim_dir,"/data/"), recursive = TRUE)
  )
  data_dir <- paste0(sim_dir,"/data/")
  
  
  # Fill NULL ---------------------------------------------------------------
  
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
  
  # Age file ----------------------------------------------------------------
  
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
  
  
  # XY File -----------------------------------------------------------------
  
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
  
  # Resistance Distance -----------------------------------------------------
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
  
  # CDPOP input ----------------------------------------------------------
  
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
  
  
  # Run CDPOP ---------------------------------------------------------------
  print("Running CDPOP...")
  
  system(paste("python", CDPOP.py, data_dir, "CDPOP_inputs.csv", sim_name))
  
  
  # Import Results ----------------------------------------------------------
  
  fi <- file.info(list.files(path = sim_dir,
                             pattern = "grid",
                             recursive = T,
                             full.names = T))
  
  ## Get latest simulation results
  newest_sim <- dirname(rownames(fi)[which.max(fi$mtime)])
  
  grid_dir <- list.files(path = newest_sim,
                         pattern = "grid",
                         recursive = T,
                         full.names = T)
  
  read.grid <- function(grid,
                        pops = NULL){
    suppressWarnings(
      cdpop_out <- read_csv(grid,
                            col_types = cols(Subpopulation = col_skip(),
                                             #XCOORD = col_skip(), YCOORD = col_skip(),
                                             sex = col_skip(), age = col_skip(),
                                             infection = col_skip(), DisperseCDist = col_skip(),
                                             hindex = col_skip()))
    )
    #geogr <- read_csv(grid)[which(cdpop_out$ID != "OPEN"),c("XCOORD", "YCOORD")]
    occ_pop <- which(cdpop_out$ID != "OPEN")
    
    if(!is.null(pops)) {
      return(occ_pop)
    } else {
      
      cd_df <- as.data.frame(cdpop_out[occ_pop,c(-1,-2,-3)])
      cd_df[,ncol(cd_df)] <- gsub(",","",cd_df[,ncol(cd_df)])
      cd_df <- apply(as.matrix(cd_df),2,as.numeric)
      
      fakedf <- data.frame(matrix(rep(paste(paste0("A", rep_len(0:(alleles-1),length.out = nrow(cd_df))),
                                            paste0("A", rep_len(0:(alleles-1),length.out = nrow(cd_df))),
                                            sep="/"), loci), ncol=loci))
      
      colnames(fakedf) <- paste0("L",1:loci)
      ncode <- 1
      gi <- adegenet::df2genind(fakedf, ploidy=2, sep="/", type="codom")
      gi@tab <- cd_df
      colnames(gi@tab) <- paste(names(gi$all.names), unlist(gi$all.names), sep=".")
      gi@other$xy <- cdpop_out[occ_pop, c(1,2)]
      gi@tab <- apply(gi@tab, 2, as.integer)
      return(gi)
    }
  }
  
  grid_list <- lapply(grid_dir, read.grid)
  pop_list <- lapply(grid_dir, read.grid, pops = TRUE)
  
  gens <- basename(grid_dir) %>% sub('.csv', '', .) %>% # <
    sub('grid', '',.) %>% as.numeric()
  
  grid_list <- grid_list[order(gens)]
  pop_list <- pop_list[order(gens)]
  
  names(pop_list) <- names(grid_list) <- paste0('gen_', sort(gens))
  
  
  # Wrap-up -----------------------------------------------------------------
  
  out <- list(grid_list = grid_list,
              pop_list = pop_list)
  return(out)
  
}

# PCA dist -------------------------------------------------------

pca_dist <- function(gi,
                     n_axes = 64){
  a_tab <- adegenet::tab(gi)
  pc <- prcomp(a_tab)
  pc_dist <- as.matrix(dist(pc$x[,1:n_axes]))
  return(pc_dist)
}

# Random Samples ----------------------------------------------------------

## Randomly select populations and individuals from within populations

gi_samp <- function(gi,
                    n_ind = 100) {
  ind_samp <- sort(sample(1:nInd(gi), n_ind))
  gi_s <- gi[ind_samp]
  
  out <- list(genind = gi_s,
              pop_samp = ind_samp)
}
