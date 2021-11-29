# 
# catraster <- raster("Data/all5.asc", crs = "EPSG:3044")
# catraster <- reclassify(catraster, cbind(c(3,4),c(4,3))) 
# 
# # PURPOSE: Landscape is a bit too big and it would make sense to have more than
# # 1 individual per cell
# # Let us try to keep linear features intact but have a coarser resolution
# 
# # Create garbage raster with 120m res
# garbage_120m <- aggregate(catraster, fact=3, fun=max, na.rm=TRUE) #max
# # Resample original using ngb
# set.seed(1) #ngb is random so final result may differ between runs
# catraster_resamp <- resample(catraster, garbage_120m, method = "ngb") 
# 
# catbrick <- layerize(catraster, falseNA=FALSE)
# catbrick <- catbrick*c(1,2,3,4,5) 
# catbrick_agg <- aggregate(catbrick, fact=3, fun=max, na.rm=TRUE)
# 
# catbrick_agg$X0[catbrick_agg$X0==0] <- NA
# catbrick_agg$X1[catbrick_agg$X1==0] <- NA
# catbrick_agg$X2[catbrick_agg$X2==0] <- NA
# catbrick_agg$X3[catbrick_agg$X3==0] <- NA
# catbrick_agg$X4[catbrick_agg$X4==0] <- NA
# catbrick_agg <- catbrick_agg-1
# plot(merge(catraster_resamp, catbrick_agg$X3))
# lol <- mosaic(catbrick_agg$X3, catbrick_agg$X4, fun=max)
# 
# catraster_coarser <- merge(lol,catraster_resamp)
# 
# ################################################################################
# ####### Checking step: Is the DYI downsampling  meaningful for analysis? #######
# # plot(catraster)
# # plot(catraster_coarser)
# ################################################################################
#      
# catraster_coarser_cropped <- crop(catraster_coarser,
#                                   extent(geosites)+c(-360,360,-360,360))
# catraster_coarser_cropped <- reclassify(catraster_coarser_cropped,
#                                         cbind(c(3,4),c(4,3)))
# plot(catraster_coarser_cropped)
# plot(geosites, add=TRUE, pch = 19)
# plot(catraster)
# plot(geosites, add=TRUE, pch = 19)
# 

# plot(catraster_SA_coarser_cropped, colNA="red")
# plot(geosites, add=TRUE, pch = 19)
# plot(catraster)
# plot(geosites, add=TRUE, pch = 19)
# 
# sum(getValues(!is.na(catraster_SA_coarser_cropped))) #1887




# Use a narrow smoother on it and get a reasonable number of individuals on it
# lizgridsmooth <- raster.gaussian.smooth(lizgrid, sigma=1, n=5)*200
# ################################################################################
# ########### Checking step: Is the smoothing  meaningful for analysis? ##########
# # plot(lizgridsmooth)
# # sum(na.omit(getValues(lizgridsmooth)))


# lizgen.genind@all.names <- lapply(lizgen.genind@all.names,
#                                   FUN=function(x){paste0("A", 0:(length(x)-1))})
# temp_list <- list()
# for (i in 1:length(lizgen.genind@all.names)){
#   temp_list[[i]] <- paste(names(lizgen.genind$all.names)[i],
#                           lizgen.genind$all.names[[i]],
#                           sep=".")
# }
