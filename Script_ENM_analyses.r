library(ade4)
library(ape)
library(blockCV)
library(diagram)
library(dismo)
library(gbm)
library(geosphere)
library(ggplot2)
library(lubridate)
library(ncdf4)
library(ncf)
library(picante)
library(pgirmess)
library(phytools)
library(RColorBrewer)
library(raster)
library(rgdal)
library(rgeos)
library(seqinr)
library(sp)
library(vioplot)

# 1. Preparation of the shapefile and different climatic rasters
# 2. Loading and ploting all the Bombus records by continent
# 3. Defining the background areas for each study area
# 4. BRT analyses with standard or spatial cross-validation
# 5. Estimation of spatial sorting bias (SSB; for spatial autocorrelation)
# 6. Comparison of the response curves for each climatic variable
# 7. Analyses of the relative influence of each climatic variable
# 8. Plotting the different climatic variables and their past projections
# 9. BRT projections based on past and present climatic variables
# 10. Post hoc analyses to compare past and present climatic niches

directory = "Bombus_obs_160620"; savingPlots = FALSE
timePeriods = c("Last GLM","Mid-Holocene","Present time")
timePeriodIDs = c("Past_GLM","Past_Holo","Present_t0")
species = read.csv("Bombus_species.csv", header=T)
continents = shapefile("Continents_shapefile/Continents.shp")
mask = raster("Mask_raster_file.nc4"); minYear = 1950; maxYear = 2000
envVariableNames1 = read.csv("Bioclimatic_vars.csv", head=T)[,"name"]

# 1. Preparation of the shapefile and different climatic rasters

	# 1.1. Preparation of the European shapefile that will be used as a mask

europe1 = subset(continents, continents$CONTINENT=="Europe")
europe1 = crop(europe1, extent(-12,29,35,72)); polygons = list(); c = 0
for (i in 1:length(europe1@polygons[[1]]@Polygons))
	{
		if (europe1@polygons[[1]]@Polygons[[i]]@area > 1)
			{
				c = c+1; polygons[[c]] = europe1@polygons[[1]]@Polygons[[i]]
			}
	}
pols = Polygons(polygons, 1); pols_list = list(); pols_list[[1]] = pols
europe2 = SpatialPolygons(pols_list); europe3 = gSimplify(europe2, 0.1)

	# 1.2. Loading all the different bioclimatic rasters for t0

envVariables = list(); envVariableNames2 = envVariableNames1
for (i in 1:length(envVariableNames1))
	{
		if (i < 10) envVariableNames2[i] = paste0("bio0",i)
		if (i >= 10) envVariableNames2[i] = paste0("bio",i)
		envVariable = raster(paste0("WorldClim_1_rasters/Present_t0/bio",i,".bil"))
		names(envVariable) = envVariableNames2[i]; envVariables[[i]] = envVariable
	}
for (i in 1:length(envVariableNames1))
	{
		envVariables[[i]][is.na(mask[])] = NA
		envVariables[[i]] = crop(envVariables[[i]], europe2, snap="out")
		envVariables[[i]] = mask(envVariables[[i]], europe2)
	}
nullRaster = envVariables[[1]]
nullRaster[!is.na(nullRaster[])] = 1
names(nullRaster) = "nullRaster"

# 2. Loading and ploting all the Bombus records by continent

observations_list = list()
for (i in 1:dim(species)[1])
	{
		tab1 = read.csv(paste0(directory,"/",species[i,1],".csv"), header=T)
		tab2 = tab1[which((tab1[,"year"]>=minYear)&(tab1[,"year"]<=maxYear)),c("longitude","latitude")]
		observations_list[[i]] = tab2[which(!is.na(raster::extract(nullRaster,tab2))),]
	}
if (savingPlots == TRUE)
	{
		pdf(paste0("All_the_figures_&_SI/Bombus_observations_NEW.pdf"), width=14, height=14)
		par(mfrow=c(6,8), oma=c(0,0,0,0), mar=c(0,1.5,1.5,0), lwd=0.4, col="gray30")
		for (i in 1:dim(species)[1])
			{
				plot(nullRaster, col="gray90", ann=F, legend=F, axes=F, box=F)
				plot(europe3, add=T, border="gray50", lwd=1.5)
				points(observations_list[[i]], col="gray30", pch=3, cex=0.3, lwd=0.3)
				mtext(paste0("B. ",species[i,1]), side=3, line=-2, at=0, cex=0.75, col="gray30")
			}
		dev.off()
	}

# 3. Defining the background areas for each study area

species = read.csv("Bombus_species.csv", header=T, sep=",")
names(nullRaster) = "nullRaster"; allObservationsOnTheContinent = c()
for (i in 1:dim(species)[1])
	{
		allObservationsOnTheContinent = rbind(allObservationsOnTheContinent, observations_list[[i]])
	}
backgroundCells = unique(raster::extract(nullRaster, allObservationsOnTheContinent, cellnumbers=T))
background = nullRaster; background[!(1:length(background[]))%in%backgroundCells] = NA

# 4. BRT analyses with standard or spatial cross-validation

samplingPtsMinDist = function(observations, minDist=500, nberOfPoints=5)
	{
		indices = rep(NA, nberOfPoints)
		selection_list = list(1:nrow(observations)) 
  		indices[1] = sample(1:dim(observations)[1], 1)
		dists = list(spDistsN1(as.matrix(observations), as.matrix(observations[indices[1],]), longlat=T))
		for (i in 2:nberOfPoints)
			{
    			selection = which(dists[[(i-1)]] > minDist)
    			if (length(selection) == 0)
    				{
    					stop("Restarts the function with a smaller minimum distance")
					}
    			selection_list[[i]] = selection
    			test = table(unlist(selection_list))
    			indices_minDist = as.numeric(names(which(test==i)))
    			indices[i] = sample(indices_minDist, 1)   
				dists[[i]] = spDistsN1(as.matrix(observations), as.matrix(observations[indices[i],]), longlat=T)
			}
		return(indices)
	}
foldSelection = function(observations, selectedPoints)
	{
		fold_selection = sapply(1:nrow(observations), function(i) which.min(spDistsN1(as.matrix(selectedPoints), as.matrix(observations[i,]), longlat=T)))
		return(fold_selection)
	}

newAnalyses = FALSE; spatialCrossValidation1 = FALSE; spatialCrossValidation2 = TRUE; datas = list()
occurrence_data_summary = matrix(nrow=dim(species)[1], ncol=2); row.names(occurrence_data_summary) = species[,1]
colnames(occurrence_data_summary) = c("n_1950_2000","n_1950_2000_filtered")
if (newAnalyses == TRUE) { for (i in 1:dim(species)[1]) {
		rasters_stack = stack(envVariables); observations = observations_list[[i]]
		probas = values(background)[!is.na(values(background))]; n = 1000
		if (n > sum(!is.na(values(background)))) n = sum(!is.na(values(background)))
		pseudo_absences = xyFromCell(background, sample(which(!is.na(values(background))), n, prob=probas, replace=F))
		presences = cbind(observations, rep(1,dim(observations)[1]), raster::extract(rasters_stack, observations))
		absences = cbind(pseudo_absences, rep(0,dim(pseudo_absences)[1]), raster::extract(rasters_stack, pseudo_absences))
		colnames(absences)[1] = "longitude"; colnames(absences)[2] = "latitude"; colnames(absences)[3] = "response"
		colnames(presences)[3] = "response"; data = rbind(presences,absences); data_to_discard = c()
		for (j in 1:length(rasters_stack@layers))
			{
				data_to_discard = c(data_to_discard, which(is.na(raster::extract(rasters_stack[[j]],data[,1:2]))))
			}
		data_to_discard = data_to_discard[order(data_to_discard)]
		data = data[which(!c(1:dim(data)[1])%in%data_to_discard),]
		occurrence_data_summary[i,1] = sum(data[,"response"])
		cellIDs = cellFromXY(rasters_stack[[1]], data[,1:2]); buffer = c()
		for (j in 1:length(unique(cellIDs)))
			{	# Keeping only one presence or pseudo-absence point per raster cell (priority = presence points):
				if (sum(cellIDs==unique(cellIDs)[j]) > 1)
					{
						tmp = data[which(cellIDs==unique(cellIDs)[j]),]
						if (sum(tmp[,"response"]==1) == 0)
							{
								buffer = rbind(buffer, tmp[sample(1:dim(tmp)[1],1),])
							}
						if (sum(tmp[,"response"]==1) == 1)
							{
								buffer = rbind(buffer, tmp[which(tmp[,"response"]==1),])
							}
						if (sum(tmp[,"response"]==1) >= 2)
							{
								indices = which(tmp[,"response"]==1)
								buffer = rbind(buffer, tmp[sample(indices,1),])
							}
					}	else	{
						buffer = rbind(buffer, data[which(cellIDs==unique(cellIDs)[j]),])
					}
			}
		data = buffer; datas[[i]] = data
		occurrence_data_summary[i,2] = sum(data[,"response"]) }
		if (!file.exists(paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_CCV_SCV_AUCs.csv"))) {
		plottingCorrelogram = FALSE
		if (plottingCorrelogram == TRUE)
			{
				correlogram = ncf::correlog(data[,"longitude"], data[,"latitude"], data[,"response"], na.rm=T, increment=10, resamp=0, latlon=T)
				dev.new(width=4.5, height=3); par(mar=c(2.2,2.2,1.5,1.5))
				plot(correlogram$mean.of.class[-1], correlogram$correlation[-1], ann=F, axes=F, lwd=0.2, cex=0.5, col=NA, ylim=c(-0.4,1.0), xlim=c(0,8500))
				abline(h=0, lwd=0.5, col="red", lty=2)
				points(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, cex=0.35, col="gray30")
				lines(correlogram$mean.of.class[-1], correlogram$correlation[-1], lwd=0.2, col="gray30")
				axis(side=1, pos=-0.4, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,-0.05,0), at=seq(0,9000,1000))
				axis(side=2, pos=0, lwd.tick=0.2, cex.axis=0.6, lwd=0.2, tck=-0.015, col.axis="gray30", mgp=c(0,0.18,0), at=seq(-0.4,1,0.2))
				title(xlab="distance (km2)", cex.lab=0.7, mgp=c(0.3,0,0), col.lab="gray30")
				title(ylab="correlation", cex.lab=0.7, mgp=c(0.4,0,0), col.lab="gray30")
			}
		theRanges = c(700,700)*1000 # distance in meters
		nberOfReplicates = 10 # one replicate = one folds partition
		gbm.x = names(rasters_stack)
		gbm.y = colnames(data)[3]
		offset = NULL
		tree.complexity = 5 # "tc" = number of nodes in the trees
		learning.rate = 0.005 # "lr" = contribution of each tree to the growing model
		bag.fraction = 0.80 # proportion of data used to train a given tree
		site.weights = rep(1, dim(data)[1])
		var.monotone = rep(0, length(gbm.x))
		n.folds = 5
		prev.stratify = TRUE
		family = "bernoulli"
		n.trees = 10 # initial number of trees
		step.size = 5 # interval at which the predictive deviance is computed and logged
					  # (at each interval, the folds are successively used as test data set
					  # nd the remaining folds as training data sets to compute the deviance)
		max.trees = 10000 # maximum number of trees that will be considered
		tolerance.method = "auto"
		tolerance = 0.001
		plot.main = TRUE
		plot.folds = FALSE
		verbose = TRUE
		silent = FALSE
		keep.fold.models = FALSE
		keep.fold.vector = FALSE
		keep.fold.fit = FALSE
		showingFoldsPlot = FALSE
		brt_model_ccvs = list() # classic cross-validations (CCVs)
		brt_model_scv1 = list() # spatial cross-validations 1 (SCV1)
		brt_model_scv2 = list() # spatial cross-validations 2 (SCV2)
		if ((spatialCrossValidation1 == TRUE)&(spatialCrossValidation2 == TRUE))	
			{
				AUCs = matrix(nrow=nberOfReplicates, ncol=3); colnames(AUCs) = c("CCV_AUC","SCV1_AUC","SCV2_AUC")
			}	else	{
				if (sum(c(spatialCrossValidation1,spatialCrossValidation2)) == 1)
					{
						AUCs = matrix(nrow=nberOfReplicates, ncol=2); colnames(AUCs) = c("CCV_AUC","SCV_AUC")
					}
				if (sum(c(spatialCrossValidation1,spatialCrossValidation2)) == 0)
					{
						AUCs = matrix(nrow=nberOfReplicates, ncol=1); colnames(AUCs) = c("AUC")
					}
			}
		for (j in 1:nberOfReplicates)
			{
				# BRT with classic (standard) cross-validation (CCV):
				pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_CCV_replicate_",j,".pdf"))
				n.trees = 10; learning.rate = 0.005; step.size = 5; fold.vector = NULL; worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								brt_model_ccvs[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
									var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
									verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
							},	error = function(cond) {
							},	warning = function(cond) {
							},	finally = {
							})
						if (length(brt_model_ccvs) == j) worked = TRUE
					}
				dev.off()
				AUCs[j,1] = brt_model_ccvs[[j]]$cv.statistics$discrimination.mean # Mean test AUC (from the AUCs computed on each fold tested as test data in the CCV)
				object = brt_model_ccvs[[j]]; df = as.data.frame(rasters_stack)
				not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
				n.trees = brt_model_ccvs[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
				prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
				rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
				if (spatialCrossValidation1 == TRUE)
					{			
						# BRT with spatial (geographic) cross-validation (SCV) based on the folds generation of Dhingra, Artois et al. (2016, eLife):
						folds_with_similar_sizes = FALSE; c = 0
						while (folds_with_similar_sizes == FALSE) # while loop to select a partition where the x folds gather at least
							{									  #  proportion = (1/(x+1)) of the total number of presence points
								data_presence = data[which(data[,3]==1),]; c = c+1; # print(c)
								fivePoints = samplingPtsMinDist(data_presence[,1:2], minDist=200, nberOfPoints=n.folds)
								fold.vector = foldSelection(data[,1:2], selectedPoints=data_presence[fivePoints,1:2])
								fold.vector_presences = fold.vector[which(data[,3]==1)]
								counts = hist(fold.vector_presences, plot=F)$counts
								props = counts[which(counts > 0)]/sum(counts); print(round(props,2))
								if (min(props) > (1/(n.folds*2))) folds_with_similar_sizes = TRUE
							}
						if (showingFoldsPlot == TRUE)
							{
								par(mar=c(0,0,0,0), oma=c(0.0,3.6,0.0,0.0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
								cols = c("olivedrab3","tan3","steelblue3","orange1","tomato2","mediumseagreen")[fold.vector]
								plot(background, col="gray90", useRaster=T, colNA=NA, box=F, axes=F, legend=F)
								pchs = c(16,3)[data[,3]+1]; cexs = c(0.25,0.5)[data[,3]+1]
								points(data[,1:2], col=cols, pch=pchs, cex=cexs, lwd=0.7)
							}
						if (spatialCrossValidation2 == TRUE) pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_SCV1_replicate_",j,".pdf"))
						if (spatialCrossValidation2 == FALSE) pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_SCV_replicate_",j,".pdf"))
						n.trees = 10; learning.rate = 0.0001; step.size = 2; worked = FALSE
						while (worked == FALSE)
							{
								trycatch = tryCatch(
									{
										brt_model_scv1[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
											var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
											verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit)
									},	error = function(cond) {
									},	warning = function(cond) {
									},	finally = {
									})
								if (length(brt_model_scv1) == j) worked = TRUE
							}
						dev.off()
						if (spatialCrossValidation2 == TRUE) AUCs[j,"SCV1_AUC"] = brt_model_scv1[[j]]$cv.statistics$discrimination.mean
						if (spatialCrossValidation2 == FALSE) AUCs[j,"SCV_AUC"] = brt_model_scv1[[j]]$cv.statistics$discrimination.mean
							# Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)		
						object = brt_model_scv1[[j]]; df = as.data.frame(rasters_stack)
						not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
						n.trees = brt_model_scv1[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
					}
				if (spatialCrossValidation2 == TRUE)
					{
						# BRT with spatial (geographic) cross-validation (SCV) based on the blocks generation of Valavi et al. (2019, MEE):
						spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,3:dim(data)[2]], proj4string=crs(nullRaster))
						worked = FALSE
						while (worked == FALSE)
							{
								trycatch = tryCatch(
									{
										myblocks = NULL
										myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRaster, k=n.folds, theRange=theRanges[1], selection="random")
									},	error = function(cond) {
									},	finally = {
									})
								if (!is.null(myblocks)) worked = TRUE
							}
						fold.vector = myblocks$foldID
						pdf(file=paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_SCV2_replicate_",j,".pdf"))
						n.trees = 10; learning.rate = 0.005; step.size = 5
						brt_model_scv2[[j]] = gbm.step(data, gbm.x, gbm.y, offset, fold.vector, tree.complexity, learning.rate, bag.fraction, site.weights,
							var.monotone, n.folds, prev.stratify, family, n.trees, step.size, max.trees, tolerance.method, tolerance, plot.main, plot.folds,
							verbose, silent, keep.fold.models, keep.fold.vector, keep.fold.fit); # summary(brt_model_scv) # gbm.plot(brt_model_scv, plot.layout=c(4,4))
						dev.off()
						if (spatialCrossValidation1 == TRUE) AUCs[j,"SCV2_AUC"] = brt_model_scv2[[j]]$cv.statistics$discrimination.mean
						if (spatialCrossValidation1 == FALSE) AUCs[j,"SCV_AUC"] = brt_model_scv2[[j]]$cv.statistics$discrimination.mean
							# Mean test AUC (from the AUCs computed on each fold tested as test data in the SCV)
						object = brt_model_scv2[[j]]; df = as.data.frame(rasters_stack)
						not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
						n.trees = brt_model_scv2[[j]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(object, newdata, n.trees, type, single.tree)
						rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction
					}
			}
		saveRDS(brt_model_ccvs, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_models_CCV.rds"))
		if (spatialCrossValidation1 == TRUE)	
			{
				if (spatialCrossValidation2 == TRUE)	
					{
						saveRDS(brt_model_scv1, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_models_SCV1.rds"))
					}	else		{
						saveRDS(brt_model_scv1, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_models_SCV.rds"))
					}
			}
		if (spatialCrossValidation2 == TRUE) saveRDS(brt_model_scv2, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_models_SCV2.rds"))
		write.csv(AUCs, paste0("BRT_projection_files/BRT_models/B_",species[i,1],"_CCV_SCV_AUCs.csv"), row.names=F, quote=F)
	}}}
if (!file.exists(paste0("Occurrence_data.csv")))
	{
		write.csv(occurrence_data_summary, "Occurrence_data.csv", quote=F)
	}
AUC_values = matrix(nrow=dim(species)[1], ncol=2); row.names(AUC_values) = species[,"species"]
colnames(AUC_values) = c("CCV","SCV")
for (i in 1:dim(species)[1])
	{
		tab = read.csv(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_CCV_SCV_AUCs.csv"), head=T)
		for (j in 1:dim(tab)[2])
			{
				AUC_values[i,j] = paste0(round(mean(tab[,j]),3)," (",round(sd(tab[,j]),3),")")
			}
	}
if (!file.exists(paste0("All_AUC_values.csv")))
	{
		write.csv(AUC_values, "All_AUC_values.csv", quote=F)
	}
AUC_values = read.csv("All_AUC_values.csv", head=T)

# 5. Estimation of spatial sorting bias (SSB; for spatial autocorrelation)

	# SSB = Dp/Da (Hijsmans 2012, Ecology), where:
		# Dp = mean distance between testing presence sites and nearest training-presence sites
		# Da = mean distance between testing absence sites and nearest training-presence sites
		# --> SSB = 1 suggests there is no spatial sorting bias
		# --> SSB = 0 suggests extreme spatial sorting bias

newAnalyses = FALSE
if (newAnalyses == TRUE)
	{
		n.folds = 5; theRanges = c(500,500)*1000
		SSBs = matrix(ncol=n.folds, nrow=dim(species)[1])
		row.names(SSBs) = species[,"species"]
		colnames(SSBs) = c("fold1","fold2","fold3","fold4","fold5")
		for (i in 1:dim(species)[1])
			{
				data = datas[[i]]
				fold.vector = rep(NA, dim(data)[1])
				for (j in 1:dim(data)[1]) fold.vector[j] = sample(1:n.folds,1)
				for (j in 1:n.folds)
					{
						p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
						a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
						reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
						SSB = ssb(p, a, reference); SSBs[i,j] = SSB[1,"p"]/SSB[1,"a"]
					}
			}
		write.csv(round(SSBs,2), paste0("BRT_projection_files/BRT_models/SSB_estimates_for_the_CCV.csv"), quote=F)
		SSBs = matrix(ncol=n.folds, nrow=dim(species)[1])
		row.names(SSBs) = species[,"species"]
		colnames(SSBs) = c("fold1","fold2","fold3","fold4","fold5")
		for (i in 1:dim(species)[1])
			{
				folds_with_similar_sizes = FALSE; c = 0
				spdf = SpatialPointsDataFrame(data[c("longitude","latitude")], data[,3:dim(data)[2]], proj4string=crs(nullRaster))
				worked = FALSE
				while (worked == FALSE)
					{
						trycatch = tryCatch(
							{
								myblocks = NULL
								myblocks = spatialBlock(spdf, species="response", rasterLayer=nullRaster, k=n.folds, theRange=theRanges[1], selection="random")
							},	error = function(cond) {
							},	finally = {
							})
						if (!is.null(myblocks)) worked = TRUE
					}
				for (j in 1:n.folds)
					{
						fold.vector = myblocks$foldID
						p = data[which((data[,"response"]==1)&(fold.vector!=j)),1:2]
						a = data[which((data[,"response"]==0)&(fold.vector!=j)),1:2]
						reference = data[which((data[,"response"]==1)&(fold.vector==j)),1:2]
						if (dim(reference)[1]>0)
							{
								SSB = ssb(p, a, reference); SSBs[i,j] = SSB[1,"p"]/SSB[1,"a"]
							}
					}
			}
		write.csv(round(SSBs,2), paste0("BRT_projection_files/BRT_models/SSB_estimates_for_the_SCV.csv"), quote=F)
	}

# 6. Comparison of the response curves for each climatic variable

newAnalyses = FALSE
if (newAnalyses == TRUE)
	{
		selectedModel = "SCV2"
		envVariableValues_list = list()
		for (i in 1:dim(species)[1])
			{
				data = datas[[i]]; data = data[which(data[,"response"]==1),]
				envVariableValues = matrix(nrow=3, ncol=length(envVariables))
				row.names(envVariableValues) = c("median","minV","maxV")
				colnames(envVariableValues) = envVariableNames2
				for (j in 1:length(envVariables))
					{
						minV = min(data[,envVariableNames2[j]], na.rm=T)
						maxV = max(data[,envVariableNames2[j]], na.rm=T)
						medianV = median(data[,envVariableNames2[j]], na.rm=T)
						envVariableValues[,j] = cbind(medianV, minV, maxV)
					}
				envVariableValues_list[[i]] = envVariableValues
			}
		pdf(paste0("All_the_figures_&_SI/All_response_curves_NEW.pdf"), width=7.5, height=8)
		par(mfrow=c(5,4), oma=c(1,1,1,1), mar=c(2,1.3,0.5,0.5), lwd=0.2, col="gray30")
		for (i in 1:length(envVariables))
			{
				predictions = list(); dfs = list()
				for (j in 1:dim(species)[1])
					{
						valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[j]]["maxV",i]-envVariableValues_list[[j]]["minV",i])/100
						df = data.frame(matrix(nrow=length(seq(envVariableValues_list[[j]]["minV",i],envVariableValues_list[[j]]["maxV",i],valuesInterval)),
											   ncol=length(envVariables))); colnames(df) = envVariableNames2
						for (k in 1:length(envVariables))
							{
								valuesInterval = 0.1; valuesInterval = (envVariableValues_list[[j]]["maxV",k]-envVariableValues_list[[j]]["minV",k])/100
								if (i == k) df[,envVariableNames2[k]] = seq(envVariableValues_list[[j]]["minV",k],envVariableValues_list[[j]]["maxV",k],valuesInterval)
								if (i != k) df[,envVariableNames2[k]] = rep(envVariableValues_list[[j]]["median",k],dim(df)[1])
							}
						dfs[[j]] = df
						brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[j,"species"],"_models_",selectedModel,".rds"))
						AUC_values = read.csv(paste0("BRT_projection_files/BRT_models/B_",species[j,"species"],"_CCV_SCV_AUCs.csv"))[,paste0(gsub("2","",selectedModel),"_AUC")]
						index = which(AUC_values==max(AUC_values))[1]
						n.trees = brt_model[[index]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
						prediction = predict.gbm(brt_model[[index]], newdata=df, n.trees, type, single.tree)
						if (j == 1)
							{
								minX = min(df[,envVariableNames2[i]]); maxX = max(df[,envVariableNames2[i]])
								minY = min(prediction); maxY = max(prediction)
							}	else	{
								if (minX > min(df[,envVariableNames2[i]])) minX = min(df[,envVariableNames2[i]])
								if (maxX < max(df[,envVariableNames2[i]])) maxX = max(df[,envVariableNames2[i]])
								if (minY > min(prediction)) minY = min(prediction)
								if (maxY < max(prediction)) maxY = max(prediction)
							}
						predictions[[j]] = prediction
					}
				for (j in 1:dim(species)[1])
					{
						if (j == 1)
							{
								plot(dfs[[j]][,envVariableNames2[i]],predictions[[j]],col="tomato2",ann=F,axes=F,lwd=0.2,type="l",xlim=c(minX,maxX),ylim=c(minY,maxY))
							}	else	{
								lines(dfs[[j]][,envVariableNames2[i]],predictions[[j]],col="tomato2",lwd=0.2)
							}
					}
				axis(side=1, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.07,0))
				axis(side=2, lwd.tick=0.2, cex.axis=0.6, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.30,0))
				title(ylab="predicted values", cex.lab=0.8, mgp=c(1.3,0,0), col.lab="gray30")
				title(xlab=envVariableNames2[i], cex.lab=0.8, mgp=c(0.9,0,0), col.lab="gray30")
				box(lwd=0.2, col="gray30")
			}
		dev.off()
	}

# 7. Analyses of the relative influence of each climatic variable

selectedModel = "SCV2"
if (!file.exists("BRT_RI_estimates.csv"))
	{
		relativeInfluences = matrix(0, nrow=dim(species)[1], ncol=length(envVariables))
		row.names(relativeInfluences) = species[,"species"]
		for (i in 1:dim(species)[1])
			{
				brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_models_",selectedModel,".rds"))
				for (j in 1:length(brt_model))
					{
						if ((i == 1)&(j == 1)) envVariableNames2 = rep(NA, length(envVariables))
						for (k in 1:length(envVariables))
							{
								if ((i == 1)&(j == 1)) envVariableNames2[k] = names(envVariables[[k]])
								relativeInfluences[i,k] = relativeInfluences[i,k] + summary(brt_model[[j]])[names(envVariables[[k]]),"rel.inf"]
							}
					}
				if (i == 1) colnames(relativeInfluences) = envVariableNames2
				relativeInfluences[i,] = relativeInfluences[i,]/length(brt_model)
			}
		write.table(round(relativeInfluences,2), "BRT_RI_estimates.csv", quote=F, sep=",")
	}
meanRelativeInfluences = matrix(0, nrow=1, ncol=length(envVariables)); c = 0
relativeInfluences = read.csv("BRT_RI_estimates.csv", header=T)
for (i in 1:dim(relativeInfluences)[1])
	{
		if (i == 1) colnames(meanRelativeInfluences) = colnames(relativeInfluences)
		meanRelativeInfluences = meanRelativeInfluences + relativeInfluences[i,]; c = c+1
	}
meanRelativeInfluences = meanRelativeInfluences/c; row.names(meanRelativeInfluences) = NULL
# print(t(round(meanRelativeInfluences,1)))

# 8. Plotting the different climatic variables and their past projections

envVariables_t0 = list()
for (i in 1:length(envVariableNames1))
	{
		envVariable = raster(paste0("WorldClim_1_rasters/Present_t0/bio",i,".bil"))
		envVariable[is.na(mask[])] = NA; envVariable = crop(envVariable, extent(-13,29,35,72))
		names(envVariable) = envVariableNames2[i]; envVariables_t0[[i]] = envVariable
	}
envVariables_list = list(); envVariables_list[[1]] = envVariables_t0
periods = c("t0","mid-Holocene","LGM"); directories = c("Present_t0","Past_Holo","Past_LGM")
period_names = c("Present times","Mid-Holocene","LGM"); model_names = c()
for (i in 1: length(models)) model_names[i] = unlist(strsplit(models[i],"_"))[1]
for (i in 2:length(periods))
	{
		buffer_list1 = list()
		models = list.files(paste0("WorldClim_1_rasters/",directories[i]))
		models = models[which(!grepl(".pdf",models))]
		for (j in 1:length(models))
			{
				files = list.files(paste0("WorldClim_1_rasters/",directories[i],"/",models[j]))
				buffer_list2 = list()
				for (k in 1:length(envVariableNames2))
					{
						envVariable = raster(paste0("WorldClim_1_rasters/",directories[i],"/",models[j],"/",files[k]))
						names(envVariable) = envVariableNames2[k]; buffer_list2[[k]] = envVariable
					}
				for (k in 1:length(envVariableNames2))
					{
						buffer_list2[[k]] = crop(buffer_list2[[k]], extent(envVariables_list[[1]][[1]]), snap="out")
					}
				buffer_list1[[j]] = buffer_list2
			}
		envVariables_list[[i]] = buffer_list1
	}
if (savingPlots)
	{
		for (i in 1:length(envVariables_list[[1]])) # variables
			{
				if (i <= 11) colour_scale = colorRampPalette(brewer.pal(9,"YlOrRd"))(150)[1:101] # temperature
				if (i >= 12) colour_scale = colorRampPalette(brewer.pal(9,"YlGnBu"))(101) # precipitation
				minV = Inf; maxV = -Inf
				for (j in 1:length(envVariables_list[[2]])) # models
					{
						for (k in length(envVariables_list):1) # periods
							{
								if (k != 1)
									{
										rast = envVariables_list[[k]][[j]][[i]]
									}	else	{
										rast = envVariables_list[[k]][[i]]
									}
								if (minV > min(rast[],na.rm=T)) minV = min(rast[],na.rm=T)
								if (maxV < max(rast[],na.rm=T)) maxV = max(rast[],na.rm=T)
							}
					}
				pdf(paste0("All_the_figures_&_SI/Maps_variable_",envVariableNames2[i],".pdf"), width=7.5, height=10.5)
				par(mfrow=c(3,3), oma=c(1.0,2.0,1.0,0.0), mar=c(0.0,0.0,0.0,0.2), lwd=0.2, col="gray30")
				for (j in 1:length(envVariables_list[[2]])) # models
					{
						for (k in length(envVariables_list):1) # periods
							{
								if ((j != 1)&(k == 1))
									{
										if ((j == 2)&(k == 1))
											{
												rast = envVariables_list[[1]][[1]]; rast[] = NA
												plot(rast, col=NA, ann=F, legend=F, axes=F, box=F)
												rastLegend = raster(t(as.matrix(c(minV,maxV))))
												plot(rastLegend, col=colour_scale, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.00,0.02,0.037,0.963), adj=3,
													 axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.6, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0,
													 mgp=c(0,0.6,0)), alpha=1, side=3, horizontal=F)
											 }	else	{
											 	plot.new()
											 }
									}	else		{
										if (k != 1)
											{
												rast = envVariables_list[[k]][[j]][[i]]
											}	else	{
												rast = envVariables_list[[k]][[i]]
											}
										index1 = round(((min(rast[],na.rm=T)-minV)/(maxV-minV))*100)+1
										index2 = round(((max(rast[],na.rm=T)-minV)/(maxV-minV))*100)+1
										plot(rast, col=colour_scale[index1:index2], ann=F, legend=F, axes=F, box=F)
										rect(-12, 67, 5, 71.5, lwd=0.2, col="white", border=NA)
										mtext(period_names[k], at=-10, line=-3, adj=0, cex=0.65, col="gray30")
										mtext(paste0("(",model_names[j],")"), at=-10, line=-4, adj=0, cex=0.6, col="gray30")
										rect(-13, 35, 29, 72, lwd=0.2, col=NA, border="gray30")
									}
							}
					}
				dev.off()
			}
	}

# 9. BRT projections based on past and present climatic variables

if (!file.exists("BRT_projection_files/All_projs.rds"))
	{
		projections1 = list(); selectedModel = "SCV2"
		for (i in 1:dim(species)[1])
			{
				projections2 = list()
				brt_model = readRDS(paste0("BRT_projection_files/BRT_models/B_",species[i,"species"],"_models_",selectedModel,".rds"))
				speciesName = as.character(species[i,"species"])
				for (j in 1:length(envVariables_list[[2]])) # models
					{
						projections3 = list()
						for (k in length(envVariables_list):1) # periods
							{
								if ((k == 1)&(j > 1))
									{
									}	else		{
										if (k == 1)
											{
												rasters_stack = stack(envVariables_list[[k]])
											}	else	{
												rasters_stack = stack(envVariables_list[[k]][[j]])
											}
										replicates = list()
										for (l in 1:length(brt_model))
											{
												df = as.data.frame(rasters_stack); not_NA = which(!is.na(rowMeans(df))); newdata = df[not_NA,]
												n.trees = brt_model[[l]]$gbm.call$best.trees; type = "response"; single.tree = FALSE
												prediction = predict.gbm(brt_model[[l]], newdata, n.trees, type, single.tree)
												rast = rasters_stack[[1]]; rast[!is.na(rast[])] = prediction; replicates[[l]] = rast
											}
										rasts = stack(replicates); projections3[[k]] = mean(rasts)
									}
							}
						projections2[[j]] = projections3
					}
				saveRDS(projections2, paste0("BRT_projection_files/Projections/",species[i,"species"],".rds"))				
				projections1[[i]] = projections2
			}
		saveRDS(projections1, "BRT_projection_files/All_projs.rds")
	}	else	{
		projections1 = readRDS("BRT_projection_files/All_projs.rds")
	}
if (savingPlots)
	{
		minV = 0; maxV = 1
		for (i in 1:length(projections1)) # species
			{
				colour_scale = c(rep("#E5E5E5",10),rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120))[21:110])
				pdf(paste0("All_the_figures_&_SI/All_BRT_projection_maps/",species[i,"species"],".pdf"), width=7.5, height=10.5)
				par(mfrow=c(3,3), oma=c(1.0,2.0,1.0,0.0), mar=c(0.0,0.0,0.0,0.2), lwd=0.2, col="gray30")
				for (j in 1:length(projections1[[i]][[2]])) # models
					{
						for (k in length(projections1[[i]]):1) # periods
							{
								if ((j != 1)&(k == 1))
									{
										if ((j == 2)&(k == 1))
											{
												rast = projections1[[i]][[1]][[1]]; rast[] = NA
												plot(rast, col=NA, ann=F, legend=F, axes=F, box=F)
												rastLegend = raster(t(as.matrix(c(minV,maxV))))
												plot(rastLegend, col=colour_scale, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.00,0.02,0.037,0.963), adj=3,
													 axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.6, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0,
													 mgp=c(0,0.6,0)), alpha=1, side=3, horizontal=F)
											 }	else	{
											 	plot.new()
											 }
									}	else		{
										rast = projections1[[i]][[j]][[k]]
										index1 = round(((min(rast[],na.rm=T)-minV)/(maxV-minV))*100)+1
										index2 = round(((max(rast[],na.rm=T)-minV)/(maxV-minV))*100)+1
										plot(rast, col=colour_scale[index1:index2], ann=F, legend=F, axes=F, box=F)
										rect(-12, 67, 5, 71.5, lwd=0.2, col="white", border=NA)
										mtext(period_names[k], at=-10, line=-3, adj=0, cex=0.65, col="gray30")
										mtext(paste0("(",model_names[j],")"), at=-10, line=-4, adj=0, cex=0.6, col="gray30")
										rect(-13, 35, 29, 72, lwd=0.2, col=NA, border="gray30")
									}
							}
					}
				dev.off()
			}
	}

# 10. Post hoc analyses to compare past and present climatic niches

cutOff = 0.10
rasters_CSI_list = list() # CSI = climatic suitability index
rasters_SRI_list = list() # SRI = species richness index
for (j in 1:length(models))
	{
		rasters_CSI = list(); rasters_SRI = list()
		bufferRasters = list(); counterRasters = list()
		for (i in 1:dim(species)[1])
			{
				bufferRasters[[i]] = projections1[[i]][[1]][[1]]
				c1 = projections1[[i]][[1]][[1]]; c1[c1[]>cutOff] = 1
				counterRasters[[i]] = c1
			}
		rasters_CSI[[1]] = mean(stack(bufferRasters)); rasters_SRI[[1]] = sum(stack(counterRasters))
		bufferRasters = list(); counterRasters = list()
		for (i in 1:dim(species)[1])
			{
				bufferRasters[[i]] = projections1[[i]][[j]][[3]]
				c1 = projections1[[i]][[j]][[3]]; c1[c1[]>cutOff] = 1
				counterRasters[[i]] = c1
			}
		rasters_CSI[[2]] = mean(stack(bufferRasters)); rasters_SRI[[2]] = sum(stack(counterRasters))
		rasters_CSI_list[[j]] = rasters_CSI; rasters_SRI_list[[j]] = rasters_SRI
	}
if (savingPlots)
	{
		colour_scale_1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120))[11:110]
		colour_scale_2 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120))[11:110]
		minV1 = Inf; maxV1 = -Inf; minV2 = Inf; maxV2 = -Inf
		for (i in 1:length(rasters_CSI_list))
			{
				for (j in 1:length(rasters_CSI_list[[i]]))
					{
						if (minV1 > min(rasters_CSI_list[[i]][[j]][],na.rm=T)) minV1 = min(rasters_CSI_list[[i]][[j]][],na.rm=T)
						if (maxV1 < max(rasters_CSI_list[[i]][[j]][],na.rm=T)) maxV1 = max(rasters_CSI_list[[i]][[j]][],na.rm=T)
					}
			}
		for (i in 1:length(rasters_SRI_list))
			{
				for (j in 1:length(rasters_SRI_list[[i]]))
					{
						if (minV2 > min(rasters_SRI_list[[i]][[j]][],na.rm=T)) minV2 = min(rasters_SRI_list[[i]][[j]][],na.rm=T)
						if (maxV2 < max(rasters_SRI_list[[i]][[j]][],na.rm=T)) maxV2 = max(rasters_SRI_list[[i]][[j]][],na.rm=T)
					}
			}
		pdf(paste0("All_the_figures_&_SI/CSI_maps_LGM_&_t0_NEW.pdf"), width=7.5, height=10.5)
		par(mfrow=c(3,3), oma=c(1.0,2.0,1.0,0.0), mar=c(0.0,0.0,0.0,0.2), lwd=0.2, col="gray30")
		for (j in length(rasters_CSI_list[[i]]):1)
			{
				for (i in 1:length(rasters_CSI_list))
					{
						rast = rasters_CSI_list[[i]][[j]]
						index1 = round(((min(rast[],na.rm=T)-minV1)/(maxV1-minV1))*100)+1
						index2 = round(((max(rast[],na.rm=T)-minV1)/(maxV1-minV1))*100)+1
						plot(rast, col=colour_scale_1[index1:index2], ann=F, legend=F, axes=F, box=F)
						rect(-12, 67, 5, 71.5, lwd=0.2, col="white", border=NA)
						mtext("Climatic suitability index", at=-10, line=-2.8, adj=0, cex=0.55, col="gray30")
						mtext(paste0("(",model_names[i],")"), at=-10, line=-3.8, adj=0, cex=0.55, col="gray30")
						rect(-13, 35, 29, 72, lwd=0.2, col=NA, border="gray30")
					}
			}
		rast = projections1[[i]][[1]][[1]]; rast[] = NA
		plot(rast, col=NA, ann=F, legend=F, axes=F, box=F)
		rastLegend = raster(t(as.matrix(c(minV1,maxV1))))
		plot(rastLegend, col=colour_scale_1, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.00,0.02,0.037,0.963), adj=3,
			 axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.6, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0,
			 mgp=c(0,0.6,0)), alpha=1, side=3, horizontal=F)
		dev.off()
		pdf(paste0("All_the_figures_&_SI/SRI_maps_LGM_&_t0_NEW.pdf"), width=7.5, height=10.5)
		par(mfrow=c(3,3), oma=c(1.0,2.0,1.0,0.0), mar=c(0.0,0.0,0.0,0.2), lwd=0.2, col="gray30")
		for (j in length(rasters_SRI_list[[i]]):1)
			{
				for (i in 1:length(rasters_SRI_list))
					{
						rast = rasters_SRI_list[[i]][[j]]
						index1 = round(((min(rast[],na.rm=T)-minV2)/(maxV2-minV2))*100)+1
						index2 = round(((max(rast[],na.rm=T)-minV2)/(maxV2-minV2))*100)+1
						plot(rast, col=colour_scale_2[index1:index2], ann=F, legend=F, axes=F, box=F)
						rect(-12, 67, 5, 71.5, lwd=0.2, col="white", border=NA)
						mtext("Species richness index", at=-10, line=-2.8, adj=0, cex=0.55, col="gray30")
						mtext(paste0("(",model_names[i],")"), at=-10, line=-3.8, adj=0, cex=0.55, col="gray30")
						rect(-13, 35, 29, 72, lwd=0.2, col=NA, border="gray30")
					}
			}		
		rast = projections1[[i]][[1]][[1]]; rast[] = NA
		plot(rast, col=NA, ann=F, legend=F, axes=F, box=F)
		rastLegend = raster(t(as.matrix(c(minV2,maxV2))))
		plot(rastLegend, col=colour_scale_2, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.00,0.02,0.037,0.963), adj=3,
			 axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.6, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0,
			 mgp=c(0,0.6,0)), alpha=1, side=3, horizontal=F)
		dev.off()
	}
if (savingPlots)
	{
		pdf(paste0("All_the_figures_&_SI/CSI_SRI_differences_NEW.pdf"), width=7.5, height=10.5)
		par(mfrow=c(3,3), oma=c(1.0,2.0,1.0,0.0), mar=c(0.0,0.0,0.0,0.2), lwd=0.2, col="gray30")
		colour_scale_1 = colorRampPalette(brewer.pal(11,"BrBG"))(120)[11:110]
		colour_scale_2 = colorRampPalette(brewer.pal(11,"PuOr"))(120)[11:110]
		colour_scale_1 = colorRampPalette(brewer.pal(11,"RdYlGn"))(120)[11:110]
		colour_scale_1 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120))[11:110]
		colour_scale_2 = rev(colorRampPalette(brewer.pal(11,"RdYlBu"))(120))[11:110]
		minV1 = Inf; maxV1 = -Inf; minV2 = Inf; maxV2 = -Inf
		for (i in 1:length(rasters_CSI_list))
			{
				if (minV1 > min(rasters_CSI_list[[i]][[1]][]-rasters_CSI_list[[i]][[2]][],na.rm=T)) minV1 = min(rasters_CSI_list[[i]][[1]][]-rasters_CSI_list[[i]][[2]][],na.rm=T)
				if (maxV1 < max(rasters_CSI_list[[i]][[1]][]-rasters_CSI_list[[i]][[2]][],na.rm=T)) maxV1 = max(rasters_CSI_list[[i]][[1]][]-rasters_CSI_list[[i]][[2]][],na.rm=T)
			}
		for (i in 1:length(rasters_SRI_list))
			{
				if (minV2 > min(rasters_SRI_list[[i]][[1]][]-rasters_SRI_list[[i]][[2]][],na.rm=T)) minV2 = min(rasters_SRI_list[[i]][[1]][]-rasters_SRI_list[[i]][[2]][],na.rm=T)
				if (maxV2 < max(rasters_SRI_list[[i]][[1]][]-rasters_SRI_list[[i]][[2]][],na.rm=T)) maxV2 = max(rasters_SRI_list[[i]][[1]][]-rasters_SRI_list[[i]][[2]][],na.rm=T)
			}
		if (abs(minV1) < abs(maxV1)) minV1 = -maxV1
		if (abs(maxV1) < abs(minV1)) maxV1 = -minV1
		if (abs(minV2) < abs(maxV2)) minV2 = -maxV2
		if (abs(maxV2) < abs(minV2)) maxV2 = -minV2
		for (i in 1:length(rasters_CSI_list))
			{
				rast = rasters_CSI_list[[i]][[1]]-rasters_CSI_list[[i]][[2]]
				index1 = round(((min(rast[],na.rm=T)-minV1)/(maxV1-minV1))*100)+1
				index2 = round(((max(rast[],na.rm=T)-minV1)/(maxV1-minV1))*100)+1
				plot(rast, col=colour_scale_1[index1:index2], ann=F, legend=F, axes=F, box=F)
				rect(-12, 67, 5, 71.5, lwd=0.2, col="white", border=NA)
				mtext("Climatic suitability index", at=-10, line=-3, adj=0, cex=0.65, col="gray30")
				mtext(paste0("(",model_names[i],")"), at=-10, line=-4, adj=0, cex=0.55, col="gray30")
				rect(-13, 35, 29, 72, lwd=0.2, col=NA, border="gray30")
			}
		for (i in 1:length(rasters_SRI_list))
			{
				rast = rasters_SRI_list[[i]][[1]]-rasters_SRI_list[[i]][[2]]
				index1 = round(((min(rast[],na.rm=T)-minV2)/(maxV2-minV2))*100)+1
				index2 = round(((max(rast[],na.rm=T)-minV2)/(maxV2-minV2))*100)+1
				plot(rast, col=colour_scale_2[index1:index2], ann=F, legend=F, axes=F, box=F)
				rect(-12, 67, 5, 71.5, lwd=0.2, col="white", border=NA)
				mtext("Species richness index", at=-10, line=-3, adj=0, cex=0.65, col="gray30")
				mtext(paste0("(",model_names[i],")"), at=-10, line=-4, adj=0, cex=0.55, col="gray30")
				rect(-13, 35, 29, 72, lwd=0.2, col=NA, border="gray30")
			}		
		rast = projections1[[i]][[1]][[1]]; rast[] = NA
		plot(rast, col=NA, ann=F, legend=F, axes=F, box=F)
		rastLegend = raster(t(as.matrix(c(minV1,maxV1))))
		plot(rastLegend, col=colour_scale_1, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.00,0.02,0.037,0.963), adj=3,
			 axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.6, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0,
			 mgp=c(0,0.6,0)), alpha=1, side=3, horizontal=F)
		plot(rast, col=NA, ann=F, legend=F, axes=F, box=F)
		rastLegend = raster(t(as.matrix(c(minV2,maxV2))))
		plot(rastLegend, col=colour_scale_2, legend.only=T, add=T, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.00,0.02,0.037,0.963), adj=3,
			 axis.args=list(cex.axis=1.0, lwd=0, lwd.tick=0.6, col.tick="gray30", tck=-0.8, col="gray30", col.axis="gray30", line=0,
			 mgp=c(0,0.6,0)), alpha=1, side=3, horizontal=F)
		dev.off()
	}

