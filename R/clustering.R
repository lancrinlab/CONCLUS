.plotDistanceGraphWithEpsilon <- function(tSNEData, minNeighbours=5,
		epsilon=1.2){
	# similar function as plotDistanceGraph,
	# but with already known epsilon value
	#
	
	dbscan::kNNdistplot(tSNEData, k=minNeighbours)
	abline(h=epsilon, lty=2)
	
}


.plotTestClustering <- function(tSNEData, minNeighbours=5,
		epsilon=1.2){
	# plots test DBSCAN on one of the pictures
	# for being ensured that clustering will
	# probably work successful
	
	dbscanResults <- fpc::dbscan(tSNEData, eps=epsilon, MinPts=minNeighbours)
	factoextra::fviz_cluster(dbscanResults, tSNEData, ellipse=TRUE, geom="point",
			legend="bottom")
}

testClustering <- function(sceObject, dataDirectory, experimentName,
		dbscanEpsilon=1.4,
		minPts=5,
		perplexities = c(30), PCs = c(4),
		randomSeed = 42,
		width=7, height=7, onefile=FALSE, ...){
	
	initialisePath(dataDirectory)
	dir.create(file.path(dataDirectory, "test_clustering"), showWarnings = F)
	
	message("Generating TSNE.")
	#1. Generating 2D tSNE plots
	tSNE <- .getTSNEresults(expr = Biobase::exprs(sceObject), cores=1,
			perplexities = perplexities,
			PCs = PCs,
			randomSeed = randomSeed)
	
	message("Saving results.")
	#picture checker
	pdf(file.path(dataDirectory, "test_clustering", "test_tSNE.pdf"),
			width=width, height=height, onefile=onefile, ...)
	print(tSNE)
	dev.off()
	
	
	#2. Clustering with dbscan
	# choosing for the best epsilon
	pdf(file.path(dataDirectory, "test_clustering", "distance_graph.pdf"),
			width=width, height=height, onefile=onefile, ...)
	.plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
			minNeighbours = minPts)
	dev.off()
	
	pdf(file.path(dataDirectory, "test_clustering", "test_clustering.pdf"),
			width=width, height=height, onefile=onefile, ...)
	print(.plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
					minNeighbours = minPts))
	dev.off()
	
	message("Pictures of test clustering were exported.")
	return(list(tSNE,
					.plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
							minNeighbours = minPts),
					.plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
							minNeighbours = minPts)))
	
}



.mkOutlierScoreDf <- function(mat){
	
	outlierScoreDf <- as.data.frame(colnames(mat))
	colnames(outlierScoreDf) <- "cellName"
	outlierScoreDf <- dplyr::mutate(outlierScoreDf, outlierScore=NA)
	
	for(i in 1:ncol(mat)){
		vec <- mat[,i]
		outlierScoreDf$outlierScore[i] <- length(vec[vec == 0])
	}
	
	outlierScoreDf$outlierScorePer <- outlierScoreDf$outlierScore / nrow(mat)
	return(outlierScoreDf)
}


.excludeOutliers <- function(dbscanMatrix, sceObject, threshold=0.3){
	# exclude outliers based on DBSCAN clustering
	# outliers are the cells which cannot be assigned
	# to any final cluster
	
	outlierInfo <- .mkOutlierScoreDf(dbscanMatrix)
	colData <- SummarizedExperiment::colData(sceObject)[SummarizedExperiment::colData(sceObject)$cellName
					%in% colnames(dbscanMatrix),]
	
	if(is.vector(colData)){
		colData <- S4Vectors::DataFrame(cellName=colData, row.names=colData)
	}
	
	# Case if coldata has the outliers scores already
	colData$outlierScorePer <- NULL
	colData$outlierScore <- NULL
	
	colData <- merge(colData, outlierInfo,
			by.x="cellName", by.y="cellName",
			all.x=TRUE, all.y=TRUE, sort=FALSE)
	rownames(colData) <- colData$cellName
	
	numberOfCellsBefore <- dim(colData)[1]
	sceObject <- sceObject[,colData$outlierScorePer < threshold]
	colData <- colData[colData$outlierScorePer < threshold,]
	numberOfCellsAfter <- dim(colData)[1]
	
	dbscanMatrix <- dbscanMatrix[,outlierInfo$outlierScorePer < threshold]
	SummarizedExperiment::colData(sceObject) <- colData
	
	message(paste(numberOfCellsBefore - numberOfCellsAfter,
					"outliers were excluded from the SingleCellExperiment object.\n"))
	
	return(list(sceObject, dbscanMatrix))
}



.runClustering <- function(tSNEResults, # for deleteOutliers = FALSE
		sceObject, dataDirectory,
		experimentName, epsilon=c(1.3, 1.4, 1.5),
		minPoints=c(3, 4), k=0, deepSplit=4,
		clusteringMethod = "ward.D2",
		cores=14,
		deleteOutliers = TRUE,
		PCs=c(4, 6, 8, 10, 20, 40, 50),
		perplexities=c(30,40), # for deleteOutliers = TRUE
		randomSeed = 42){
	# It combines all the clustering parts. Takes tSNE coordinates and gives
	# results of final clustering: sceObject and cell correlation matrix
	
	# run dbscan
	message("Running dbscan using ", cores, " cores.")
	dbscanResults <- runDBSCAN(tSNEResults, sceObject, dataDirectory,
			experimentName, epsilon=epsilon,
			minPoints=minPoints, cores=cores)
	if(deleteOutliers){
		# filter cluster outliers
		message("Excluding clustering outliers.")
		filteredResults <- .excludeOutliers(dbscanResults, sceObject)
		sceObjectFiltered <- filteredResults[[1]]
		#dbscanResultsFiltered <- filteredResults[[2]]
		message("Getting TSNE coordinates for the filtered sceObject.")
		tSNEResultsFiltered <- generateTSNECoordinates(sceObjectFiltered,
				dataDirectory,
				experimentName, PCs=PCs,
				perplexities=perplexities,
				randomSeed=randomSeed,
				cores = cores)
		
		dbscanResultsFiltered <- runDBSCAN(tSNEResultsFiltered, sceObjectFiltered,
				dataDirectory,
				experimentName, epsilon=epsilon,
				minPoints=minPoints, cores=cores)
	}else{
		sceObjectFiltered = sceObject
		dbscanResultsFiltered = dbscanResults
	}
	
	# assign cells to cluster
	clusteringResults <- clusterCellsInternal(dbscanResultsFiltered, sceObjectFiltered,
			clusterNumber=k, deepSplit=deepSplit,
			clusteringMethod=clusteringMethod,
			cores=cores)
	sceObjectFiltered <- clusteringResults[[1]]
	cellsSimilarityMatrix <- clusteringResults[[2]]
	
	return(list(sceObjectFiltered, cellsSimilarityMatrix))
}



### This function calculates how many time a cell were not assigned to
### any clusters by dbscan. ###
### it returns a data frame ###


.mkSimMat <- function(mat, cores=14){
	
	myCluster <- parallel::makeCluster(cores, # number of cores to use
			type = "PSOCK") # type of cluster
	doParallel::registerDoParallel(myCluster)
	
	simMats <- foreach::foreach(i=1:nrow(mat)) %dopar% {
		simMat <- matrix(0, ncol=ncol(mat), nrow=ncol(mat))
		colnames(simMat) <- colnames(mat)
		rownames(simMat) <- colnames(mat)
		
		vec <- unique(mat[i,])
		for(j in vec[vec!=0]){
			selCol <- colnames(mat)[mat[i,] == j]
			simMat[rownames(simMat) %in% selCol,
					colnames(simMat) %in% selCol] <-
					simMat[rownames(simMat) %in% selCol,
							colnames(simMat) %in% selCol] + 1
		}
		rm(cl, selCol, i, j)
		return(simMat)
	}
	parallel::stopCluster(myCluster)
	
	simMat <- matrix(0, ncol=ncol(mat), nrow=ncol(mat))
	colnames(simMat) <- colnames(mat)
	rownames(simMat) <- colnames(mat)
	
	for(i in 1:nrow(mat)){
		simMat <- simMat + simMats[[i]]
	}
	
	rm(simMats)
	simMat <- simMat / nrow(mat)
	stopifnot(isSymmetric(simMat))
	
	return(simMat)
}



clusterCellsInternal <- function(dbscanMatrix, sceObject, clusterNumber=0,
		deepSplit = 4, cores=1,
		clusteringMethod = "ward.D2") {
	# 
	
	cellsSimilarityMatrix <- .mkSimMat(dbscanMatrix, cores=cores)
	
	distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
	clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
	
	if(clusterNumber == 0){
		message(paste0("Assigning cells to clusters. DeepSplit = ", deepSplit))
		clusters <- unname(cutreeDynamic(clusteringTree,
						distM=as.matrix(distanceMatrix),
						verbose=0, deepSplit = deepSplit))
	} else {
		message(paste0("Assigning cells to ", clusterNumber, " clusters."))
		clusters <- cutree(clusteringTree, k=clusterNumber)
	}
	
	SummarizedExperiment::colData(sceObject)$clusters <- factor(clusters)
	
	return(list(sceObject, cellsSimilarityMatrix))
}




### This function returns a matrix with "protocells" representing clusters ###
### values show how much two "protocells" are similar ###
### 1 if clusters are very similar, 0 if very different ###

.mkSimMed <- function(simMat, clusters){
	
	clusMed <- matrix(ncol=length(unique(clusters)), nrow=nrow(simMat))
	clusterNames <- levels(clusters)
	
	for(i in 1:ncol(clusMed)){
		clusMed[,i] <- matrixStats::rowMedians(simMat[,clusters == clusterNames[i]])
	}
	
	clusMed <- t(clusMed)
	
	simMed <- matrix(ncol=length(unique(clusters)),
			nrow=length(unique(clusters)))
	
	for(i in 1:ncol(simMed)){
		simMed[,i] <- matrixStats::rowMedians(clusMed[,clusters == clusterNames[i]])
	}
	
	# colnames(simMed) = 1:length(unique(clusters))
	# rownames(simMed) = 1:length(unique(clusters))
	
	colnames(simMed) <- clusterNames
	rownames(simMed) <- clusterNames
	
	return(simMed)
}

calculateClustersSimilarity <- function(cellsSimilarityMatrix, sceObject,
		clusteringMethod = "ward.D2"){
	
	# Calculating cluster similarity for plotting picture
	# and ranking genes result is the square matrix with
	# dimension equal to number of cluster. Numbers in matrix
	# are similarity between cluster.
	
	clusters <- SummarizedExperiment::colData(sceObject)$clusters
	clustersNumber <- length(unique(clusters))
	clustersNames <- levels(clusters)
	
	# Calculating matrix
	clustersSimilarityMatrix <-
			.mkSimMed(simMat=as.matrix(cellsSimilarityMatrix), clusters=clusters)
	
	# Plotting matrix
	distanceMatrix <- as.dist(sqrt((1-clustersSimilarityMatrix)/2))
	clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
	
	clustersSimOrdered <- data.frame(clusterNames = clustersNames,
			clusterIndexes = 1:clustersNumber)
	rownames(clustersSimOrdered) <- clustersSimOrdered$clusterIndexes
	clustersSimOrdered <- clustersSimOrdered[clusteringTree$order,]
	
	return(list(clustersSimilarityMatrix, clustersSimOrdered$clusterNames))
	
}



addClusteringManually <- function(fileName, sceObject, dataDirectory,
		experimentName, columnName = "clusters"){
	
	tableData <- read.table(file.path(dataDirectory, "output_tables",
					paste0(experimentName,"_", fileName)), sep="\t")
	
	if(all(rownames(SummarizedExperiment::colData(sceObject)) %in% rownames(tableData))){
		if(ncol(tableData) == 1){
			SummarizedExperiment::colData(sceObject)$clusters <-
					factor(tableData[rownames(SummarizedExperiment::colData(sceObject)),])
		} else {
			SummarizedExperiment::colData(sceObject)$clusters <-
					factor(tableData[rownames(SummarizedExperiment::colData(sceObject)),][,columnName])
		}
		
		return(sceObject)
	} else {
		message("Rownames in colData are not equal to rownames in table.
						Returning SCE object with cells intersecting with clusters_table.")
		sceObject <- sceObject[ ,colnames(sceObject) %in%
						intersect(colnames(sceObject),
								rownames(tableData))]
		tableData$randomColumn <- NA
		tableData <- tableData[rownames(tableData) %in%
						intersect(colnames(sceObject),
								rownames(tableData)), ]
		SummarizedExperiment::colData(sceObject)$clusters <-
				factor(tableData[rownames(SummarizedExperiment::colData(sceObject)),][,columnName])
		return(sceObject)
	}
}
