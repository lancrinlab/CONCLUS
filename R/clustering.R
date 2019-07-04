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

#' To check one iteration of clustering before running full workflow CONCLUS.
#' 
#' This function generates a single clustering iteration of CONCLUS to check whether
#' chosen parameters for dbscan are suitable for your data.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param dbscanEpsilon a parameter of fpc::dbscan() function.
#' @param minPts a parameter of fpc::dbscan() function.
#' @param PCs a vector of PCs for plotting.
#' @param perplexities vector of perplexities (t-SNE parameter).
#' @param randomSeed random seed for reproducibility.
#' @param width plot width.
#' @param height plot height.
#' @param ... other pdf() arguments.
#'
#' @return t-SNE results, a distance graph plot, a t-SNE plot colored by test clustering solution.
#' @export

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



#' DBSCAN clustering on t-SNE results.
#' 
#' This function provides consensus DBSCAN clustering based on the results of t-SNE. 
#' You can tune algorithm parameters in options to get the number of clusters you want.
#'
#' @param tSNEResults the result of conclus::generateTSNECoordinates() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory of a given CONCLUS run (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param epsilon a parameter of fpc::dbscan() function.
#' @param minPoints a parameter of fpc::dbscan() function.
#' @param k preferred number of clusters. Alternative to deepSplit.
#' @param PCs a vector of first principal components.
#' For example, to take ranges 1:5 and 1:10 write c(5, 10).
#' @param perplexities a vector of perplexity for t-SNE.
#' @param randomSeed random seed for reproducibility.
#' @param deepSplit intuitive level of clustering depth. Options are 1, 2, 3, 4.
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param deleteOutliers Whether cells which were often defined as outliers by dbscan must be deleted.
#' It will require recalculating of the similarity matrix of cells. Default is FALSE.
#' Usually those cells appear in an "outlier" cluster and can be easier distinguished and deleted later
#' if necessary.
#' 
#'
#' @return A list containing filtered from outliers SingleCellExperiment object and cells similarity matrix.
#' @export
runClustering <- function(tSNEResults, # for deleteOutliers = FALSE
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



#' Cluster cells and get similarity matrix of cells.
#' 
#' The function returns consensus clusters by using hierarchical clustering on the similarity matrix of cells.
#' It provides two options: to specify an exact number of clusters (with clusterNumber parameter)
#' or to select the depth of splitting (deepSplit parameter).
#' 
#' @param dbscanMatrix an output matrix of conclus::runDBSCAN() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param clusterNumber a parameter, specifying the exact number of cluster.
#' @param deepSplit a parameter, specifying how deep we will split the clustering tree. It takes integers from 1 to 4.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param clusteringMethod a clustering methods passed to hclust() function.
#'
#' @return A SingleCellExperiment object with modified/created "clusters" column in the colData, and cells similarity matrix.
#' @export
clusterCellsInternal <- function(dbscanMatrix, sceObject, clusterNumber=0,
		deepSplit, cores=14,
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
