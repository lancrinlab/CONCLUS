
runCONCLUS <- function(dataDirectory, experimentName, columnsMetaData,
		species,
		colorPalette="default",
		statePalette="default",
		clusteringMethod="ward.D2",
		epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4), k=0,
		PCs=c(4, 6, 8, 10, 20, 40, 50),
		perplexities=c(30,40),
		randomSeed = 42,
		deepSplit=4, preClustered = F,
		orderClusters = FALSE,
		cores=1,
		plotPDFcellSim = TRUE,
		deleteOutliers = TRUE,
		tSNEalreadyGenerated = FALSE,
		tSNEresExp = ""){
	
	initialisePath(dataDirectory)
	
	sceObject <- normaliseCountMatrix(countMatrix, species = species, 
			colData = columnsMetaData)
	
	# Generating 2D tSNE plots
	if(!tSNEalreadyGenerated){
		tSNEResults <- generateTSNECoordinates(sceObject, dataDirectory,
				experimentName, PCs=PCs,
				perplexities=perplexities,
				randomSeed = randomSeed,
				cores = cores)
	}else{
		tSNEResults <- readRDS(file.path(dataDirectory, "output_tables",
						paste0(tSNEresExp,"_tSNEResults.rds")))
	}
	
	if(preClustered){
		# Running dbscan
		message("Running dbscan using ", cores, " cores.")
		dbscanResults <- runDBSCAN(tSNEResults, sceObject, dataDirectory,
				experimentName, epsilon=epsilon,
				minPoints=minPoints,
				cores=cores)
		
		# assigning cells to clusters
		message("Calculating cells similarity matrix.")
		cellsSimilarityMatrix <- clusterCellsInternal(dbscanResults, sceObject, clusterNumber=k,
				deepSplit=deepSplit, cores=cores,
				clusteringMethod=clusteringMethod)[[2]]
		sceObjectFiltered <- sceObject
	} else {
		# Running clustering
		clusteringResults <- runClustering(tSNEResults, sceObject, dataDirectory,
				experimentName,
				epsilon=epsilon, minPoints=minPoints,
				k=k, deepSplit=deepSplit,
				cores=cores,
				clusteringMethod=clusteringMethod,
				deleteOutliers = deleteOutliers,
				PCs=PCs,
				perplexities=perplexities,
				randomSeed = randomSeed)
		sceObjectFiltered <- clusteringResults[[1]]
		cellsSimilarityMatrix <- clusteringResults[[2]]
	}
	
	print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters,
					dnn=list("Cells distribution by clusters")))
	
	clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
	colorPalette <- .choosePalette(colorPalette, clustersNumber)
	
	# Plotting cluster stablility and 2D visualisations
	plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, dataDirectory,
			experimentName, colorPalette,
			orderClusters = orderClusters,
			statePalette = statePalette,
			clusteringMethod = clusteringMethod,
			plotPDF = plotPDFcellSim)
	
	plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
			PCs=PCs, perplexities=perplexities, colorPalette,
			columnName = "clusters",
			tSNEresExp = tSNEresExp)
	plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
			PCs=PCs, perplexities=perplexities, colorPalette,
			columnName = "noColor",
			tSNEresExp = tSNEresExp)
	if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) %in% "state")){
		plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
				PCs=PCs, perplexities=perplexities, statePalette,
				columnName = "state",
				tSNEresExp = tSNEresExp)
	}
	
	# Calculating cluster similarity and marker genes
	clustersSimilarityMatrix <-
			calculateClustersSimilarity(cellsSimilarityMatrix,
					sceObject = sceObjectFiltered,
					clusteringMethod = clusteringMethod)[[1]]
	
	# Exporting key matrices
	exportMatrix(cellsSimilarityMatrix, dataDirectory, experimentName,
			"cellsSimilarityMatrix")
	exportMatrix(clustersSimilarityMatrix, dataDirectory, experimentName,
			"clustersSimilarityMatrix")
	
	plotClustersSimilarity(clustersSimilarityMatrix,
			sceObject = sceObjectFiltered,
			dataDirectory = dataDirectory,
			experimentName = experimentName,
			colorPalette = colorPalette,
			statePalette = statePalette,
			clusteringMethod = clusteringMethod)
	
	rankGenes(sceObjectFiltered, clustersSimilarityMatrix, dataDirectory,
			experimentName)
	
	return(sceObjectFiltered)
}
