#' Run CONCLUS in one click
#'
#' This function performs core CONCLUS workflow. It generates PCA and t-SNE coordinates, 
#' runs DBSCAN, calculates similarity matrices of cells and clusters, assigns cells to clusters,
#' searches for positive markers for each cluster. The function saves plots and tables into dataDirectory.
#'
#' @param sceObject a SingleCellExperiment object with your data.
#' @param dataDirectory CONCLUS will create this directory if it doesn't exist and store there all output files.
#' @param experimentName most of output file names of CONCLUS are hardcoded.
#' experimentName will stay at the beginning of each output file name to
#' distinguish different runs easily.
#' @param colorPalette a vector of colors for clusters.
#' @param statePalette a vector of colors for states.
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param epsilon a parameter of fpc::dbscan() function.
#' @param minPoints a parameter of fpc::dbscan() function.
#' @param k preferred number of clusters. Alternative to deepSplit. A parameter of cutree() function.
#' @param PCs a vector of first principal components.
#' For example, to take ranges 1:5 and 1:10 write c(5, 10).
#' @param perplexities a vector of perplexity for t-SNE.
#' @param randomSeed random seed for reproducibility.
#' @param deepSplit intuitive level of clustering depth. Options are 1, 2, 3, 4.
#' @param preClustered if TRUE, it will not change the column clusters after the run.
#' However, it will anyway run DBSCAN to calculate similarity matrices.
#' @param orderClusters can be either FALSE (default) of "name".
#' If "name", clusters in the similarity matrix of cells will be ordered by name.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param plotPDFcellSim if FALSE, the similarity matrix of cells will be saved in png format.
#' FALSE is recommended for count matrices with more than 2500 cells due to large pdf file size.
#' @param deleteOutliers whether cells which were often defined as outliers by dbscan must be deleted.
#' It will require recalculating of the similarity matrix of cells. Default is FALSE.
#' Usually those cells form a separate "outlier" cluster and can be easier distinguished and deleted later
#' if necessary.
#' @param tSNEalreadyGenerated if you already ran CONCLUS ones and have t-SNE coordinated saved
#' You can set TRUE to run the function faster since it will skip the generation of t-SNE coordinates and use the stored ones. 
#' Option TRUE requires t-SNE coordinates to be located in your 'dataDirectory/tsnes' directory.
#' @param tSNEresExp experimentName of t-SNE coordinates which you want to use.
#' This argument allows copying and pasting t-SNE coordinates between different CONCLUS runs without renaming the files.
#'
#' @keywords CONCLUS
#' @export
#' @return A SingleCellExperiment object.

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
