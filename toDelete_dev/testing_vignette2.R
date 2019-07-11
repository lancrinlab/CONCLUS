library(foreach)
library(ggplot2)
library(pheatmap)
library(zoo)
library(dynamicTreeCut)
library(factoextra) 
library(digest)
library(RColorBrewer)
library(doParallel)
library(BiocParallel)
library(scran) 
library(monocle)
library(KEGGREST)
library(AnnotationDbi)
library(SingleCellExperiment)
library(Cairo)
library(curl)

source("/home/descostes/git/CONCLUS/R/clustering.R")
source("/home/descostes/git/CONCLUS/R/dbScan.R")
source("/home/descostes/git/CONCLUS/R/export.R")
source("/home/descostes/git/CONCLUS/R/markers.R")
source("/home/descostes/git/CONCLUS/R/normalization.R")
source("/home/descostes/git/CONCLUS/R/plot.R")
source("/home/descostes/git/CONCLUS/R/runCONCLUS.R")
source("/home/descostes/git/CONCLUS/R/sharedInternals.R")
source("/home/descostes/git/CONCLUS/R/tSNE.R")


dataDirectory <- "/home/descostes/Documents/analysis/lancrin/conclus/testing/"
experimentName <- "Bergiers"
load("/home/descostes/Documents/analysis/lancrin/conclus/colData.Rdat")
load("/home/descostes/Documents/analysis/lancrin/conclus/countMatrix.Rdat")


## Quick start

		
outputDirectory <- dataDirectory
experimentName <- "Bergiers"

countMatrix <- read.delim(file.path(system.file("extdata", package = "conclus"),
				"Bergiers_counts_matrix_filtered.tsv"), 
		stringsAsFactors = FALSE)
columnsMetaData <- read.delim(file.path(system.file("extdata", package = "conclus"), 
				"Bergiers_colData_filtered.tsv"))

sceObjectCONCLUS <- runCONCLUS(outputDirectory, 
		experimentName,
		columnsMetaData,
		species = "mmu", 
		plotPDFcellSim = TRUE, # FALSE for > 2500 cells
		k = 10, cores = 8,
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"),
		deleteOutliers = FALSE)

exportClusteringResults(sceObjectCONCLUS, outputDirectory, experimentName, 
		"clusters_table.tsv")

## Test clustering

sceObject <- normaliseCountMatrix(countMatrix, species = "mmu", 
		colData = columnsMetaData)
p <- testClustering(sceObject, outputDirectory, experimentName)
p[[1]]
p[[3]]


# CONCLUS step by step


## Normalization of the counts matrix

sceObject <- normaliseCountMatrix(countMatrix, species = "mmu", 
		colData = columnsMetaData)

dim(sceObject)
SingleCellExperiment::counts(sceObject)[1:5,1:5]
Biobase::exprs(sceObject)[1:5,1:5]
coldataSCE <- as.data.frame(SummarizedExperiment::colData(sceObject))
head(coldataSCE)
rowdataSCE <- as.data.frame(SummarizedExperiment:::rowData(sceObject))
head(rowdataSCE)


## Generation of t-SNE coordinates

initialisePath(outputDirectory)
PCs=c(4, 6, 8, 10, 20, 40, 50)
perplexities=c(30, 40)
randomSeed = 42
tSNEResults <- generateTSNECoordinates(sceObject, outputDirectory, 
		experimentName, PCs=PCs, 
		perplexities=perplexities,
		randomSeed = randomSeed, cores=8)
ncol(tSNEResults)
head(tSNEResults[1,3][[1]])

## Clustering with DB-SCAN

epsilon=c(1.3, 1.4, 1.5)
minPoints=c(3, 4)
cores=14
message("Running dbscan using ", cores, " cores.")
dbscanResults <- runDBSCAN(tSNEResults, sceObject, outputDirectory, 
		experimentName, epsilon=epsilon, 
		minPoints=minPoints,
		cores=cores)

## Cell and cluster similarity matrix calculation

clusteringMethod="ward.D2"
k=10 # parameter for cutree
deepSplit = 0 # 0 to avoid cutreeDynamic, 1 to 4 to use it
message("Calculating cells similarity matrix.")
clusteringResults <- clusterCellsInternal(dbscanResults, sceObject, clusterNumber=k, 
		deepSplit=deepSplit, cores=cores,
		clusteringMethod=clusteringMethod)
sceObjectFiltered <- clusteringResults[[1]]
cellsSimilarityMatrix <- clusteringResults[[2]]

print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters, 
				dnn=list("Cells distribuion by clusters")))

clustersSimilarityMatrix <- 
		calculateClustersSimilarity(cellsSimilarityMatrix, 
				sceObject = sceObjectFiltered,
				clusteringMethod = "ward.D2")[[1]]


## Plotting


### t-SNE colored by clusters or conditions

tSNEclusters <- plotClusteredTSNE(sceObjectFiltered, outputDirectory, experimentName,
		PCs=PCs, perplexities=perplexities, colorPalette = "default",
		columnName = "clusters", returnPlot = TRUE)
tSNEnoColor <- plotClusteredTSNE(sceObjectFiltered, outputDirectory, experimentName,
		PCs=PCs, perplexities=perplexities, colorPalette = "default",
		columnName = "noColor", returnPlot = TRUE)

if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) %in% "state")){
	tSNEstate <- plotClusteredTSNE(sceObjectFiltered, outputDirectory, experimentName,
			PCs=PCs, perplexities=perplexities, colorPalette = "default",
			columnName = "state", returnPlot = TRUE)
}

tSNEclusters[[5]]
tSNEnoColor[[5]]
tSNEstate[[5]]


### Cell similarity heatmap

colorPalette="default"
statePalette="default"
plotPDFcellSim = TRUE
orderClusters = FALSE
clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
colorPalette <- choosePalette(colorPalette, clustersNumber)

plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, outputDirectory,
		experimentName, colorPalette, 
		orderClusters = orderClusters, 
		statePalette = statePalette, 
		clusteringMethod = clusteringMethod,
		plotPDF = plotPDFcellSim,
		returnPlot = TRUE)

### Cluster similarity heatmap
  
plotClustersSimilarity(clustersSimilarityMatrix, 
		sceObjectFiltered,
		dataDirectory = outputDirectory, 
		experimentName = experimentName, 
		colorPalette = colorPalette,
		statePalette = statePalette,
		clusteringMethod = clusteringMethod,
		returnPlot = TRUE)

## Results export

exportMatrix(cellsSimilarityMatrix, outputDirectory, experimentName, 
		"cellsSimilarityMatrix")
exportMatrix(clustersSimilarityMatrix, outputDirectory, experimentName, 
		"clustersSimilarityMatrix")
exportData(sceObjectFiltered, outputDirectory, experimentName)

## Marker genes identification

conclus::rankGenes(sceObjectFiltered, clustersSimilarityMatrix, outputDirectory, 
		experimentName)
rankedGenesClus5 <- read.delim(file.path(outputDirectory, "marker_genes",
				"Bergiers_cluster_5_genes.tsv"),
		stringsAsFactors = FALSE)
head(rankedGenesClus5, n = 10)


# Plot a heatmap with positive marker genes

genesNumber <- 10
markersClusters <- conclus::getMarkerGenes(outputDirectory, sceObjectCONCLUS, 
		experimentName = experimentName,
		genesNumber = genesNumber)

orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- T  # F to show normalized counts

conclus::plotCellHeatmap(markersClusters, sceObjectCONCLUS, outputDirectory, 
		experimentName, 
		paste0("clusters",
				length(levels(SummarizedExperiment::colData(sceObjectCONCLUS)$clusters)),
				"_meanCentered",meanCentered,
				"_orderClusters",orderClusters,
				"_orderGenes",orderGenes,"_top",
				genesNumber, "markersPerCluster"), 
		meanCentered = meanCentered, 
		colorPalette = RColorBrewer::brewer.pal(10, "Paired"),
		orderClusters = orderClusters,
		orderGenes = orderGenes,
		fontsize_row = 4,
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"),
		color = colorRampPalette(c("#023b84","#4b97fc", 
						"#FEE395", 
						"#F4794E", "#D73027",
						"#a31008","#7a0f09"))(100),
		returnPlot = TRUE,
		width = 7.5, height = 6.5)

orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- F  # F to show normalized counts
conclus::plotCellHeatmap(markersClusters, sceObjectCONCLUS, outputDirectory, 
		experimentName, 
		paste0("clusters",
				length(levels(SummarizedExperiment::colData(sceObjectCONCLUS)$clusters)),
				"_meanCentered",meanCentered,
				"_orderClusters",orderClusters,
				"_orderGenes",orderGenes,"_top",
				genesNumber, "markersPerCluster"), 
		meanCentered = meanCentered, 
		colorPalette = RColorBrewer::brewer.pal(10, "Paired"),
		orderClusters = orderClusters,
		orderGenes = orderGenes,
		fontsize_row = 4,
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"),
		color = colorRampPalette(c("#023b84","#4b97fc", 
						"#FEE395", 
						"#F4794E", "#D73027",
						"#a31008","#7a0f09"))(100),
		returnPlot = TRUE)


# Plot t-SNE colored by expression of a selected gene

plotGeneExpression("Cd34", experimentName, dataDirectory, sceObjectCONCLUS, tSNEpicture = 10, returnPlot=T)
plotGeneExpression("Ctsc", experimentName, dataDirectory, sceObjectCONCLUS, tSNEpicture = 10, returnPlot=T)
plotGeneExpression("H19", experimentName, dataDirectory, sceObjectCONCLUS, tSNEpicture = 10, returnPlot=T)
plotGeneExpression("Postn", experimentName, dataDirectory, sceObjectCONCLUS, tSNEpicture = 10, returnPlot=T)

# Collect publicly available info about marker genes

result <- getGenesInfo(markersClusters, groupBy = "clusters",
		getUniprot = FALSE) # please change to getUniprot = TRUE

outputDir <- file.path(outputDirectory, "/marker_genes/getGenesInfo")
dir.create(outputDir, showWarnings=F)
write.table(result, file = file.path(outputDir, 
				"Bergiers_markersClusters_top10_clusters9_genesInfo.csv"),
		quote = FALSE, sep = ";", row.names = FALSE)

saveMarkersLists(experimentName, outputDirectory)

saveGenesInfo(outputDirectory, sep = ";", header = TRUE, 
		startFromFile = 1, getUniprot = FALSE) # please change to getUniprot = TRUE

# Supervised clustering 

exportClusteringResults(sceObjectCONCLUS, dataDirectory, experimentName, "clusters_table.tsv")
clustersTable <- read.delim(file.path(dataDirectory, "output_tables", paste0(experimentName, "_clusters_table.tsv")), stringsAsFactors = FALSE) 
clustersTable$clusters[clustersTable$clusters == "3"] = "2"
clustersTable$clusters[clustersTable$clusters == "4"] = "2"
clustersTable$clusters[clustersTable$clusters == "9"] = "5"
write.table(clustersTable, file.path(dataDirectory, "output_tables", paste0(experimentName, "_clusters_table_manual.tsv")), quote = FALSE, sep = "\t")


sceObjectCONCLUS <- addClusteringManually(fileName = "clusters_table_manual.tsv", 
		dataDirectory = dataDirectory, 
		experimentName = experimentName,
		sceObject = sceObjectCONCLUS, 
		columnName = "clusters")
		

sceObjectCONCLUS <- runCONCLUS(dataDirectory, experimentName, 
		statePalette= c("bisque", "cadetblue2", "coral1", "cornflowerblue"),
		preClustered = TRUE, manualClusteringObject = sceObjectCONCLUS, cores = 8) 




genesNumber <- 10
markersClusters <- getMarkerGenes(dataDirectory, sceObjectCONCLUS, 
		experimentName = experimentName,
		genesNumber = genesNumber)

orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- F  
plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
		experimentName, 
		paste0("clusters",
				length(levels(colData(sceObjectCONCLUS)$clusters)),
				"_meanCentered",meanCentered,
				"_orderClusters",orderClusters,
				"_orderGenes",orderGenes,"_top",
				genesNumber, "markersPerCluster"), 
		meanCentered = meanCentered, 
		colorPalette = brewer.pal(10, "Paired"),
		orderClusters = orderClusters,
		orderGenes = orderGenes,
		fontsize_row = 5,
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"),
		color = colorRampPalette(c("#023b84","#4b97fc", 
						"#FEE395", 
						"#F4794E", "#D73027",
						"#a31008","#7a0f09"))(100),
		returnPlot = TRUE)

meanCentered <- T  # F to show normalized counts,
plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
		experimentName, 
		paste0("clusters",
				length(levels(colData(sceObjectCONCLUS)$clusters)),
				"_meanCentered",meanCentered,"_orderClusters",orderClusters,
				"_orderGenes",orderGenes,"_top",
				genesNumber, "markersPerCluster"), 
		meanCentered = meanCentered, 
		colorPalette = brewer.pal(10, "Paired"),
		orderClusters = orderClusters,
		orderGenes = orderGenes,
		fontsize_row = 5,
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"),
		color = colorRampPalette(c("#023b84","#4b97fc", 
						"#FEE395", 
						"#F4794E", "#D73027",
						"#a31008","#7a0f09"))(100),
		returnPlot = TRUE)


plotGeneExpression("Cd34", experimentName, dataDirectory, sceObjectCONCLUS,
		tSNEpicture = 10, returnPlot = TRUE)
plotGeneExpression("Apoe", experimentName, dataDirectory, sceObjectCONCLUS,
		tSNEpicture = 10, returnPlot = TRUE)
plotGeneExpression("Ctsc", experimentName, dataDirectory, sceObjectCONCLUS,
		tSNEpicture = 10, returnPlot = TRUE)
plotGeneExpression("Clec4n", experimentName, dataDirectory, sceObjectCONCLUS,
		tSNEpicture = 10, returnPlot = TRUE)




	