
# Quick start

library(conclus)

outputDirectory <- "/home/descostes/Documents/analysis/lancrin/conclus/vignette_images/"
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
		k = 10, cores = 1,
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"),
		deleteOutliers = FALSE)

exportClusteringResults(sceObjectCONCLUS, outputDirectory, experimentName, 
		"clusters_table.tsv")

## Test clustering

sceObject <- conclus::normaliseCountMatrix(countMatrix, species = "mmu", 
		colData = columnsMetaData)
p <- conclus::testClustering(sceObject, outputDirectory, experimentName)

ggsave(filename = paste0(outputDirectory, "fig1_TestClustering_tSNEnoColor.svg"), plot = p[[1]], device = svg)
ggsave(filename = paste0(outputDirectory, "fig1_TestClustering_tSNEdbScan.svg"), plot = p[[3]], device = svg)


# CONCLUS step by step

## Normalization of the counts matrix

sceObject <- conclus::normaliseCountMatrix(countMatrix, species = "mmu", 
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
# default parameters, can be selected by a user
PCs=c(4, 6, 8, 10, 20, 40, 50)
perplexities=c(30, 40)
randomSeed = 42
tSNEResults <- generateTSNECoordinates(sceObject, outputDirectory, 
		experimentName, PCs=PCs, 
		perplexities=perplexities,
		randomSeed = randomSeed)

ncol(tSNEResults)
# the third matrix of t-SNE coordinates with PC = 8 and perplixities = 30
# it is saved as "tsnes/Bergiers_tsne_coordinates_3_8PCs_30perp.tsv"
head(tSNEResults[1,3][[1]])


## Clustering with DB-SCAN

epsilon=c(1.3, 1.4, 1.5)
minPoints=c(3, 4)
cores=1
dbscanResults <- conclus::runDBSCAN(tSNEResults, sceObject, outputDirectory, 
		experimentName, epsilon=epsilon, 
		minPoints=minPoints,
		cores=cores)

## Cell and cluster similarity matrix calculation

clusteringMethod="ward.D2"
k=10 # parameter for cutree
deepSplit = 0 # 0 to avoid cutreeDynamic, 1 to 4 to use it
message("Calculating cells similarity matrix.")
clusteringResults <- conclus::clusterCellsInternal(dbscanResults, sceObject, clusterNumber=k, 
		deepSplit=deepSplit, cores=cores,
		clusteringMethod=clusteringMethod)
sceObjectFiltered <- clusteringResults[[1]]
cellsSimilarityMatrix <- clusteringResults[[2]]

print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters, 
				dnn=list("Cells distribuion by clusters")))

clustersSimilarityMatrix <- 
		conclus::calculateClustersSimilarity(cellsSimilarityMatrix, 
				sceObject = sceObjectFiltered,
				clusteringMethod = "ward.D2")[[1]]


## Plotting


### t-SNE colored by clusters or conditions

tSNEclusters <- conclus::plotClusteredTSNE(sceObjectFiltered, outputDirectory, 
		experimentName, PCs=PCs, perplexities=perplexities, 
		colorPalette = "default", columnName = "clusters", 
		returnPlot = TRUE)

tSNEnoColor <- conclus::plotClusteredTSNE(sceObjectFiltered, outputDirectory, 
		experimentName, PCs=PCs, perplexities=perplexities, 
		colorPalette = "default", columnName = "noColor", 
		returnPlot = TRUE)

if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) 
				%in% "state")){
	tSNEstate <- conclus::plotClusteredTSNE(sceObjectFiltered, outputDirectory, 
			experimentName, PCs=PCs, perplexities=perplexities, 
			colorPalette = "default", columnName = "state", 
			returnPlot = TRUE)
}

ggsave(filename = paste0(outputDirectory, "fig3_tSNEColor_clusters.svg"), plot = tSNEclusters[[5]], device = svg)
ggsave(filename = paste0(outputDirectory, "fig4_tSNEnoColor.svg"), plot = tSNEnoColor[[5]], device = svg)
ggsave(filename = paste0(outputDirectory, "fig5_tSNEState.svg"), plot = tSNEstate[[5]], device = svg)





