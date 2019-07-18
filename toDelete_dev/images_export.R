
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


png(paste0(outputDirectory, "fig1_KNN.png"), width = 700, height=700, res=100)
p <- conclus::testClustering(sceObject, outputDirectory, experimentName)
dev.off()

png(paste0(outputDirectory, "fig2_tSNEnoColor.png"), width = 700, height=700, res=100, type = "cairo-png")
p[[1]]
dev.off()

png(paste0(outputDirectory, "fig3_tSNEdbScan.png"), width = 700, height=700, res=100, type = "cairo-png")
p[[3]]
dev.off()


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


### Cell similarity heatmap


colorPalette="default"
statePalette="default"
plotPDFcellSim = TRUE
orderClusters = FALSE
clustersNumber <- length(unique(
				SummarizedExperiment::colData(sceObjectFiltered)$clusters))
colorPalette <- conclus::choosePalette(colorPalette, clustersNumber)

png(paste0(outputDirectory, "fig6_plotCellSimilarity.png"), width = 700, height=700, res = 100)
conclus::plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, outputDirectory,
		experimentName, colorPalette, 
		orderClusters = orderClusters, 
		statePalette = statePalette, 
		clusteringMethod = clusteringMethod,
		plotPDF = plotPDFcellSim,
		returnPlot = TRUE)
dev.off()

### Cluster similarity heatmap

png(paste0(outputDirectory, "fig7_plotClustersSimilarity.png"), width = 700, height=700)
conclus::plotClustersSimilarity(clustersSimilarityMatrix, 
		sceObjectFiltered,
		dataDirectory = outputDirectory, 
		experimentName = experimentName, 
		colorPalette = colorPalette,
		statePalette = statePalette,
		clusteringMethod = clusteringMethod,
		returnPlot = TRUE)
dev.off()

## Results export

conclus::exportMatrix(cellsSimilarityMatrix, outputDirectory, experimentName, 
		"cellsSimilarityMatrix")
conclus::exportMatrix(clustersSimilarityMatrix, outputDirectory, experimentName, 
		"clustersSimilarityMatrix")
conclus::exportData(sceObjectFiltered, outputDirectory, experimentName)


## Marker genes identification

conclus::rankGenes(sceObjectFiltered, clustersSimilarityMatrix, outputDirectory, 
		experimentName)
rankedGenesClus5 <- read.delim(file.path(outputDirectory, "marker_genes",
				"Bergiers_cluster_5_genes.tsv"),
		stringsAsFactors = FALSE)
head(rankedGenesClus5, n = 10)


# Plot a heatmap with positive marker genes

genesNumber <- 10
markersClusters <- conclus::getMarkerGenes(outputDirectory, sceObjectFiltered, 
		experimentName = experimentName,
		genesNumber = genesNumber)

orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- T  # F to show normalized counts

png(paste0(outputDirectory, "fig8_positivemarkers.png"), width = 700, height=700, res = 150)
conclus::plotCellHeatmap(markersClusters, sceObjectFiltered, outputDirectory, 
		experimentName, 
		paste0("clusters",
				length(levels(
								SummarizedExperiment::colData(sceObjectFiltered)$clusters)),
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
dev.off()



orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- F  # F to show normalized counts

png(paste0(outputDirectory, "fig9_positivemarkers-unnorm.png"), width = 700, height=700, res = 150)
conclus::plotCellHeatmap(markersClusters, sceObjectFiltered, outputDirectory, 
		experimentName, 
		paste0("clusters",
				length(levels(SummarizedExperiment::colData(sceObjectFiltered)$clusters)),
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
dev.off()



# Plot t-SNE colored by expression of a selected gene


ccl3 <- plotGeneExpression("Ccl3", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)

ggsave(filename = paste0(outputDirectory, "fig10_tSNESccl3.svg"), plot = ccl3, device = svg)


plotGeneExpression("Gp9", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)

plotGeneExpression("Fn1", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)

plotGeneExpression("Alox5ap", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)
		```





