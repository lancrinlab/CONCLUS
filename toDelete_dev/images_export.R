
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

png(paste0(outputDirectory, "fig4_tSNEColor_clusters.png"), width = 700, height=700, res=100, type = "cairo-png")
tSNEclusters[[5]]
dev.off()

png(paste0(outputDirectory, "fig5_tSNEnoColor.png"), width = 700, height=700, res=100, type = "cairo-png")
tSNEnoColor[[5]]
dev.off()

png(paste0(outputDirectory, "fig6_tSNEState.png"), width = 700, height=700, res=100, type = "cairo-png")
tSNEstate[[5]]
dev.off()



### Cell similarity heatmap


colorPalette="default"
statePalette="default"
plotPDFcellSim = TRUE
orderClusters = FALSE
clustersNumber <- length(unique(
				SummarizedExperiment::colData(sceObjectFiltered)$clusters))
colorPalette <- conclus::choosePalette(colorPalette, clustersNumber)

png(paste0(outputDirectory, "fig7_plotCellSimilarity.png"), width = 700, height=700, res = 100)
conclus::plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, outputDirectory,
		experimentName, colorPalette, 
		orderClusters = orderClusters, 
		statePalette = statePalette, 
		clusteringMethod = clusteringMethod,
		plotPDF = plotPDFcellSim,
		returnPlot = TRUE)
dev.off()

### Cluster similarity heatmap

png(paste0(outputDirectory, "fig8_plotClustersSimilarity.png"), width = 700, height=700, res=100)
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

png(paste0(outputDirectory, "fig9_positivemarkers.png"), width = 700, height=700, res = 160)
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

png(paste0(outputDirectory, "fig10_positivemarkers-unnorm.png"), width = 700, height=700, res = 160)
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

png(paste0(outputDirectory, "fig11_positivemarkers-ccl3.png"), width = 700, height=700, res = 100)
plotGeneExpression("Ccl3", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)
dev.off()


png(paste0(outputDirectory, "fig12_positivemarkers-Gp9.png"), width = 700, height=700, res = 100)
plotGeneExpression("Gp9", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)
dev.off()

png(paste0(outputDirectory, "fig13_positivemarkers-Fn1.png"), width = 700, height=700, res = 100)
plotGeneExpression("Fn1", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)
dev.off()

png(paste0(outputDirectory, "fig14_positivemarkers-Alox5ap.png"), width = 700, height=700, res = 100)
plotGeneExpression("Alox5ap", experimentName, outputDirectory, sceObjectFiltered, 
		tSNEpicture = 10, returnPlot=T)
dev.off()


# Collect publicly available info about marker genes

## Collect information for the top 10 markers for each cluster

result <- getGenesInfo(markersClusters, groupBy = "clusters",
		getUniprot = TRUE)

outputDir <- file.path(outputDirectory, "/marker_genes/getGenesInfo")
dir.create(outputDir, showWarnings=F)
write.table(result, file = file.path(outputDir, 
				"Bergiers_markersClusters_top10_clusters9_genesInfo.csv"),
		quote = FALSE, sep = ";", row.names = FALSE)


## Collect information for the top X markers for each cluster

saveMarkersLists(experimentName, outputDirectory)
saveGenesInfo(outputDirectory, sep = ";", header = TRUE, 
		startFromFile = 10, getUniprot = TRUE)


# Supervised clustering 

exportClusteringResults(sceObjectFiltered, outputDirectory, experimentName, 
		"clusters_table.tsv")
clustersTable <- read.delim(file.path(outputDirectory, "output_tables", 
				paste0(experimentName, "_clusters_table.tsv")), stringsAsFactors = FALSE) 

clustersTable$clusters[clustersTable$clusters == "3"] = "2"
clustersTable$clusters[clustersTable$clusters == "9"] = "8"
clustersTable$clusters[clustersTable$clusters == "10"] = "8"

write.table(clustersTable, file.path(outputDirectory, "output_tables", 
				paste0(experimentName, "_clusters_table_manual.tsv")), 
		quote = FALSE, sep = "\t")



sceObjectFiltered <- addClusteringManually(fileName = "clusters_table_manual.tsv", 
		dataDirectory = outputDirectory, 
		experimentName = experimentName,
		sceObject = sceObjectFiltered, 
		columnName = "clusters")

# Redo the analysis with manual clustering
sceObjectFiltered2 <- runCONCLUS(outputDirectory, experimentName, 
		statePalette= c("bisque", "cadetblue2", "coral1", "cornflowerblue"),
		preClustered = TRUE, manualClusteringObject = sceObjectFiltered, k=7, cores=10)


dataDirectory = outputDirectory 
statePalette= c("bisque", "cadetblue2", "coral1", "cornflowerblue")
preClustered = TRUE
manualClusteringObject = sceObjectFiltered
cores=10
columnsMetaData = NA
species = NA
colorPalette="default"
clusteringMethod="ward.D2"
epsilon=c(1.3, 1.4, 1.5)
minPoints=c(3, 4)
k=7
PCs=c(4, 6, 8, 10, 20, 40, 50)
perplexities=c(30,40)
randomSeed = 42
deepSplit=4
orderClusters = FALSE
plotPDFcellSim = TRUE
deleteOutliers = TRUE
tSNEalreadyGenerated = FALSE
tSNEresExp = ""




meanCentered <- T  # F to show normalized counts, 
plotCellHeatmap(markersClusters, sceObjectFiltered, outputDirectory, 
		experimentName, 
		paste0("clusters",
				length(levels(
								SummarizedExperiment::colData(sceObjectFiltered)$clusters)),
				"_meanCentered",meanCentered,"_orderClusters",orderClusters,
				"_orderGenes",orderGenes,"_top",
				genesNumber, "markersPerCluster"), 
		meanCentered = meanCentered, 
		colorPalette = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
				"#FB9A99", "#E31A1C", "#FDBF6F"),
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


