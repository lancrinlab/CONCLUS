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


## 3 Genes and cells filtering, normalization


sceObject <- normaliseCountMatrix(countMatrix, species = "mmu", 
		colData = colData)
dim(sceObject)
SingleCellExperiment::counts(sceObject)[1:5,1:5]
Biobase::exprs(sceObject)[1:5,1:5]
coldataSCE <- as.data.frame(SummarizedExperiment::colData(sceObject))
head(coldataSCE)
rowdataSCE <- as.data.frame(SummarizedExperiment:::rowData(sceObject))
head(rowdataSCE)


## 4 Test clustering

p <- testClustering(sceObject, dataDirectory, experimentName)
p[[1]]
p[[3]]


## 5.1 Generate t-SNE coordinates

# !! make this initialise path optional
initialisePath(dataDirectory)
# default parameters, can be selected by a user
PCs=c(4, 6, 8, 10, 20, 40, 50)
perplexities=c(30, 40)
randomSeed = 42

tSNEResults <- generateTSNECoordinates(sceObject, dataDirectory, 
		experimentName, PCs=PCs, 
		perplexities=perplexities,
		randomSeed = randomSeed)


## 5.2 Option 1: upload your clustering solution to CONCLUS

SummarizedExperiment::colData(sceObject)$clusters = factor(c(rep(1, 100), rep(2, 200), rep(3, (ncol(sceObject)-300) ) ))
table(SummarizedExperiment::colData(sceObject)$clusters)
epsilon=c(1.3, 1.4, 1.5)
minPoints=c(3, 4)
cores=14
message("Running dbscan using ", cores, " cores.")
dbscanResults <- runDBSCAN(tSNEResults, sceObject, dataDirectory, 
		experimentName, epsilon=epsilon, 
		minPoints=minPoints,
		cores=cores)
dim(dbscanResults)
dbscanResults[1:7, 1:10]

## 5.3 Get cells similarity matrix

clusteringMethod="ward.D2"
k=10 # parameter for cutree
message("Calculating cells similarity matrix.")
cellsSimilarityMatrix <- clusterCellsInternal(dbscanResults, sceObject, clusterNumber=k, 
		deepSplit=deepSplit, cores=cores,
		clusteringMethod=clusteringMethod)[[2]]
sceObjectFiltered <- sceObject

print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters, 
				dnn=list("Cells distribuion by clusters")))

colorPalette="default"
statePalette="default"
plotPDFcellSim = TRUE
orderClusters = FALSE
clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
colorPalette <- choosePalette(colorPalette, clustersNumber)

plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, dataDirectory,
		experimentName, colorPalette, 
		orderClusters = orderClusters, 
		statePalette = statePalette, 
		clusteringMethod = clusteringMethod,
		plotPDF = plotPDFcellSim,
		returnPlot = TRUE)


## 5.4 Option 2: receive CONCLUS clustering based on dbscan

deepSplit = 0 # 0 to avoid cutreeDynamic, 1 to 4 to use it
deleteOutliers = FALSE
epsilon=c(1.3, 1.4, 1.5)
minPoints=c(3, 4)
cores=14
clusteringMethod="ward.D2"
k=10 # split the dendrogram with cutree function into 10 groups
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

print(table(SummarizedExperiment::colData(sceObjectFiltered)$clusters, 
				dnn=list("Cells distribuion by clusters")))

colorPalette="default"
statePalette="default"
plotPDFcellSim = TRUE
orderClusters = FALSE
clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
colorPalette <- choosePalette(colorPalette, clustersNumber)


plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, dataDirectory,
		experimentName, colorPalette, 
		orderClusters = orderClusters, 
		statePalette = statePalette, 
		clusteringMethod = clusteringMethod,
		plotPDF = plotPDFcellSim,
		returnPlot = TRUE)

## 5.5 Plot t-SNE colored by clusters or conditions

tSNEclusters <- plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
		PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
		columnName = "clusters", returnPlot = TRUE)
tSNEnoColor <- plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
		PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
		columnName = "noColor", returnPlot = TRUE)
if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) %in% "state")){
	tSNEstate <- plotClusteredTSNE(sceObjectFiltered, dataDirectory, experimentName,
			PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
			columnName = "state", returnPlot = TRUE)
}

tSNEclusters[[5]]
tSNEnoColor[[5]]
tSNEstate[[5]]

## 5.6 Calculate cluster similarity matrix

clustersSimilarityMatrix <- 
		calculateClustersSimilarity(cellsSimilarityMatrix, 
				sceObject = sceObjectFiltered,
				clusteringMethod = "ward.D2")[[1]]

plotClustersSimilarity(clustersSimilarityMatrix, 
		sceObjectFiltered,
		dataDirectory = dataDirectory, 
		experimentName = experimentName, 
		colorPalette = colorPalette,
		statePalette = statePalette,
		clusteringMethod = clusteringMethod,
		returnPlot = TRUE)

## 5.7 Identify marker genes for each cell cluster

rankGenes(sceObjectFiltered, clustersSimilarityMatrix, dataDirectory, 
		experimentName)

rankedGenesClus5 <- read.delim(file.path(dataDirectory, "marker_genes",
				"Bergiers_cluster_5_genes.tsv"),
		stringsAsFactors = FALSE)
head(rankedGenesClus5, n = 10)


## 5.8 Export key matrices

exportMatrix(cellsSimilarityMatrix, dataDirectory, experimentName, 
		"cellsSimilarityMatrix")
exportMatrix(clustersSimilarityMatrix, dataDirectory, experimentName, 
		"clustersSimilarityMatrix")


## 6 runCONCLUS in one function

sceObjectCONCLUS <- runCONCLUS(sceObject, dataDirectory, experimentName, 
		plotPDFcellSim = TRUE, # FALSE for > 2500 cells
		k = 10,
		cores = 14, # 14 for servers, 1 for PC
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"),
		deleteOutliers = FALSE  # TRUE takes more time
)

conclus::exportClusteringResults(sceObjectCONCLUS, dataDirectory, experimentName, 
		"clusters_table.tsv")

## 7 Plot a heatmap with positive marker genes

genesNumber <- 10
markersClusters <- getMarkerGenes(dataDirectory, sceObjectCONCLUS, 
		experimentName = experimentName,
		genesNumber = genesNumber)
orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- T  # F to show normalized counts
plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
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
plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
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

## 8 Give names to the clusters

clustersTable <- read.delim(file.path(dataDirectory, "output_tables",                                                    
				paste0(experimentName, "_clusters_table.tsv")), 
		stringsAsFactors = FALSE)
clustersTable$clusters[clustersTable$clusters == "9"] <- "newCluster"
clustersTable$clusters[clustersTable$clusters == "10"] <- "newCluster"
write.table(clustersTable, file.path(dataDirectory, "output_tables",                                                    
				paste0(experimentName, "_clusters_table.tsv")), 
		quote = FALSE, sep = "\t")


sceObjectCONCLUS <- addClusteringManually(fileName = "clusters_table.tsv", 
		dataDirectory = dataDirectory, 
		experimentName = experimentName,
		sceObject = sceObjectCONCLUS, 
		columnName = "clusters")


sceObjectCONCLUS <- runCONCLUS(sceObjectCONCLUS, dataDirectory, experimentName, 
		preClustered = TRUE,
		tSNEalreadyGenerated = TRUE, # to use t-SNE coords from dataDirectory/tsnes
		tSNEresExp = experimentName,
		cores = 14, # 14 for servers, 1 for PC
		statePalette = c("bisque", "cadetblue2", 
				"coral1", "cornflowerblue"))

genesNumber <- 10
markersClusters <- getMarkerGenes(dataDirectory, sceObjectCONCLUS, 
		experimentName = experimentName,
		genesNumber = genesNumber)
orderClusters <- T # F will apply hierarchical clustering to all cells
orderGenes <- T    # F will apply hierarchical clustering to all genes
meanCentered <- T  # F to show normalized counts
plotCellHeatmap(markersClusters, sceObjectCONCLUS, dataDirectory, 
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
						"#a31008", "#921912"))(100),
		returnPlot = TRUE,
		width = 7, height = 5.5)


plotGeneExpression("Ccl3", experimentName, dataDirectory, 
		sceObject = sceObjectCONCLUS,
		tSNEpicture = 10, returnPlot = TRUE)


conclus::plotGeneExpression("Lama1", experimentName, dataDirectory, 
		sceObject = sceObjectCONCLUS,
		tSNEpicture = 10, returnPlot = TRUE)

tSNEstate[[10]]


## 9 Collect publicly available info about marker genes

result <- getGenesInfo(markersClusters, groupBy = "clusters",
		getUniprot = TRUE) # please change to getUniprot = TRUE

outputDir <- file.path(dataDirectory, "/marker_genes/getGenesInfo")
dir.create(outputDir, showWarnings=F)
write.table(result, file = file.path(outputDir, 
				"Bergiers_markersClusters_top10_clusters9_genesInfo.csv"),
		quote = FALSE, sep = ";", row.names = FALSE)








































		






#------------ testing ncbi retrieval

genes = markersClusters 
databaseDir = system.file("extdata", package = "conclus") 
groupBy = "clusters"
orderGenes = "initial"
getUniprot = TRUE
silent = FALSE
coresGenes = 20


NCBI <- unname(unlist( foreach::foreach(NCBIid=result$Entrez.Gene.ID) %dopar% getNCBIentry(NCBIid) ))


for(NCBIid in result$Entrez.Gene.ID){
	
	print(NCBIid)
	data <- NULL
	if(!is.na(NCBIid)){
		url <- paste0("https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=full_report&list_uids=",
				NCBIid)
		webpage <- read_html(url)
		data_html <- html_nodes(webpage,'#summaryDl dd:nth-child(20)')
		data <- html_text(data_html)
		data <- gsub(" See more", "", data)
		data <- gsub("\n          human\n          all\n", "", data)
		data <- gsub(";", ",", data)
		
		rm(url, data_html, webpage)
	}
	#if(!S4Vectors::isEmpty(data)){
	#	return(data)
	#}else{
	#	return(NA)
	#}
}







