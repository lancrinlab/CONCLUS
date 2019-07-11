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

sceObject <- conclus::normaliseCountMatrix(countMatrix, species = "mmu", 
		colData = columnsMetaData)
p <- conclus::testClustering(sceObject, outputDirectory, experimentName)
p[[1]]
p[[3]]


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

tSNEclusters <- conclus::plotClusteredTSNE(sceObjectFiltered, outputDirectory, experimentName,
		PCs=PCs, perplexities=perplexities, colorPalette = "default",
		columnName = "clusters", returnPlot = TRUE)
tSNEnoColor <- conclus::plotClusteredTSNE(sceObjectFiltered, outputDirectory, experimentName,
		PCs=PCs, perplexities=perplexities, colorPalette = "default",
		columnName = "noColor", returnPlot = TRUE)

if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) %in% "state")){
	tSNEstate <- conclus::plotClusteredTSNE(sceObjectFiltered, outputDirectory, experimentName,
			PCs=PCs, perplexities=perplexities, colorPalette = "default",
			columnName = "state", returnPlot = TRUE)
}

tSNEclusters[[5]]
tSNEnoColor[[5]]
tSNEstate[[5]]


