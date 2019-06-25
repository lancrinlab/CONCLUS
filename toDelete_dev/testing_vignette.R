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

source("/home/descostes/git/CONCLUS/R/visualisation_and_clustering_functions.R")


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


## Test clustering

p <- testClustering(sceObject, dataDirectory, experimentName)
p[[1]]
p[[3]]






































































#------------------------------ test code for testClustering

dbscanEpsilon=1.4
minPts=5
perplexities = c(30)
PCs = c(4)
randomSeed = 42
width=7
height=7
onefile=FALSE


tSNE <- getTSNEresults(expr = Biobase::exprs(sceObject), cores=1,
		perplexities = perplexities,
		PCs = PCs,
		randomSeed = randomSeed)



expressionMatrix = Biobase::exprs(sceObject)
cores=1
PCs=PCs
perplexities=perplexities
randomSeed=42

PCAData <- prcomp(t(expressionMatrix))$x
myCluster <- parallel::makeCluster(cores, # number of cores to use
		type = "PSOCK") # type of cluster
doParallel::registerDoParallel(myCluster)





#------------------------------------------------