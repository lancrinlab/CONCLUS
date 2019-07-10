\name{plotCellSimilarity}
\alias{plotCellSimilarity}
\title{Generate a cell similarity heatmap.}
\description{
This function plots similarity matrix as a heatmap, so one can see similarity between parts of different clusters.
}
\usage{
plotCellSimilarity(sceObject, cellsSimilarityMatrix, dataDirectory,
  experimentName, colorPalette = "default", statePalette = "default",
  clusteringMethod = "ward.D2", orderClusters = FALSE,
  plotPDF = TRUE, returnPlot = FALSE, width = 7, height = 6, ...)
}
\arguments{

\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{cellsSimilarityMatrix}{Output matrix from conclus::clusterCellsInternal() function.}

\item{dataDirectory}{Output directory for CONCLUS results.}

\item{experimentName}{Prefix used for output files.}

\item{colorPalette}{A vector of colors for clusters. Default = "default", see details.}

\item{statePalette}{a vector of colors for states or conditions.See details.}

\item{clusteringMethod}{Clustering method passed to hclust() function. See ?hclust for a list of method. Default = "ward.D2"}

\item{orderClusters}{If True, clusters in the similarity matrix of cells will be ordered by name. Default = FALSE}

\item{plotPDF}{If TRUE export to pdf, if FALSE export to png.}

\item{returnPlot}{Boolean indicating if the pHeatmap object should  be returned by the function. Default = FALSE.}

\item{width}{Width of the plot.}

\item{height}{Height of the plot.}

\item{...}{Additional parameters to pass to pdf() or png(), and pheatmap() functions.}
}
\value{
A ggplot object if returnPlot = TRUE. Saves the pdf in "dataDirectory/pictures" folder.
}
\details{
   colorPalette/statePalette -- A vector of colors for clusters/states or 'default' value. If 'default' is selected, the number of clusters is limited to 16. 
If an error message is thrown, re-run the function with your own color vector.
}
\seealso{
  	\code{\link{clusterCellsInternal}}
}
\examples{
 countMatrix <- read.delim(file.path(system.file("extdata", package = "conclus"),
                                    "Bergiers_counts_matrix_filtered.tsv"), 
                          stringsAsFactors = FALSE)
  
  columnsMetaData <- read.delim(file.path(system.file("extdata", package = "conclus"), 
                                "Bergiers_colData_filtered.tsv"))
  head(columnsMetaData)
  
  sceObject <- conclus::normaliseCountMatrix(countMatrix, species = "mmu", 
                                  colData = columnsMetaData)
                                  
  outputDirectory <- "MyOutputDirectory"
  experimentName <- "Bergiers"
  
  initialisePath(outputDirectory)
  # default parameters, can be selected by a user
  tSNEResults <- generateTSNECoordinates(sceObject, outputDirectory, 
                                          experimentName, PCs= c(4, 6, 8, 10, 20, 40, 50), 
                                          perplexities= c(30, 40),
                                          randomSeed = 42)
                                          
  epsilon=c(1.3, 1.4, 1.5)
  minPoints=c(3, 4)
  cores=14
  message("Running dbscan using ", cores, " cores.")
  dbscanResults <- conclus::runDBSCAN(tSNEResults, sceObject, outputDirectory, 
                           experimentName, epsilon=epsilon, 
                           minPoints=minPoints,
                           cores=cores)
  
  cellsSimilarityMatrix <- conclus::clusterCellsInternal(dbscanResults, sceObject, clusterNumber=10, 
                                  deepSplit= 0, cores=cores,
                                  clusteringMethod= "ward.D2")[[2]]
                                  
  

  clustersNumber <- length(unique(SummarizedExperiment::colData(sceObjectFiltered)$clusters))
  colorPalette <- conclus::choosePalette("default", clustersNumber)


  conclus::plotCellSimilarity(sceObjectFiltered, cellsSimilarityMatrix, outputDirectory,
                 experimentName, colorPalette = colorPalette, orderClusters = FALSE, 
                 statePalette = "default", clusteringMethod = "ward.D2",
                 plotPDF = TRUE,
                 returnPlot = TRUE)
}
