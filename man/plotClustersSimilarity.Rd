\name{plotClustersSimilarity}
\alias{plotClustersSimilarity}
\title{Generate a cluster similarity heatmap.}
\description{
Plot a similarity heatmap pulling the cells per cluster.
}
\usage{
plotClustersSimilarity(clustersSimilarityMatrix, sceObject, dataDirectory,
  experimentName, colorPalette = "default", statePalette = "default", 
  clusteringMethod = "ward.D2", returnPlot = FALSE, width = 7, height = 5.5, 
  ...)
}
\arguments{
\item{clustersSimilarityMatrix}{Result of conclus::calculateClustersSimilarity() function.}

\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{dataDirectory}{Output directory for CONCLUS results.}

\item{experimentName}{Prefix used for output files.}

\item{colorPalette}{A vector of colors for clusters. Default = "default", see details.}

\item{statePalette}{A vector of colors for states or conditions. Default = "default", See details.}

\item{clusteringMethod}{Clustering method passed to hclust() function. See ?hclust for a list of method. Default = "ward.D2".}

\item{returnPlot}{Boolean indicating if the pHeatmap object should  be returned by the function. Default = FALSE.}

\item{width}{Width of the plot. Default = 7.}

\item{height}{Height of the plot. Default = 5.5.}

\item{...}{Additional parameters to pass to pdf() and pheatmap() functions.}
}
\value{
A ggplot object if returnPlot = TRUE. It saves the pdf in "dataDirectory/pictures" folder.
}
\details{
   colorPalette/statePalette -- A vector of colors for clusters/states or 'default' value. If 'default' is selected, the number of clusters is limited to 16. 
If an error message is thrown, re-run the function with your own color vector.
}
\seealso{
  	\code{\link{calculateClustersSimilarity}}
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
                                  
  sceObjectFiltered <- sceObject
  
  clustersSimilarityMatrix <- 
      conclus::calculateClustersSimilarity(cellsSimilarityMatrix, 
          sceObject = sceObjectFiltered,
          clusteringMethod = "ward.D2")[[1]]
          
  conclus::plotClustersSimilarity(clustersSimilarityMatrix, 
                       sceObjectFiltered,
                       dataDirectory = outputDirectory, 
                       experimentName = experimentName, 
                       colorPalette = colorPalette,
                       statePalette = statePalette,
                       clusteringMethod = clusteringMethod,
                       returnPlot = TRUE)
}
