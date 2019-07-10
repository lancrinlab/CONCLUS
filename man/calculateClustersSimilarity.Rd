\name{calculateClustersSimilarity}
\alias{calculateClustersSimilarity}
\title{Pool cell similarity information into clusters similarity.}
\description{
Having computed cells similarity, pools information into clusters.
}
\usage{
calculateClustersSimilarity(cellsSimilarityMatrix, sceObject,
  clusteringMethod = "ward.D2")
}
\arguments{
\item{cellsSimilarityMatrix}{A similarity matrix given by conclus::clusterCellsInternal() function.}

\item{sceObject}{a SingleCellExperiment object.}

\item{clusteringMethod}{Clustering method passed to hclust() function. See ?hclust for a list of methods. Default = "ward.D2"}
}
\value{
A list containing the cluster similarity matrix and cluster names (order).
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
                                  
  sceObjectFiltered <- sceObject
  
  clustersSimilarityMatrix <- 
      conclus::calculateClustersSimilarity(cellsSimilarityMatrix, 
          sceObject = sceObjectFiltered,
          clusteringMethod = "ward.D2")[[1]]
}

