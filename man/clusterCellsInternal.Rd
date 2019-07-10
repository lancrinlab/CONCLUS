\name{clusterCellsInternal}
\alias{clusterCellsInternal}
\title{Cluster cells and get similarity matrix of cells.}
\description{
Returns consensus clusters by using hierarchical clustering on the similarity matrix of cells.
}
\usage{
clusterCellsInternal(dbscanMatrix, sceObject, clusterNumber = 0,
  deepSplit = 4, cores = 1, clusteringMethod = "ward.D2")
}
\arguments{
\item{dbscanMatrix}{An output matrix of conclus::runDBSCAN() function.}

\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{clusterNumber}{Exact number of cluster. Default = 0}

\item{deepSplit}{Intuitive level of clustering depth. Options are 1, 2, 3, 4. Default = 4}

\item{cores}{Maximum number of jobs that CONCLUS can run in parallel. Default = 1}

\item{clusteringMethod}{Clustering method passed to hclust() function. See ?hclust for a list of method. Default = "ward.D2".}
}
\value{
A SingleCellExperiment object with modified/created "clusters" column in the colData, and cells similarity matrix.
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
}
\seealso{
 	\code{\link{runDBSCAN}}
}
