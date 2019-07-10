\name{runDBSCAN}
\alias{runDBSCAN}
\title{Run clustering iterations with selected parameters using DBSCAN.}
\usage{
runDBSCAN(tSNEResults, sceObject, dataDirectory, experimentName,
  cores = 1, epsilon = c(1.3, 1.4, 1.5), minPoints = c(3, 4))
}
\description{
This function returns a matrix of clustering iterations of DBSCAN.
}
\arguments{
\item{tSNEResults}{Results of conclus::generateTSNECoordinates() function.}

\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{dataDirectory}{Output directory for CONCLUS results.}

\item{experimentName}{Prefix used for output files.}

\item{cores}{Maximum number of jobs that CONCLUS can run in parallel. Default = 1}

\item{epsilon}{Reachability distance parameter of fpc::dbscan() function. See Ester et al. (1996) for more details. Default = c(1.3, 1.4, 1.5)}

\item{minPoints}{Reachability minimum no. of points parameter of fpc::dbscan() function. See Ester et al. (1996) for more details. Default = c(3, 4)}
}
\value{
A matrix of DBSCAN results.
}
\details{

  Following the calculation of t-SNE coordinates, DBSCAN is run with a range of epsilon and MinPoints values which will yield a total of 84 clustering solutions (PCs x perplexities x MinPoints x epsilon). *minPoints* is the minimum 
cluster size which you assume to be meaningful for your experiment and *epsilon* is the radius around the cell where the algorithm will try to find *minPoints* dots. Optimal *epsilon* must lay one the knee of the k-NN function as 
shown in the "test_clustering/distance_graph.pdf"
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
}
