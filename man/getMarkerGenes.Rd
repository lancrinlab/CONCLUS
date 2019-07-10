\name{getMarkerGenes}
\alias{getMarkerGenes}
\title{Get top N marker genes from each cluster.}
\description{
This function reads results of conclus::rankGenes() from "dataDirectory/marker_genes" and selects top N markers for each cluster.
}
\usage{
getMarkerGenes(dataDirectory, sceObject, experimentName, genesNumber = 10,
removeDuplicates = TRUE)
}
\arguments{
\item{dataDirectory}{Output directory used to generate CONCLUS results.}

\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{experimentName}{Prefix used to generate output files.}

\item{genesNumber}{Top N number of genes to get from one cluster. Default = 10}

\item{removeDuplicates}{Boolean indicating if duplicated genes should be deleted or not.}
}
\value{
A data frame where the first columns are marker genes ("geneName") and 
the second column is the associated group ("clusters").
}
\seealso{

 \code{\link{runCONCLUS}}
}
\examples{

  outputDirectory <- "./YourOutputDirectory" ## Replace by the path to the directory in which results should be written
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

  genesNumber <- 10
  markersClusters <- conclus::getMarkerGenes(outputDirectory, sceObjectCONCLUS, 
                                  experimentName = experimentName,
                                  genesNumber = genesNumber)
  
}

