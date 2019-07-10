\name{generateTSNECoordinates}
\alias{generateTSNECoordinates}
\title{Generate and save t-SNE coordinates with selected parameters.}
\description{
The function generates several t-SNE coordinates based on given perplexity and ranges of PCs. 
Final number of t-SNE plots is length(PCs)*length(perplexities).
}
\usage{
generateTSNECoordinates(sceObject, dataDirectory, experimentName,
  randomSeed = 42, cores = 1, PCs = c(4, 6, 8, 10, 20, 40, 50),
  perplexities = c(30, 40))
}
\arguments{
\item{sceObject}{A SingleCellExperiment object containing experiment data. See ?conclus::normaliseCountMatrix.}

\item{dataDirectory}{Output directory for CONCLUS results (supposed to be the same for one experiment during the workflow).}

\item{experimentName}{Prefix used for output files.}

\item{randomSeed}{Random seed for reproducibility.}

\item{cores}{Maximum number of jobs that CONCLUS can run in parallel. Default = 1}

\item{PCs}{Vector of first principal components. For example, to take ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)}

\item{perplexities}{A vector of perplexity (t-SNE parameter). See details. Default = c(30, 40)}
}
\value{
An object with t-SNE results (coordinates for each plot). Also writes coordinates in "dataDirectory/tsnes" subfolder.
}
\details{

  Generates an object of fourteen (by default) tables with tSNE coordinates. Fourteen because it will vary seven values of principal components *PCs=c(4, 6, 8, 10, 20, 40, 50)* and two values of perplexity *perplexities=c(30, 40)* in 
all possible combinations. The chosen values of PCs and perplexities can be changed if necessary. We found that this combination works well for sc-RNA-seq datasets with 400-2000 cells. If you have 4000-9000 cells and expect more than 15 
clusters, we recommend you to take more first PCs and higher perplexity, for example, *PCs=c(8, 10, 20, 40, 50, 80, 100)* and *perplexities=c(200, 240)*. For details about perplexities parameter see ‘?Rtsne’.
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
}

\seealso{
    \code{\link{normaliseCountMatrix}}
}
