\name{normaliseCountMatrix}
\alias{normaliseCountMatrix}
\title{normaliseCountMatrix}
\description{
Create a SingleCellExperiment object and perform normalization.
}
\usage{
normaliseCountMatrix(countMatrix, species,
  sizes = c(20, 40, 60, 80, 100), rowData = NULL, colData = NULL,
  alreadyCellFiltered = FALSE, runQuickCluster = TRUE,
  databaseDir = TRUE)
}
\arguments{
\item{countMatrix}{Matrix with non-normalised gene expression.}

\item{species}{Currently limited to human and mouse. Possible values are either 'mmu' or 'human'.}

\item{sizes}{A numeric vector of pool sizes, i.e., number of cells per pool. See ?scran::computeSumFactors() for more details.}

\item{rowData}{A data frame with Ensembl IDs for each row of the count matrix}

\item{colData}{A data frame with information about cells. See details.}

\item{alreadyCellFiltered}{if TRUE, cells quality check and filtering will not be applied. 
However, the function may delete some cells if they have negative size factors after scran::computeSumFactors.}

\item{runQuickCluster}{if scran::quickCluster() function must be applied.
Usually, it allows to improve normalization for medium-size count matrices. 
However, it is not recommended for datasets with less than 200 cells and
may take too long for datasets with more than 10000 cells.}

\item{databaseDir}{Boolean indicating if the default database for mouse should be used. See details. Default = TRUE}
}
\value{
A SingleCellExperiment object with normalized gene expression, colData, and rowData.
}
\details{

  This function uses the normalization method of the scater package. For more details about the normalization used see 
?scater::normalize. The size factors used in the normalization are calculated with scran::computeSumFactors.
 
  colData -- Dataframe containing three columns: cellName, state, and cellBarcode.
  
  databaseDir -- If TRUE, uses the path to the annotation database provided with CONCLUS called "Mmus_gene_database_secretedMol.tsv". 
This database only handles MusMusculus 'mmu'. The function will work also without the database (FALSE) but slower because it will 
retrieve genes info from biomaRt. 

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
}
\seealso{
  \code{\link[scater]{normalize}}
  \code{\link[scran]{computeSumFactors}}
  \code{\link[scran]{quickCluster}}
}
