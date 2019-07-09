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
  databaseDir = system.file("extdata", package = "conclus"))
}
\arguments{
\item{countMatrix}{Matrix with non-normalised gene expression.}

\item{species}{Currently limited to human and mouse. Possible values are either 'mmu' or 'human'.}

\item{sizes}{A numeric vector of pool sizes, i.e., number of cells per pool. See ?scran::computeSumFactors() for more details.}

\item{rowData}{a data frame with information about genes}

\item{colData}{a data frame with information about cells}

\item{alreadyCellFiltered}{if TRUE, cells quality check and filtering will not be applied. 
However, the function may delete some cells if they have negative size factors after scran::computeSumFactors.}

\item{runQuickCluster}{if scran::quickCluster() function must be applied.
Usually, it allows to improve normalization for medium-size count matrices. 
However, it is not recommended for datasets with less than 200 cells and
may take too long for datasets with more than 10000 cells.}

\item{databaseDir}{a path to annotation database provided with CONCLUS called 
"Mmus_gene_database_secretedMol.tsv" (only for MusMusculus 'mmu').
The function will work also without the database but slower because it will retrieve genes info from biomaRt.}
}
\value{
A SingleCellExperiment object with normalized gene expression, colData, and rowData.
}
