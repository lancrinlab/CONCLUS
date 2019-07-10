\name{exportClusteringResults}
\alias{exportClusteringResults}
\title{exportClusteringResults}
\description{
The function saves clustering results into a table. Row names are cell names in the same order as in the sceObject.
}
\usage{
exportClusteringResults(sceObject, dataDirectory, experimentName, fileName)
}
\arguments{
\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{dataDirectory}{output directory (supposed to be the same for one experiment during the workflow).}

\item{experimentName}{Prefix used for output files.}

\item{fileName}{Name of the cluster table file.}
}
