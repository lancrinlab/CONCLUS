\name{addClusteringManually}
\alias{addClusteringManually}
\title{add modified Clustering to CONCLUS object}
\description{
The function replaces the content of the column "clusters" in the colData(sceObject) 
with the clustering provided in the user table.
The function will return the sceObject with cells which intersect with the cells from the input table.
}
\usage{
addClusteringManually(fileName, sceObject, dataDirectory, experimentName,
  columnName = "clusters")
}
\arguments{
\item{fileName}{A file with the clustering solution (for example, from previous CONCLUS runs).}

\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{dataDirectory}{Output directory (supposed to be the same for one experiment during the workflow).}

\item{experimentName}{Prefix used in runCONCLUS}

\item{columnName}{Name of the column with the clusters. Default = "clusters"}
}
\value{
A SingleCellExperiment object with the created/renewed column "clusters" in the colData(sceObject).
}
