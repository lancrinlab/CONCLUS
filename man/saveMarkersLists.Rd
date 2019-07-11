\name{saveMarkersLists}
\alias{saveMarkersLists}
\title{Save top N marker genes for each cluster}
\description{
The function takes the output files of conclus::rankGenes(), extracts top N markers and saves
them into the first "geneName" column of the output table. The second column "clusters" contains the 
name of the corresponding cluster.
}
\usage{
saveMarkersLists(experimentName, dataDirectory,
  inputDir = file.path(dataDirectory, "marker_genes"),
  outputDir = file.path(dataDirectory,paste0("marker_genes/markers_lists")), 
  pattern = "genes.tsv", Ntop = 100)
}
\arguments{
\item{experimentName}{Name of the experiment which appears at the beginning of the file name 
(supposed to be the same for one experiment during the workflow).}

\item{dataDirectory}{Experiment directory (supposed to be the same for one experiment during the workflow).}

\item{inputDir}{Input directory, usually "marker_genes" created automatically after conclus::runCONCLUS(). Default = file.path(dataDirectory, "marker_genes")}

\item{outputDir}{Output directory. Default = file.path(dataDirectory,paste0("marker_genes/markers_lists"))}

\item{pattern}{A pattern of the input file names to take. Default = "genes.tsv"}

\item{Ntop}{Number of top markers to take from each cluster. Default = 100}
}
\value{
It saves files into the outputDir/marker_genes/markers_lists. The number of files is equal to the number of clusters.
}
