\name{plotClusteredTSNE}
\alias{plotClusteredTSNE}
\title{Plot t-SNE.}
\description{
Plot t-SNE generated with different PCs and perplexities. It can also use a coloring scheme by clusters or states.
}
\usage{
plotClusteredTSNE(sceObject, dataDirectory, experimentName,
  tSNEresExp = "", colorPalette = "default", PCs = c(4, 6, 8, 10, 20,
  40, 50), perplexities = c(30, 40), columnName = "clusters",
  returnPlot = FALSE, width = 6, height = 5, ...)
}
\arguments{
\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{dataDirectory}{Output directory for CONCLUS results.}

\item{experimentName}{Prefix used for output files.}

\item{tSNEresExp}{If t-SNE coordinates were generated in a different CONCLUS run, you can use them without renaming the files.
Path to the tsnes coordinates folder.}

\item{colorPalette}{A vector of colors for clusters. Default = "default", see details.}

\item{PCs}{A vector of first principal components.
For example, to take ranges 1:5 and 1:10 write c(5, 10). Default = c(4, 6, 8, 10, 20, 40, 50)}

\item{perplexities}{Numeric scalar defining the perplexity parameter, see ‘?Rtsne’ for more details. Default = c(30, 40)}

\item{columnName}{Name of the column to plot on t-SNE dimensions. Possible values are clusters, noColor, or state.}

\item{returnPlot}{Boolean indicating if the pHeatmap object should  be returned by the function. Default = FALSE.}

\item{width}{Width of the plot. Default = 6.}

\item{height}{Height of the plot. Default = 5.}

\item{...}{Additional parameters to pass to pdf() functions.}
}
\value{
A ggplot object if returnPlot = TRUE. It saves the pdf to "dataDirectory/pictures/tSNE_pictures".
}
\details{
   colorPalette -- A vector of colors for clusters/states or 'default' value. If 'default' is selected, the number of clusters is limited to 16. 
If an error message is thrown, re-run the function with your own color vector.
}
\examples{

  outputDirectory <- "./YourOutputDirectory" ## Replace by the path to the directory in which results should be written
  experimentName <- "Bergiers"
  countMatrix <- read.delim(file.path(system.file("extdata", package = "conclus"),
                                    "Bergiers_counts_matrix_filtered.tsv"), 
                          stringsAsFactors = FALSE)
  columnsMetaData <- read.delim(file.path(system.file("extdata", package = "conclus"), 
                                "Bergiers_colData_filtered.tsv"))
  initialisePath(outputDirectory)
  # default parameters, can be selected by a user
  PCs=c(4, 6, 8, 10, 20, 40, 50)
  perplexities=c(30, 40)
  randomSeed = 42
  epsilon=c(1.3, 1.4, 1.5)
  minPoints=c(3, 4)
  cores=1
  clusteringMethod="ward.D2"
  k=10 # parameter for cutree
  deepSplit = 0 # 0 to avoid cutreeDynamic, 1 to 4 to use it
  
  sceObject <- conclus::normaliseCountMatrix(countMatrix, species = "mmu", 
                                  colData = columnsMetaData)
  
  tSNEResults <- generateTSNECoordinates(sceObject, outputDirectory, 
                                          experimentName, PCs=PCs, 
                                          perplexities=perplexities,
                                          randomSeed = randomSeed)
  
  

  dbscanResults <- conclus::runDBSCAN(tSNEResults, sceObject, outputDirectory, 
                           experimentName, epsilon=epsilon, 
                           minPoints=minPoints,
                           cores=cores)
  
  cellsSimilarityMatrix <- conclus::clusterCellsInternal(dbscanResults, sceObject, 
                               clusterNumber=k, deepSplit=deepSplit, cores=cores,
                                  clusteringMethod=clusteringMethod)[[2]]
  
  tSNEclusters <- conclus::plotClusteredTSNE(sceObject, outputDirectory, experimentName,
                    PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
                    columnName = "clusters", returnPlot = TRUE)
  tSNEnoColor <- conclus::plotClusteredTSNE(sceObject, outputDirectory, experimentName,
                PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
                columnName = "noColor", returnPlot = TRUE)
                
  if(any(colnames(SummarizedExperiment::colData(sceObjectFiltered)) %in% "state")){
    tSNEstate <- conclus::plotClusteredTSNE(sceObject, outputDirectory, experimentName,
                    PCs=PCs, perplexities=perplexities, colorPalette = colorPalette,
                    columnName = "state", returnPlot = TRUE)

  tSNEclusters[[5]]
  tSNEnoColor[[5]]
  tSNEstate[[5]]
}

