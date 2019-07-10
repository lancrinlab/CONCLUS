\name{plotGeneExpression}
\alias{plotGeneExpression}
\title{plotGeneExpression}
\description{
The function saves a t-SNE plot colored by expression of a given gene. 
Warning: filename with t-SNE results is hardcoded, so please don't rename the output file.
}
\usage{
plotGeneExpression(geneName, experimentName, dataDirectory, sceObject,
  graphsDirectory = "pictures", tSNEpicture = 1,
  commentName = "", palette = c("grey", "red", "#7a0f09", "black"),
  returnPlot = FALSE, savePlot = TRUE, alpha = 1, limits = NA,
  pointSize = 1, width = 6, height = 5, ...)
}
\arguments{
\item{geneName}{Name of the gene to highlight on the t-SNE plot.}

\item{experimentName}{Prefix used for output files (supposed to be the same than the one used with runCONCLUS).}

\item{dataDirectory}{Output directory of a given CONCLUS run (supposed to be the same than the one used with runCONCLUS).}

\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{graphsDirectory}{Name of the output subdirectory. Default is "dataDirectory/pictures".}

\item{tSNEpicture}{Number of picture that you want to use for plotting. 
Please check "dataDirectory/tsnes" or "dataDirectory/pictures/tSNE_pictures/clusters" to get the number which corresponds to the number of files, it is usually from 1 to 14. Default = 1}

\item{commentName}{Comment that you want to specify in the filename.}

\item{palette}{Color palette for the legend.}

\item{returnPlot}{Boolean, should the function return a ggplot object or not. Default = FALSE}

\item{savePlot}{Boolean, should the function export the plot to pdf or not. Default = TRUE}

\item{alpha}{Opacity of the points of the plot. Default = 1}

\item{limits}{Range of the gene expression shown in the legend. See details.}

\item{pointSize}{Size of the point. Default = 1}

\item{width}{Width of the plot. Default = 6.}

\item{height}{Height of the plot. Default = 5.}

\item{...}{other parameters of the pdf() function.}
}
\value{
A ggplot object of the plot if needed.
}
\details{

  limits -- This option allows generating t-SNE plots with equal color scale to 
  compare the expression of different genes. By default, limits are the range
  of expression of a selected gene.
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
  
  plotGeneExpression("Cd34", experimentName, dataDirectory, sceObjectCONCLUS, tSNEpicture = 10)

}

