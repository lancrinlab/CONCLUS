\name{choosePalette}
\alias{choosePalette}
\title{Choose palette for a plot.}
\description{
Function used for choosing the palette for clusters.
}
\usage{
choosePalette(colorPalette, clustersNumber)
}
\arguments{
\item{colorPalette}{A vector of colors for clusters or "default", see details.}

\item{clustersNumber}{Number of clusters in the output palette.}
}
\value{
Color palette with the number of colors equal to the clusterNumber parameter.
}
\details{

  colorPalette -- A vector of colors for clusters/states or 'default' value. If 'default' is selected, the number of clusters is limited to 16. 
If an error message is thrown, re-run the function with your own color vector.
}