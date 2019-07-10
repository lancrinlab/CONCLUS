\name{exportMatrix}
\alias{exportMatrix}
\title{Export matrix to a file.}
\description{
The function allows you to export a matrix to a .csv file with a hard-coded filename (according to experimentName) 
in the "dataDirectory/output_tables" directory for further analysis.
}
\usage{
exportMatrix(matrix, dataDirectory, experimentName, name)
}
\arguments{
\item{matrix}{Your matrix (e.g., expression matrix)}

\item{dataDirectory}{CONCLUS output directory for a given experiment (supposed to be the same for one experiment during the workflow).}

\item{experimentName}{Prefix used for output files.}

\item{name}{Name of the file. Will be placed after the experimentName header.}
}
