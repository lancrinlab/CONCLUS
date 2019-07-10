\name{exportData}
\alias{exportData}
\title{Export normalized expression matrix, columns and row data as the full workspace to the *output_tables* subfolder.}
\description{
The function allows you to export the expression matrix, columns and row data as the full workspace in the "dataDirectory/output_tables" directory for further analysis.
}
\usage{
exportData(sceObject, dataDirectory, experimentName)
}
\arguments{
\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{dataDirectory}{CONCLUS output directory for a given experiment (supposed to be the same for one experiment during the workflow).}

\item{experimentName}{Prefix used for output files.}
}
