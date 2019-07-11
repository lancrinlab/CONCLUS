\name{saveGenesInfo}
\alias{saveGenesInfo}
\title{Save gene information into a table or tables for multiple inputs.}
\description{
This function runs conclus::getGenesInfo() function for all tables into the inputDir 
and saves the result into the outputDir.
}
\usage{
saveGenesInfo(dataDirectory = "", inputDir = "", outputDir = "",
  pattern = "", databaseDir = system.file("extdata", package =
  "conclus"), sep = ";", header = TRUE, startFromFile = 1,
  groupBy = "clusters", orderGenes = "initial", getUniprot = TRUE,
  silent = FALSE, cores = 1)
}
\arguments{
\item{dataDirectory}{A directory with CONCLUS output. You can specify either 
dataDirectory, then inputDir and outputDir will be hardcoded, or inputDir and outputDir only.
The first is recommended during running CONCLUS workflow when the second option
is comfortable when you created input tables with genes manually. Default = ""}

\item{inputDir}{Input directory containing text files. These files can be obtained by 
applying conclus::saveMarkersLists() function or created manually. Each file must be a 
data frame with the first column called "geneName" containing gene symbols and (or) ENSEMBL IDs. Default = ""}

\item{outputDir}{Output directory. Default = ""}

\item{pattern}{A pattern of file names to take. Default = ""}

\item{databaseDir}{A path to the database "Mmus_gene_database_secretedMol.tsv". It is provided with the conclus package. Default = system.file("extdata", package =
  "conclus")}

\item{sep}{A parameter of read.delim() function. Default = ";"}

\item{header}{Whether or not your input files have a header. Default = TRUE}

\item{startFromFile}{Number of the input file to start with. The function approaches files one by one.
It uses web scraping method to collect publicly available info from MGI, NCBI and UniProt websites.
Sometimes, if the Internet connection is not reliable, the function can drop. 
In this case, it is comfortable to start from the failed file and not to redo the previous ones. Default = 1}

\item{groupBy}{A column in the input table used for grouping the genes in the output tables. Default = "clusters"}

\item{orderGenes}{If "initial" then the order of genes will not be changed. Other option is "alphabetical". Default = "initial"}

\item{getUniprot}{Boolean, whether to get information from UniProt or not. Default=TRUE}

\item{silent}{Whether to show messages from intermediate steps or not. Default = FALSE}

\item{cores}{Maximum number of jobs that the function can run in parallel. Default = 1}
}
\value{
It saves text files either in the 'dataDirectory/marker_genes/saveGenesInfo' or outputDir 
depending on whether you specify dataDirectory or (inpitDir and outputDir) explicitly.
}
