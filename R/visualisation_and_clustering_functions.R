

# getKEGGGenes
# Extracting genes from KEGG database by pathway ID.
# Returns only the genes that found in expression matrix.

# Please do not export this function because it does not work now.
.getKEGGGenes <- function(pathwayID, sceObject, species="mmu"){
	#
	
	keggOutput <- keggGet(paste0(species, pathwayID))[[1]]$GENE
	keggOutput <- keggOutput[seq(2, length(keggOutput), by = 2)]
	geneList <- unique(unname(sapply(keggOutput,
							function(x) strsplit(x, ";")[[1]][1])))
	filteredGeneList <- geneList[geneList %in% rownames(sceObject)]
	
	message(paste0("Taking ", length(filteredGeneList), " genes out of ",
					length(geneList), " from KEGG pathway: ", species, pathwayID,
					" ", keggGet(paste0(species, pathwayID))[[1]]$NAME))
	
	return(filteredGeneList)
	
}






















#' Save gene information into a table or tables for multiple inputs.
#'
#' This function runs conclus::getGenesInfo() function for all tables into the inputDir 
#' and saves the result into the outputDir.
#'
#' @param dataDirectory a directory with CONCLUS output. You can specify either 
#' dataDirectory, then inputDir and outputDir will be hardcoded, or inputDir and outputDir only.
#' The first is recommended during running CONCLUS workflow when the second option
#' is comfortable when you created input tables with genes manually.
#' @param inputDir input directory containing text files. These files can be obtained by 
#' applying conclus::saveMarkersLists() function or created manually. Each file must be a 
#' data frame with the first column called "geneName" containing gene symbols and (or) ENSEMBL IDs.
#' @param pattern a pattern of file names to take.
#' @param outputDir output directory.
#' @param databaseDir a path to the database "Mmus_gene_database_secretedMol.tsv". It is provided with the conclus package.
#' @param sep a parameter of read.delim() function.
#' @param header whether or not your input files have a header.
#' @param startFromFile number of the input file to start with. The function approaches files one by one.
#' It uses web scraping method to collect publicly available info from MGI, NCBI and UniProt websites.
#' Sometimes, if the Internet connection is not reliable, the function can drop. 
#' In this case, it is comfortable to start from the failed file and not to redo the previous ones.
#' @param groupBy a column in the input table used for grouping the genes in the output tables.
#' @param orderGenes if "initial" then the order of genes will not be changed.
#' @param getUniprot boolean, whether to get information from UniProt or not. Default is TRUE.
#' Sometimes, the connection to the website is not reliable. 
#' If you tried a couple of times and it failed, select FALSE. 
#' @param silent whether to show messages from intermediate steps or not.
#' @param coresGenes maximum number of jobs that the function can run in parallel.
#'
#' @return It saves text files either in the 'dataDirectory/marker_genes/saveGenesInfo' or outputDir 
#' depending on whether you specify dataDirectory or (inpitDir and outputDir) explicitly.
#' @export
saveGenesInfo <- function(dataDirectory = "",
                          inputDir = "", 
                          outputDir = "", 
                          pattern = "", 
                          databaseDir = system.file("extdata", package = "conclus"),
                          sep = ";", header = TRUE,
                          startFromFile = 1, #outputFormat = c("csv", "xlsx"),
                          groupBy = "clusters", # getGenesInfo params
                          orderGenes = "initial",
                          getUniprot = TRUE,
                          silent = FALSE, coresGenes = 20){

    if(dataDirectory != ""){
        inputDir = file.path(dataDirectory, "/marker_genes/markers_lists")
        outputDir = file.path(dataDirectory, "/marker_genes/saveGenesInfo")
        pattern = "markers.csv"
        dir.create(outputDir, showWarnings=F)
    }

    filesList <- list.files(inputDir, pattern = pattern)
    filePrefix <- do.call(rbind, strsplit(filesList, "[.]"))[,1]

    for(file in startFromFile:length(filesList)) {
        cat(file, "\n")
        genes <- read.delim(file.path(inputDir, filesList[file]),
                            sep = sep, header = header,
                            stringsAsFactors = FALSE)

        result <- getGenesInfo(genes, databaseDir,
                               groupBy = groupBy,
                               orderGenes = orderGenes,
                               getUniprot = getUniprot,
                               silent = silent, coresGenes = coresGenes)

        message("Writing the output file number ", file, "\n")

        #if(any(outputFormat %in% "csv")){
            write.table(result, file = file.path(outputDir,
                                                 paste0(filePrefix[file],
                                                        "_genesInfo.csv")),
                        quote = FALSE, sep = ";", row.names = FALSE)
        #}else if(any(outputFormat %in% "xlsx")){
            #write.xlsx2(result, file = file.path(outputDir,
            #                                     paste0(filePrefix[file],
            #                                            "_genesInfo.xlsx")),
            #            row.names=FALSE)
        #}
    }
}

#' Save top N marker genes for each cluster into a format suitable for conclus::saveGenesInfo() function.
#' 
#' The function takes the output files of conclus::rankGenes(), extracts top N markers and saves
#' them into the first "geneName" column of the output table. The second column "clusters" contains the 
#' name of the corresponding cluster.
#'
#' @param experimentName name of the experiment which appears at the beginning of the file name 
#' (supposed to be the same for one experiment during the workflow).
#' @param dataDirectory experiment directory (supposed to be the same for one experiment during the workflow).
#' @param inputDir input directory, usually "marker_genes" created automatically after conclus::runCONCLUS().
#' @param outputDir output directory.
#' @param pattern a pattern of the input file names to take.
#' @param Ntop number of top markers to take from each cluster.
#'
#' @return It saves files into the outputDir. The number of files is equal to the number of clusters.
#' @export
saveMarkersLists <- function(experimentName, dataDirectory,
                             inputDir = file.path(dataDirectory, "marker_genes"),
                             outputDir = file.path(dataDirectory,
                                                   paste0("marker_genes/markers_lists")),
                             pattern = "genes.tsv", Ntop = 100){

    dir.create(outputDir, showWarnings=F)

    filesList <- list.files(outputDir, pattern = "_markers.csv")
    deletedFiles <- sapply(filesList, function(fileName) 
				file.remove(file.path(outputDir, fileName)))
	
    fnames <- list.files(inputDir, pattern = pattern)
    for(i in 1:length(fnames)){
        tmp <- read.delim(file.path(inputDir, fnames[i]), stringsAsFactors = FALSE)
        markerList <- as.data.frame(tmp$Gene[1:Ntop])
        outputName <- gsub(paste0("_", pattern), "", fnames[i])
        clusterName <- gsub(paste0(experimentName, "_"), "", outputName)
        colnames(markerList) <- "geneName"
        markerList$clusters <- clusterName
        write.table(markerList, file.path(outputDir,
                                          paste0(outputName, "_markers.csv")),
                    row.names = FALSE, quote = FALSE, sep = ";")
    }
}
