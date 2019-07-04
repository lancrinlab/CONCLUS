

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

















# from 2_mk_coldata_light.R
# I rewrote it with the new style
# This function creates colData or add columns mtGenes, genesNum,
# codGenes, genesSum, codSum, mtPer, codPer, sumMtPer, sumCodPer to the
# existing colData.
addCellsInfo <- function(countMatrix, rowData, colData = NULL){
    message("Adding cell info for cells filtering.")
    coldata <- data.frame(cellName = colnames(countMatrix),
                          stringsAsFactors = FALSE)

    ### add info about all genes in a cell ###
    coldata <- dplyr::mutate(coldata, genesNum = NA, genesSum = NA, oneUMI = NA)
    coldata$genesSum <- colSums(countMatrix)
    for(i in 1:ncol(countMatrix)){
        vec <- countMatrix[,i]
        coldata$genesNum[coldata$cellName == colnames(countMatrix)[i]] <-
            length(vec[vec > 0])
        coldata$oneUMI[coldata$cellName == colnames(countMatrix)[i]] <-
            length(vec[vec == 1])
    }
    rm(vec)
    coldata <- dplyr::mutate(coldata,
                            oneUMIper = 100 * coldata$oneUMI / coldata$genesNum)

    ### add info about mitochondrial and protein-coding genes ###
    coldata <- dplyr::mutate(coldata,
                             mtGenes = NA, mtSum = NA,
                             codGenes = NA, codSum = NA)
    for(i in 1:ncol(countMatrix)){
        mt <- countMatrix[rownames(countMatrix) %in%
                            rowData$nameInCountMatrix[rowData$chromosome_name ==
                                                                        "MT"],i]
        coldata$mtGenes[coldata$cellName == colnames(countMatrix)[i]] <-
            length(mt[mt > 0])
        coldata$mtSum[coldata$cellName == colnames(countMatrix)[i]]   <- sum(mt)

        cod <- countMatrix[rownames(countMatrix) %in%
                    rowData$nameInCountMatrix[rowData$gene_biotype ==
                                                            "protein_coding"],i]
        coldata$codGenes[coldata$cellName == colnames(countMatrix)[i]] <-
            length(cod[cod > 0])
        coldata$codSum[coldata$cellName == colnames(countMatrix)[i]] <- sum(cod)
    }

    coldata <- dplyr::mutate(coldata,
                    mtPer      = 100 * coldata$mtGenes  / coldata$genesNum,
                    codPer     = 100 * coldata$codGenes / coldata$genesNum,
                    sumMtPer  = 100 * coldata$mtSum    / coldata$genesSum,
                    sumCodPer = 100 * coldata$codSum   / coldata$genesSum)
    rm(mt, cod)

    if(!is.null(colData)){
        colData$cellName <- rownames(colData)
        coldata <- merge(colData, coldata,
                         by.x = "cellName", by.y = "cellName",
                         all.x = FALSE, all.y = TRUE, sort = FALSE)
    }

    rownames(coldata) <- coldata$cellName
    coldata <- coldata[colnames(countMatrix),]
    stopifnot(all(rownames(coldata) == colnames(countMatrix)))

    return(coldata)
}


filterGenes <- function(countMatrix, rowData){

    # internal function, filters genes which are more than in 10 cells and less than (all-10) cells

    selRows <- ((rowSums(countMatrix[,] >= 1)) > 10)
    countMatrix <- countMatrix[selRows,]
    rowData <- rowData[rowData$nameInCountMatrix %in% rownames(countMatrix),]

    return(list(countMatrix, rowData))
}


#' Collect genes information to one table.
#'
#' The function takes a data frame containing gene symbols and (or) ENSEMBL IDs and returns
#' a data frame with such information as gene name, feature type, chromosome,
#' gene IDs in different annotations, knockout information from MGI, a summary from NCBI 
#' and UniProt, and whether or not a gene belongs to GO terms containing proteins on the cell surface or 
#' involved in secretion.
#'
#' @param genes a data frame with the first column called "geneName" containing gene symbols and (or) ENSEMBL IDs.
#' Other columns are optional. For example, the second column could be "clusters" with the name of the cluster 
#' for which the gene is a marker.
#' @param databaseDir a path to the database provided with CONCLUS called "Mmus_gene_database_secretedMol.tsv".
#' @param groupBy a column in the input table used for grouping the genes in the output tables.
#' This option is useful if a table contains genes from different clusters.
#' @param orderGenes if "initial" then the order of genes will not be changed.
#' @param getUniprot boolean, whether to get information from UniProt or not. Default is TRUE.
#' Sometimes, the connection to the website is not reliable. 
#' If you tried a couple of times and it failed, select FALSE. 
#' @param silent whether to show messages from intermediate steps or not.
#' @param coresGenes maximum number of jobs that the function can run in parallel.
#'
#' @return Returns a data frame.
#' @export
getGenesInfo <- function(genes, databaseDir = system.file("extdata", package = "conclus"), 
                         groupBy = "clusters",
                         orderGenes = "initial",
                         getUniprot = TRUE,
                         silent = FALSE, coresGenes = 20){

    # MGI
    getMGIentry <- function(MGIid){
        library('rvest')
        library(S4Vectors)
        data <- NULL
        if(!is.na(MGIid)){
            url <- paste0("http://www.informatics.jax.org/marker/",
                          MGIid)
            webpage <- read_html(url)
            data_html <- html_nodes(webpage,'#mpMarkerClip')
            data <- html_text(data_html)
            data <- gsub("\t", "", data, perl = TRUE)
            data <- gsub("\n", "", data, perl = TRUE)
            data <- gsub(";", ",", data)

            rm(url, data_html, webpage)
        }
        if(!S4Vectors::isEmpty(data)){
            return(data)
        }else{
            return(NA)
        }
    }

    #NCBI
    getNCBIentry <- function(NCBIid){
        library('rvest')
        library(S4Vectors)
        data <- NULL
        if(!is.na(NCBIid)){
            url <- paste0("https://www.ncbi.nlm.nih.gov/gene?cmd=Retrieve&dopt=full_report&list_uids=",
                          NCBIid)
            webpage <- read_html(url)
            data_html <- html_nodes(webpage,'#summaryDl dd:nth-child(20)')
            data <- html_text(data_html)
            data <- gsub(" See more", "", data)
            data <- gsub("\n          human\n          all\n", "", data)
            data <- gsub(";", ",", data)

            rm(url, data_html, webpage)
        }
        if(!S4Vectors::isEmpty(data)){
            return(data)
        }else{
            return(NA)
        }
    }

    # Uniprot
    getUniprotEntry <- function(UniprotID){
        library('rvest')
        library(S4Vectors)
        data <- NULL
        if(!is.na(UniprotID)){
            url <- paste0("https://www.uniprot.org/uniprot/",
                          UniprotID, "#function")
            webpage <- read_html(url)
            data_html <- html_nodes(webpage,'h2+ .annotation > span')
            data <- html_text(data_html)
            data <- gsub("* Publication.*", " Publications", data)
            data <- gsub('<p><a href="/manual/evidences#ECO:0000250">More...</a></p>',
                         "", data)
            data <- gsub('\r\n .* <p>Manually curated information .* toiUniProtKB:.* (.*_.*).',
                         "", data)
            data <- gsub('\r\n .* <p>Manual validated information.*',
                         "", data)
            data <- gsub('UniRule annotation', "", data)
            data <- gsub(";", ",", data)

            rm(url, data_html, webpage)
        }
        if(!S4Vectors::isEmpty(data)){
            return(data)
        }else{
            return(NA)
        }
    }

    library(DataCombine)
    library(doParallel)
    #database <- read.delim(file.path(databaseDir, "Mmus_gene_database.tsv"),
    #                       stringsAsFactors = FALSE)
    database <- read.delim(file.path(databaseDir, "Mmus_gene_database_secretedMol.tsv"),
                           stringsAsFactors = FALSE)
    database <- database[!duplicated(database$Symbol),]

    genesEnsembl <- genes[grepl("ENSMUSG", genes$geneName),]
    genesSymbol <- genes[!grepl("ENSMUSG", genes$geneName),]

    resultEnsembl <- merge(genesEnsembl, database,
                           by.x = "geneName", by.y = "Ensembl",
                           all.x = TRUE, all.y = FALSE, sort = FALSE)
    resultEnsembl$Ensembl <- resultEnsembl$geneName

    resultSymbol <- merge(genesSymbol, database,
                          by.x = "geneName", by.y = "Symbol",
                          all.x = TRUE, all.y = FALSE, sort = FALSE)
    resultSymbol$Symbol <- resultSymbol$geneName

    colnames(resultEnsembl)
    colnamesOrder = c("geneName", "clusters", "Ensembl", "Symbol", "Name",
                      "Feature.Type", "MGI.Gene.Marker.ID", "Entrez.Gene.ID",
                      "Uniprot.ID", "chromosome_name", "go_id", "name_1006")
    resultEnsembl <- resultEnsembl[,colnamesOrder]
    resultSymbol <- resultSymbol[,colnamesOrder]

    result <- rbind(resultSymbol, resultEnsembl)
    #colnames(result)[11:12] <- c("CellSurface.GOid", "CellSurface.GOname")

    rownames(result) <- result$geneName
    result <- result[genes$geneName,]

    if(orderGenes == "alphabetical"){
        result <- result[order(result$geneName),]
    }

    myCluster <- parallel::makeCluster(coresGenes, # number of cores to use
                             type = "PSOCK") # type of cluster
    doParallel::registerDoParallel(myCluster)

    # MGI vector
    if(!silent){
        message("Collecting knockout phenotype information from MGI.")
    }

    MGI <- unname(unlist(foreach::foreach(MGIid=result$MGI.Gene.Marker.ID) %dopar% getMGIentry(MGIid)))

    #NCBI vector
    if(!silent){
        message("Retrieving info from NCBI.")
    }

    NCBI <- unname(unlist( foreach::foreach(NCBIid=result$Entrez.Gene.ID) %dopar% getNCBIentry(NCBIid) ))

    if(getUniprot){
        # Uniprot
        if(!silent){
            message("Getting summary from Uniprot.")
        }
        UniprotFunction <- unname(unlist( foreach::foreach(UniprotID=result$Uniprot.ID) %dopar% getUniprotEntry(UniprotID) ))
    }

    parallel::stopCluster(myCluster)

    #MGI <- data.frame(MGI, stringsAsFactors = FALSE)
    #NCBI <- data.frame(NCBI, stringsAsFactors = FALSE)
    #UniprotFunction <- data.frame(UniprotFunction, stringsAsFactors = FALSE)

    if(getUniprot){
        result <- cbind(result, MGI, NCBI, UniprotFunction)
    }else{
        result <- cbind(result, MGI, NCBI)
    }

    result$MGI <- as.character(result$MGI)
    result$NCBI <- as.character(result$NCBI)
    if(getUniprot){
        result$UniprotFunction <- as.character(result$UniprotFunction)
    }

    rownames(result) <- c(1:nrow(result))

    # inserting space for comments
    if(any(colnames(result) %in% groupBy) &
       (orderGenes == "initial") &
       (length(unique(result$clusters)) > 1)){
        resultFinal <- result
        groupingTable <- table(resultFinal[,groupBy])
        groupingTable <- groupingTable[unique(resultFinal$clusters)]
        resultFinal <- InsertRow(resultFinal, c("For notes:",
                                                rep("", (ncol(result))-1)),
                                 RowNum = 1)

        RowNum <- groupingTable[1] + 1
        for(i in 1:(length(groupingTable)-1)){
            resultFinal <- InsertRow(resultFinal, rep("", ncol(result)),
                                     RowNum = (RowNum + 1))
            resultFinal <- InsertRow(resultFinal, c("For notes:",
                                                    rep("", (ncol(result))-1)),
                                     RowNum = (RowNum + 2))
            RowNum <- RowNum + 2 + groupingTable[i+1]
        }
        result <- resultFinal
        rm(resultFinal)
    }

    rm(resultEnsembl, resultSymbol, database, colnamesOrder)

    #result <- result[,c("geneName", "clusters", "Name",
    #                    "Feature.Type",  "CellSurface.GOid",
    #                    "CellSurface.GOname", "MGI", "NCBI",
    #                    "UniprotFunction", "chromosome_name",
    #                    "Symbol", "Ensembl", "MGI.Gene.Marker.ID",
    #                    "Entrez.Gene.ID", "Uniprot.ID")]
    if(getUniprot){
        result <- result[,c("geneName", "clusters", "Name",
                        "Feature.Type",  "go_id",
                        "name_1006", "MGI", "NCBI",
                        "UniprotFunction", "chromosome_name",
                        "Symbol", "Ensembl", "MGI.Gene.Marker.ID",
                        "Entrez.Gene.ID", "Uniprot.ID")]
    }else{
        result <- result[,c("geneName", "clusters", "Name",
                        "Feature.Type",  "go_id",
                        "name_1006", "MGI", "NCBI",
                        "chromosome_name", # no "UniprotFunction"
                        "Symbol", "Ensembl", "MGI.Gene.Marker.ID",
                        "Entrez.Gene.ID", "Uniprot.ID")]
    }
    return(result)
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
