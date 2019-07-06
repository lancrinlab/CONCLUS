

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
