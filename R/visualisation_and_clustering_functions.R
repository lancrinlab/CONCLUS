

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
