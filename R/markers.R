#' Rank marker genes by statistical significance.
#'
#' This function searches marker genes for each cluster. It saves tables in the "dataDirectory/marker_genes" directory,
#' one table per cluster.
#' 
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param clustersSimilarityMatrix matrix, result of conclus::calculateClustersSimilarity() function.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @export
#' @param column name of the column with a clustering result.

rankGenes <- function(sceObject, clustersSimilarityMatrix, dataDirectory,
		experimentName, column="clusters"){
	
	markerGenesDirectory <- "marker_genes"
	expr <- Biobase::exprs(sceObject)
	colData <- SummarizedExperiment::colData(sceObject)
	simMed <- clustersSimilarityMatrix
	outputDir <- file.path(dataDirectory, markerGenesDirectory)
	
	message("Ranking marker genes for each cluster.")
	
	filesList <- list.files(outputDir, pattern = "_genes.tsv")
	deletedFiles <- sapply(filesList, function(fileName)
				file.remove(file.path(outputDir, fileName)))
	
	stopifnot(all(colnames(expr) == rownames(colData)))
	groups <- unique(colData[,column])
	simMed = simMed + 0.05
	for(i in 1:length(groups)){
		message(paste("Working on the file", i))
		tTestPval <- data.frame(Gene=rownames(expr))
		otherGroups <- groups[groups!=groups[i]]
		
		for(k in 1:length(otherGroups)){
			tTestPval[,paste0("vs_", otherGroups[k])] <- NA
			x <- expr[,colData[,c(column)] == groups[i]]
			y <- expr[,colData[,c(column)] == otherGroups[k]]
			t <- (rowMeans(x) - rowMeans(y)) /
					sqrt(apply(expr, 1, var)*(1/ncol(x) + 1/ncol(y)))
			df <- ncol(x) + ncol(y) - 2
			tTestPval[, paste0("vs_", otherGroups[k])] <-
					pt(t, df, lower.tail=FALSE)
		}
		
		tTestFDR <- data.frame(Gene=tTestPval$Gene)
		for(l in 1:length(otherGroups)){
			tTestFDR[,paste0("vs_", otherGroups[l])] <-
					p.adjust(tTestPval[,paste0("vs_", otherGroups[l])],
							method="fdr")
		}
		
		submat <- as.matrix(tTestFDR[,2:(length(otherGroups)+1)])
		tTestFDR$mean_log10_fdr <- rowMeans(log10(submat+1e-300))
		tTestFDR$n_05 <- apply(submat, 1, function(x)
					length(x[!is.na(x) & x < 0.05]))
		
		weights <- simMed[i, otherGroups]
		tTestFDR$score <- apply(submat, 1, function(x)
					sum(-log10(x+1e-300) * weights) / ncol(submat))
		
		tTestFDR <- tTestFDR[order(tTestFDR$score,decreasing=TRUE),]
		
		write.table(tTestFDR, file.path(outputDir,
						paste0(experimentName, "_cluster_",
								groups[i], "_genes.tsv")),
				col.names=TRUE, row.names=FALSE, quote=FALSE,
				sep="\t")
	}
	
	rm(tTestFDR, tTestPval, i, k, l, x, y, t, df, otherGroups, submat,
			weights, groups)
}




#' Get top N marker genes from each cluster. 
#' 
#' This function reads results of conclus::rankGenes() from "dataDirectory/marker_genes" and selects top N markers for each cluster.
#' 
#' @param dataDirectory output directory for a run of CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param genesNumber top N number of genes to get from one cluster.
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param removeDuplicates boolean, if duplicated genes must be deleted or not.
#'
#' @return A data frame where the first columns are marker genes ("geneName") and 
#' the second column is the groups ("clusters").
#' @export
getMarkerGenes <- function(dataDirectory, sceObject, genesNumber=14,
		experimentName, removeDuplicates = TRUE){
	
	markerGenesDirectory <- "marker_genes"
	numberOfClusters <- length(unique(SummarizedExperiment::colData(sceObject)$clusters))
	dir = file.path(dataDirectory, markerGenesDirectory)
	nTop = genesNumber
	clusters = unique(SummarizedExperiment::colData(sceObject)$clusters)
	
	markersClusters <- as.data.frame(matrix(ncol = 2,
					nrow = nTop*numberOfClusters))
	colnames(markersClusters) = c("geneName", "clusters")
	
	(fns <- list.files(dir, pattern = "_genes.tsv", full.names = FALSE))
	if(length(fns) != numberOfClusters){
		message(paste("Something wrong with number of files.
								It is supposed to be equal to number of clusters:", numberOfClusters))
		message(paste("Returning the marker genes from
								first", clusters, "clusters."))
		runUntil = numberOfClusters
	}else{
		runUntil = length(fns)
	}
	for(i in 1:runUntil){
		tmpAll <- read.delim(file.path(dir, paste(experimentName,
								"cluster", clusters[i], "genes.tsv", sep="_")),
				stringsAsFactors = FALSE)
		markersClusters$clusters[(nTop*(i-1)+1):(nTop*i)] <- as.character(clusters[i])
		markersClusters$geneName[(nTop*(i-1)+1):(nTop*i)] <- tmpAll$Gene[1:nTop]
	}
	if(removeDuplicates){
		markersClusters <- markersClusters[!duplicated(markersClusters$geneName),]
	}
	return(markersClusters)
}
