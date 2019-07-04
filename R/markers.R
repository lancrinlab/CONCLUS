# rankGenesInternal() saves n files in your outputDir, n=number of groups
rankGenesInternal <- function(expr, colData, column, simMed, outputDir,
		experimentName){
	
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
	# 
	
	markerGenesDirectory <- "marker_genes"
	rankGenesInternal(Biobase::exprs(sceObject), SummarizedExperiment::colData(sceObject), column,
			clustersSimilarityMatrix,
			file.path(dataDirectory, markerGenesDirectory), experimentName)
	
}
