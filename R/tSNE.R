#' Generate and save t-SNE coordinates with selected parameters.
#'
#' The function generates several t-SNE coordinates based on given perplexity and ranges of PCs. 
#' Final number of t-SNE plots is length(PCs)*length(perplexities)
#' It writes coordinates in "dataDirectory/tsnes" subfolder.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param randomSeed random seed for reproducibility.
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param PCs a vector of first principal components.
#' For example, to take ranges 1:5 and 1:10 write c(5, 10).
#' @param perplexities a vector of perplexity (t-SNE parameter).
#'
#' @return An object with t-SNE results (coordinates for each plot).
#' @export
generateTSNECoordinates <- function(sceObject, dataDirectory, experimentName,
		randomSeed=42, cores=14,
		PCs=c(4, 6, 8, 10, 20, 40, 50),
		perplexities=c(30,40)){
	
	tSNEDirectory <- "tsnes"
	message(paste0("Running TSNEs using ", cores, " cores."))
	TSNEres <- .getTSNEresults(Biobase::exprs(sceObject), cores=cores, PCs=PCs,
			perplexities=perplexities, randomSeed=randomSeed)
	
	PCA <- rep(PCs, length(perplexities))
	perp <- rep(perplexities, each=length(PCs))
	
	outputDir <- file.path(dataDirectory, tSNEDirectory)
	filesList <- list.files(outputDir, pattern = "_tsne_coordinates_")
	deletedFiles <- sapply(filesList, function(fileName) 
				file.remove(file.path(outputDir, fileName)))
	
	for (i in 1:(length(PCs)*length(perplexities))){
		write.table(TSNEres[1, i][[1]],
				file=file.path(dataDirectory, tSNEDirectory,
						paste0(experimentName,'_tsne_coordinates_', i, "_" ,
								PCA[i], "PCs_", perp[i], "perp.tsv")),
				quote=FALSE, sep='\t')
	}
	saveRDS(TSNEres, file = file.path(dataDirectory, "output_tables",
					paste0(experimentName,"_tSNEResults.rds")))
	rm(tSNEDirectory, PCA, perp)
	return(TSNEres)
}
