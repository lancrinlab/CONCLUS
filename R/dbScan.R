### This function calculates dbscan for all t-SNE from TSNEtables with all
### combinations of paramenters from epsilon and minPoints
### it does not set random seed. It allows to vary this parameter automatically
### it returns a matrix where columns are iterations
### number of iterations is equal to
### ncol is (TSNEtables)*length(epsilon)*length(epsilon)

.mkDbscan <- function(TSNEtables, cores = 14, epsilon = c(1.2, 1.5, 1.8),
		minPoints = c(15, 20)){
	myCluster <- parallel::makeCluster(cores, # number of cores to use
			type = "PSOCK") # type of cluster
	doParallel::registerDoParallel(myCluster)
	dbscanResults <- foreach::foreach(i=rep(rep(1:ncol(TSNEtables),
									each=length(minPoints)),
							length(epsilon)),
					eps=rep(epsilon,
							each=ncol(TSNEtables)*length(minPoints)),
					MinPts=rep(minPoints,
							ncol(TSNEtables)*length(epsilon)),
					.combine='cbind') %dopar% {
				fpc::dbscan(TSNEtables[1,i][[1]], eps=eps,
						MinPts=MinPts)$cluster
			}
	parallel::stopCluster(myCluster)
	return(dbscanResults)
}

#' Run clustering iterations with selected parameters using DBSCAN.
#'
#' This function returns a matrix of clustering iterations of DBSCAN.
#'
#' @param tSNEResults results of conclus::generateTSNECoordinates() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames 
#' (supposed to be the same for one experiment during the workflow).
#' @param cores maximum number of jobs that CONCLUS can run in parallel.
#' @param epsilon a fpc::dbscan() parameter.
#' @param minPoints a fpc::dbscan() parameter.
#'
#' @return A matrix of DBSCAN results.
#' @export
runDBSCAN <- function(tSNEResults, sceObject, dataDirectory, experimentName,
		cores=14, epsilon=c(1.3, 1.4, 1.5), minPoints=c(3, 4)){
	
	outputDataDirectory <- "output_tables"
	# taking only cells from the sceObject
	for(i in 1:ncol(tSNEResults)){
		tmp <- tSNEResults[1,i][[1]]
		tSNEResults[1,i][[1]] <- tmp[colnames(sceObject),]
	}
	dbscanResults <- .mkDbscan(tSNEResults, cores = cores, epsilon = epsilon,
			minPoints = minPoints)
	dbscanResults <- t(dbscanResults)
	colnames(dbscanResults) <- SummarizedExperiment::colData(sceObject)$cellName
	write.table(dbscanResults, file.path(dataDirectory, outputDataDirectory,
					paste0(experimentName, "_dbscan_results.tsv")), quote = FALSE,
			sep = "\t")
	rm(tmp)
	return(dbscanResults)
}
