# Do not export this function.
.normaliseCountMatrix <- function(countMatrix,
		species,
		method="default",
		sizes=c(20,40,60,80,100),
		rowData=NULL,
		colData=NULL,
		alreadyCellFiltered = FALSE,
		runQuickCluster = TRUE,
		databaseDir = TRUE){
	# Does normalisation of count matrix with.
	# There are 2 possible methods: "default" or "census"
	# The function returns SCE object with normalised count matrix
	if(method == "default"){
		rowData <- annotateGenes(countMatrix, species = species,
				rowData = rowData, databaseDir = databaseDir)
		colData <- addCellsInfo(countMatrix, rowData = rowData,
				colData = colData)
		if(!alreadyCellFiltered){
			filterCellsResult <- filterCells(countMatrix, colData)
			countMatrix <- filterCellsResult[[1]]
			colData <- filterCellsResult[[2]]
			rm(filterCellsResult)
		}
		filterGenesResult <- filterGenes(countMatrix, rowData)
		countMatrix <- filterGenesResult[[1]]
		rowData <- filterGenesResult[[2]]
		rm(filterGenesResult)
		
		stopifnot(all(rownames(countMatrix) == rownames(rowData)))
		stopifnot(all(colnames(countMatrix) == rownames(colData)))
		
		sce <-
				SingleCellExperiment::SingleCellExperiment(assays = list(counts = as.matrix(countMatrix)),
						colData=colData, rowData=rowData)
		
		# normalization
		message("Running normalization. It can take a while, depends on the
						number of cells.")
		if(runQuickCluster){
			cl <- tryCatch(scran::quickCluster(sce), error=function(e) NULL)
		}else{
			cl <- NULL
		}
		
		# compute sizeFactors which will be used for normalization
		sceNorm <- scran::computeSumFactors(sce, sizes = sizes, clusters = cl)
		
		message("summary(sizeFactors(sceObject)):")
		print(summary(SingleCellExperiment::sizeFactors(sceNorm)))
		if(length(SingleCellExperiment::sizeFactors(sceNorm)[SingleCellExperiment::sizeFactors(sceNorm) <= 0]) > 0){
			message("Cells with negative sizeFactors will be deleted before the
							downstream analysis.")
		}
		sceNorm <- sceNorm[, SingleCellExperiment::sizeFactors(sceNorm) > 0]
		sceNorm <- scater::normalize(sceNorm)
		rm(sce)
		
		return(sceNorm)
		
	}else if(method == "census"){
		message("Method 'census' is currently unavailable. Please select 'default'.")
		message("Unmodified count matrix returned.")
		return(countMatrix)
		#    sceObject <- normalize_dataset(as.matrix(countMatrix))
		#    SummarizedExperiment::colData(sceObject)$cellName = rownames(SummarizedExperiment::colData(sceObject))
		#    return(sceObject)
	}else{
		message("Wrong method. Unmodified count matrix returned.")
		return(countMatrix)
	}
}

#' normaliseCountMatrix
#'
#' Create a SingleCellExperiment object and perform normalization. The same as conclus::normalizeCountMatrix.
#'
#' @param countMatrix a matrix with non-normalised gene expression.
#' @param species either 'mmu' or 'human'.
#' @param method a method of clustering: available option is "default" using scran and scater.
#' @param sizes a vector of size factors from scran::computeSumFactors() function.
#' @param rowData a data frame with information about genes
#' @param colData a data frame with information about cells
#' @param alreadyCellFiltered if TRUE, cells quality check and filtering will not be applied. 
#' However, the function may delete some cells if they have negative size factors after scran::computeSumFactors.
#' @param runQuickCluster if scran::quickCluster() function must be applied.
#' Usually, it allows to improve normalization for medium-size count matrices. 
#' However, it is not recommended for datasets with less than 200 cells and
#' may take too long for datasets with more than 10000 cells.
#' @param databaseDir a path to annotation database provided with CONCLUS called 
#' "Mmus_gene_database_secretedMol.tsv" (only for MusMusculus 'mmu').
#' The function will work also without the database but slower because it will retrieve genes info from biomaRt.
#'
#' @return A SingleCellExperiment object with normalized gene expression, colData, and rowData.
#' @export

normaliseCountMatrix <- function(countMatrix,
		species,
		method="default",
		sizes=c(20,40,60,80,100),
		rowData=NULL,
		colData=NULL,
		alreadyCellFiltered = FALSE,
		runQuickCluster = TRUE,
		databaseDir = TRUE # FALSE for not using the database but download from biomaRt
){
	
	.normaliseCountMatrix(countMatrix,
			species,
			method=method,
			sizes=sizes,
			rowData=rowData,
			colData=colData,
			alreadyCellFiltered = alreadyCellFiltered,
			runQuickCluster = runQuickCluster,
			databaseDir = databaseDir)
}
# deleted from the description of the normaliseCountMatrix() function:
# #' @param method a method of clustering: "default" (using scran and scater) or "census" (using Census from Monocle).

#' @export
normalizeCountMatrix <- function(countMatrix,
		species,
		method="default",
		sizes=c(20,40,60,80,100),
		rowData=NULL,
		colData=NULL,
		alreadyCellFiltered = FALSE,
		runQuickCluster = TRUE,
		databaseDir = ""){
	
	.normaliseCountMatrix(countMatrix,
			species,
			method="default",
			sizes=c(20,40,60,80,100),
			rowData=NULL,
			colData=NULL,
			alreadyCellFiltered = FALSE,
			runQuickCluster = TRUE,
			databaseDir = "")
}
