# This function assumes that rownames(countMatrix) are either ENSEMBL IDs or
# or SYMBOLs. It will return a rowData with the same rownames as in countMatrix
# but genes which are not ENSEMBL IDs or SYMBOLs will not receive the annotation.
# Both types are possible in one matrix but not in one gene.
# Genes like Lrrc8a_ENSMUSG00000007476 will not be annotated.
# Iformation about cell surface localization will be used in the Shine application.
# It is useful to find cell surface markers. But this part of the function
# takes some time, so if you need this info, you can disable this option by
# cellSurface = FALSE to speed up the calculations.
.annotateGenes <- function(countMatrix, species = "mmu",
		genomeAnnot, ensemblPattern, rowData = NULL,
		databaseDir = TRUE){
	
	if(!databaseDir){
		
		if(missing(species) & (missing(genomeAnnot) | missing(ensemblPattern))){
			message("Species is either not selected or not equal to 'mmu' or 'human'.
							Please, select among default species or use options genomeAnnot and
							ensemblPattern. See example in the help page.")
			return(NULL)
		}else if(!missing(genomeAnnot) & !missing(ensemblPattern)){
			species = "manual"
		}else if(species == "mmu"){
			suppressMessages(library(org.Mm.eg.db, warn.conflicts = F))
			genomeAnnot <- org.Mm.eg.db
			ensemblPattern <- "ENSMUSG"
		}else if(species == "human"){
			suppressMessages(library(org.Hs.eg.db, warn.conflicts = F))
			genomeAnnot <- org.Hs.eg.db
			ensemblPattern <- "ENSG"
		}
		ensemblGenes <- rownames(countMatrix)[grep(ensemblPattern,
						rownames(countMatrix))]
		symbolGenes <- rownames(countMatrix)[!grepl(ensemblPattern,
						rownames(countMatrix))]
		message(paste0("Annotating ",length(ensemblGenes), " genes containing ",
						ensemblPattern, " pattern."))
		# if none of ENSEMBL genes are in database we fill their rowData rows with NA
		if(length(intersect(AnnotationDbi::keys(genomeAnnot,
								keytype = "ENSEMBL"),
						ensemblGenes)) == 0){
			rowdataEnsembl <- data.frame(ENSEMBL = ensemblGenes,
					SYMBOL = NA,
					GENENAME = NA)
		}else{
			rowdataEnsembl <- AnnotationDbi::select(genomeAnnot, keys=ensemblGenes,
					keytype="ENSEMBL",
					columns=c("SYMBOL",
							"GENENAME"),
					multiVals="first")
			rowdataEnsembl <- rowdataEnsembl[!duplicated(rowdataEnsembl$ENSEMBL),]
		}
		rowdataEnsembl$nameInCountMatrix <- ensemblGenes
		message("Annotating rest ", length(symbolGenes), " genes
						considering them as SYMBOLs.")
		rowdataSymbol <- AnnotationDbi::select(genomeAnnot, keys=symbolGenes,
				keytype="SYMBOL",
				columns=c("ENSEMBL",
						"GENENAME"),
				multiVals="first")
		rowdataSymbol <- rowdataSymbol[!duplicated(rowdataSymbol$SYMBOL),]
		rowdataSymbol$nameInCountMatrix <- symbolGenes
		rowdata <- base::rbind(rowdataSymbol, rowdataEnsembl)
		rm(rowdataEnsembl, rowdataSymbol)
		
		# sometimes several ensembls give one symbol
		(mult_sym <- rowdata$SYMBOL[!is.na(rowdata$SYMBOL) &
									duplicated(rowdata$SYMBOL)])
		
		# we decided don't combine such ensembls, but leave them unique with
		#"symbol_ensembl" ###
		(rowdata$SYMBOL[rowdata$SYMBOL %in% mult_sym] <-
					paste(rowdata$SYMBOL[rowdata$SYMBOL %in% mult_sym],
							rowdata$ENSEMBL[rowdata$SYMBOL %in% mult_sym],
							sep = "_"))
		rm(mult_sym)
		
		ensembl <- useMart(biomart = "ensembl", dataset="mmusculus_gene_ensembl")
		message("Retrieving information about genes from biomaRt.
						It can take up to five minutes, depends on Internet connection.")
		res <- getBM(attributes=c("ensembl_gene_id", "go_id", "name_1006",
						"chromosome_name", "gene_biotype"),
				mart=ensembl)
		tmp <- res[!duplicated(res$ensembl_gene_id),]
		rowdata <- merge(rowdata, tmp[c("ensembl_gene_id",
								"chromosome_name", "gene_biotype")],
				by.x = "ENSEMBL", by.y = "ensembl_gene_id",
				all.x = TRUE, all.y = FALSE, sort = FALSE)
		rowdata_GO <- merge(rowdata, res[c("ensembl_gene_id",
								"go_id", "name_1006")],
				by.x = "ENSEMBL", by.y = "ensembl_gene_id",
				all.x = TRUE, all.y = FALSE, sort = FALSE)
		rowdata_GO <- rowdata_GO[!is.na(rowdata_GO$name_1006) &
						((rowdata_GO$name_1006 == "cell surface") |
							(rowdata_GO$name_1006=="cell surface receptor signaling pathway")),]
		rowdata_GO$name_1006[duplicated(rowdata_GO$ENSEMBL)] <-
				"cell surface receptor signaling pathway"
		rowdata_GO <- rowdata_GO[!duplicated(rowdata_GO$ENSEMBL),]
		rowdata <- merge(rowdata, rowdata_GO[c("nameInCountMatrix",
								"go_id", "name_1006")],
				by.x = "nameInCountMatrix", by.y = "nameInCountMatrix",
				all.x = TRUE, all.y = TRUE, sort = FALSE)
		rm(tmp, ensembl, res, rowdata_GO)
	}else{
		# only for mouse
		databaseDir = system.file("extdata", package = "conclus")
		ensemblPattern <- "ENSMUSG"
		database <- read.delim(file.path(databaseDir,
						"Mmus_gene_database_secretedMol.tsv"),
				stringsAsFactors = FALSE)
		database <- database[!duplicated(database$Symbol),]
		
		ensemblGenes <- rownames(countMatrix)[grep(ensemblPattern,
						rownames(countMatrix))]
		ensemblGenesInternal <- gsub(paste0(".*_", ensemblPattern),
				ensemblPattern, ensemblGenes)
		symbolGenes <- rownames(countMatrix)[!grepl(ensemblPattern,
						rownames(countMatrix))]
		
		rowdataEnsembl <- data.frame(ensemblGenesInternal = ensemblGenesInternal,
				nameInCountMatrix = ensemblGenes)
		rowdataSymbol <- data.frame(nameInCountMatrix = symbolGenes)
		
		message(paste0("Annotating ",length(ensemblGenes), " genes containing ",
						ensemblPattern, " pattern."))
		
		rowdataEnsembl <- merge(rowdataEnsembl, database,
				by.x = "ensemblGenesInternal", by.y = "Ensembl",
				all.x = TRUE, all.y = FALSE, sort = FALSE)
		rowdataEnsembl <- rowdataEnsembl[,-1]
		
		rowdataEnsembl$Ensembl <- rowdataEnsembl$nameInCountMatrix
		rowdataSymbol$Symbol <- rowdataSymbol$nameInCountMatrix
		
		message("Annotating rest ", length(symbolGenes), " genes
						considering them as SYMBOLs.")
		
		rowdataSymbol <- merge(rowdataSymbol, database,
				by.x = "nameInCountMatrix", by.y = "Symbol",
				all.x = TRUE, all.y = FALSE, sort = FALSE)
		
		rowdata <- base::rbind(rowdataSymbol, rowdataEnsembl)
		colnames(rowdata)[colnames(rowdata) == "Ensembl"] <- "ENSEMBL"
		colnames(rowdata)[colnames(rowdata) == "Symbol"] <- "SYMBOL"
		colnames(rowdata)[colnames(rowdata) == "Name"] <- "GENENAME"
		colnames(rowdata)[colnames(rowdata) == "Feature.Type"] <- "gene_biotype"
	}
	
	if(!is.null(rowData)){
		rowData$nameInCountMatrix <- rownames(rowData)
		rowdata <- merge(rowData, rowdata,
				by.x = "nameInCountMatrix", by.y = "nameInCountMatrix",
				all.x = TRUE, all.y = TRUE, sort = FALSE)
	}
	
	rownames(rowdata) <- rowdata$nameInCountMatrix
	rowdata <- rowdata[rownames(countMatrix),]
	rowdata$SYMBOL[(S4Vectors::isEmpty(rowdata$SYMBOL)) | (rowdata$SYMBOL == "")] <- NA
	stopifnot(all(rownames(rowdata) == rownames(countMatrix)))
	
	return(rowdata)
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
		databaseDir = TRUE){
	# Does normalisation of count matrix with.
	# There are 2 possible methods: "default" or "census"
	# The function returns SCE object with normalised count matrix
	if(method == "default"){
		rowData <- .annotateGenes(countMatrix, species = species,
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
