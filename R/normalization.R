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




.filterCells <- function(countMatrix, colData, genesSumThr = 100,
		MoreBetter = c("genesNum", "sumCodPer", "genesSum"),
		MoreWorse = c("sumMtPer")){
	message("Running filterCells.")
	countMatrix <- countMatrix[,colSums(countMatrix) > genesSumThr]
	colData <- colData[colnames(countMatrix),]
	mb <- MoreBetter
	mw <- MoreWorse
	
	reportTable <- data.frame(matrix(NA, ncol = length(mb)+length(mw),
					nrow = nrow(colData)))
	colnames(reportTable) <- c(mb, mw)
	reportTable <- cbind(cellName = colData$cellName, reportTable)
	rownames(reportTable) <- reportTable$cellName
	
	stopifnot(all(colData$cellName==reportTable$cellName))
	
	for(j in 1:length(mb)){
		quan <- quantile(colData[,colnames(colData) == mb[j]])
		threshold <- 2.5*quan[2] - 1.5*quan[4]
		if(threshold < 0){
			threshold <- (quan[1]+quan[2]) / 2
		}
		reportTable[colData[,colnames(colData)==mb[j]] >= as.numeric(threshold),
				colnames(reportTable)==mb[j]] = 1
		reportTable[colData[,colnames(colData)==mb[j]] < as.numeric(threshold),
				colnames(reportTable)==mb[j]] = 0
	}
	for(j in 1:length(mw)){
		quan <- quantile(colData[,colnames(colData) == mw[j]])
		threshold <- 2.5*quan[4] - 1.5*quan[2]
		if(threshold > quan[5]){
			threshold <- (quan[3]+quan[4]) / 2
		}
		reportTable[colData[,colnames(colData)==mw[j]] <= as.numeric(threshold),
				colnames(reportTable)==mw[j]] = 1
		reportTable[colData[,colnames(colData)==mw[j]] > as.numeric(threshold),
				colnames(reportTable)==mw[j]] = 0
	}
	
	### add columns with filtering score and verdict ###
	reportTable <- dplyr::mutate(reportTable, score = NA)
	reportTable$score <- rowSums(reportTable[,colnames(reportTable) %in%
							c(mb,mw)])
	reportTable <- dplyr::mutate(reportTable, filterPassed = NA)
	reportTable$filterPassed[reportTable$score >= length(mb)+length(mw)] <- 1
	reportTable$filterPassed[reportTable$score < length(mb)+length(mw)] <- 0
	
	### add filtering verdict to colData ###
	colData <- dplyr::mutate(colData, filterPassed = NA)
	colData$filterPassed[colData$cellName %in%
					reportTable$cellName[reportTable$filterPassed == 1]] <- 1
	colData$filterPassed[colData$cellName %in%
					reportTable$cellName[reportTable$filterPassed == 0]] <- 0
	
	reportTable <- reportTable[order(reportTable$score, decreasing  =  FALSE), ]
	
	rm(threshold, j, mb, mw)
	
	rownames(colData) <- colData$cellName
	colData <- colData[colnames(countMatrix), ]
	stopifnot(all(rownames(colData) == colnames(countMatrix)))
	
	countMatrix <- countMatrix[,colData$filterPassed == 1]
	colData <- colData[colData$filterPassed == 1,]
	
	return(list(countMatrix, colData))
}


# from 2_mk_coldata_light.R
# I rewrote it with the new style
# This function creates colData or add columns mtGenes, genesNum,
# codGenes, genesSum, codSum, mtPer, codPer, sumMtPer, sumCodPer to the
# existing colData.
.addCellsInfo <- function(countMatrix, rowData, colData = NULL){
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
		colData <- .addCellsInfo(countMatrix, rowData = rowData,
				colData = colData)
		if(!alreadyCellFiltered){
			filterCellsResult <- .filterCells(countMatrix, colData)
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
