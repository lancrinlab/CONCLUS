.generateAnnotationColors <- function(colData, colorPaletteParameter,
		statePalette){
	
	clusters <- levels(colData$clusters)
	states <- unique(colData$state)
	clusterNumber <- length(unique(colData$clusters))
	
	colorAnnotationClusters <- .choosePalette(colorPaletteParameter, clusterNumber)
	#colorAnnotationState <- chooseStatePalette(length(states))
	colorAnnotationState <- .choosePalette(statePalette, length(states))
	names(colorAnnotationState) <- states
	names(colorAnnotationClusters) <- clusters
	
	return(list(state=colorAnnotationState, clusters=colorAnnotationClusters))
}



.plotCellSimilarity <- function(sceObject, cellsSimilarityMatrix, dataDirectory,
		experimentName, colorPalette="default",
		statePalette="default", clusteringMethod="ward.D2",
		orderClusters = FALSE,
		plotPDF = TRUE,
		returnPlot = FALSE,
		width=7, height=6, onefile=FALSE, #pdf
		family, title, fonts, version,
		paper, encoding, bg, fg, pointsize,
		pagecentre, colormodel,
		useDingbats, useKerning,
		fillOddEven, compress,
		color = colorRampPalette(rev( #pheatmap
						RColorBrewer::brewer.pal(n = 7,
								name = "RdYlBu")))(100),
		kmeans_k = NA, breaks = NA,
		border_color = "grey60",
		cellwidth = NA, cellheight = NA,
		scale = "none",
		cluster_rows = FALSE, #not default
		cluster_cols = FALSE, #not default
		clustering_distance_rows = "euclidean",
		clustering_distance_cols = "euclidean",
		clustering_method = "complete",
		clustering_callback = identity2,
		cutree_rows = NA, cutree_cols = NA,
		treeheight_row = ifelse((
							class(cluster_rows) == "hclust") ||
						cluster_rows, 50, 0),
		treeheight_col = ifelse((
							class(cluster_cols) == "hclust") ||
						cluster_cols, 50, 0),
		legend = TRUE,
		legend_breaks = NA,
		legend_labels = NA,
		annotation_row = NA,
		annotation_col = NA,
		annotation = NA,
		annotation_colors = NA,
		annotation_legend = TRUE,
		annotation_names_row = TRUE,
		annotation_names_col = TRUE,
		drop_levels = TRUE,
		show_rownames = FALSE, #not default
		show_colnames = FALSE, #not default
		fontsize = 7.5, #not default (10)
		fontsize_row = 0.03,
		fontsize_col = fontsize,
		display_numbers = F,
		number_format = "%.2f",
		number_color = "grey30",
		fontsize_number = 0.8 * fontsize,
		gaps_row = NULL, gaps_col = NULL,
		labels_row = NULL, labels_col = NULL,
		filename = NA,
		widthHeatmap = NA, heightHeatmap = NA, #edited
		silent = FALSE,
		main = paste0("Cells similarity matrix ",
				ncol(cellsSimilarityMatrix),
				" columns, ",
				nrow(cellsSimilarityMatrix),
				" rows."),
		widthPNG = 800, heightPNG = 750 #png
){
	# plots cells correlation matrix gained form
	# clusterCellsInternal() function as the result of DBSCAN
	
	graphsDirectory <- "pictures"
	colData <- SummarizedExperiment::colData(sceObject)
	clustersNumber <- length(unique(colData$clusters))
	
	if(orderClusters == "name"){
		# Ordering expressionMatrixrix
		newOrder <- unname(unlist(sapply(levels(colData$clusters),
								function(cluster)
									.orderCellsInCluster(cluster,
											colData,
											expressionMatrix,
											clusteringMethod=clusteringMethod))))
		
		cellsSimilarityMatrix <- cellsSimilarityMatrix[newOrder, newOrder]
		cluster_cols <- FALSE
		cluster_rows <- FALSE
	}else if(orderClusters == FALSE){
		distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
		clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
		cluster_cols <- clusteringTree
		cluster_rows <- clusteringTree
	}else{
		message("Unknown option of orderClusters. Options are 'FALSE' or 'name'.
						Using the default version 'FALSE'.")
		distanceMatrix <- as.dist(sqrt((1-cellsSimilarityMatrix)/2))
		clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
		cluster_cols <- clusteringTree
		cluster_rows <- clusteringTree
	}
	
	annotationColors <- .generateAnnotationColors(colData, colorPalette,
			statePalette)
	columnsToPlot <- switch(is.null(colData$state) + 1, c("clusters", "state"),
			c("clusters"))
	
	pheatmapObject <- pheatmap::pheatmap(cellsSimilarityMatrix,
			show_colnames=show_colnames,
			show_rownames=show_rownames,
			annotation_col=as.data.frame(colData[columnsToPlot]),
			annotation_colors=annotationColors,
			fontsize_row=fontsize_row,
			cluster_cols=cluster_cols,
			cluster_rows=cluster_rows,
			fontsize=fontsize,
			main = main,
			color = color, kmeans_k = kmeans_k, # not changed by default
			breaks = breaks,
			border_color = border_color,
			cellwidth = cellwidth, cellheight = cellheight, scale = scale,
			clustering_distance_rows = clustering_distance_rows,
			clustering_distance_cols = clustering_distance_cols,
			clustering_method = clustering_method,
			clustering_callback = clustering_callback,
			treeheight_row = treeheight_row,
			treeheight_col = treeheight_col, legend = legend,
			legend_breaks = legend_breaks,
			legend_labels = legend_labels, annotation_row = annotation_row,
			annotation = annotation, annotation_legend = annotation_legend,
			annotation_names_row = annotation_names_row,
			annotation_names_col = annotation_names_col,
			drop_levels = drop_levels,
			fontsize_col = fontsize_col,
			display_numbers = display_numbers,
			number_format = number_format,
			number_color = number_color,
			fontsize_number = fontsize_numbere, gaps_row = gaps_row,
			gaps_col = gaps_col, labels_row = labels_row,
			labels_col = labels_col, filename = filename, width = widthHeatmap,
			height = heightHeatmap, silent = silent)
	
	if(plotPDF){
		pdf(file=file.path(dataDirectory, graphsDirectory, paste(experimentName,
								"cells_correlation", clustersNumber,"clusters.pdf", sep="_")),
				width=width, height=height, onefile=onefile, # not changed by default
				family=family, title=title, fonts=fonts, version=version,
				paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize,
				pagecentre=pagecentre, colormodel=colormodel,
				useDingbats=useDingbats, useKerning=useKerning,
				fillOddEven=fillOddEven, compress=compress)
	}else{
		message("Plot type is not pdf. Saving in png.")
		png(filename=file.path(dataDirectory, graphsDirectory,
						paste(experimentName,
								"cells_correlation", clustersNumber,
								"clusters.png", sep="_")),
				width = widthPNG, height = heightPNG, type = "cairo")
		#message("Plot type is not pdf. Saving in svg.")
		#svg(file.path(dataDirectory, graphsDirectory,
		#              paste(experimentName,
		#                    "cells_correlation", clustersNumber,
		#                    "clusters.svg", sep="_")),width=8,height=8)
	}
	grid::grid.newpage()
	grid::grid.draw(pheatmapObject$gtable)
	dev.off()
	
	if(returnPlot){
		return(pheatmapObject)
	}
}

#' Save a cells similarity matrix.
#' 
#' This function plots similarity matrix as a heatmap, so one can see similarity between parts of different clusters.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param cellsSimilarityMatrix an output matrix from the conclus::clusterCellsInternal() function.
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7"). 
#' @param statePalette "default" or a vector of colors for the column "state" in the colData, for example c("yellow", "#CC79A7"). 
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param orderClusters boolean, order clusters or not.
#' @param plotPDF if TRUE export to pdf, if FALSE export to png. 
#' FALSE is recommended for datasets with more than 2500 cells due to large pdf file size.
#' @param returnPlot boolean, return plot or not. Default if FALSE.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters of pdf(), pheatmap() and png() functions.
#'
#' @return A ggplot object or nothing (depends on the returnPlot parameter).
#' It saves the pdf in "dataDirectory/pictures" folder.
#' @export
plotCellSimilarity <- function(sceObject, cellsSimilarityMatrix, dataDirectory,
		experimentName, colorPalette="default",
		statePalette="default", clusteringMethod="ward.D2",
		orderClusters = FALSE,
		plotPDF = TRUE,
		returnPlot = FALSE,
		width=7, height=6, ...){
	
	.plotCellSimilarity(sceObject, cellsSimilarityMatrix, dataDirectory,
			experimentName, colorPalette=colorPalette,
			statePalette=statePalette, clusteringMethod=clusteringMethod,
			orderClusters = orderClusters,
			plotPDF = plotPDF,
			returnPlot = returnPlot,
			width=width, height=height, ...)
}
	

#' Plot t-SNE. Addtionally, it can highlight clusters or states.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param tSNEresExp if t-SNE coordinates were generated in a different CONCLUS run, you can use them without renaming the files.
#' Please copy tsnes folder from the source run to the current one and write that experimentName in the tSNEresExp argument.
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7").
#' @param PCs vector of PCs (will be specified in filenames).
#' @param perplexities vector of perplexities (will be specified in filenames).
#' @param columnName name of the column to plot on t-SNE dimensions.
#' @param returnPlot boolean, return plot or not.
#' @param width plot width.
#' @param height plot height.
#' @param ... other arguments of the pdf() function.
#'
#' @return A ggplot object or nothing (depends on the returnPlot parameter).
#' @export

plotClusteredTSNE <- function(sceObject, dataDirectory, experimentName,
		tSNEresExp = "",
		colorPalette = "default",
		PCs=c(4, 6, 8, 10, 20, 40, 50),
		perplexities=c(30, 40),
		columnName="clusters",
		returnPlot = FALSE,
		width=6, height=5, onefile=FALSE, ...){
	# plots picture based on t-SNE coordinates from
	# generateTSNECoordinates() and clustering results
	# from clusterCellsInternal() or runClustering()
	
	tSNEDirectory <- "tsnes"
	graphsDirectory <- "pictures"
	graphsTSNEDirectory <- "tSNE_pictures"
	
	if(tSNEresExp == ""){
		tSNEresExp <- experimentName
	}
	
	### Plot all precalculated pSNEs to show your clusters ###
	
	if(columnName == "noColor"){
		numberElements <- NULL
	}else{
		numberElements <- length(unique(SummarizedExperiment::colData(sceObject)[,columnName]))
		colorPalette <- .choosePalette(colorPalette, numberElements)
	}
	
	outputDir <- file.path(dataDirectory, graphsDirectory, graphsTSNEDirectory,
			paste("tSNE", numberElements, columnName, sep="_"))
	dir.create(outputDir, showWarnings = F)
	
	filesList <- list.files(outputDir, pattern = "_tsne_coordinates_")
	deletedFiles <- sapply(filesList, function(fileName) 
				file.remove(file.path(outputDir, fileName)))
	
	PCA <- rep(PCs, length(perplexities))
	perp <- rep(perplexities, each=length(PCs))
	
	tSNEplots <- rep(list(NA),(length(PCs)*length(perplexities)))
	
	for(i in 1:(length(PCs)*length(perplexities))){
		
		coordinatesName <- paste0(tSNEresExp, '_tsne_coordinates_', i, "_",
				PCA[i], "PCs_", perp[i], "perp")
		
		#fns <- list.files(file.path(dataDirectory, tSNEDirectory), full.names = TRUE)
		
		TSNEres <- read.delim(file.path(dataDirectory, tSNEDirectory,
						paste0(coordinatesName, ".tsv")),
				stringsAsFactors = FALSE)
		#TSNEres <- read.delim(fns[grepl(paste0(PCA[i], "PCs_", perp[i], "perp.tsv"), fns)],
		#                            stringsAsFactors = FALSE)
		
		TSNEres <- TSNEres[rownames(TSNEres) %in% SummarizedExperiment::colData(sceObject)$cellName, ]
		
		if(columnName != "noColor"){
			TSNEres[columnName] <- factor(SummarizedExperiment::colData(sceObject)[,columnName])
		}
		
		pdf(file.path(dataDirectory, graphsDirectory, graphsTSNEDirectory,
						paste("tSNE", numberElements, columnName, sep="_"),
						paste0(coordinatesName, ".pdf")),
				width=width, height=height, onefile=onefile, ...)
		if(columnName == "noColor"){
			tmp <- ggplot2::ggplot(TSNEres, aes_string(x=names(TSNEres)[1],
									y=names(TSNEres)[2])) +
					geom_point(size=I(1)) + theme_bw()
		}else{
			tmp <- ggplot2::ggplot(TSNEres, aes_string(x=names(TSNEres)[1],
									y=names(TSNEres)[2],
									color=columnName)) +
					geom_point(size=I(1)) +
					scale_color_manual(values=colorPalette) + theme_bw()
		}
		print(tmp)
		dev.off()
		tSNEplots[[i]] <- tmp
	}
	if(returnPlot){
		return(tSNEplots)
	}
	rm(PCA, perp)
}


.orderCellsInCluster <- function(cluster, colData, mtx,
		clusteringMethod="ward.D2"){
	# Order cells according to clustering results
	# Uses for ordering matrix to further plot it with pheatmap()
	
	cells <- colData[colData$clusters == cluster, ]$cellName
	if(length(cells) > 2){
		tree <- hclust(dist(t(mtx[, cells])), method=clusteringMethod)
		return(cells[tree$order])
	}else{
		return(cells)
	}
	
}

.orderGenesInCluster <- function(cluster, markersClusters, mtx,
		clusteringMethod="ward.D2"){
	# Order cells according to clustering results
	# Uses for ordering matrix to further plot it with pheatmap()
	
	genes <- markersClusters[markersClusters$clusters == cluster, ]$geneName
	if(length(genes) > 2){
		tree <- hclust(dist(mtx[genes, ]), method=clusteringMethod)
		return(genes[tree$order])
	}else{
		return(genes)
	}
	
	
}

.plotCellHeatmap <- function(markersClusters, sceObject, dataDirectory,
		experimentName,
		fileName, meanCentered=TRUE, colorPalette="default",
		statePalette="default", clusteringMethod="ward.D2",
		orderClusters = FALSE, #FALSE, TRUE, name, similarity
		orderGenes = FALSE, # FALSE, TRUE (will be ordered the same as clusters)
		returnPlot = FALSE,
		saveHeatmapTable = FALSE,
		width=10, height=8.5, onefile=FALSE, #pdf
		family, title, fonts, version,
		paper, encoding, bg, fg, pointsize,
		pagecentre, colormodel,
		useDingbats, useKerning,
		fillOddEven, compress,
		color = colorRampPalette(c("#023b84","#4b97fc",
						"#c9d9ef","#FEE395", #not default
						"#F4794E", "#D73027",#pheatmap
						"#a31008","#7a0f09"))(100),
		kmeans_k = NA, breaks = NA,
		border_color = "grey60",
		cellwidth = NA, cellheight = NA,
		scale = "none", cluster_rows = TRUE,
		cluster_cols = FALSE, # not original default
		clustering_distance_rows = "euclidean",
		clustering_distance_cols = "euclidean",
		clustering_method = "complete",
		clustering_callback = identity2,
		cutree_rows = NA, cutree_cols = NA,
		treeheight_row = ifelse((
							class(cluster_rows) == "hclust") ||
						cluster_rows, 50, 0),
		treeheight_col = ifelse((
							class(cluster_cols) == "hclust") ||
						cluster_cols, 50, 0),
		legend = TRUE,
		legend_breaks = NA,
		legend_labels = NA,
		annotation_row = NA,
		annotation_col = NA,
		annotation = NA,
		annotation_colors = NA,
		annotation_legend = TRUE,
		annotation_names_row = TRUE,
		annotation_names_col = TRUE,
		drop_levels = TRUE,
		show_rownames = TRUE,
		show_colnames = FALSE, # not original default (T)
		main = NA,
		fontsize = 7.5, # not original default (10)
		fontsize_row = 8, # not original default
		fontsize_col = fontsize,
		display_numbers = F,
		number_format = "%.2f",
		number_color = "grey30",
		fontsize_number = 0.8 * fontsize,
		gaps_row = NULL, gaps_col = NULL,
		labels_row = NULL, labels_col = NULL,
		filename = NA, widthHeatmap = NA,
		heightHeatmap = NA, silent = FALSE){
	# plots correlation between cells between clusters
	
	colData <- SummarizedExperiment::colData(sceObject)
	expressionMatrix <- Biobase::exprs(sceObject)[rownames(Biobase::exprs(sceObject)) %in%
					markersClusters$geneName,]
	
	if(meanCentered == TRUE){
		meanRows <- rowSums(expressionMatrix)/ncol(expressionMatrix)
		expressionMatrix <- expressionMatrix - meanRows
	}
	
	if(orderClusters == FALSE & orderGenes == TRUE){
		message("Genes cannot be ordered without clusters.
						Returning heatmap with orderCluster = similarity")
		orderClusters <- "similarity"
	}
	if(orderClusters == "name"){
		# Ordering expressionMatrixrix
		newOrder <- unname(unlist(sapply(levels(colData$clusters),
								function(cluster)
									.orderCellsInCluster(cluster,
											colData,
											expressionMatrix,
											clusteringMethod=clusteringMethod))))
		
		expressionMatrix <- expressionMatrix[, newOrder]
		cluster_cols <- FALSE
		if(orderGenes){
			newOrder <- unname(unlist(sapply(levels(colData$clusters),
									function(cluster)
										.orderGenesInCluster(cluster,
												markersClusters,
												expressionMatrix,
												clusteringMethod=clusteringMethod))))
			expressionMatrix <- expressionMatrix[newOrder, ]
			cluster_rows <- FALSE
		}
	}else if((orderClusters == TRUE) | (orderClusters == "similarity")){
		cellsSimilarityMatrix <- read.delim(file.path(dataDirectory,
						"output_tables",
						paste0(experimentName,
								"_cellsSimilarityMatrix.csv")),
				stringsAsFactors = FALSE,
				header = TRUE,
				sep = ",")
		
		clustersSimOrdered <- calculateClustersSimilarity(cellsSimilarityMatrix,
				sceObject,
				clusteringMethod)[[2]]
		
		newOrder <- unname(unlist(sapply(clustersSimOrdered,
								function(cluster)
									.orderCellsInCluster(cluster,
											colData,
											expressionMatrix,
											clusteringMethod=clusteringMethod))))
		
		expressionMatrix <- expressionMatrix[, newOrder]
		cluster_cols <- FALSE
		if(orderGenes){
			newOrder <- unname(unlist(sapply(clustersSimOrdered,
									function(cluster)
										.orderGenesInCluster(cluster,
												markersClusters,
												expressionMatrix,
												clusteringMethod=clusteringMethod))))
			expressionMatrix <- expressionMatrix[newOrder, ]
			cluster_rows <- FALSE
		}
	}else if(orderClusters == FALSE){
		distanceMatrix <- dist(t(expressionMatrix))
		cluster_cols <- hclust(distanceMatrix, method="ward.D2")
	}else{
		message("Unknown option of orderClusters. Options are 'TRUE' ('similarity'),
						'FALSE', 'name'. Using the default version 'FALSE'.")
		distanceMatrix <- dist(t(expressionMatrix))
		cluster_cols <- hclust(distanceMatrix, method="ward.D2")
	}
	
	if(!orderGenes){
		cluster_rows <- hclust(dist(expressionMatrix), method="ward.D2")
	}
	
	annotationColors <- .generateAnnotationColors(colData, colorPalette,
			statePalette)
	columnsToPlot <- switch(is.null(colData$state) + 1, c("clusters", "state"),
			c("clusters"))
	
	if(is.null(colData$clusters)){
		annCol <- switch(is.null(colData$state) + 1,
				as.data.frame(colData["state"]), NA)
		annColors <- switch(is.null(colData$state) + 1,
				annotationColors[names(annotationColors) == "state"],
				NA)
	}else{
		annCol <- as.data.frame(colData[columnsToPlot])
		annColors <- annotationColors
	}
	
	pheatmapObject <- pheatmap::pheatmap(expressionMatrix,
			show_colnames=show_colnames,
			show_rownames=show_rownames,
			annotation_col=annCol,
			annotation_colors=annColors,
			fontsize_row=fontsize_row,
			cluster_cols=cluster_cols,
			main=main,
			cluster_rows=cluster_rows,
			color=color,
			fontsize=fontsize, kmeans_k = kmeans_k, # not changed by default
			breaks = breaks,
			border_color = border_color,
			cellwidth = cellwidth, cellheight = cellheight, scale = scale,
			clustering_distance_rows = clustering_distance_rows,
			clustering_distance_cols = clustering_distance_cols,
			clustering_method = clustering_method,
			clustering_callback = clustering_callback,
			treeheight_row = treeheight_row,
			treeheight_col = treeheight_col, legend = legend,
			legend_breaks = legend_breaks,
			legend_labels = legend_labels,
			annotation = annotation,
			annotation_row = annotation_row,
			annotation_legend = annotation_legend,
			annotation_names_row = annotation_names_row,
			annotation_names_col = annotation_names_col,
			drop_levels = drop_levels,
			fontsize_col = fontsize_col,
			display_numbers = display_numbers,
			number_format = number_format,
			number_color = number_color,
			fontsize_number = fontsize_numbere, gaps_row = gaps_row,
			gaps_col = gaps_col, labels_row = labels_row,
			labels_col = labels_col, filename = filename, width = widthHeatmap,
			height = heightHeatmap, silent = silent)
	
	pdf(file.path(dataDirectory, "pictures",
					paste0(experimentName, "_", fileName, ".pdf")),
			width=width, height=height, onefile=onefile, # not changed by default
			family=family, title=title, fonts=fonts, version=version,
			paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize,
			pagecentre=pagecentre, colormodel=colormodel,
			useDingbats=useDingbats, useKerning=useKerning,
			fillOddEven=fillOddEven, compress=compress)
	grid::grid.newpage()
	grid::grid.draw(pheatmapObject$gtable)
	dev.off()
	
	if(saveHeatmapTable){
		exportMatrix(expressionMatrix, dataDirectory, experimentName, fileName)
	}
	
	if(returnPlot){
		return(pheatmapObject)
	}
}

#' Save markers heatmap.
#' 
#' This function plots heatmap with marker genes on rows and clustered cells on columns. 
#'
#' @param markersClusters a data frame where the first column is "geneName" containing genes names from sceObject, 
#' and the second column is corresponding "clusters". All names from that column must come from the column "clusters" in the colData(sceObject).
#' The data frame can be obtained from conclus::getMarkerGenes() function or created manually.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory of a given CONCLUS run (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param fileName name of the ouput file
#' @param meanCentered boolean, should mean centering be applied to the expression data or not.
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7").
#' @param statePalette "default" or a vector of colors for the column "state" in the colData, for example c("yellow", "#CC79A7").
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param orderClusters boolean, should the heatmap be structured by clusters.
#' @param orderGenes boolean, should the heatmap be structured by genes.
#' @param returnPlot boolean, whether to return a ggplot object with the plot or not.
#' @param saveHeatmapTable boolean, whether to save the expression matrix used for heatmap into a .csv file or not.
#' The file will be saved into 'dataDirectory/output_tables' with the same name as the .pdf plot.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters from pdf() and pheatmap() functions.
#'
#' @return A ggplot object of the plot if needed. The function saves pdf in "dataDirectiry/pictures" folder.
#' @export
plotCellHeatmap <- function(markersClusters, sceObject, dataDirectory,
		experimentName,
		fileName, meanCentered=TRUE, colorPalette="default",
		statePalette="default", clusteringMethod="ward.D2",
		orderClusters = FALSE, #FALSE, TRUE, name, similarity
		orderGenes = FALSE, # FALSE, TRUE (will be ordered the same as clusters)
		returnPlot = FALSE,
		saveHeatmapTable = FALSE,
		width=10, height=8.5, ...){
	
	.plotCellHeatmap(markersClusters, sceObject, dataDirectory,
			experimentName,
			fileName, meanCentered=meanCentered, colorPalette=colorPalette,
			statePalette=statePalette, clusteringMethod=clusteringMethod,
			orderClusters = orderClusters, #FALSE, TRUE, name, similarity
			orderGenes = orderGenes, # FALSE, TRUE (will be ordered the same as clusters)
			returnPlot = returnPlot,
			saveHeatmapTable = saveHeatmapTable,
			width=width, height=height, ...)
}




#' plotGeneExpression
#'
#' The function saves a t-SNE plot colored by expression of a given gene. 
#' Warning: filename with t-SNE results is hardcoded, so please don't rename the output file.
#'
#' @param geneName name of the gene you want to plot.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which appears in filenames (supposed to be the same for one experiment during the workflow).
#' @param graphsDirectory name of the subdirectory where to put graphs. Default is "dataDirectory/pictures".
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param tSNEpicture number of the picture you want to use for plotting. 
#' Please check "dataDirectory/tsnes" or "dataDirectory/pictures/tSNE_pictures/clusters" to get the number, it is usually from 1 to 14.
#' @param commentName comment you want to specify in the filename.
#' @param palette color palette for the legend.
#' @param returnPlot boolean, should the function return a ggplot object or not.
#' @param savePlot boolean, should the function export the plot to pdf or not.
#' @param alpha opacity of the points of the plot.
#' @param limits range of the gene expression shown in the legend.
#' This option allows generating t-SNE plots with equal color
#' scale to compare the expression of different genes. By default, limits are the range
#' of expression of a selected gene.
#' @param pointSize size of the point.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters of the pdf() function.
#'
#' @return A ggplot object of the plot if needed.
#' @export

plotGeneExpression <- function(geneName, experimentName, dataDirectory,
		graphsDirectory = "pictures",
		sceObject, tSNEpicture=1,
		commentName = "", palette = c("grey","red",
				"#7a0f09",
				"black"),
		returnPlot = FALSE,
		savePlot = TRUE,
		alpha = 1, limits = NA,
		pointSize = 1,
		width=6, height=5, ...){
	
	experimentName <- experimentName
	dataDirectory <- dataDirectory
	tSNEDirectory <- "tsnes"
	
	### Plot all precalculated t-SNEs to show your clusters ###
	
	clustersNumber <- length(unique(SummarizedExperiment::colData(sceObject)$clusters))
	
	coordsName <- list.files(file.path(dataDirectory, tSNEDirectory),
			pattern = paste0(experimentName,'_tsne_coordinates_',
					tSNEpicture, "_"))
	
	tSNECoords <- read.delim(file.path(dataDirectory, tSNEDirectory, coordsName),
			stringsAsFactors=FALSE)
	
	tSNECoords <- tSNECoords[SummarizedExperiment::colData(sceObject)$cellName, ]
	
	if(!geneName %in% rownames(Biobase::exprs(sceObject))){
		print("Gene is not found in expression matrix")
	}
	
	stopifnot(all(rownames(tSNECoords) == colnames(sceObject)))
	tSNECoords$expression <- Biobase::exprs(sceObject)[geneName, ]
	
	if(length(limits) == 1){
		limits <- c(min(tSNECoords$expression),
				max(tSNECoords$expression))
	}
	
	if(savePlot){
		pdf(file.path(dataDirectory, graphsDirectory, paste0(paste(experimentName,
										"tSNE", clustersNumber, "clusters" , geneName, commentName,
										"tSNEpicture", tSNEpicture, "_alpha", alpha,
										sep="_"), ".pdf")),
				width=width, height=height, ...)
	}
	tmp <- ggplot2::ggplot(tSNECoords, aes(x=tSNECoords[,1],
							y=tSNECoords[,2], color=expression)) +
			geom_point(size=I(pointSize), alpha = alpha) + theme_bw() +
			scale_colour_gradientn(colours=alpha(colorRampPalette(palette)(100), 0.8),
					limits = limits) +
			ggtitle(geneName)
	#brewer.pal(9, "OrRd")[0:9]
	print(tmp)
	
	if(savePlot){
		dev.off()
	}
	
	if(returnPlot){
		return(tmp)
	}
}




.plotClustersSimilarity <- function(clustersSimilarityMatrix, sceObject,
		dataDirectory,
		experimentName, colorPalette,
		statePalette,
		clusteringMethod,
		returnPlot = FALSE,
		width=7, height=5.5, onefile=FALSE, #pdf
		family, title, fonts, version,
		paper, encoding, bg, fg, pointsize,
		pagecentre, colormodel,
		useDingbats, useKerning,
		fillOddEven, compress,
		color = colorRampPalette(rev( #pheatmap
						RColorBrewer::brewer.pal(n = 7,
								name = "RdYlBu")))(100),
		kmeans_k = NA, breaks = NA,
		border_color = "grey60",
		cellwidth = NA, cellheight = NA,
		scale = "none", cluster_rows = TRUE,
		cluster_cols = TRUE,
		clustering_distance_rows = "euclidean",
		clustering_distance_cols = "euclidean",
		clustering_method = "complete",
		clustering_callback = identity2,
		cutree_rows = NA, cutree_cols = NA,
		treeheight_row = ifelse((
							class(cluster_rows) == "hclust") ||
						cluster_rows, 50, 0),
		treeheight_col = ifelse((
							class(cluster_cols) == "hclust") ||
						cluster_cols, 50, 0),
		legend = TRUE,
		legend_breaks = NA,
		legend_labels = NA,
		annotation_row = NA,
		annotation_col = NA,
		annotation = NA,
		annotation_colors = NA,
		annotation_legend = TRUE,
		annotation_names_row = TRUE,
		annotation_names_col = TRUE,
		drop_levels = TRUE,
		show_rownames = TRUE,
		show_colnames = TRUE,
		fontsize = 7.5, #not default (10)
		fontsize_row = fontsize,
		fontsize_col = fontsize,
		display_numbers = F,
		number_format = "%.2f",
		number_color = "grey30",
		fontsize_number = 0.8 * fontsize,
		gaps_row = NULL, gaps_col = NULL,
		labels_row = NULL, labels_col = NULL,
		filename = NA, widthHeatmap = NA,
		heightHeatmap = NA, silent = FALSE,
		main =
				paste0("Clusters similarity matrix ",
						ncol(clustersSimilarityMatrix),
						" columns, ", nrow(clustersSimilarityMatrix),
						" rows.")){
	
	clusters <- SummarizedExperiment::colData(sceObject)$clusters
	clustersNumber <- length(unique(clusters))
	clustersNames <- levels(clusters)
	dataDirectory <- dataDirectory
	experimentName <- experimentName
	graphsDirectory <- "pictures"
	
	# Plotting matrix
	distanceMatrix <- as.dist(sqrt((1-clustersSimilarityMatrix)/2))
	clusteringTree <- hclust(distanceMatrix, method=clusteringMethod)
	
	colDataSimilarity <- data.frame(clusters = clustersNames)
	rownames(colDataSimilarity) <- colDataSimilarity$clusters
	
	annotationColors <- .generateAnnotationColors(SummarizedExperiment::colData(sceObject),
			colorPalette,
			statePalette)
	
	pheatmapObject <- pheatmap::pheatmap(clustersSimilarityMatrix,
			show_colnames=show_colnames,
			show_rownames=show_rownames,
			annotation_col=colDataSimilarity,
			annotation_colors=annotationColors,
			cluster_cols=clusteringTree,
			cluster_rows=clusteringTree,
			fontsize=fontsize,
			main = main,
			color = color, kmeans_k = kmeans_k, # not changed by default
			breaks = breaks,
			border_color = border_color,
			cellwidth = cellwidth, cellheight = cellheight, scale = scale,
			clustering_distance_rows = clustering_distance_rows,
			clustering_distance_cols = clustering_distance_cols,
			clustering_method = clustering_method,
			clustering_callback = clustering_callback,
			treeheight_row = treeheight_row,
			treeheight_col = treeheight_col, legend = legend,
			legend_breaks = legend_breaks,
			legend_labels = legend_labels, annotation_row = annotation_row,
			annotation = annotation, annotation_legend = annotation_legend,
			annotation_names_row = annotation_names_row,
			annotation_names_col = annotation_names_col,
			drop_levels = drop_levels,
			fontsize_row = fontsize_row, fontsize_col = fontsize_col,
			display_numbers = display_numbers,
			number_format = number_format,
			number_color = number_color,
			fontsize_number = fontsize_numbere, gaps_row = gaps_row,
			gaps_col = gaps_col, labels_row = labels_row,
			labels_col = labels_col, filename = filename, width = widthHeatmap,
			height = heightHeatmap, silent = silent)
	
	message("\nSaving a heatmap with the clusters similarity matrix.")
	pdf(file.path(dataDirectory, graphsDirectory,
					paste(experimentName,"clusters_similarity", clustersNumber,
							"clusters.pdf", sep="_")),
			width=width, height=height, onefile=onefile, # not changed by default
			family=family, title=title, fonts=fonts, version=version,
			paper=paper, encoding=encoding, bg=bg, fg=fg, pointsize=pointsize,
			pagecentre=pagecentre, colormodel=colormodel,
			useDingbats=useDingbats, useKerning=useKerning,
			fillOddEven=fillOddEven, compress=compress)
	
	grid::grid.newpage()
	grid::grid.draw(pheatmapObject$gtable)
	dev.off()
	
	if(returnPlot){
		return(pheatmapObject)
	}
}

#' Save a similarity cluster matrix.
#'
#' @param clustersSimilarityMatrix a matrix, result of conclus::calculateClustersSimilarity() function.
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory for CONCLUS (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param colorPalette "default" or a vector of colors for the column "clusters" in the colData, for example c("yellow", "#CC79A7").
#' @param statePalette "default" or a vector of colors for the column "state" in the colData, for example c("yellow", "#CC79A7").
#' @param clusteringMethod a clustering methods passed to hclust() function.
#' @param returnPlot boolean, return plot or not.
#' @param width plot width.
#' @param height plot height.
#' @param ... other parameters of pdf() and pheatmap() functions.
#'
#' @return A ggplot object or nothing (depends on returnPlot parameter). It saves the pdf in "dataDirectory/pictures" folder.
#' @export
plotClustersSimilarity <- function(clustersSimilarityMatrix, sceObject,
		dataDirectory,
		experimentName, colorPalette,
		statePalette,
		clusteringMethod,
		returnPlot = FALSE,
		width=7, height=5.5, ...) {
	
	.plotClustersSimilarity(clustersSimilarityMatrix, sceObject,
			dataDirectory, 
			experimentName, colorPalette,
			statePalette,
			clusteringMethod,
			returnPlot = returnPlot,
			width=width, height=height, ...)
}
