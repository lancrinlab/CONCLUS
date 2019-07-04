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
									orderCellsInCluster(cluster,
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
	
	annotationColors <- generateAnnotationColors(colData, colorPalette,
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
