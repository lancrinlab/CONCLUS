.plotDistanceGraphWithEpsilon <- function(tSNEData, minNeighbours=5,
		epsilon=1.2){
	# similar function as plotDistanceGraph,
	# but with already known epsilon value
	#
	
	dbscan::kNNdistplot(tSNEData, k=minNeighbours)
	abline(h=epsilon, lty=2)
	
}


.plotTestClustering <- function(tSNEData, minNeighbours=5,
		epsilon=1.2){
	# plots test DBSCAN on one of the pictures
	# for being ensured that clustering will
	# probably work successful
	
	dbscanResults <- fpc::dbscan(tSNEData, eps=epsilon, MinPts=minNeighbours)
	factoextra::fviz_cluster(dbscanResults, tSNEData, ellipse=TRUE, geom="point",
			legend="bottom")
}

#' To check one iteration of clustering before running full workflow CONCLUS.
#' 
#' This function generates a single clustering iteration of CONCLUS to check whether
#' chosen parameters for dbscan are suitable for your data.
#'
#' @param sceObject a SingleCellExperiment object with your experiment.
#' @param dataDirectory output directory (supposed to be the same for one experiment during the workflow).
#' @param experimentName name of the experiment which will appear in filenames (supposed to be the same for one experiment during the workflow).
#' @param dbscanEpsilon a parameter of fpc::dbscan() function.
#' @param minPts a parameter of fpc::dbscan() function.
#' @param PCs a vector of PCs for plotting.
#' @param perplexities vector of perplexities (t-SNE parameter).
#' @param randomSeed random seed for reproducibility.
#' @param width plot width.
#' @param height plot height.
#' @param ... other pdf() arguments.
#'
#' @return t-SNE results, a distance graph plot, a t-SNE plot colored by test clustering solution.
#' @export

testClustering <- function(sceObject, dataDirectory, experimentName,
		dbscanEpsilon=1.4,
		minPts=5,
		perplexities = c(30), PCs = c(4),
		randomSeed = 42,
		width=7, height=7, onefile=FALSE, ...){
	
	initialisePath(dataDirectory)
	dir.create(file.path(dataDirectory, "test_clustering"), showWarnings = F)
	
	message("Generating TSNE.")
	#1. Generating 2D tSNE plots
	tSNE <- .getTSNEresults(expr = Biobase::exprs(sceObject), cores=1,
			perplexities = perplexities,
			PCs = PCs,
			randomSeed = randomSeed)
	
	message("Saving results.")
	#picture checker
	pdf(file.path(dataDirectory, "test_clustering", "test_tSNE.pdf"),
			width=width, height=height, onefile=onefile, ...)
	print(tSNE)
	dev.off()
	
	
	#2. Clustering with dbscan
	# choosing for the best epsilon
	pdf(file.path(dataDirectory, "test_clustering", "distance_graph.pdf"),
			width=width, height=height, onefile=onefile, ...)
	.plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
			minNeighbours = minPts)
	dev.off()
	
	pdf(file.path(dataDirectory, "test_clustering", "test_clustering.pdf"),
			width=width, height=height, onefile=onefile, ...)
	print(.plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
					minNeighbours = minPts))
	dev.off()
	
	message("Pictures of test clustering were exported.")
	return(list(tSNE,
					.plotDistanceGraphWithEpsilon(tSNE$data, epsilon=dbscanEpsilon,
							minNeighbours = minPts),
					.plotTestClustering(tSNE$data, epsilon=dbscanEpsilon,
							minNeighbours = minPts)))
	
}
