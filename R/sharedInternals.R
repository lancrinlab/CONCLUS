### this function calculates PCA and then tSNE with PCs and perplexities ###
### it returns a list of pSNE = PCA+tSNE results ###
### to get XY coordinates call psne_res[1,i][[1]] ###
### i = [1:length(PCs)*length(perplexities)] is a number of iteration ###

.getTSNEresults <- function(expressionMatrix, cores=1,
		PCs=c(4, 6, 8, 10, 20, 40, 50),
		perplexities=c(30, 40), randomSeed=42){
	PCAData <- prcomp(t(expressionMatrix))$x
	myCluster <- parallel::makeCluster(cores, # number of cores to use
			type = "PSOCK") # type of cluster
	doParallel::registerDoParallel(myCluster)
	tSNECoordinates <- foreach::foreach(PCA=rep(PCs, length(perplexities)),
					perp=rep(perplexities, each=length(PCs)),
					.combine='cbind') %dopar% {
				library(SingleCellExperiment)
				tmp <- scater::runTSNE(
						SingleCellExperiment(assays=list(
										logcounts=t(PCAData[,1:PCA]))),
						scale_features=FALSE, perplexity=perp,
						rand_seed=randomSeed, theme_size=13, 
						return_SCESet=FALSE)
				
				scater::plotTSNE(tmp)
			}
	parallel::stopCluster(myCluster)
	message(paste("Calculated", length(PCs)*length(perplexities),
					"2D-tSNE plots."))
	return(tSNECoordinates)
}



#' Depending if the number of clusters is more than 12 or not, one of two built-in palettes will be applied.
#' If you give your vector of colors, the function will not change them.
#' If the number of clusters is more than 26, it will copy colors to get the needed length of the palette.
#'
#' @param colorPalette Either "default" or a vector of colors, for example c("yellow", "#CC79A7"). 
#' @param clustersNumber number of clusters in the output palette.
#'
#' @return Color palette with the number of colors equal to the clusterNumber parameter.
#' @export
.choosePalette <- function(colorPalette, clustersNumber){
	
	colorPalette26 <- c( "yellow", "darkgoldenrod1", "coral1", "deeppink",
			"indianred", "coral4", "darkblue", "darkmagenta",
			"darkcyan", "mediumorchid", "plum2", "gray73",
			"cadetblue3", "khaki",
			"darkolivegreen1", "lightgreen", "limegreen",
			"darkolivegreen4", "green", "#CC79A7", "violetred3",
			"brown3", "darkslategray1", "gray51", "slateblue2",
			"blue")
	
	pickDefaultPalette <- function(clustersNumber, colorPalette26){
		if(clustersNumber < 13) return(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
							"#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
							"#CAB2D6", "#6A3D9A", "#FFFF99",
							"#B15928")[1:clustersNumber])
		return(rep(colorPalette26,
						round(clustersNumber/length(colorPalette26))+1)[1:clustersNumber])
	}
	
	if(length(colorPalette) == 1){
		if(colorPalette == "default"){
			return(pickDefaultPalette(clustersNumber, colorPalette26))
		}
	}
	
	if(clustersNumber > length(colorPalette)){
		message("The number of clusters is greater than the number of colors.
						Using default CONCLUS palette instead.")
		return(pickDefaultPalette(clustersNumber, colorPalette26))
	}
	
	return(colorPalette[1:clustersNumber])
}


#' Create all needed directories for CONCLUS output.
#'
#' @param dataDirectory output directory for a given CONCLUS run (supposed to be the same for one experiment during the workflow).
#' @export
initialisePath <- function(dataDirectory){
	# creates directories for further writing of results.
	# names of directories are hardcoded.
	# no idea if it is good or bad.
	
	graphsDirectory <- "pictures"
	markerGenesDirectory <- "marker_genes"
	tSNEDirectory <- "tsnes"
	outputDataDirectory <- "output_tables"
	tSNEPicturesDirectory <- "tSNE_pictures"
	
	
	dir.create(dataDirectory, showWarnings=F)
	dir.create(file.path(dataDirectory, graphsDirectory), showWarnings=F)
	dir.create(file.path(dataDirectory, graphsDirectory, tSNEPicturesDirectory),
			showWarnings=F)
	dir.create(file.path(dataDirectory, markerGenesDirectory), showWarnings=F)
	dir.create(file.path(dataDirectory, tSNEDirectory), showWarnings=F)
	dir.create(file.path(dataDirectory, outputDataDirectory), showWarnings=F)
	
}
