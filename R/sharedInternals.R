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
### this function calculates PCA and then UMAP with PCs and nearest neighbors ###
### it returns a list of pUMAP = PCA+UMAP results ###
### to get XY coordinates call pumap_res[1,i][[1]] ###
### i = [1:length(PCs)*length(neighbors)] is a number of iteration ###

.getUMAPresults <- function(expressionMatrix, cores=1,
                            PCs=c(4, 6, 8, 10, 20, 40, 50),
                            neighbors=c(10, 15),Metric="euclidean",
                            Spread=1,Min_dist=0.15){
  
  PCAData <- prcomp(t(expressionMatrix))$x
  myCluster <- parallel::makeCluster(cores, # number of cores to use
                                     type = "PSOCK") # type of cluster
  doParallel::registerDoParallel(myCluster)
  UMAPcoordinates <- foreach::foreach(PCA=rep(PCs, length(neighbors)),
                                      perp=rep(neighbors, each=length(PCs)),
                                      .combine='cbind') %dopar% {
                                        library(ggplot2)
                                        tmp <- uwot::umap(X=(PCAData[,1:PCA]),
                                          n_neighbors=perp,metric = Metric,
                                          spread = Spread, min_dist = Min_dist)
                                        tmp <- data.frame(tmp)
                                        colnames(tmp) <- c("UMAP_1","UMAP_2")
                                        rownames(tmp) <- rownames(PCAData)
                                        tmp <- ggplot(data = tmp,aes(x=UMAP_1,y=UMAP_2))
                                        return(tmp)
                                      
                                      }
  parallel::stopCluster(myCluster)
  message(paste("Calculated", length(PCs)*length(neighbors),
                "UMAP plots."))
  return(UMAPcoordinates)
}

choosePalette <- function(colorPalette, clustersNumber){
	
	colorPalette26 <- c( "yellow", "darkgoldenrod1", "coral1", "deeppink",
			"indianred", "coral4", "darkblue", "darkmagenta",
			"darkcyan", "mediumorchid", "plum2", "gray73",
			"cadetblue3", "khaki",
			"darkolivegreen1", "lightgreen", "limegreen",
			"darkolivegreen4", "green", "#CC79A7", "violetred3",
			"brown3", "darkslategray1", "gray51", "slateblue2",
			"blue")
	
	pickDefaultPalette <- function(clustersNumber, colorPalette26){
		if(clustersNumber < 13){
			return(c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C",
							"#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
							"#CAB2D6", "#6A3D9A", "#FFFF99",
							"#B15928")[1:clustersNumber])
		}else{
			return(colorPalette26[seq_len(length(clustersNumber))])
		} 
	}
	
	if(colorPalette == "default" && clustersNumber > 26)
			stop("The default option is limited to 26 colors, please provide your own color vector.")
	
	if(colorPalette != "default" && clustersNumber > length(colorPalette)){
		stop("The number of clusters is greater than the number of given colors.")
	}
	
	if(colorPalette == "default")
			return(pickDefaultPalette(clustersNumber))
	
	return(colorPalette)
}


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
