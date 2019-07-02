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
