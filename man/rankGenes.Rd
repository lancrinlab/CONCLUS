\name{rankGenes}
\alias{rankGenes}
\title{Rank marker genes by statistical significance.}
\description{
This function searches marker genes for each cluster. It saves tables in the "dataDirectory/marker_genes" directory,
one table per cluster.
}
\usage{
rankGenes(sceObject, clustersSimilarityMatrix, dataDirectory,
  experimentName, column = "clusters")
}
\arguments{
\item{sceObject}{A SingleCellExperiment object with your experiment.}

\item{clustersSimilarityMatrix}{Result of conclus::calculateClustersSimilarity() function.}

\item{dataDirectory}{Output directory for CONCLUS results}

\item{experimentName}{Prefix used for output files.}

\item{column}{Name of the column with a clustering result.}
}
\details{

  To understand the nature of the consensus clusters identified by CONCLUS, it is essential to identify genes which could be classified as marker genes for each cluster. To this aim, each gene should be "associated" to a particular 
cluster. This association is performed by looking at upregulated genes in a particular cluster compared to the others (multiple comparisons). The function *rankGenes* performs multiple comparisons of all genes from the 
sceObjectFiltered and rank them according to a score reflecting a FDR power. It saves one file per cluster in the folder *marker_genes*.

For each table corresponding to a particular consensus cluster,  the first column is a gene name. The following columns represent adjusted p-values (FDR) of a one-tailed T-test between the considered cluster and all others. 
For example, If one considers cluster 5 (“Bergiers_cluster_5_genes.tsv”), the null hypothesis is that each gene, on each row, are not upregulated.  The FDR of the hypothesis that a gene was not upregulated in the cluster 5 versus 
cluster 1, is stored in the column “vs_1” (FDR < 0.05 enables to reject the null hypothesis. In other words, if FDR < 0.05, the gene is upregulated in cluster 5 vs cluster 1). 

Top genes with significant FDR in most of the comparisons can be assumed as positive markers of a cluster. The column mean_log10_fdr is the mean power of FDR in all comparisons; the column n_05 is the number of comparisons in which the 
gene was significantly upregulated. The score for marker genes is the average power of FDR among all comparisons for a cluster multiplied to weights taken from the clustersSimilarityMatrix + 0.05. Taking into account both FDRs of all 
comparisons and clustersSimilarityMatrix allows us to keep the balance between highlighting markers for individual clusters and their “families” which makes the final heatmap as informative as possible. 

Note: Adding 0.05 to the clustersSimilarityMatrix in calculating the score helps avoiding the following problem: in case you have a cluster very different from all others, it will have the value 1 on the diagonal and 0 similarities to 
all others groups in the clustersSimilarityMatrix. So all weights for that cluster will be zeros meaning that the score would also be zero and genes will be ordered in alphabetical order in the corresponding marker genes list file.

For a cluster k and a gene G, a scoreG was defined in the following way:
 
scoreG= ∑i(−log10(fdrk,i + ϵ)∗weightk,i) / nClusters−1, 

Where 

1. fdrk,i is an adjusted p-value obtained by comparing expression of G in cluster k versus expression of G in cluster i 
2. weightk,i is a similarity between these two groups taken from the  element in the clustersSimilarityMatrix 
3. nClusters is a number of consensus clusters given to the rankGenes() 
4. ϵ = 10−300 is a small number which does not influence the ranking and added to avoid an error when fdr is equal to zero 
5. k = [1,…,nClusters] 
6. I = ([1,…,nClusters]exceptfor[k])
}
