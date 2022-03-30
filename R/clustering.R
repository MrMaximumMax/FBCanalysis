#Function to cluster Earth Mover's Distance Matrix

#' Cluster Earth Mover's Distance Matrix
#'
#' Cluster Earth Mover's Distance Square Matrix data and record cluster assignments
#' for involved Patient_IDs for a specified clustering technique and number of
#' clusters.
#'
#' @param matrix Object of type matrix storing Earth Mover's Distances for patient time series data distribution pairs
#' @param method Clustering method (hierarchical, kmeans, diana, fanny, som, modelbased, sota, pam, clara)
#' @param nclust Number of clusters (if not specified, user will be asked in the terminal)
#' @param plotclust TRUE/FALSE if clustering data should be visualized (TRUE by default)
#'
#' @return Object of type list storing cluster data and clustering assignments for the Patient_IDs from the Earth Mover's Distance matrix
#'
#' @import stats
#' @import FCPS
#' @import cluster
#'
#' @details Hierarchical clustering describes a general agglomerative hierarchical
#' clustering approach in which the optimum value of an objective function is used
#' to choose which pair of clusters should merge at each step (see \link{hclust} for
#' more details).
#'
#' K-means clustering is a vector quantization technique that divides
#' a set of n observations into k groups, with each observation belonging to the cluster
#' with the closest mean or centroid (see \link{kmeans} for more details).
#'
#' The divisive analytic method (DIANA) constructs a hierarchical clustering structure,
#' starting with a single huge cluster containing all n observations. Clusters are further
#' split until each has a single observation. At each step, the cluster with the
#' largest diameter is selected, where the diameter of a cluster is defined as the
#' biggest difference between any two of its observations (see \link{diana} for
#' more details).
#'
#' While partitioning around medoids (PAM) is comparable to k-means, it is considered
#' more robust since it allows for the use of dissimilarities other than euclidean
#' distance. As with k-means, the number of clusters is determined in advance, and
#' an initial set of cluster centers is required to begin the process (see \link{pam}
#' for more details).
#'
#' Clustering large applications (CLARA) is a system that involves sampling to apply
#' PAM to a sequence of sub-datasets. When the number of observations is big, this
#' results in shorter run times. It is substantially faster than other partitioning
#' algorithms such as PAM at handling huge datasets. Internally, this is performed by
#' examining fixed-size sub-datasets, which results in linear rather than quadratic
#' time and storage requirements for n (see \link{clara} for more details).
#'
#' FANNY explains a fuzzy analysis clustering method. Each observation is distributed
#' through-out the numerous groups in a fuzzy clustering (see \link{fanny} for more
#' details).
#'
#' Self-organizing maps (SOM) are a widespread unsupervised learning technique used
#' by com- putational biologists and academics in machine learning. SOM is a neural
#' network-based system that is well-known for its ability to map and display
#' two-dimensional data (see \link{SOMclustering} for more details).
#'
#' Modelbased clusterinng fits the data to a statistical model composed of a finite
#' mixture of Gaussian distributions. Each mixture component represents a cluster,
#' and the maximum like- lihood method, or estimation maximum (EM), is used to
#' determine the mixture components and group memberships (see \link{ModelBasedClustering}
#' for more details).
#'
#' SOTA, self-organizing tree algorithm, denotes an unsupervised network with a
#' hierarchical and divisive binary tree topology. It is a fast approach, which
#' makes it suitable for clustering a large number of elements. It combines the
#' advantages of hierarchical clustering with those of SOMs. The algorithm chooses
#' the most diverse node and separates it into two nodes referred as cells (see
#' \link{sota} for more details).
#'
#' @references Joe H Ward Jr. Hierarchical grouping to optimize an objective function.
#' Journal of the American statistical association, 58(301):236–244, 1963.
#'
#' Stephane Tuffery. Data mining and statistics for decision making. John Wiley
#' & Sons, 2011.
#'
#' Fionn Murtagh. A survey of recent advances in hierarchical clustering algorithms.
#' The computer journal, 26(4):354–359, 1983.
#'
#' Fionn Murtagh. Clustering in massive data sets. In Handbook of massive data sets,
#' pages 501–543. Springer, 2002.
#'
#' Greg Hamerly and Charles Elkan. Alternatives to the k-means algorithm that find
#' better clusterings. In Proceedings of the eleventh international conference on
#' Information and knowledge management, pages 600–607, 2002.
#'
#' R Wehrens and J Kruisselbrink. kohonen: Supervised and unsupervised self-organising
#' maps r package version 3.0. 10, 2019.
#'
#' Leonard Kaufman and Peter J Rousseeuw. Finding groups in data: an introduction
#' to cluster analysis. John Wiley & Sons, 2009.
#'
#' Mark Van der Laan, Katherine Pollard, and Jennifer Bryan. A new partitioning
#' around medoids algorithm. Journal of Statistical Computation and Simulation,
#' 73(8):575–584, 2003.
#'
#' Chris Fraley and Adrian E Raftery. Model-based clustering, discriminant analysis,
#' and density estimation. Journal of the American statistical Association,
#' 97(458):611–631, 2002.
#'
#' Javier Herrero, Alfonso Valencia, and Joaquın Dopazo. A hierarchical unsupervised
#' growing neural network for clustering gene expression patterns.
#' Bioinformatics, 17(2):126–136, 2001.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' matrix <- emd_matrix(list, "FEV1")
#' clustering <- clust_matrix(matrix, method = "hierarchical", nclust = 2)
#'
#' @export
clust_matrix <- function(matrix, method, nclust, plotclust) {

  #Create a list to store output of clustering in a standardized format, reusable
  #for further functionalities, regardless with clustering method was applied
  newlist <- list()
  #In case user has not specified number of cluster
  if (missing(nclust)) {
    #Ask user how many cluster should be used and m
    nclust <- as.numeric(readline("How many clusters should be defined?: "))
  }

  if (method == "hierarchical") {
    #Cluster EMD data hierarchically according to Ward's Minimum Variance method
    hc <- hclust(as.dist(matrix),method = "ward.D2")

    #In case user has not specified "plotclust" or "plotclust = TRUE"
    if (missing(plotclust) || plotclust == TRUE) {
      #Design for node style in dendrogram
      nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
      #Plot dendrogram for applied hierarchical clustering
      plot(as.dendrogram(hc), type = "rectangle", ylab = "Height", nodePar = nodePar,
           leaflab = "none", main = "Dendrogram, Ward's Min. Var. Method")
      #Outline number of clusters with visualized rectangles
      rect.hclust(hc, k = nclust)
    }
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- cutree(hc , k = nclust)
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Hierarchical (Ward's Minimum Variance)"
    newlist

  } else if (method == "kmeans") {
    #Cluster EMD data according to kmeans algorithm
    clust <- kmeans(matrix, nclust)

    #In case user has not specified "plotclust" or "plotclust = TRUE"
    if (missing(plotclust) || plotclust == TRUE) {
      #Visualize two-dimensionally result of kmeans clustering
      clusplot(matrix, clust$cluster, main = "kmeans Clustering")
    }

    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$cluster
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Kmeans"
    newlist

  } else if (method == "diana") {

    #Check what user has specified on plotclust
    if (missing(plotclust) || plotclust == TRUE) {
      condition <- TRUE
    } else {
      condition <- FALSE
    }

    #Apply diana clustering on EMD data
    clust <- DivisiveAnalysisClustering(matrix, ClusterNo = nclust, PlotTree = condition)
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$Cls
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Diana"
    newlist

  } else if (method == "fanny") {

    #Apply fanny clustering on EMD data
    clust <- FannyClustering(matrix, ClusterNo = nclust)

    #In case user has not specified "plotclust" or "plotclust = TRUE"
    if (missing(plotclust) || plotclust == TRUE) {
      #Visualize two-dimensionally result of fanny clustering
      clusplot(fanny(matrix, nclust), main = "fanny clustering")
    }
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$Cls
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Fanny"
    newlist

  } else if (method == "som") {

    #Apply Som clustering on EMD data
    #Currently no cluster visualization available
    clust <- SOMclustering(matrix, ClusterNo = nclust, PlotIt = FALSE)
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$Cls
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Som"
    newlist

  } else if (method == "modelbased") {

    #Apply Som clustering on EMD data
    #Currently no cluster visualization available
    clust <- ModelBasedClustering(matrix, ClusterNo = nclust, PlotIt = FALSE)
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$Cls
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Modelbased"
    newlist

  } else if (method == "sota") {

    #Apply Sota clustering on EMD data
    #Currently no cluster visualization available
    clust <- SOTAclustering(matrix, ClusterNo = nclust, PlotIt = FALSE)
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$Cls
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Sota"
    newlist

  } else if (method == "pam") {

    #Apply pam clustering on EMD data
    clust <- PAMclustering(matrix, ClusterNo = nclust)
    #In case user has not specified "plotclust" or "plotclust = TRUE"
    if (missing(plotclust) || plotclust == TRUE) {
      #Visualize two-dimensionally result of pam clustering
      clusplot(pam(matrix, nclust), main = "pam clustering")
    }
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$Cls
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Pam"
    newlist

  } else if (method == "clara") {

    #Apply clara clustering on EMD data
    clust <- cluster::clara(matrix, nclust)
    #In case user has not specified "plotclust" or "plotclust = TRUE"
    if (missing(plotclust) || plotclust == TRUE) {
      #Visualize two-dimensionally result of clara clustering
      clusplot(cluster::clara(matrix, nclust), main = "clara clustering")
    }
    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$cluster
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Clara"
    newlist

  } else {
    stop("Incorrect input of clustering method")
  }
}
