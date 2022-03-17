#Function to cluster Earth Mover's Distance Matrix

#' Cluster Earth Mover Distance Square Matrix data
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
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be pulled from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo_and_testfiles/ts_demofiles1)
#' #Sampling frequency is supposed to be daily
#' matrix <- emd_matrix(list, "FEV1")
#' clustering <- clust_matrix(matrix, method = "hierarchical", nclust = 2)
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
