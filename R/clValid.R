#Functions for Cluster validation measure analysis

#' Initialize Cluster Validation Measure Analysis
#'
#' Initialize Cluster Validation Measure Analysis in the context of Fluctuation
#' Based Clustering (FBC) analysis. The call initialized a interactive console
#' workflow where the user may indicate the clustering techniques of interest
#' as well as the cluster numbers of interest.
#'
#' @return Object of type list storing cluster method(s) and number of cluster range of interest (to be used for function: \link{clValid_flow})
#'
#' @details See \link{clust_matrix} for more details on possible clustering techniques
#'
#' @references Guy Brock, Vasyl Pihur, Susmita Datta, and Somnath Datta. clvalid:
#' An r package for cluster validation. Journal of Statistical Software, 25:1–22, 2008.
#'
#' @examples
#' init_clValid()
#'
#' @export
init_clValid <- function() {

  #List the methods to present them to the user
  methods <- as.data.frame(c("hierarchical", "kmeans", "diana", "fanny", "som", "modelbased", "sota",
                             "pam", "clara"))
  #Corresponding numbers to methods so that user can choose easily which method(s) to apply
  num <- as.data.frame(c(1:9))
  #Combine methods and corresponding numbers
  methods <- cbind(num, methods)
  colnames(methods) <- c("Number", "Clustering method:")
  #Present methods to user
  cat("The following methods are available for analysis", "\n\n")
  print(methods, row.names = FALSE)
  cat("\n\n")
  #Ask user which method(s) to apply
  numbers <- as.character(readline("Indicate which methods to use (e.g. 1,2,8): "))
  numbers <- strsplit(numbers,",")
  numbers <- numbers[[1]]
  numbers <- as.numeric(numbers)
  #New list to store the parameter for clValid analysis (to recall them easily later)
  newList <- list()
  #Filter for methods by indicated number
  new <- subset(methods, Number %in% numbers)
  new <- as.character(new[,2])
  #Add the chosen clustering methods to list
  newList[["cl_methods"]] <- new
  #Ask user which cluster number to apply
  cat("\n\n")
  cl_num <- readline("Which cluster numbers do you want to include (e.g. 2:5)?: ")
  #Split by ":"
  cl_num <- strsplit(cl_num,":")[1]
  #Extract min number of clusters
  min_cl_num <- cl_num[[1]][1]
  #Extract max number of clusters
  max_cl_num <- cl_num[[1]][2]
  #Store min and max number of clusters in list
  newList[["min_cl_num"]] <- as.numeric(min_cl_num)
  newList[["max_cl_num"]] <- as.numeric(max_cl_num)
  newList
}

#' Cluster validation measure analysis workflow
#'
#' Interactive console workflow to calculate and evaluate cluster validation measures
#' which have been determined previously by the call \link{init_clValid}.
#'
#' @param matrix Earth Mover's Distance Matrix for processed patient time series data (also see functions: \link{emd_matrix}, \link{patient_list})
#' @param par Object of type list storing clustering methods and cluster range of interest; initialized via function: \link{init_clValid}
#'
#' @return Object of type list storing chosen clustering method and number of clusters (can be then used for function \link{clust_matrix})
#'
#' @import clValid
#' @import RankAggreg
#'
#' @details The call guides through an interactive workflow and generates cluster
#' evaluation measures, stores and lists, visualizes corresponding plots and lets
#' the user decide which technique is the prefered one. Once the user has chosen
#' his favourite, the flow continues to the function \link{clust_matrix} and generates
#' the respective clustering output. The internal cluster validation methods utilize
#' just the dataset and the clustering partition as input and evaluates the
#' clustering’s quality by using intrinsic information included in the data.
#'
#' The call calculates Connectivity, Silhouette width and Dunn index. Connectivity
#' describes the connectness to neighbors of particular clustering partition and
#' should be minimized. Silhouette width defines the average silhouette value for
#' each observation and should be maximized. The Dunnn index is a definition for
#' Ratio of shortest distance between non-cluster observation and greatest intra-cluster
#' distance and should be maximized likewise.
#'
#' Furthermore, cluster stability measures are available, namely Average proportion
#' of non-overlap (APN), Average distance (AD), Average distance between means (ADM)
#' and Figure of merit (FOM). APN is the average proportion of observations that are
#' not clustered using complete and leaky data. AD defines the average distance in
#' observations for both complete and leaky data. ADM deals with the average distance
#' between cluster centers in complete and leaky data. FOM is a measure for average intra-
#' cluster variance in leaky data. All measures should be minimized. Furthermore,
#' Rank Aggregation may be performed. It approaches to provide a generic and flexible
#' framework for objectively integrating several ordered lists in a suitable and
#' efficient way. The used technique for evaluating clustering the rank is by a
#' Cross-entropy approach, which is incorporating Spearman’s footrule distance measure.
#' In the end a recommendation for the best fitting clustering model is given.
#'
#' @references Guy Brock, Vasyl Pihur, Susmita Datta, and Somnath Datta. clvalid:
#' An r package for cluster validation. Journal of Statistical Software, 25:1–22, 2008.
#'
#' Julia Handl, Joshua Knowles, and Douglas B Kell. Computational cluster validation
#' in post-genomic data analysis. Bioinformatics, 21(15):3201–3212, 2005.
#'
#' Peter J Rousseeuw. Silhouettes: a graphical aid to the interpretation and
#' validation of cluster analysis. Journal of computational and applied mathematics,
#' 20:53–65, 1987.
#'
#' Joseph C Dunn. Well-separated clusters and optimal fuzzy partitions. Journal of
#' cybernetics, 4(1):95–104, 1974.
#'
#' Vasyl Pihur, Susmita Datta, and Somnath Datta. Weighted rank aggregation of
#' cluster validation measures: a monte carlo cross-entropy approach. Bioinformatics,
#' 23(13):1607–1615, 2007.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' distmat <- emd_matrix(list, "PEF", maxIter = 5000)
#' parameters <- init_clValid()
#' output <- clValid_flow(distmat, parameters)
#'
#' @export
clValid_flow <- function(matrix, par) {

  #Output list to store the finally chosen clustering combination
  output_list = list()
  #First question for user (internal cluster validation)
  question <- as.character(readline("Do you want to perform internal cluster validation (y/n)?: "))
  #If yes
  if (question == "y") {
    #clValid's function to receive measures on internal cluster validation on
    #chosen clustering combinations
    int <- clValid::clValid(matrix, as.numeric(par$min_cl_num):as.numeric(par$max_cl_num), clMethods = par$cl_methods,
                            validation = "internal")
    #Print summary to console so that user receives result on measures
    clValid::summary(int)
    #Show plots in development on stability measures with varying clustering
    #methods and number of clusters
    clValid::plot(int)
  }
  #Second question for user (cluster stability validation)
  question <- as.character(readline("Do you want to perform cluster stability validation (y/n)?: "))
  #If yes
  if (question == "y") {
    #clValid's function to receive measures on cluster stability validation on
    #chosen clustering combinations
    stab <- clValid::clValid(matrix, par$min_cl_num:par$max_cl_num, clMethods = par$cl_methods,
                             validation = "stability")
    #Print summary to console so that user receives result on measures
    clValid::summary(stab)
    #Show plots in development on stability measures with varying clustering
    #methods and number of clusters
    clValid::plot(stab)
  }
  #Third question for user (Rank Aggregation)
  question_r <- as.character(readline("Do you want to perform rank aggregation? (y/n): "))
  #If yes
  if (question_r == "y") {
    #Store both results from internal and cluster stability validation in list
    result <- clValid::clValid(matrix, par$min_cl_num:par$max_cl_num, clMethods = par$cl_methods,
                               validation = c("internal","stability"))
    #Calculate rank weights on super list
    res <- clValid::getRanksWeights(result)
    cat("\n\n", "Ranks of cluster measures: ", "\n\n")
    print(res$ranks[,1:3])
    cat("\n")
    #Fourth question (Cross Entropy Search Algorithm)
    question_r2 <- as.character(readline("Let Cross-Entropy search for ideal method? (y/n): "))
    #If yes
    if (question_r2 == "y") {
      #Apply RankAggreg for Cross Entropy search on super list with Spearman's footrule
      if(require("RankAggreg")) {
        CEWS <- RankAggreg(x=res$ranks, k=5, weights=res$weights, seed=123, verbose=FALSE)
        cat("\n\n")
        #Print ideal clustering combination from super list
        print(CEWS)
        cat("\n\n")
      }
    }
  }
    #Empty line
    cat("\n")
    #Present in a matrix the previously analysed clustering methods
    present = matrix(0, nrow = 2, ncol = length(par$cl_methods))
    #Assign numbers to the methods so that the user can choose method by number
    present[1,] <- c(1:length(par$cl_methods))
    present[2,] <- par$cl_methods
    #Print clustering methods to user
    print(present)
    cat("\n\n")
    #Choose number for clustering combination
    new_question <- as.integer(readline("Which of your analyzed methods should be used (please indicate number)?: "))
    #Add clustering method to output list
    output_list[["method"]] <- par$cl_methods[new_question]
    #Choose number of cluster
    next_question <- as.integer(readline("How many clusters should be applied?: "))
    #Add number of clusters to output list
    output_list[["clust_num"]] <- as.numeric(next_question)
    final_output <- clust_matrix(matrix, method = output_list$method, nclust = output_list$clust_num)
}

