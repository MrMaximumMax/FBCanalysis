#Functions for Cluster validation measure analysis

#' Initialize Cluster Validation Measure Analysis in the context of Fluctuation Based Clustering (FBC) analysis
#'
#' @return Object of type list storing cluster method(s) and number of cluster range of interest (to be used for function: \link{clValid_flow})
#'
#' @examples
#' init_clValid()
#'
#' @export
init_clValid <- function() {

  #List the methods to present them to the user
  methods <- as.data.frame(c("hierarchical", "kmeans", "diana", "fanny", "som", "modelbased", "sota",
                             "pam", "clara", "agnes"))
  #Corresponding numbers to methods so that user can choose easily which method(s) to apply
  num <- as.data.frame(c(1:10))
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

#' Interactive console workflow to calculate and evaluate cluster validation measures
#'
#' @param matrix Earth Mover's Distance Matrix for processed patient time series data (also see functions: \link{emd_matrix}, \link{patient_list})
#' @param par Object of type list storing clustering methods and cluster range of interest; initialized via function: \link{init_clValid}
#'
#' @return Object of type list storing chosen clustering method and number of clusters (can be then used for function \link{clust_matrix})
#'
#' @import clValid
#' @import RankAggreg
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

