#Generate, process and visualize Earth Mover's Distance data

#' Generate a Earth Mover's Distance Matrix out of time series data list entries
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list(path)})
#' @param parameter Parameter of interest to determine Earth Mover's Distances between distributions
#' @param maxIter Maximum of iterations to calculate Earth Mover's Distance (default: 500)
#'
#' @return Earth Mover's Distance Square Matrix of type matrix
#'
#' @import emdist
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be pulled from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo_and_testfiles/ts_demofiles1)
#' #Sampling frequency is supposed to be daily
#' matrix <- emd_matrix(list, "FEV1")
#'
#' @export
emd_matrix <- function (plist, parameter, maxIter) {

  #Determine number of patients in list
  N <- length(plist)
  #Make a NxN matrix, filled with 0's, as basis for distance matrix
  #For loops can go to each position then and calculate Earth Mover's
  #distance for each patient-parameter-distribution-pair
  distmat <- matrix(0, nrow = N, ncol = N);
  #Initializes progress bar as feedback on progress of whole calculation procedure
  pb <- txtProgressBar(min = 0, max = N, style = 3)

  #Start first for-loop (from 1 to number of patients)
  for (i in c(1:N)) {
    #Update progress bar at each new i
    setTxtProgressBar(pb, i)
    #Take data from first patient out of list
    distr_a <- plist[[i]]
    #Filter patient data for specified parameter and normalize data
    norm_a <- as.matrix(as.data.frame(table((distr_a[,parameter] - min(distr_a[,parameter]))/(max(distr_a[,parameter]) - min(distr_a[,parameter])))))
    #Start second for-loop (from 1 to number of patients)
    for (j in c(1:N)) {
      #Check that i and j are not equal to remain diagonal 0's on distance matrix
      if (j!=i) {
        #Take data from second patient out of list
        distr_b <- plist[[j]]
        #Filter patient data for specified parameter and normalize data
        norm_b <- as.matrix(as.data.frame(table((distr_b[,parameter] - min(distr_b[,parameter]))/(max(distr_b[,parameter]) - min(distr_b[,parameter])))))
        #In case max. Iterations not specified, default number of iterations on
        #calculating EMD is set to 500
        #In some cases not sufficient (leads to a warning message)
        if (missing(maxIter)) {
          #Calculate EMD on patient data distribution pair
          distmat[i,j] <- emd(norm_a, norm_b, max.iter = 500);
        } else {
          #Calculate EMD on patient data distribution pair with specified max. iterations
          distmat[i,j] <- emd(norm_a, norm_b, max.iter = maxIter)
        }
      }
    }
  }
  #Rename columns and rows from matrix with Patient_ID's to recognize pairs easily
  colnames(distmat) <- names(plist)
  rownames(distmat) <- names(plist)
  distmat
}

#' Visualize an Earth Mover's Distance Square Matrix as a heatmap
#'
#' @param input Earth Mover's Distance Matrix or list storing patient time series data (also see function: \link{patient_list})
#' @param parameter In case list is input, the parameter of interest from time series data list
#'
#' @return Visualized Earth Mover's Distance Matrix as a heatmap
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be pulled from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo_and_testfiles/ts_demofiles1)
#' #Sampling frequency is supposed to be daily
#' matrix <- emd_matrix(list, "FEV1")
#' emd_heatmap(matrix)
#'
#' @export
emd_heatmap <- function(input, parameter) {

  #In case input is either of type "matrix" or "double" and no parameter is specified
  if (missing(parameter) & typeof(input) == "matrix" || typeof(input) == "double") {
    #Apply heatmap on input, disable any further functionalities (e.g dendrogram),
    #remain order of matrix data entries and scale font size of legend
    heatmap(input, Colv = NA, Rowv = NA, scale = "column", cexRow = 0.7, cexCol = 0.7)

    #In case input is list containing patient data
  } else if (typeof(input) == "list") {
    #Calculate EMD matrix on data from for specified parameter and then apply heatmap
    #on EMD matrix, disable any further functionalities (e.g. dendrogram), remain
    #order of matrix data entries and scale font size of legend
    heatmap((emd_matrix(input, parameter)), Colv = NA, Rowv = NA, scale = "column",
            cexRow = 0.7, cexCol = 0.7)
    #In case type of input is neither "matrix", "double" nor "list", stop and give
    #feedback to user
  } else {
    stop("Incorrect input. Either indicater a distance matrix or list of time series dara")
  }
}

#' Determine pair of maximum fluctuation difference in a list storing time series data
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest from time series data list
#'
#' @return Console output with Patient_ID pair, corresponding Earth Mover's Distance and visualized boxplot of both time series data distributions
#'
#' @import emdist
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be pulled from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo_and_testfiles/ts_demofiles1)
#' #Sampling frequency is supposed to be daily
#' matrix <- emd_matrix(list, "FEV1")
#' max_fluc(list, "PEF")
#'
#' @export
max_fluc <- function(plist, parameter) {

  #Calculate EMD matrix out of specified list and parameter
  distmat <- emd_matrix(plist, parameter)
  #Determine matrix position where highest EMD was find in distmat
  max_pair <- as.vector(which(distmat==max(distmat), arr.ind = TRUE))
  #Take colname for first and second index from max_pair (= Patient-pair)
  pat1 <- colnames(distmat)[max_pair[1]]
  pat2 <- colnames(distmat)[max_pair[2]]
  #Take out the value for the determined patient-pair
  value <- distmat[max_pair[1],max_pair[2]]
  #Output for user
  cat("\n\n","Maximum fluctuation difference between: ", pat1, " and ", pat2, "\n",
      "Earth Mover's Distance = ", value)
  #Store patient names in vector
  patients <- c(pat1, pat2)
  #Visualize boxplot for patients and specified parameter
  patient_boxplot(plist, patients, parameter)
}
