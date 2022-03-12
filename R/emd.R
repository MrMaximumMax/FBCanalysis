#Generate, process and visualize Earth Mover's Distance data
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
