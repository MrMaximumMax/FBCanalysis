#Functions for Jaccard index determination/random data removal

#' Remove random data from time series data list
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param removal Amount of data removal (0 = 0%, 1 = 100%)
#'
#' @return Object of type list storing patient time series data with indicated amount of data removed randomly
#'
#' @import dplyr
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be retrieved from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo/phys)
#' #Sampling frequency is supposed to be daily
#' list_rm <- rnd_dat_rm(testlist, 0.95)
#'
#' @export
rnd_dat_rm <- function(plist, removal) {

  #Make one time series data frame out of all list elements
  df <- do.call(rbind, plist)
  #Determine number of time series data entries and number of entries to be removed
  n <- round(nrow(df)*removal, digits = 0)
  #For loop, for how many entries to remove
  for (i in 1:n) {
    #Select a random row from 1:n_row
    random_row <- sample(1:nrow(df), 1)
    #Remove the current random row
    df <- df[-random_row,]
  }
  #Make a new list
  newlist = list()
  #Determine the Patient_ID's from the time series data frame
  patnames <- names(table(df[,"Patient_ID"]))
  #For each Patient_ID
  for (i in 1:length(patnames)) {
    #Filter the time series data frame by current Patient_ID
    new_df <- dplyr::filter(df, Patient_ID %in% patnames[i])
    #Add the current new data frame to the list
    newlist[[patnames[i]]] <- new_df
  }
  newlist
}

#' Simulate random data removal from time series data list and determine Jaccard index via Cognate Cluster approach
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param removal Amount of random data removal to determine Jaccard index
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param n_clust Number of clusters (also see function: \link{clust_matrix})
#' @param Iter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix}); default is 5,000 for this function
#'
#' @return Object of type matrix storing received Jaccard indices for indicated amount of random data removal for all clusters
#'
#' @import tibble
#' @import dplyr
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be retrieved from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo/phys)
#' #Sampling frequency is supposed to be daily
#' output <- sim_jaccard_cognate(list, "PEF", 0.05, 10, "hierarchical", 2, 1000)
#'
#' @export
sim_jaccard_cognate <- function(plist, parameter, removal, n_simu, method, n_clust, Iter) {
  #Simulate random data removal and Jaccard index determination by Cognate Cluster Approach

  #In case, user did not specify maximum for EMD calculations, the default value is
  #set to a high number (5,000)
  if (missing(Iter)) {
    Iter <- 5000
  }
  #Make a gold standard EMD matrix for complete data without data removal
  distmat_complete <- emd_matrix(plist, parameter, maxIter = Iter)
  #Create a summary matrix to fill up with Jaccard indices for each simulation run
  #and each cluster
  summary <- matrix(0, nrow = n_simu, ncol = n_clust)
  #Cluster the gold standard EMD matrix
  dat_complete <- clust_matrix(distmat_complete, method, nclust = n_clust, plotclust = FALSE)
  #Take the cluster assignments by patient out of cluster result list and order them
  #by Patient_ID (later used to calculate Jaccard index)
  compare <- as.data.frame(dat_complete$Cls)
  compare <- tibble::rownames_to_column(compare,"Patient_ID")
  compare <- compare[order(compare[,"Patient_ID"]),]
  #For each simulation run
  for (i in 1:n_simu) {
    #Remove randomly data from time series data list according to specified removal amount
    list_rm <- rnd_dat_rm(plist, removal)
    #Calculate a new EMD matrix after random data removal
    distmat_rm <- emd_matrix(list_rm, parameter, maxIter = Iter)
    #Cluster random data removal EMD matrix according to specified parameters
    dat_rm <- clust_matrix(distmat_rm, method, nclust = n_clust, plotclust = FALSE)
    #Take the cluster assignments by patient out of cluster result list and order them
    #by Patient_ID
    compare2 <- as.data.frame(dat_rm$Cls)
    compare2 <- tibble::rownames_to_column(compare2,"Patient_ID")
    compare2 <- compare2[order(compare2[,"Patient_ID"]),]
    compare2 <- compare2[,2]
    #Combine both cluster assignment to compare assignments from gold standard
    #and random data removal
    combined <- cbind(compare, compare2)
    colnames(combined) <- c("Patient_ID","all_dat","data_removal")
    #Calculate for each gold standard cluster Jaccard index with every random
    #data removal cluster to determine cognate cluster
    for (m in 1:n_clust) {
      #filter for current gold standard cluster
      overlap1 <- dplyr::filter(combined, all_dat == m)
      #Make a vector to store all calculated Jaccard indices; Highest Jaccard
      #index then represents the cognate cluster
      jaccard_vector <- vector()
      #For-loop to filter for each data removal cluster
      for (n in 1:n_clust) {
        #Filter for data removal cluster
        overlap2 <- dplyr::filter(combined, data_removal == n)
        #Find intersection between current gold standard cluster and current
        #data removal cluster
        common <- nrow(intersect(overlap1, overlap2))
        #Quantify all; Both current gold standard and data removal cluster
        all <- nrow(union(overlap1, overlap2))
        #Calculate current Jaccard index
        jaccard <- common/all
        #Add Jaccard index to vector
        jaccard_vector[n] <- jaccard
      }
      #Add highest entry from Jaccard vector to summary-matrix (this represents
      #the cognate cluster value)
      summary[i,m] <- max(jaccard_vector)
    }
  }
  summary
}

#' Simulate random data removal from time series data list and determine Jaccard index via Earth Mover's Distance approach
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param removal Amount of random data removal to determine Jaccard index
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param n_clust Number of clusters (also see function: \link{clust_matrix})
#' @param Iter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix}); default is 5,000 for this function
#'
#' @return Object of type matrix storing received Jaccard indices for indicated amount of random data removal for all clusters
#'
#' @import tibble
#' @import dplyr
#' @import emdist
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be retrieved from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo/phys)
#' #Sampling frequency is supposed to be daily
#' output <- sim_jaccard_emd(list, "PEF", 0.05, 10, "hierarchical", 2, 1000)
#'
#' @export
sim_jaccard_emd <- function(plist, parameter, removal, n_simu, method, n_clust, Iter) {

  #Simulate random data removal and Jaccard index determination by EMD Approach
  #In case, user did not specify maximum for EMD calculations, the default value is
  #set to a high number (5,000)
  if (missing(Iter)) {
    Iter <-  5000
  }
  #Make a gold standard EMD matrix for complete data without data removal
  distmat_complete <- emd_matrix(plist, parameter, maxIter = Iter)
  #Cluster the gold standard EMD matrix
  dat_complete <- clust_matrix(distmat_complete, method, nclust = n_clust, plotclust = FALSE)
  #Take the cluster assignments by patient out of cluster result list and order them
  #by Patient_ID (later used to calculate Jaccard index)
  compare <- as.data.frame(dat_complete$Cls)
  compare <- tibble::rownames_to_column(compare,"Patient_ID")
  compare <- compare[order(compare[,"Patient_ID"]),]
  #Make a masterlist to store EMD's from each new cluster after random data removal
  #and gold standard clusters
  masterlist <- list()
  #Create a summary matrix to fill up with Jaccard indices for each simulation run
  #and each cluster
  summary <- matrix(nrow = n_simu, ncol = n_clust)
  #For each cluster
  for (t in 1:n_clust) {
    #Take the Patient_ID's from current cluster according to clustering data output
    patnames <- dplyr::filter(compare, dat_complete$Cls == t)[,1]
    #Make data frame out of time series data list (to be used later to extract distributions)
    distribution <- do.call(rbind,plist)
    #List to store normalized distributions of each Patient_ID for gold standard
    dat <- list()
    #For each Patient_ID in time series data list
    for (u in 1:length(patnames)) {
      #Filter for current Patient_ID to receive current distribution
      distr_a <- dplyr::filter(distribution, Patient_ID %in% patnames[u])
      #Normalize current distribution
      norm_a <- ((distr_a[,parameter] - min(distr_a[,parameter]))/(max(distr_a[,parameter]) - min(distr_a[,parameter])))
      #Add current distribution to list
      dat[[patnames[u]]] <- norm_a
    }
    #Add for each current cluster in masterlist the normalized time series values
    #of all cluster members
    masterlist[[t]] <- do.call(c,dat)
  }
  #Initialize a progress bar
  N <- length(n_simu)
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  #For each simulation run
  for (i in 1:n_simu) {
    #Update progress bar
    setTxtProgressBar(pb, i)
    #Randomly remove data from time series data list according to specified
    #removal amount
    list_rm <- rnd_dat_rm(plist, removal)
    #data frame with each Patient_ID to later store the EMDs to every gold
    #standard cluster
    assignments <- as.data.frame(names(list_rm))
    #For every Patient_ID in list after random data removal
    for (q in 1:length(list_rm)) {
      #New vector to store EMD for current Patient_ID to each gold standard cluster
      distances <- vector()
      #For every normalized distribution from gold standard clusters
      for(r in 1:length(masterlist)) {
        #Take current random data removal distribution and normalize
        distr_b <- list_rm[[q]]
        norm_b <- as.matrix(as.data.frame(table((distr_b[,parameter] - min(distr_b[,parameter]))/(max(distr_b[,parameter]) - min(distr_b[,parameter])))))
        #Calculate EMD between current random data removal distribution and
        #current gold standard distribution to distances vector
        distances[r] <- emd(norm_b, as.matrix(as.data.frame(table(masterlist[[r]]))))
      }
      #Add cluster with lowest EMD to gold standard from distances vector to
      #current Patient_ID
      assignments[q,2] <- which(distances == min(distances))
    }
    #Order the assignments data frame and bind with gold standard cluster assignments
    compare2 <- assignments[order(assignments[,1]),]
    combined <- cbind(compare, compare2)
    combined[,3] <- NULL
    colnames(combined) <- c("Patient_ID","all_dat","data_removal")
    #For each cluster calculate Jaccard index from current simulation
    for (m in 1:n_clust) {
      #Filter for cluster in gold standard
      overlap1 <- dplyr::filter(combined, all_dat == m)
      #Filter for cluster after random data removal
      overlap2 <- dplyr::filter(combined, data_removal == m)
      #Determine intersection between gold standard and random data removal
      common <- nrow(intersect(overlap1, overlap2))
      #Determine unification quantity
      all <- nrow(union(overlap1, overlap2))
      #Calculate Jaccard index
      jaccard <- common/all
      #Store jaccard index for current simulation and current cluster in summary
      summary[i,m] <- jaccard
    }
  }
  summary
}

#' Simulate random data removal from time series data list and determine Jaccard index via Cognate Cluster approach for multiple random data removal steps
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param clust_num Cluster of interest
#' @param n_clust Number of clusters
#' @param range Range to simulate random data removal (e.g. c(0.1,0.2,0.5,0.7,0.8))
#'
#' @return Object of type list storing Jaccard indices for each indicated random data removal step and visualized results in a boxplot
#'
#' @import tibble
#' @import dplyr
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be retrieved from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo/phys)
#' #Sampling frequency is supposed to be daily
#' output <- jaccard_run_cognate(list,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
#'
#' @export
jaccard_run_cognate <- function(plist, parameter, n_simu, method, clust_num, n_clust, range) {

  #List to store Jaccard indices for each simulation
  jaccard_list = list()
  #Determines how many times the function sim_jaccard_cognate must be called for
  #each random data removal step
  n <- length(range)
  #For each step in range
  for (i in 1:n) {
    #Simulate random data removal and Jaccard index determination by cognate cluster
    #approach n_simu times and add Jaccard indices for each cluster to list
    jaccard_list[[as.character(range[i])]] <- sim_jaccard_cognate(plist, parameter, range[i], n_simu = n_simu, method = method, n_clust = n_clust)
  }
  #Plotlist to easily visualize boxplots for each step
  plotlist = list()
  #For each step in range
  for (j in 1:n) {
    #Add Jaccard indices for current step and cluster of interest to list
    #This list is a preparation for boxplots
    plotlist[[as.character(range[j])]] <- jaccard_list[[as.character(range[j])]][,as.numeric(clust_num)]
  }
  #Visualize the results via boxplot
  boxplot(plotlist, main = paste("Random Data Removal: ",parameter, " Cluster Number: ",clust_num),
          names = range, ylab = "Jaccard index", xlab = "Data removal", ylim = c(0,1))
  jaccard_list
}

#' Simulate random data removal from time series data list and determine Jaccard index via Earth Mover's Distance approach for multiple random data removal steps
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param clust_num Cluster of interest
#' @param n_clust Number of clusters
#' @param range Range to simulate random data removal (e.g. c(0.1,0.2,0.5,0.7,0.8))
#'
#' @return Object of type list storing Jaccard indices for each indicated random data removal step and visualized results in a boxplot
#'
#' @import tibble
#' @import dplyr
#' @import emdist
#'
#' @examples
#' list <- patient_list('.../ts_demofiles1') #Just folder; files can be retrieved from GitHub demo files
#' #(https://github.com/MrMaximumMax/FBCanalysis/tree/master/demo/phys)
#' #Sampling frequency is supposed to be daily
#' output <- jaccard_run_emd(list,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
#'
#' @export
jaccard_run_emd <- function(plist, parameter, n_simu, method, clust_num, n_clust, range) {

  #List to store Jaccard indices for each simulation
  jaccard_list = list()
  #Determines how many times the function sim_jaccard_cognate must be called for
  #each random data removal step
  n <- length(range)
  #For each step in range
  for (i in 1:n) {
    #Simulate random data removal and Jaccard index determination by cognate cluster
    #approach n_simu times and add Jaccard indices for each cluster to list
    jaccard_list[[as.character(range[i])]] <- sim_jaccard_emd(plist, parameter, range[i], n_simu = n_simu, method = method, n_clust = n_clust)
  }
  #Plotlist to easily visualize boxplots for each step
  plotlist = list()
  #For each step in range
  for (j in 1:n) {
    #Add Jaccard indices for current step and cluster of interest to list
    #This list is a preparation for boxplots
    plotlist[[as.character(range[j])]] <- jaccard_list[[as.character(range[j])]][,as.numeric(clust_num)]
  }
  #Visualize the results via boxplot
  boxplot(plotlist, main = paste("Random Data Removal: ",parameter, " Cluster Number: ",clust_num),
          names = range, ylab = "Jaccard index", xlab = "Data removal", ylim = c(0,1))
  jaccard_list
}
