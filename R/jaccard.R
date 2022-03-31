#Functions for Jaccard index determination/random data removal

#' Randomly remove data from time series data list
#'
#' Remove a specific amount of data randomly from a time series data list.
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param removal Amount of data removal (0 = 0%, 1 = 100%)
#'
#' @return Object of type list storing patient time series data with indicated amount of data removed randomly
#'
#' @import dplyr
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
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

#' Simulate random data removal via Global Cognate Cluster cluster assignment approach
#'
#' Simulate random data removal for a removal amount with indicated
#' number of simulations from time series data list and determine Jaccard
#' index for all cluster via Global Cognate Cluster cluster assignment approach
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param removal Amount of random data removal to determine Jaccard index
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param n_clust Number of clusters (also see function: \link{clust_matrix})
#' @param maxIter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#'
#' @details The cognate cluster approach works in the manner that first a Gold
#' Standard cluster is determined meaning the cluster assignments without any
#' data removal. Subsequently, random data is removed from the original, complete
#' data and clustering is performed again on the leaky data. The cluster determined
#' to be cognate to the Gold Standard cluster is the one with the highest overlap
#' in cluster members, meaning hte cluster with highest acheived Jaccard index.
#' Afterwards, the Jaccard indices are calculated, comparing cluster members
#' with complete and leaky data, for each cluster.
#'
#' @references Anja Jochmann, Luca Artusio, Hoda Sharifian, Angela Jamalzadeh,
#' Louise J Fleming, Andrew Bush, Urs Frey, and Edgar Delgado-Eckert.
#' Fluctuation-based clustering reveals phenotypes of patients with different
#' asthma severity. ERJ open research, 6(2), 2020.
#'
#' Edgar Delgado-Eckert, Oliver Fuchs, Nitin Kumar, Juha Pekkanen, Jean-Charles
#' Dalphin, Josef Riedler, Roger Lauener, Michael Kabesch, Maciej Kupczyk,
#' Sven-Erik Dahlen, et al. Functional phenotypes determined by fluctuation-based
#' clustering of lung function measurements in healthy and asthmatic cohort participants.
#' Thorax, 73(2):107–115, 2018.
#'
#' @return Object of type matrix storing received Jaccard indices for indicated amount of
#' random data removal for all clusters
#'
#' @import tibble
#' @import dplyr
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' output <- sim_jaccard_global(list, "PEF", 0.05, 10, "hierarchical", 2, 1000)
#'
#' @export
sim_jaccard_global <- function(plist, parameter, removal, n_simu, method, n_clust, maxIter, normalize) {
  #Simulate random data removal and Jaccard index determination by Global Cognate Cluster Approach

  #In case, user did not specify maximum for EMD calculations, the default value is
  #set to a high number (5,000)
  if (missing(maxIter)) {
    maxIter <- 5000
  }
  #If data must be normalized or not
  if (missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  #Make a gold standard EMD matrix for complete data without data removal
  distmat_complete <- emd_matrix(plist, parameter, maxIter = maxIter, normalize = normal)
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
    distmat_rm <- emd_matrix(list_rm, parameter, maxIter = maxIter, normalize = normal)
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

#' Simulate random data removal via Earth Mover's Distance Cognate Cluster cluster
#' assignment approach
#'
#' Simulate random data removal for a removal amount with indicated
#' number of simulations from time series data list and determine Jaccard
#' index for all clusters via Earth Mover's distance cognate cluster assignment
#' approach.
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param removal Amount of random data removal to determine Jaccard index
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param n_clust Number of clusters (also see function: \link{clust_matrix})
#' @param maxIter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#'
#' @return Object of type matrix storing received Jaccard indices for indicated amount
#' of random data removal for all clusters
#'
#' @details This method represents a novel approach and potential complementary
#' method to \link{sim_jaccard_global}. First, clustering is performed on complete
#' data without removal, serving as Gold Standard clusters. For every Gold Standard
#' cluster then, all time series data from all patients is z-normalized and then
#' assumed to be as one Gold Standard distribution. Subsequently, random data is
#' removed form the time series data. Each leaky data distribution is then compared
#' via Earth Mover's Distance to each Gold Standard Distribution. The Gold Standard
#' cluster distribution to which the observed leaky distribution exhibits the
#' lowest Earth Mover's Distance gets the assignment. This process is repeated until
#' every leaky time series data distribution is assigned to a cluster. Afterwards,
#' the Jaccard indices are calculated, comparing cluster members with complete
#' and leaky data, for each cluster.
#'
#' @references Yossi Rubner, Carlo Tomasi, and Leonidas J Guibas. A metric for
#' distributions with applications to image databases. In Sixth International
#' Conference on Computer Vision (IEEE Cat. No. 98CH36271), pages 59–66. IEEE, 1998.
#'
#' @import tibble
#' @import dplyr
#' @import emdist
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' output <- sim_jaccard_emd(list, "PEF", 0.05, 10, "hierarchical", 2, 100)
#'
#' @export
sim_jaccard_emd <- function(plist, parameter, removal, n_simu, method, n_clust, maxIter, normalize) {

  #Simulate random data removal and Jaccard index determination by EMD Approach
  #In case, user did not specify maximum for EMD calculations, the default value is
  if (missing(maxIter)) {
    maxIter <- 5000
  }
  #Indicates if data needs to be normalized or not
  if (missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  #Make a gold standard EMD matrix for complete data without data removal
  distmat_complete <- emd_matrix(plist, parameter, maxIter = maxIter, normalize = normal)
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
      if (missing(normalize) || normalize == TRUE) {
      norm_a <- (((distr_a[,parameter] - min(distr_a[,parameter]))/(max(distr_a[,parameter]) - min(distr_a[,parameter]))))
      } else {
        norm_a <- distr_a[,parameter]
      }
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
        if (missing(normalize) || normalize == TRUE) {
        norm_b <- as.matrix(as.data.frame(table((distr_b[,parameter] - min(distr_b[,parameter]))/(max(distr_b[,parameter]) - min(distr_b[,parameter])))))
        } else {
          norm_b <- as.matrix(as.data.frame(table(distr_b[,parameter])))
        }
        #Calculate EMD between current random data removal distribution and
        #current gold standard distribution to distances vector
        distances[r] <- emd(norm_b, as.matrix(as.data.frame(table(masterlist[[r]]))), maxIter = maxIter)
      }
      #Add cluster with lowest EMD to gold standard from distances vector to
      #current Patient_ID
      #Check up; EMD might be 0 if data removal did not affect it at all
      if (min(distances) == 0) {
        assignments[q,2] <- compare[q,2]
      } else {
      assignments[q,2] <- which(distances == min(distances))
      }
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
      common <- length(intersect(overlap1[,1], overlap2[,1]))
      #Determine unification quantity
      all <- length(union(overlap1[,1], overlap2[,1]))
      #Calculate Jaccard index
      jaccard <- common/all
      #Store jaccard index for current simulation and current cluster in summary
      summary[i,m] <- jaccard
    }
  }
  summary
}

#' Simulate random data removal range via Global Cognate Cluster cluster assignment
#' approach
#'
#' Simulate amount of random data removal from time series data list and determine
#' Jaccard index via Global Cognate Cluster approach for multiple random data removal
#' steps for a specific cluster of interest.
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param clust_num Cluster of interest
#' @param n_clust Number of clusters
#' @param range Range to simulate random data removal (e.g. c(0.1,0.2,0.5,0.7,0.8))
#' @param maxIter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#'
#' @return Object of type list storing Jaccard indices for each indicated random data removal step and visualized results in a boxplot
#'
#' @import tibble
#' @import dplyr
#'
#' @details See \link{sim_jaccard_global} for more detailed approach on Jaccard
#' index determination. The difference in this function is that now only one cluster
#' is observed für multiple amounts of random data removal where for each data
#' removal step defined the resulting Jaccard indices are stored in a list object.
#' Furthermore, a boxplot visualization is generated, in the style of recent
#' publications.
#'
#' @references Anja Jochmann, Luca Artusio, Hoda Sharifian, Angela Jamalzadeh,
#' Louise J Fleming, Andrew Bush, Urs Frey, and Edgar Delgado-Eckert.
#' Fluctuation-based clustering reveals phenotypes of patients with different
#' asthma severity. ERJ open research, 6(2), 2020.
#'
#' Edgar Delgado-Eckert, Oliver Fuchs, Nitin Kumar, Juha Pekkanen, Jean-Charles
#' Dalphin, Josef Riedler, Roger Lauener, Michael Kabesch, Maciej Kupczyk,
#' Sven-Erik Dahlen, et al. Functional phenotypes determined by fluctuation-based
#' clustering of lung function measurements in healthy and asthmatic cohort participants.
#' Thorax, 73(2):107–115, 2018.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' output <- jaccard_run_global(list,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
#'
#'
#' @export
jaccard_run_global <- function(plist, parameter, n_simu, method, clust_num, n_clust, range, maxIter, normalize) {

  if(missing(maxIter)) {
    Iter <- 5000
  } else {
    Iter <- maxIter
  }
  if(missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  #List to store Jaccard indices for each simulation
  jaccard_list = list()
  #Determines how many times the function sim_jaccard_global must be called for
  #each random data removal step
  n <- length(range)
  #For each step in range
  for (i in 1:n) {
    #Simulate random data removal and Jaccard index determination by cognate cluster
    #approach n_simu times and add Jaccard indices for each cluster to list
    jaccard_list[[as.character(range[i])]] <- sim_jaccard_global(plist, parameter,
                                                                  range[i], n_simu = n_simu,
                                                                  method = method, n_clust = n_clust,
                                                                  maxIter = Iter, normalize = normal)
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

#' Simulate random data removal range via Earth Mover's Distance Cognate Cluster
#' cluster assignment approach
#'
#' Simulate amount of random data removal from time series data list and determine
#' Jaccard index via Earth Mover's Distance approach for multiple random data
#' removal steps for a specific cluster of interest.
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param clust_num Cluster of interest
#' @param n_clust Number of clusters
#' @param range Range to simulate random data removal (e.g. c(0.1,0.2,0.5,0.7,0.8))
#' @param maxIter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#'
#' @return Object of type list storing Jaccard indices for each indicated random data
#' removal step and visualized results in a boxplot
#'
#' @import tibble
#' @import dplyr
#' @import emdist
#'
#' @details See \link{sim_jaccard_emd} for more detailed approach on Jaccard
#' index determination. The difference in this function is that now only one cluster
#' is observed für multiple amounts of random data removal where for each data
#' removal step defined the resulting Jaccard indices are stored in a list object.
#' Furthermore, a boxplot visualization is generated, in the style of recent
#' publications.
#'
#' @references Yossi Rubner, Carlo Tomasi, and Leonidas J Guibas. A metric for
#' distributions with applications to image databases. In Sixth International
#' Conference on Computer Vision (IEEE Cat. No. 98CH36271), pages 59–66. IEEE, 1998.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' output <- jaccard_run_emd(list,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
#'
#'
#' @export
jaccard_run_emd <- function(plist, parameter, n_simu, method, clust_num, n_clust, range, maxIter, normalize) {

  if(missing(maxIter)) {
    maxIter <- 5000
  }
  if(missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  #List to store Jaccard indices for each simulation
  jaccard_list = list()
  #Determines how many times the function sim_jaccard_global must be called for
  #each random data removal step
  n <- length(range)
  #For each step in range
  for (i in 1:n) {
    #Simulate random data removal and Jaccard index determination by cognate cluster
    #approach n_simu times and add Jaccard indices for each cluster to list
    jaccard_list[[as.character(range[i])]] <- sim_jaccard_emd(plist, parameter, range[i],
                                                              n_simu = n_simu, method = method,
                                                              n_clust = n_clust, maxIter = maxIter,
                                                              normalize = normal)
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

#' Determine reappearance frequency
#'
#' Determine the reappearance frequency for clustered elements in perturbed data
#' for specific amount of random data removal
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param removal Amount of random data removal to determine Jaccard index
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param n_clust Number of clusters (also see function: \link{clust_matrix})
#' @param Iter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#'
#' @return Vector of length of n_simu where reappearance frequency is stored for each simulation run
#'
#' @details To begin, each participant's collection of measured values is randomly
#' depleted of a certain proportion of measurements. The clustering technique is
#' then performed using the perturbed data, and the resulting clusters are compared
#' to the original clusters created with the unperturbed gold standard.
#' This technique is done iteratively in order to provide statistics indicating
#' the original clusters' stability following the elimination of random data.
#' These stability statistics are calculated using two cluster similarity metrics:
#' Jaccard's index, a measure of global similarity that quantifies the extent to
#' which the original and modified clusters overlap. Additionally, it is utilized
#' to decide which cluster is considered the cognate cluster. Then it was transformed
#' into a local measure, meaning the frequency with which each member of the original
#' clusters reappeared between iterations.
#'
#' @references Edgar Delgado-Eckert, Oliver Fuchs, Nitin Kumar, Juha Pekkanen,
#' Jean-Charles Dalphin, Josef Riedler, Roger Lauener, Michael Kabesch, Maciej Kupczyk,
#' Sven-Erik Dahlen, et al. Functional phenotypes determined by fluctuation-based
#' clustering of lung function measurements in healthy and asthmatic cohort participants.
#' Thorax, 73(2):107–115, 2018.
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' #' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' output <- reap_freq(list,"PEF",0.1,10,"hierarchical,2,5000)
#'
reap_freq <- function(plist, parameter, removal, n_simu, method, n_clust, Iter, normalize) {
  #Simulate random data removal and Jaccard index determination by Cognate Cluster Approach

  #In case, user did not specify maximum for EMD calculations, the default value is
  #set to a high number (5,000)
  if(missing(Iter)) {
    Iter <- 5000
  }
  if(missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  master_vector <- vector(length = n_simu)
  #Make a gold standard EMD matrix for complete data without data removal
  distmat_complete <- emd_matrix(plist, parameter, maxIter = Iter, normalize = normal)
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
    reap_vector <- vector(length = length(plist))
    #Remove randomly data from time series data list according to specified removal amount
    list_rm <- rnd_dat_rm(plist, removal)
    #Calculate a new EMD matrix after random data removal
    distmat_rm <- emd_matrix(list_rm, parameter, maxIter = Iter, normalize = normal)
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
    patnames <- names(plist)
    for (m in 1:length(patnames)) {
      current_gold_cluster <- dplyr::filter(combined, Patient_ID == patnames[m])[,2]
      overlap1 <- dplyr::filter(combined, all_dat == current_gold_cluster)[,1]
      #Make a vector to store all calculated Jaccard indices; Highest Jaccard
      #index then represents the cognate cluster
      jaccard_vector <- vector()
      #For-loop to filter for each data removal cluster
      for (n in 1:n_clust) {
        #Filter for data removal cluster
        overlap2 <- dplyr::filter(combined, data_removal == n)[,1]
        #Find intersection between current gold standard cluster and current
        #data removal cluster
        common <- length(intersect(overlap1, overlap2))
        #Quantify all; Both current gold standard and data removal cluster
        all <- length((union(overlap1, overlap2)))

        #Calculate current Jaccard index
        jaccard <- common/all
        #Add Jaccard index to vector
        jaccard_vector[n] <- jaccard
      }
        cognate_cluster <- which.max(jaccard_vector)
        members <- dplyr::filter(combined, data_removal == cognate_cluster)
        members <- members[,1]
        if (patnames[m] %in% members) {
          reap_vector[m] <- 1
        } else {
          reap_vector[m] <- 0
        }
    }
    current_freq <- sum(reap_vector)/(length(reap_vector))
    master_vector[i] <- current_freq
  }
  master_vector
}

#' Determine reappearance frequency for range
#'
#' Determine the reappearance frequency for clustered elements in perturbed data
#' for specified amounts of random data removal
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param n_clust Number of clusters
#' @param range Range to simulate random data removal (e.g. c(0.1,0.2,0.5,0.7,0.8))
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#' @param Iter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#'
#' @return List storing vectors with reappearance frequencies for each simulation run
#' for each removal amount indicated as well as illustrated results in a boxplot
#'
#' @details The procedure on determining the reappearance frequenncy is described
#' in \link{reap_freq}.
#'
#' @references Edgar Delgado-Eckert, Oliver Fuchs, Nitin Kumar, Juha Pekkanen,
#' Jean-Charles Dalphin, Josef Riedler, Roger Lauener, Michael Kabesch, Maciej Kupczyk,
#' Sven-Erik Dahlen, et al. Functional phenotypes determined by fluctuation-based
#' clustering of lung function measurements in healthy and asthmatic cohort participants.
#' Thorax, 73(2):107–115, 2018.
#'
#' @import dplyr
#'
#' @export
#'
#' @examples
#' #' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' output <- reap_freq(list,"PEF",,10,"hierarchical,2,c(0.01,0.05,0.1,0.2))
#'
reap_freq_run <- function(plist, parameter, n_simu, method, n_clust, range, Iter, normalize) {

  if(missing(Iter)) {
    Iter <- 5000
  }
  if(missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  #List to store frequencies for each removal amount
  freq_list = list()
  #Determines how many times the function sim_jaccard_global must be called for
  #each random data removal step
  n <- length(range)
  #For each step in range
  for (i in 1:n) {
    #Simulate random data removal and Jaccard index determination by cognate cluster
    #approach n_simu times and add Jaccard indices for each cluster to list
    freq_list[[as.character(range[i])]] <- reap_freq(plist, parameter, range[i], n_simu = n_simu,
                                                     method = method, n_clust = n_clust, Iter = Iter,
                                                     normalize = normal)
  }
  #Visualize the results via boxplot
  boxplot(freq_list, main = paste("Reappearance frequency: ", parameter),
          names = range, ylab = "Reappearance frequency", xlab = "Data removal", ylim = c(0,1))
  freq_list
}

#' Simulate random data removal via alternative Earth Mover's Distance
#' Cognate Cluster cluster assignment approach
#'
#' Simulate random data removal for a removal amount with indicated
#' number of simulations from time series data list and determine Jaccard
#' index for all clusters via Earth Mover's distance cognate cluster assignment
#' approach.
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param removal Amount of random data removal to determine Jaccard index
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param n_clust Number of clusters (also see function: \link{clust_matrix})
#' @param maxIter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#'
#' @return Object of type matrix storing received Jaccard indices for indicated amount
#' of random data removal for all clusters
#'
#' @details This method represents a novel approach and potential complementary
#' method to \link{sim_jaccard_global} and alternative to \link{sim_jaccard_emd}.
#' First, clustering is performed on complete data without removal, serving as
#' Gold Standard clusters. Subsequently, random data is removed form the time series
#' data. Each leaky data distribution is then compared via Earth Mover's Distance
#' to each member's distribution of each Gold Standard cluster The Gold Standard
#' cluster to which the observed leaky distribution exhibits the lowest avergae
#' Earth Mover's Distance gets the assignment. This process is repeated until
#' every leaky time series data distribution is assigned to a cluster. Afterwards,
#' the Jaccard indices are calculated, comparing cluster members with complete
#' and leaky data, for each cluster.
#'
#' @references Yossi Rubner, Carlo Tomasi, and Leonidas J Guibas. A metric for
#' distributions with applications to image databases. In Sixth International
#' Conference on Computer Vision (IEEE Cat. No. 98CH36271), pages 59–66. IEEE, 1998.
#'
#' @import tibble
#' @import dplyr
#' @import emdist
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' output <- sim_jaccard_emd_2(list, "PEF", 0.05, 10, "hierarchical", 2, 100)
#'
#' @export
sim_jaccard_emd_2 <- function(plist, parameter, removal, n_simu, method, n_clust, maxIter, normalize) {

  #Simulate random data removal and Jaccard index determination by EMD Approach
  #In case, user did not specify maximum for EMD calculations, the default value is
  if (missing(maxIter)) {
    maxIter <- 5000
  }
  #Indicates if data needs to be normalized or not
  if (missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  #Make a gold standard EMD matrix for complete data without data removal
  distmat_complete <- emd_matrix(plist, parameter, maxIter = maxIter, normalize = normal)
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
      if (missing(normalize) || normalize == TRUE) {
        norm_a <- (((distr_a[,parameter] - min(distr_a[,parameter]))/(max(distr_a[,parameter]) - min(distr_a[,parameter]))))
      } else {
        norm_a <- distr_a[,parameter]
      }
      #Add current distribution to list
      dat[[patnames[u]]] <- norm_a
    }
    #Add for each current cluster in masterlist the normalized time series values
    #of each cluster members
    masterlist[[t]] <- dat
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
    #For every Patient_ID in list after random data removal
    for (q in 1:length(list_rm)) {
      #New vector to store EMD for current Patient_ID to each gold standard cluster
      averaged_distances <- vector()
      current_pat <- names(list_rm)[q]
      #For every gold standard cluster
      for(r in 1:length(masterlist)) {
        #For every distribution in current gold standard cluster
        averaged_dist_current_clust <- vector()
        for(s in 1:length(masterlist[[r]])) {
        #Take current random data removal distribution and normalize if necessary
        distr_a <- list_rm[[current_pat]]
        if (missing(normalize) || normalize == TRUE) {
          norm_a <- as.matrix(as.data.frame(table((distr_a[,parameter] - min(distr_a[,parameter]))/(max(distr_a[,parameter]) - min(distr_a[,parameter])))))
        } else {
          norm_a <- as.matrix(as.data.frame(table(distr_a[,parameter])))
        }
        #Take current gold standard patient distribution and normalize if necessary
        distr_b <- masterlist[[r]][[s]]
        if (missing(normalize) || normalize == TRUE) {
          norm_b <- as.matrix(as.data.frame(table((distr_b - min(distr_b))/(max(distr_b) - min(distr_b)))))
        } else {
          norm_b <- as.matrix(as.data.frame(table(distr_b)))
        }
        #Calculate EMD between current random data removal distribution and
        #current gold standard distribution to distances vector
        averaged_dist_current_clust[s] <- emd(norm_b, norm_a, maxIter = maxIter)
        }
        averaged_distances[r] <- mean(averaged_dist_current_clust)
      }
      cognate_cluster <- which.min(averaged_distances)
      compare[q,3] <- cognate_cluster
    }
    colnames(compare) <- c("Patient_ID","all_dat","data_removal")
    for (m in 1:n_clust) {
      #filter for current gold standard cluster
      overlap1 <- dplyr::filter(compare, all_dat == m)
      overlap2 <- dplyr::filter(compare, data_removal == m)
      common <- length(intersect(overlap1[,1],overlap2[,1]))
      all <- length(union(overlap1[,1],overlap2[,1]))
      jaccard <- common/all
      summary[i,m] <- jaccard
    }
  }
  summary
}

#' Simulate random data removal range via alternative Earth Mover's Distance
#' Cognate Cluster cluster assignment approach
#'
#' Alternative: Simulate amount of random data removal from time series data list
#' and determine Jaccard index via Earth Mover's Distance approach for multiple
#' random data removal steps for a specific cluster of interest.
#'
#' @param plist Object of type list storing patient time series data (also see function: \link{patient_list})
#' @param parameter Parameter of interest in time series data list
#' @param n_simu Number of simulations
#' @param method Clustering method (also see function: \link{clust_matrix})
#' @param clust_num Cluster of interest
#' @param n_clust Number of clusters
#' @param range Range to simulate random data removal (e.g. c(0.1,0.2,0.5,0.7,0.8))
#' @param maxIter Maximum iterations to determine Earth Mover's Distances (also see function: \link{emd_matrix});
#' default is 5,000 for this function
#' @param normalize Indicates if parameter indicated needs to be normalized or not (TRUE by default)
#'
#' @return Object of type list storing Jaccard indices for each indicated random data
#' removal step and visualized results in a boxplot
#'
#' @import tibble
#' @import dplyr
#' @import emdist
#'
#' @details See \link{sim_jaccard_emd_2} for more detailed approach on Jaccard
#' index determination. The difference in this function is that now only one cluster
#' is observed für multiple amounts of random data removal where for each data
#' removal step defined the resulting Jaccard indices are stored in a list object.
#' Furthermore, a boxplot visualization is generated, in the style of recent
#' publications.
#'
#' @references Yossi Rubner, Carlo Tomasi, and Leonidas J Guibas. A metric for
#' distributions with applications to image databases. In Sixth International
#' Conference on Computer Vision (IEEE Cat. No. 98CH36271), pages 59–66. IEEE, 1998.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' output <- jaccard_run_emd_2(list,"PEF",10,"hierarchical",1,3,c(0.005,0.01,0.05,0.1,0.2))
#'
#'
#' @export
jaccard_run_emd_2 <- function(plist, parameter, n_simu, method, clust_num, n_clust, range, maxIter, normalize) {

  if(missing(maxIter)) {
    maxIter <- 5000
  }
  if(missing(normalize) || normalize == TRUE) {
    normal <- TRUE
  } else {
    normal <- FALSE
  }
  #List to store Jaccard indices for each simulation
  jaccard_list = list()
  #Determines how many times the function sim_jaccard_global must be called for
  #each random data removal step
  n <- length(range)
  #For each step in range
  for (i in 1:n) {
    #Simulate random data removal and Jaccard index determination by cognate cluster
    #approach n_simu times and add Jaccard indices for each cluster to list
    jaccard_list[[as.character(range[i])]] <- sim_jaccard_emd_2(plist, parameter, range[i],
                                                              n_simu = n_simu, method = method,
                                                              n_clust = n_clust, maxIter = maxIter,
                                                              normalize = normal)
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
