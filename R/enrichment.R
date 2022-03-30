#General functionalities to cluster data and add cluster assignments to ts- and extended data

#' Add extended data
#'
#' Add extended data from a csv file, match with Patient ID entries from a
#' previously generated time series data list and preprocess for further analysis.
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list})
#' @param path Path where enrichment csv file is stored
#'
#' @return Processed data as object of type data frame; Enrichment data Patient_IDs are matched with Time Series Data List Patient IDs; In case it was indicated, NA values in the enrichment data are filled up by random sampling
#'
#' @import utils
#' @import dplyr
#'
#' @details the extended csv file should have a column including the Patient ID.
#' Additionally, one specifies the list in which time series data is saved.
#' This is advantageous since the function can now do matching, i.e. determine
#' which Patient IDs occur in both the extended dataset and time series datalist.
#' So for a result, any Patient ID that appears in the extended dataset but
#' does not exist in the time series datalist will be deleted from the extended
#' dataset, as it cannot be used in any further investigation. Nonetheless,
#' Patient IDs from the time series data that do not appear in the enrichment dataset
#' will be added to the enrichment dataset, but each new parameter will be featureless,
#' so added as NA value.
#'
#' If one selects option 1 (leave missing values as NA), no further processing
#' of the input occurs. The extended data set will be added to the environment as
#' a data frame. In this situation, the NA values from the extended dataset will be
#' included in the summary indicating, for example, that a certain cluster has a given
#' percentage of missing values. This may also lead to some additional findings, such
#' as that a specific parameter considerably enriches a cluster yet many data is absent.
#' If the one selects option 2 (sample missing values), the function loops over each
#' NA entry and selects a random value from the whole distribution for the parameter
#' for which the data is missing. This cycle is repeated until the whole dataset has
#' been processed and the data will be added as a data frame to the environment.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' enr <- add_enrich(list,
#' 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv')
#'
#' @export
add_enrich <- function(plist, path) {

  #Read csv-file on indicated path; empty fields in csv-file are filled up as NA
  dat <- read.csv(path, na.strings=c("","NA"))

  #Print first 5 lines of csv and present to user
  cat("\n\n") #two empty lines as space before printing
  print(head(dat[1:5,]))
  cat("\n\n") #two empty lines after printing

  #Create a data frame with the colnmaes and assigned ascending numbers to the colnames
  feedback <- matrix(nrow = 2, ncol = length(colnames(dat)))
  feedback[1,] <- colnames(dat)
  feedback[2,] <- 1:length(dat)
  feedback <- as.data.frame(feedback)
  #Print the created feedback data frame
  print(feedback, row.names=F, col.names=F, quote=F)
  cat("\n")
  #Ask the user which number (= which column) represents the Patient_ID
  #This is of interest since later functionalities recognize the Patient_ID and
  #do not consider them as a population parameter itself
  question1 <- as.numeric(readline("Which column represents the Patient_ID (please indicate number)?: "))
  #Rename column as "Patient_ID" so that later functionalities always recognize
  #the patient column
  names(dat)[names(dat) == feedback[1,question1]] <- 'Patient_ID'

  #The enrichment data frame should be so prepared that Patient_ID's that appear in
  #the indicated csv-file but not in the patient data list are thrown out and
  #Patient_ID's that appear in the patient data list but not in the csv-file are
  #added to the data frame and every parameter is filled with NA

  #Get Patient_ID's from both previously loaded csv-file and patient time
  #series data list
  names_enr <- as.character(as.data.frame(table(dat[,'Patient_ID']))[,1])
  names_ts <- names(plist)

  #Check which patients are in enrich but not ts and throw them out
  throw <- as.vector(dplyr::setdiff(names_enr,names_ts))
  #Throw out Patient_ID's that do not appear in time series data list
  dat_new1 <- dplyr::filter(dat, !Patient_ID %in% throw)

  #Check which patients are in ts but not enrich and add them with NAs
  add <- dplyr::setdiff(names_ts,names_enr)

  #In case there are Patient_ID's to add to the enrichment data frame
  if (length(add > 0)) {
    #Make a new matrix with same number of columns and colnames than csv-file
    #Number of rows are number of Patient_ID's to be added
    dat_new2 <- matrix(data = NA, nrow = length(add), ncol = length(colnames(dat)))
    #Transform matrix to data frame
    dat_new2 <- as.data.frame(dat_new2)
    #Change the colnames to the colnames of the csv-file
    colnames(dat_new2) <- colnames(dat)
    #Add missing Patient_ID's to the Patient_ID column
    dat_new2[,"Patient_ID"] <- add
    #Rowbind the previously created data frame with the original enrichment data frame
    dat_new <- rbind(dat_new1,dat_new2)
  } else {
    #In case there are no Patient_ID's to add, the previous data frame can be taken over
    dat_new <- dat_new1
  }

  #In case there are NA entries in the data frame
  if (any(is.na(dat_new))) {
    #Ask user how NA values shoud be handled
    question <- readline("Do you want to leave missing values as NA (1) or let them sample (2)?: ")

    if (question == 2) {
      #Go through each row and column
      for (i in 1:ncol(dat_new)) {
        for (j in 1:nrow(dat_new)) {
          #Check in each cell if value is NA
          if (is.na(dat_new[j,i])) {
            #Take all non-NA values from column
            samp <- na.omit(dat_new[,i])
            #Randomly sample one value on NA position
            dat_new[j,i] <- sample(samp,1)
          }
        }
      }
    }
  }
  dat_new
}

#' Add clustering assignments to extended data
#'
#' Add clustering assignments from clustering output to pre-processed extended/enrichment
#' data frame.
#'
#' @param enrich Preprocessed enrichment data frame (also see function: \link{add_enrich})
#' @param clustdat Object of type list storing clustering data (also see function: \link{clust_matrix})
#'
#' @return Processed data frame with added column indicating cluster assignments
#'
#' @import tibble
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' clustering <- clust_matrix(matrix, method = "kmeans", nclust = 3)
#' enr <- add_enrich(list,
#' 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv')
#' enr <- add_clust2enrich(enr, clustering)
#'
#' @export
add_clust2enrich <- function(enrich, clustdat) {

  #Get vector from clustdat list where cluster assignments to patients are stored
  #(also see function clust_maxtrix())
  new <- as.data.frame(clustdat$Cls)
  new <- tibble::rownames_to_column(new)
  #Ordering by Patient_ID later facilitates assignment since data frame "new"
  #and enrichment data frame have the same order so only cbind needs to be done
  new <- new[order(new$rowname),]
  #Remove column with Patient_ID's afterwards, otherwise Patient_ID would be
  #stored twice in final data frame
  new[,"rowname"] <- NULL
  #Rename cluster column
  colnames(new) <- "Cluster"
  #Do similar as above: Order enrichment data frame by Patient_ID
  enrich <- enrich[order(enrich$Patient_ID),]
  #Combine both data frame
  enrich <- cbind(enrich,new)
}

#' Add clustering assignments to time series data
#'
#' Add clustering assignments from clustering output to time series data list
#' and store the data in a data frame.
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list})
#' @param clustdat Object of type list storing clustering data (also see function: \link{clust_matrix})
#'
#' @return Processed data frame storing time series data with added column indicating cluster assignments
#'
#' @import tibble
#' @import dplyr
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' clustering <- clust_matrix(matrix, method = "kmeans", nclust = 3)
#' ts <- add_clust2ts(list, clustering)
#'
#' @export
add_clust2ts <- function(plist, clustdat) {

  #Make one data frame out of all dataframe within time series data list
  ts_df <- do.call(rbind,plist)
  #New list to store the processed data
  datalist <- list()
  #Get vector from clustdat list where cluster assignments to patients are stored
  #(also see function clust_maxtrix())
  new <- as.data.frame(clustdat$Cls)
  #Add Patient_ID's to extra column (helps for later processes)
  new <- tibble::rownames_to_column(new)
  #Rename the columns
  colnames(new) <- c("Patient_ID","Cluster")

  #For-loop for every cluster
  for (i in 1:max(new[,"Cluster"])) {
    #Filter for current cluster number
    compare <- dplyr::filter(new, Cluster %in% i)
    #Remove column "Cluster"
    compare[,"Cluster"] <- NULL
    #Make a vector, only containing the Patient_ID's for current cluster
    compare <- as.character(compare[,1])
    #Filter time series data frame for Patient_ID's from current cluster
    new2 <- dplyr::filter(ts_df, Patient_ID %in% compare)
    #Add a new column "Cluster" to time series data frame storing cluster number
    new2[,"Cluster"] <- i
    #Add current cluster time series data to list
    datalist[[i]] <- new2
  }
  #Make one time series data frame out of all list elements
  assigned <- do.call(rbind,datalist)
}

#' Observe a specific cluster of interest
#'
#' Observe a specific cluster of interest of preporcessed extended and time
#' series data for overview and p-values.
#'
#' @param ts.dat Processed data frame storing time series data and cluster assignments (also see function: \link{add_clust2ts})
#' @param enrich Processed data frame storing enrichment data and cluster assignments (also see function: \link{add_clust2enrich})
#' @param clustno Cluster number of interest
#' @param numeric Statistical test to be performed on continuous data
#' @param categorical Statistical test to be performed on categorical data
#'
#' @return Terminal output presenting summary of time series and enrichment data with corresponding p-values
#'
#' @import arsenal
#' @import dplyr
#'
#' @details There are five techniques available to compute the p-value for continuous
#' or categorical data. In order to determine the relevant p-value inside a cluster
#' of interest, the data distribution within the cluster should be compared to the
#' data distribution outside the cluster. Prior to conducting the related probability
#' tests, one data processing step is performed, namely the construction of two data
#' distributions, one including only data included inside the cluster and another
#' comprising data from outside the cluster, for the purpose of comparing them.
#'
#' **Mann-Whitney Test** *(numeric = "wt")*
#' The Mannn-Whitney Test, sometimes referred to as the Wilcoxon rank-sum test (WRS),
#' is used to measure the significance of continuous variables within the observed
#' distribution. The WRS is used to test if the central tendency of two independent
#' samples is different. When the t-test for independent samples does not meet the
#' requirements, the WRS is used. The null hypothesis H0 states that the populations’
#' distributions are equal. H1 is the alternative hypothesis meaning that the
#' distributions are not equal. The test is consistent under the broader formulation
#' only when the following happens under H1.
#'
#' **Analysis of Variance** *(numeric = "anova")*
#' A further approach on determining significance for continuous distributions is the
#' Analysis of Variance (ANOVA). The perquisites for ANOVA are that samples are
#' sampled independently from each other. Additionally, variance homogeneity and
#' normal distribution must be given. A independent variable, consisting
#' of I categories should be given. H0 indicates that no differences in the means for
#' each I is given. It is used to compare two or more independent samples with
#' comparable or dissimilar sample sizes.
#'
#' **Kruskall-Wallis test** *(numeric = "kwt")*
#' With several samples are not normally distributed but also small sample sizes
#' and outliers, the Kruskall-Wallis test may be preferred. It is used to compare
#' two or more independent samples with comparable or dissimilar sample sizes.
#' It expands the WRS, which is only used to compare two groups. If the researcher
#' can make the assumption that all groups have an identically shaped and scaled
#' distribution, except for differences in medians, H0 is that all groups have
#' equal medians.
#'
#' **Fisher’s exact test** *(categorical = "fe")*
#' The hypergeometric test or Fisher’s exact test (FET) is used to analyze categorical
#' variables within the enriched data set. It is a statistical significance test
#' for contingency tables that is employed in the study of them. The test is helpful
#' for categorical data derived from object classification. It is used to assess
#' the importance of associations and inconsistencies between classes. The FET is
#' often used in conjunction with a 2 × 2 contingency table that represents two
#' categories for a variable, as well as assignment inside or outside of the cluster.
#' The p-value is calculated as if the table’s margins are fixed. This results in
#' a hypergeometric distribution of the numbers in the table cells under the null
#' hypothesis of independence. A hypergeometric distribution is a discrete probability
#' distribution that describes the probability of k successes, defined as random draws
#' for which the object drawn has a specified feature in n draws without replacement
#' from a finite population of size N containing exactly K objects with that feature,
#' where each draw is either successful or unsuccessful. The test is only practicable for
#' normal computations in the presence of a 2 × 2 contingency table. However, the
#' test’s idea may be extended to the situation of a m × n table in general.
#' Statistics programs provide a Monte Carlo approach for approximating the more
#' general case.
#'
#' **Chi-Square test** *(categorical = "chisq")*
#' Additionally, one may choose to do a Chi-Square test. This is a valid statistical
#' hypothesis test when the test statistic is normally distributed under the null
#' hypothesis. According to Pear- son, the difference between predicted and actual
#' frequencies in one or more categories of a contingency table is statistically
#' significant.
#'
#' @references Siegel Sidney. Nonparametric statistics for the behavioral sciences.
#' The Journal of Nervous and Mental Disease, 125(3):497, 1957.
#'
#' Kinley Larntz. Small-sample comparisons of exact levels for chi-squared
#' goodness-of-fit statistics. Journal of the American Statistical Association,
#' 73(362):253–263, 1978.
#'
#' Cyrus R Mehta and Nitin R Patel. A network algorithm for performing fisher’s
#' exact test in r× c contingency tables. Journal of the American Statistical Association,
#' 78(382):427–434, 1983.
#'
#' Aravind Subramanian, Pablo Tamayo, Vamsi K Mootha, Sayan Mukherjee, Benjamin
#' L Ebert, Michael A Gillette, Amanda Paulovich, Scott L Pomeroy, Todd R Golub,
#' Eric S Lander, et al. Gene set enrichment analysis: a knowledge-based approach
#' for interpreting genome-wide expression profiles. Proceedings of the National
#' Academy of Sciences, 102(43):15545–15550, 2005.
#'
#' William H Kruskal and W Allen Wallis. Errata: Use of ranks in one-criterion variance
#' analysis. Journal of the American statistical Association, 48(264):907–911, 1953.
#'
#' Kinley Larntz. Small-sample comparisons of exact levels for chi-squared goodness-
#' of-fit statistics. Journal of the American Statistical Association,
#' 73(362):253–263, 1978.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' clustering <- clust_matrix(matrix, method = "kmeans", nclust = 3)
#' enr <- add_enrich(list,
#' 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv')
#' enr <- add_clust2enrich(enr, clustering)
#' ts <- add_clust2ts(list, clustering)
#' enr_obs_clust(ts, enr, 1, numeric = "anova", categorical = "fe")
#'
#' @export
enr_obs_clust <- function(ts.dat, enrich, clustno, numeric, categorical) {

  #Methods only applicable on pre-processed enrichment and time series data
  #from functions: add_clust2enrich & add_clust2ts

  #Remove "Patient_ID" column from enrichment data, otherwise they would be
  #recognized as categorical variable for summary
  enrich[,"Patient_ID"] <- NULL
  #Split the data so that p-values can be calculated
  #(in cluster = 0 vs. not in cluster = 1)
  enrich_in <- dplyr::filter(enrich, Cluster == clustno)
  enrich_in[,"Cluster"] <- 0
  enrich_out <- dplyr::filter(enrich, Cluster != clustno)
  enrich_out[,"Cluster"] <- 1
  #Now bind them both together again
  enrich_new <- rbind(enrich_in,enrich_out)
  #From package "Arsenal"; Creates and calculates a summary of the data by cluster
  #By default: p-values for continuous variables are calculated by Mann-Whitney test
  #and Hypergeometric test for categorical variable
  #The tableby() function creates a list where all the data is stored
  table_enr <- tableby(Cluster ~., data = enrich_new, numeric.test = numeric, cat.test = categorical)
  #Now this list is transformed to a data frame to print easily the intermediary results
  table_enr <- as.data.frame(summary(table_enr))
  #Remove some reoccuring characters to make it readable
  table_enr[,1] <- gsub("&nbsp;","",as.character(table_enr[,1]))
  #Take out first part of summary, showing the parameters
  col1 <- as.data.frame(table_enr[,1])
  #Take out second part, showing the data for parameters inside cluster
  col2 <- as.data.frame(table_enr[,2])
  #Take out third part, showing the corresponding p-values
  col3 <- as.data.frame(table_enr[,5])
  #Bind the three parts
  table_enr <- cbind(col1,col2,col3)
  #Assign new colnames
  colnames(table_enr) <- c("Parameter"," ","p-value")

  #Count how many patients are in currently observed cluster
  #Filter for current cluster
  n_pat <- dplyr::filter(enrich, Cluster %in% clustno)
  #Length = Number of patients
  n_pat <- length(n_pat[,1])

  #find the Patient_IDs from the cluster of interest
  patnames <- dplyr::filter(ts.dat, Cluster == clustno)
  patnames <- names(table(patnames[,"Patient_ID"]))
  #find the parameters for cluster of interest
  parameters <- colnames(ts.dat)
  parameters <- gsub('Patient_ID', NA, parameters)
  parameters <- gsub('Time', NA, parameters)
  parameters <- gsub('Interpolated', NA, parameters)
  parameters <- gsub('Cluster', NA, parameters)
  parameters <- na.omit(parameters)
  #Start an analogous procedure for time series data
  #Remove both Patient_ID and Time so that only time series data, information
  #on interpolation and cluster assignment remain
  ts.dat[,"Patient_ID"] <- NULL
  ts.dat[,"Time"] <- NULL
  #Split the data so that p-values can be calculated
  #(in cluster = 0 vs. not in cluster = 1)
  ts.dat_in <- dplyr::filter(ts.dat, Cluster == clustno)
  ts.dat_in[,"Cluster"] <- 0
  ts.dat_out <- dplyr::filter(ts.dat, Cluster != clustno)
  ts.dat_out[,"Cluster"] <- 1
  #Now bind them both together again
  ts.dat_new <- rbind(ts.dat_in,ts.dat_out)
  #Make an evaluation by cluster in and out with arsenal functionality
  table_ts <- tableby(Cluster ~., data = ts.dat_new, numeric.test = numeric)
  #Similar procedure as above (line 978)
  table_ts <- as.data.frame(summary(table_ts))
  table_ts[,1] <- gsub("&nbsp;","",as.character(table_ts[,1]))
  table_ts <- table_ts[-c(7:9),]
  cola <- as.data.frame(table_ts[,1])
  colb <- as.data.frame(table_ts[,2])
  colc <- as.data.frame(table_ts[,5])
  table_ts <- cbind(cola,colb,colc)
  colnames(table_ts) <- c("Parameter", " ", "p-value")

  #for summary: number of measurements; number of interpolated data
  #Filter time series data for current cluster and determine length (= number of measurements)
  n_mes <- dplyr::filter(ts.dat, Cluster %in% clustno)
  n_val <- length(n_mes[,1])
  #Filter and determine length for both interpolated and non-interpolated data
  n_real <- length(dplyr::filter(n_mes, Interpolated == FALSE)[,1])
  n_inter <- length(dplyr::filter(n_mes, Interpolated == TRUE)[,1])
  #Determine the ratio interpolated vs. non-interpolated data
  inter_ratio <- round(100*(n_inter/(n_inter + n_real)), digits = 4)
  #Printed feedback on observation of current cluster for user with previously
  #prepared data
  cat("\n\n", "Summary Cluster No.: ", clustno, "\n\n",
      "Number of patients: ", n_pat, "\n", "Number of ts values: ", n_val,
      "\n\n", "Real measurements: ", n_real, "\n", "Interpolated data: ",
      n_inter, "(", inter_ratio, "%)", "\n\n", "Time series data:","\n")
  print(table_ts)
  cat("\n\n", "Population data:","\n")
  print(table_enr)
}

#' Simulate random sampling in extended data
#'
#' Simulate random sampling for NA entries in extended data and
#' check stability of resulting p-values for the parameters for an
#' indicated number of random sampling simulations.
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list})
#' @param path Path where enrichment csv file is stored
#' @param clustdat Object of type list storing clustering data (also see function: \link{clust_matrix})
#' @param clustno Cluster number of interest
#' @param n_sim Number of simulations
#'
#' @return Object of type list storing the received p-values for each parameter in a vector and boxplot visualizing the received p-values
#'
#' @import arsenal
#' @import dplyr
#'
#' @details It allows the sampling in NA entries to be repeated for each parameter
#' in the extended data set. The primary objective here is to validate the random
#' sampling process for missing data by running many simulations and comparing the
#' resultant p-values. An extended data frame with NA elements is saved as a
#' simulation foundation. This data frame will always serve as the foundation
#' for any subsequent simulations added. Following that, the program runs through
#' each NA item in the dataset and generates a random sample of the current
#' parameter’s distribution. After completing this step for each parameter, the
#' function generates the associated p-values as explained in \link{enr_obs_clust}.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' path <- 'https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/enrich/enrichment.csv'
#' test <- sim_sample_enr(list,path,clustering,1,100, numeric = "anova", categorical = "fe")
#' sim_sample_enr <- function(plist, path, clustdat, clustno, n_sim)
#'
#' @export
sim_sample_enr <- function(plist, path, clustdat, clustno, n_sim, numeric, categorical) {

  #Analogous start as in function: add_enrich()
  #Read csv-file on indicated path; empty fields in csv-file are filled up as NA
  dat <- read.csv(path, na.strings=c("","NA"))

  #Print first 5 lines of csv and present to user
  cat("\n\n") #two empty lines as space before printing
  print(head(dat[1:5,]))
  cat("\n\n") #two empty lines after printing

  #Create a data frame with the colnmaes and assigned ascending numbers to the colnames
  feedback <- matrix(nrow = 2, ncol = length(colnames(dat)))
  feedback[1,] <- colnames(dat)
  feedback[2,] <- 1:length(dat)
  feedback <- as.data.frame(feedback)
  #Print the created feedback data frame
  print(feedback, row.names=F, col.names=F, quote=F)
  cat("\n")
  #Ask the user which number (= which column) represents the Patient_ID
  #This is of interest since later functionalities recognize the Patient_ID and
  #do not consider them as a population parameter itself
  question1 <- as.numeric(readline("Which column represents the Patient_ID (please indicate number)?: "))
  #Rename column as "Patient_ID" so that later functionalities always recognize
  #the patient column
  names(dat)[names(dat) == feedback[1,question1]] <- 'Patient_ID'

  #The enrichment data frame should be so prepared that Patient_ID's that appear in
  #the indicated csv-file but not in the patient data list are thrown out and
  #Patient_ID's that appear in the patient data list but not in the csv-file are
  #added to the data frame and every parameter is filled with NA

  #Get Patient_ID's from both previously loaded csv-file and patient time
  #series data list
  names_enr <- as.character(as.data.frame(table(dat[,'Patient_ID']))[,1])
  names_ts <- names(plist)

  #Check which patients are in enrich but not ts and throw them out
  throw <- as.vector(dplyr::setdiff(names_enr,names_ts))
  #Throw out Patient_ID's that do not appear in time series data list
  dat_new1 <- dplyr::filter(dat, !Patient_ID %in% throw)

  #Check which patients are in ts but not enrich and add them with NAs
  add <- dplyr::setdiff(names_ts,names_enr)

  #In case there are Patient_ID's to add to the enrichment data frame
  if (length(add > 0)) {
    #Make a new matrix with same number of columns and colnames than csv-file
    #Number of rows are number of Patient_ID's to be added
    dat_new2 <- matrix(data = NA, nrow = length(add), ncol = length(colnames(dat)))
    #Transform matrix to data frame
    dat_new2 <- as.data.frame(dat_new2)
    #Change the colnames to the colnames of the csv-file
    colnames(dat_new2) <- colnames(dat)
    #Add missing Patient_ID's to the Patient_ID column
    dat_new2[,"Patient_ID"] <- add
    #Rowbind the previously created data frame with the original enrichment data frame
    dat_new <- rbind(dat_new1,dat_new2)
  } else {
    #In case there are no Patient_ID's to add, the previous data frame can be taken over
    dat_new <- dat_new1
  }
  #Add cluster assignments to enrichment data
  dat_new <- add_clust2enrich(dat_new,clustdat)
  #Make a new list to store the p-values for each parameter from enrichment data
  newlist <- list()
  #Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_sim, style = 3)
  #For loop for each simulation
  for (k in 1:n_sim) {
    #Update progress bar
    setTxtProgressBar(pb, k)
    #dat_new_add is always taken for each simulation to sample values for the NA
    #cells and then determine the respective p-values for the parameters
    dat_new_add <- dat_new
    #Go through every column and row and find cells with NA
    for (i in 1:ncol(dat_new_add)) {
      for (j in 1:nrow(dat_new_add)) {
        #In case there is NA in current cell
        if (is.na(dat_new_add[j,i])) {
          #Take all values from current column where NA is found but disregard NAs
          samp <- na.omit(dat_new_add[,i])
          #Sample randomly for current cell
          dat_new_add[j,i] <- sample(samp,1)
        }
      }
    }
    #Remove "Patient_ID" column from enrichment data, otherwise they would be
    #recognized as categorical variable for summary
    dat_new_add[,"Patient_ID"] <- NULL
    #Split the data so that p-values can be calculated
    #(in cluster = 0 vs. not in cluster = 1)
    enrich_in <- dplyr::filter(dat_new_add, Cluster == clustno)
    enrich_in[,"Cluster"] <- 0
    enrich_out <- dplyr::filter(dat_new_add, Cluster != clustno)
    enrich_out[,"Cluster"] <- 1
    #Now bind them both together again
    enrich_new <- rbind(enrich_in,enrich_out)
    #From package "Arsenal"; Creates and calculates a summary of the data by cluster
    #By default: p-values for continuous variables are calculated by Mann-Whitney test
    #and Hypergeometric test for categorical variable
    #The tableby() function creates a list where all the data is stored
    table_enr <- tableby(Cluster ~., data = enrich_new, numeric.test = numeric, cat.test = categorical)
    #Now this list is transformed to a data frame to print easily the intermediary results
    table_enr <- as.data.frame(summary(table_enr))
    #Remove some reoccuring characters to make it readable
    table_enr[,1] <- gsub("&nbsp;","",as.character(table_enr[,1]))
    #Take out first part of summary, showing the parameters
    col1 <- as.data.frame(table_enr[,1])
    #Take out second part, showing the data for parameters inside cluster
    col2 <- as.data.frame(table_enr[,2])
    #Take out third part, showing the corresponding p-values
    col3 <- as.data.frame(table_enr[,5])
    #Bind the three parts
    table_enr <- cbind(col1,col2,col3)
    #Assign new colnames
    colnames(table_enr) <- c("Parameter"," ","p-value")
    #Modify the data frame so that first column stores the parameter and the second
    #column the corresponding p-value
    table_enr[,2] <- NULL
    table_enr[ table_enr == "" ] <- NA
    table_enr <- table_enr[complete.cases(table_enr),]
    #For each parameter from the data frame
    for (l in 1:nrow(table_enr)) {
      #extract the corresponding p-value and store for k'th run of simulation
      #in list
      #Check if value = "< 0.001", replace with "0.001", otherwise it is not
      #possible to make a boxplot; Else just take the received value
      if (table_enr[l,2] == "< 0.001") {
        newlist[[table_enr[l,1]]][k] <- 0.001
      } else {
        newlist[[table_enr[l,1]]][k] <-  as.numeric(table_enr[l,2])
      }
    }
  }
  #Boxplot for sampling results
  boxplot(newlist, main = paste("Enrichment sampling / Cluster No: ", clustno,
                                "/ n_simulations = ", n_sim), ylab = "p-value")
  newlist
}
