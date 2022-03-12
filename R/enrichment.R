#General functionalities to cluster data and add cluster assignments to ts- and enrichment data
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

  } else if (method == "agnes") {

    #Apply agnes clustering on EMD data
    clust <- AgglomerativeNestingClustering(matrix, ClusterNo = nclust)
    #In case user has not specified "plotclust" or "plotclust = TRUE"
    if (missing(plotclust) || plotclust == TRUE) {
      #Visualize two-dimensionally result of clara clustering
      plot(agnes(matrix, nclust))
    }

    #Add an object of class "integer" to list, recording which Patient_ID belongs
    #to which cluster
    newlist[["Cls"]] <- clust$Cls
    #Record number of clusters
    newlist[["nclust"]] <- nclust
    #Add character to list to record which clustering method was used
    newlist[["Method"]] <- "Agnes"
    newlist

  } else {
    stop("Incorrect input of clustering method")
  }
}
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
  dat_new1 <- filter(dat, !Patient_ID %in% throw)

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
    compare <- filter(new, Cluster %in% i)
    #Remove column "Cluster"
    compare[,"Cluster"] <- NULL
    #Make a vector, only containing the Patient_ID's for current cluster
    compare <- as.character(compare[,1])
    #Filter time series data frame for Patient_ID's from current cluster
    new2 <- filter(ts_df, Patient_ID %in% compare)
    #Add a new column "Cluster" to time series data frame storing cluster number
    new2[,"Cluster"] <- i
    #Add current cluster time series data to list
    datalist[[i]] <- new2
  }
  #Make one time series data frame out of all list elements
  assigned <- do.call(rbind,datalist)
}
enr_obs_clust <- function(ts.dat, enrich, clustno) {

  #Methods only applicable on pre-processed enrichment and time series data
  #from functions: add_clust2enrich & add_clust2ts

  #Remove "Patient_ID" column from enrichment data, otherwise they would be
  #recognized as categorial variable for summary
  enrich[,"Patient_ID"] <- NULL
  #Split the data so that p-values can be calculated
  #(in cluster = 0 vs. not in cluster = 1)
  enrich_in <- filter(enrich, Cluster == clustno)
  enrich_in[,"Cluster"] <- 0
  enrich_out <- filter(enrich, Cluster != clustno)
  enrich_out[,"Cluster"] <- 1
  #Now bind them both together again
  enrich_new <- rbind(enrich_in,enrich_out)
  #From package "Arsenal"; Creates and calculates a summary of the data by cluster
  #By default: p-values for continuous variables are calculated by Mann-Whitney test
  #and Hypergeometric test for categorial variable
  #The tableby() function creates a list where all the data is stored
  table_enr <- tableby(Cluster ~., data = enrich_new)
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
  n_pat <- filter(enrich, Cluster %in% clustno)
  #Length = Number of patients
  n_pat <- length(n_pat[,1])

  #Start an analogous procedure for time series data
  #Remove both Patient_ID and Time so that only time series data, information
  #on interpolation and cluster assignment remain
  ts.dat[,"Patient_ID"] <- NULL
  ts.dat[,"Time"] <- NULL
  #Split the data so that p-values can be calculated
  #(in cluster = 0 vs. not in cluster = 1)
  ts.dat_in <- filter(ts.dat, Cluster == clustno)
  ts.dat_in[,"Cluster"] <- 0
  ts.dat_out <- filter(ts.dat, Cluster != clustno)
  ts.dat_out[,"Cluster"] <- 1
  #Now bind them both together again
  ts.dat_new <- rbind(ts.dat_in,ts.dat_out)
  table_ts <- tableby(Cluster ~., data = ts.dat_new)
  #Similar procedure as above (line 978)
  table_ts <- as.data.frame(summary(table_ts))
  table_ts[,1] <- gsub("&nbsp;","",as.character(table_ts[,1]))
  cola <- as.data.frame(table_ts[,1])
  colb <- as.data.frame(table_ts[,2])
  colc <- as.data.frame(table_ts[,5])
  table_ts <- cbind(cola,colb,colc)
  colnames(table_ts) <- c("Parameter", " ", "p-value")

  #for summary: number of measurements; number of interpolated data
  #Filter time series data for current cluster and determine length (= number of measurements)
  n_mes <- filter(ts.dat, Cluster %in% clustno)
  n_val <- length(n_mes[,1])
  #Filter and determine length for both interpolated and non-interpolated data
  n_real <- length(filter(n_mes, Interpolated == FALSE)[,1])
  n_inter <- length(filter(n_mes, Interpolated == TRUE)[,1])
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
sim_sample_enr <- function(plist, path, clustdat, clustno, n_sim) {

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
  dat_new1 <- filter(dat, !Patient_ID %in% throw)

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
    #recognized as categorial variable for summary
    dat_new_add[,"Patient_ID"] <- NULL
    #Split the data so that p-values can be calculated
    #(in cluster = 0 vs. not in cluster = 1)
    enrich_in <- filter(dat_new_add, Cluster == clustno)
    enrich_in[,"Cluster"] <- 0
    enrich_out <- filter(dat_new_add, Cluster != clustno)
    enrich_out[,"Cluster"] <- 1
    #Now bind them both together again
    enrich_new <- rbind(enrich_in,enrich_out)
    #From package "Arsenal"; Creates and calculates a summary of the data by cluster
    #By default: p-values for continuous variables are calculated by Mann-Whitney test
    #and Hypergeometric test for categorial variable
    #The tableby() function creates a list where all the data is stored
    table_enr <- tableby(Cluster ~., data = enrich_new)
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
