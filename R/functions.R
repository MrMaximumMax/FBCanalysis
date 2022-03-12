#General functions to read out, process and observe time series data
patient_list <- function (path) {

  raw.files <- tibble(filename = list.files(path))
  raw.file.paths <- raw.files  %>%
    mutate(filepath = paste0(path,"/", filename))
  raw.data <- raw.file.paths %>%
    # 'do' the function for each row in turn
    rowwise() %>%
    do(., read_csv(file=.$filepath, show_col_types = FALSE))

  raw.data <- as.data.frame(raw.data)

  cat("\n\n")
  print(head(raw.data[1:5,]))
  cat("\n\n")

  feedback1 <- colnames(raw.data)
  feedback2 <- 1:length(feedback1)
  feedback <- matrix(nrow = 2, ncol = length(feedback1))
  feedback[1,] <- feedback1
  feedback[2,] <- feedback2
  feedback <- as.data.frame(feedback)

  print(feedback, row.names=F, col.names=F, quote=F)
  cat("\n")

  question1 <- as.numeric(readline("Which column represents the Patient_ID (please indicate number)?: "))
  names(raw.data)[names(raw.data) == feedback[1,question1]] <- 'Patient_ID'
  question2 <- as.numeric(readline("Which column represents the time (please indicate number)?: "))
  names(raw.data)[names(raw.data) == feedback[1,question2]] <- 'Time'

  cat("\n\n", "Day:Month:Year", "\t", "(1)", "\n", "Month:Day:Year", "\t", "(2)", "\n",
      "Year:Month:Day", "\t", "(3)", "\n", "Year:Day:Month", "\t", "(4)", "\n",
      "Day:Month:Year:Hour:Min", "(5)", "\n", "Month:Day:Year:Hour:Min", "(6)", "\n",
      "Year:Month:Day:Hour:Min", "(7)", "\n", "Year:Day:Month:Hour:Min", "(8)", "\n\n")

  question3 <- readline("What's the time format (please indicate number)?: ")

  cat("\n\n", "Twice daily", "\t", "(1)", "\n", "Daily", "\t\t", "(2)", "\n",
      "Twice a week", "\t", "(3)", "\n", "Weekly", "\t", "(4)", "\n",
      "Twice a month", "\t", "(5)", "\n", "Monthly", "\t", "(6)", "\n",
      "Twice a quarter", "(7)", "\n", "Quarterly", "\t", "(8)", "\n\n")
  question4 <- readline("What's the sampling frequency (please indicate number)?: ")


  cat("\n\n", "Sample missing values with 5 HIGHEST measured", "\t", "(1)", "\n",
      "Sample missing values with 5 LOWEST measured", "\t", "(2)", "\n\n",
      "Interpolate missing values and apply L1 Regularization / Lasso Regression (3)", "\n",
      "Interpolate missing values and apply L2 Regulariztation / Ridge Regression (4)","\n\n")
  usecase <- readline("Select measure for filling up NA values: ")

  patnames <- names(table(raw.data[,"Patient_ID"]))
  n <- length(patnames)

  datalist <- list()

  for (i in 1:n) {

    filter <- patnames[i]
    patdat <- filter(raw.data, Patient_ID %in% filter)

    if (question3 == 1) {
      patdat$Time <- dmy(patdat$Time)
    }
    if (question3 == 2) {
      patdat$Time <- mdy(patdat$Time)
    }
    if (question3 == 3) {
      patdat$Time <- ymd(patdat$Time)
    }
    if (question3 == 4) {
      patdat$Time <- ydm(patdat$Time)
    }
    if (question3 == 5) {
      patdat$Time <- dmy_hm(patdat$Time)
    }
    if (question3 == 6) {
      patdat$Time <- mdy_hm(patdat$Time)
    }
    if (question3 == 7) {
      patdat$Time <- ymd_hm(patdat$Time)
    }
    if (question3 == 8) {
      patdat$Time <- ydm_hm(patdat$Time)
    }

    patdat <- arrange(patdat, Time)

    firstdate <- patdat[1,"Time"]
    lastdate <- patdat[max(nrow(patdat)),"Time"]

    if (question4 == 1) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "day"))
      firstdate_2 <- firstdate + 0.5
      lastdate_2 <- lastdate + 0.5
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "day"))
      newdates <- full_join(newdates,newdates_2)
    }
    if (question4 == 2) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "day"))
    }
    if (question4 == 3) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "week"))
      firstdate_2 <- firstdate + 3
      lastdate_2 <- lastdate + 3
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "week"))
      newdates <- full_join(newdates,newdates_2)
    }
    if (question4 == 4) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "week"))
    }
    if (question4 == 5) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "month"))
      firstdate_2 <- firstdate + 15
      lastdate_2 <- lastdate + 15
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "month"))
      newdates <- full_join(newdates,newdates_2)
    }
    if (question4 == 6) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "month"))
    }
    if (question4 == 7) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "quarter"))
      firstdate_2 <- firstdate + 45
      lastdate_2 <- lastdate + 45
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "quarter"))
      newdates <- full_join(newdates,newdates_2)
    }
    if (question4 == 8) {
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "quarter"))
    }

    patdat <- full_join(patdat,newdates)
    patdat[,"Patient_ID"] <- filter

    patdat1 <- patdat[complete.cases(patdat),]
    patdat1[,"Interpolated"] <- FALSE

    patdat2 <- patdat[!complete.cases(patdat),]
    patdat2[,"Interpolated"] <- TRUE

    patdat <- rbind(patdat1,patdat2)[order(patdat$Time),]

    case <- filter(patdat, Interpolated == TRUE)
    parameters <- (names(patdat))
    parameters <- gsub("Patient_ID",NA,parameters)
    parameters <- gsub("Time",NA,parameters)
    parameters <- gsub("Interpolated",NA,parameters)
    parameters <- na.omit(parameters)

    if (usecase == 1) {

      for (i in 1:length(parameters)) {

        lookup <- patdat[,parameters[i]]
        lookup <- na.omit(lookup)
        lookup <- tail(sort(lookup),5)

        for (j in 1:nrow(case)) {

          case[j,parameters[i]] <- sample(lookup, 1)

        }
      }
      rest <- filter(patdat, Interpolated == FALSE)
      patdat <- rbind(case,rest)

      datalist[[filter]] <- patdat[order(patdat$Time),]
    }

    if (usecase == 2) {

      for (i in 1:length(parameters)) {

        lookup <- patdat[,parameters[i]]
        lookup <- na.omit(lookup)
        lookup <- sort(lookup)[1:5]

        for (j in 1:nrow(case)) {

          case[j,parameters[i]] <- sample(lookup, 1)

        }
      }
      rest <- filter(patdat, Interpolated == FALSE)
      patdat <- rbind(case,rest)

      datalist[[filter]] <- patdat[order(patdat$Time),]
    }

    if (usecase == 3 || usecase == 4) {

      patdat <- patdat[order(patdat$Time),]
      sequence <- 1:nrow(patdat)
      patdat$Seq <- sequence

      patdat <- na_interpolation(patdat, option = 'spline')

      for (j in 1:length(parameters)) {

        y <- patdat[,parameters[j]]
        x <- data.matrix(patdat[,c('Time','Seq')])

        if (usecase == 3) {
          alph <- 1
        }
        if (usecase == 4) {
          alph <- 0
        }

        model <- glmnet(x, y, alpha = alph)
        cv_model <- cv.glmnet(x, y, alpha = alph)

        best_lambda <- cv_model$lambda.min
        best_model <- glmnet(x, y, alpha = alph, lambda = best_lambda)

        y_predicted <- as.data.frame(predict(model, s = best_lambda, newx = x))

        y_predicted <- tibble::rownames_to_column(y_predicted,"Seq")
        y_predicted$Seq <- as.numeric(y_predicted$Seq)
        y_predicted <- y_predicted[order(y_predicted$Seq),]
        names(y_predicted)[names(y_predicted) == "s1"] <- parameters[j]

        seq_interpol <- filter(patdat, Interpolated == TRUE)
        seq_interpol <- seq_interpol$Seq
        y_keep <- filter(y_predicted, Seq %in% seq_interpol)

        y_old <- filter(patdat, Interpolated == FALSE)
        y_old <- y_old[,names(y_old) %in% c(parameters[j],"Seq")]

        y_new <- rbind(y_keep,y_old)
        y_new <- y_new[order(y_new$Seq),]

        patdat[,parameters[j]] <- y_new[,2]
        patdat <- patdat[order(patdat$Seq),]

      }
      patdat$Seq <- NULL
      datalist[[filter]] <- patdat
    }
  }
  datalist
}
patient_ts_plot <- function(plist, Patient_ID, parameter, normalized) {

  #Take specific data frame out of list according to specified Patient_ID
  patient_df <- as.data.frame(plist[[Patient_ID]])

  #In case "normalized" is not specified or "normalized = TRUE"
  if (missing(normalized) || normalized == TRUE) {

    #Plot time on x-axis and z-normalize specified parameter from data frame on y-axis
    plot(patient_df$Time,znorm(patient_df[, which(names(patient_df) %in% parameter)]),
         xlab = "Time", ylab = parameter, main = Patient_ID, col = "blue")

    #In case "normalized = FALSE"
  } else {
    #Plot time on x-axis and non-normalized specified parameter from data frame on y-axis
    plot(patient_df$Time,patient_df[, which(names(patient_df) %in% parameter)],
         xlab = "Time", ylab = parameter, main = Patient_ID, col = "blue")
  }
}
patient_boxplot <- function(plist, patients, parameter, normalized) {

  #In case boxplot for only one patient
  if (length(patients) == 1) {
    #Take list entry from specified patient and transform to data frame
    patient_df <- plist[[patients]]
    #Filter patient data for specified parameter
    dat <- patient_df[, which(names(patient_df) %in% parameter)]

    #In case "normalized" is not specified or "normalized = TRUE"
    if (missing(normalized) || normalized == TRUE) {
      #Use written function to z-normalize data
      dat <- znorm(dat)
    }
    #Apply boxplot on z-normalized data
    boxplot(dat, main = patients,ylab = parameter)

    #In case boxplot for more than one patient
  } else {
    #Transform whole list to one data frame
    patient_df <- do.call("rbind", plist)
    #Filter data frame for specified Patient_ID's
    dat <- patient_df %>% filter(Patient_ID %in% patients)

    #In case "normalized" not specified or "normalized = TRUE"
    if (missing(normalized) || normalized == TRUE) {
      #Filter data frame for specified parameter and z-normalize
      #Apply boxplot of z-normalized data against Patient_ID's
      boxplot(znorm(dat[, which(names(dat) %in% parameter)]) ~ dat$Patient_ID,
              main = "Combined boxplots", ylab = parameter, xlab = NULL)

      #In case "normalized = FALSE"
    } else {
      #Apply boxplot of non-normalized data against Patient_ID's
      boxplot(dat[, which(names(dat) %in% parameter)] ~ dat$Patient_ID,
              main = "Combined boxplots", ylab = parameter, xlab = NULL)
    }
  }
}
patient_hist <- function(plist, Patient_ID, parameter, normalized) {

  #Take data frame out of list for specified Patient_ID
  df <- plist[[Patient_ID]]
  #Filter data for specified parameter
  dat <- df[,parameter]

  #In case "normalized" not specified or "normalized = TRUE"
  if (missing(normalized) || normalized == TRUE) {
    #z-normalize the data for specified parameter
    dat <- znorm(dat)
  }
  #Apply histogram on data
  hist(dat, xlab = parameter, main = Patient_ID)
}

#Generate Earth Mover's Distance data
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

#Functions for Jaccard/random data removal
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
    new_df <- filter(df, Patient_ID %in% patnames[i])
    #Add the current new data frame to the list
    newlist[[patnames[i]]] <- new_df
  }
  newlist
}
sim_jaccard_cognate <- function(plist, parameter, removal, n_simu, method, n_clust, Iter) {
  #Simulate random data removal and Jaccard index determination by Cognate Cluster Approach

  #In case, user did not specify maximum for EMD calcilations, the default value is
  #set to a high number (5,000)
  if (missing(Iter)) {
    Iter = 5000
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
      overlap1 <- filter(combined, all_dat == m)
      #Make a vector to store all calculated Jaccard indices; Highest Jaccard
      #index then represents the cognate cluster
      jaccard_vector <- vector()
      #For-loop to filter for each data removal cluster
      for (n in 1:n_clust) {
        #Filter for data removal cluster
        overlap2 <- filter(combined, data_removal == n)
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
sim_jaccard_emd <- function(plist, parameter, removal, n_simu, method, n_clust, Iter) {

  #Simulate random data removal and Jaccard index determination by EMD Approach
  #In case, user did not specify maximum for EMD calcilations, the default value is
  #set to a high number (5,000)
  if (missing(Iter)) {
    Iter = 5000
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
    patnames <- filter(compare, dat_complete$Cls == t)[,1]
    #Make data frame out of time series data list (to be used later to extract distributions)
    distribution <- do.call(rbind,plist)
    #List to store normalized distributions of each Patient_ID for gold standard
    dat <- list()
    #For each Patient_ID in time series data list
    for (u in 1:length(patnames)) {
      #Filter for current Patient_ID to receive current distribution
      distr_a <- filter(distribution, Patient_ID %in% patnames[u])
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
      overlap1 <- filter(combined, all_dat == m)
      #Filter for cluster after random data removal
      overlap2 <- filter(combined, data_removal == m)
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
    jaccard_list[[as.character(range[i])]] <- sim_jaccard_cognate(testlist, parameter, range[i], n_simu = n_simu, method = method, n_clust = n_clust)
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
    jaccard_list[[as.character(range[i])]] <- sim_jaccard_emd(testlist, parameter, range[i], n_simu = n_simu, method = method, n_clust = n_clust)
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

#Functions for Cluster validation measure analysis
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
      #Apply RankAggreg for Cross Entropy search on suoper list with Spearman's footrule
      if(require("RankAggreg")) {
        CEWS <- RankAggreg(x=res$ranks, k=5, weights=res$weights, seed=123, verbose=FALSE)
        cat("\n\n")
        #Print ideal clustering combination from super list
        print(CEWS)
        cat("\n\n")
      }
      #Fifth question
      new_question <- as.character(readline("Do you want to take best result (y/n)?: "))
      #If yes
      if(new_question == "y") {
        #Take result from CEWS list
        split <- CEWS$top.list[1]
        #Split the character by "-"
        split <- strsplit(split, "-")
        #Add both clustering method string and number to output list
        output_list[["method"]] <- split[[1]][1]
        output_list[["clust_num"]] <- as.numeric(split[[1]][2])
      } else {
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
        new_question <- as.integer(readline("Which of your analyzed methods should be used instead (please indicate number)?: "))
        #Add clustering method to output list
        output_list[["method"]] <- par$cl_methods[new_question]
        #Choose number of cluster
        next_question <- as.integer(readline("How many clusters should be applied?: "))
        #Add number of clusters to output list
        output_list[["clust_num"]] <- as.numeric(next_question)
      }
      #Now generate EMD matrix clustering data in standardized format to use
      #for further functions regarding random data removal and enrichment
      #analysis
      final_output <- clust_matrix(matrix, method = output_list$method, nclust = output_list$clust_num)
    }
  } else {
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
    new_question <- as.integer(readline("Which of your analyzed methods should be ussed instead (please indicate number)?: "))
    #Add clustering method to output list
    output_list[["method"]] <- par$cl_methods[new_question]
    #Choose number of cluster
    next_question <- as.integer(readline("How many clusters should be applied?: "))
    #Add number of clusters to output list
    output_list[["clust_num"]] <- as.numeric(next_question)
    final_output <- clust_matrix(matrix, method = output_list$method, nclust = output_list$clust_num)
  }
}
