#General functions to read out, process and observe time series data

#' FBCanalysis: A package for developing and evaluating biomedical time series
#' data clustering model based on Fluctuation Based Clustering (FBC)
#'
#' This R package aims to offer researchers with fast tools for clustering
#' patient time series data and confirming the distinction using additional
#' metrics such as population parameter enrichment analysis, stability after
#' random data removal, and conventional cluster stability measures. The package
#' attempts to conveniently apply computational methods and capabilities for
#' developing and evaluating unsupervised clustering models with the goal of
#' data-drivenly categorizing asthmatic patients according to their illness dynamics.
#'
#' clustering may be used within the proposed package to identify significant
#' diverse groupings in a patient population, and enrichment analysis is used to
#' examine any possible correlations with clinically relevant characteristics.
#'
#' The R package thus aims to offer researchers with fast tools for clustering
#' patient time series data and confirming the distinction using additional metrics
#' such as population parameter enrichment analysis, stability after random data
#' removal, and conventional cluster stability measures
#'
#' @section Time series data preparation and visualization functions:
#'
#'  \link{patient_list}
#'
#'  \link{patient_ts_plot}
#'
#'  \link{patient_boxplot}
#'
#'  \link{patient_hist}
#'
#' @section Earth Mover's Distance processing functions:
#'
#'  \link{emd_matrix}
#'
#'  \link{emd_heatmap}
#'
#'  \link{max_fluc}
#'
#' @section Clustering to determine heterogeneous groups:
#'
#'  \link{clust_matrix}
#'
#' @section Enrichment analysis functions:
#'
#'  \link{add_enrich}
#'
#'  \link{add_clust2enrich}
#'
#'  \link{add_clust2ts}
#'
#'  \link{enr_obs_clust}
#'
#'  \link{sim_sample_enr}
#'
#' @section Determine cluster stability upon random data removal:
#'
#'  \link{sim_jaccard_cognate}
#'
#'  \link{sim_jaccard_emd}
#'
#'  \link{jaccard_run_cognate}
#'
#'  \link{jaccard_run_emd}
#'
#' @section Cluster stability measure validation:
#'
#'  \link{init_clValid}
#'
#'  \link{clValid_flow}
#'
#' @section Helper function:
#'
#'  \link{znorm}
#'
#' @seealso \href{https://github.com/MrMaximumMax/FBCanalysis}{GitHub Repository}
#'
#' @docType package
#' @name FBCanalysis
NULL
#> NULL
#' Process patient time series data by interpolation options and store data in
#' an object of type list.
#'
#' @param path Path where csv file(s) are stored (only folder, not specific file(s))
#' @param GitHub Set TRUE when csv file comes form GitHub (FALSE by default); only in demo needed
#'
#' @return Object of type list storing patient time series data
#'
#' @import glmnet
#' @import lubridate
#' @import imputeTS
#' @import tibble
#' @import dplyr
#' @import readr
#' @import utils
#'
#' @details Prior to undertaking an analysis using one of the FBC procedures,
#' it is necessary to adequately process and prepare the relevant time series data.
#' The function then creates an interactive flow using the console in R Studio.
#' To begin, the method retrieves all csv files in the provided folder, indicating
#' that it is capable of handling multiple files. The function extracts all csv files
#' from the given directory and merges them into a single raw data frame. The user
#' then indicates which column represents Patient ID and time for adequate processing.
#' The csv files are merged, columns are selected where the Patient ID column will
#' be renamed ”Patient_ID” and the time column will be titled ”Time”. This
#' standardization approach is critical for subsequent features because it enables
#' the easy detection of time series data and the consistent computation and
#' processing of data, for example z-normalization.
#'
#' The user should also indicate in the interactive console the time formate which
#' will be standardized with the help of \link{lubridate}. This is crucial because
#' the technique can now filter the raw data by Patient ID, extract the start and
#' end timestamps for each Patient ID, and then align the data if any records are
#' missing while maintaining the indicated sample frequency.
#'
#' The user may choose between seven approaches: L1 Regularization/Least absolute
#' shrinkage and selection operator (LASSO) Regression, L2 Regularization/Ridge
#' Regression, Elastic Net Regularization, Linear interpolation, Cubic C2 interpolation
#' or, according to recent articles, fill in missing values using the highest or
#' lowest quartile of measurements in the given time series data distribution.
#'
#' The Regression and Regularization techniques generate adequate polynomials for
#' each possible degree n-1 (where n is the total number of data points). Afterwards,
#' cross-validation (from \link{glmnet}) is applied to determine the lambda value
#' for the lowest MSE of the model. Afterwards, the model with polynomial degree
#' for the lowest MSE is chosen and the missing data is interpolated with the
#' regularized model.
#'
#' It is also possible to apply a simple linear interpolation in between missing
#' time series data points. It may be the easiest option to employ straight lines
#' between neighboring points (also see \link{na_interpolation}). Nevertheless,
#' these basic spline polynomials may be notoriously inexact. Cubic spline polynomials
#' mostly provide better results.
#'
#' Another option for the user is to apply the interpolation by using a cubic C spline.
#' It implies that the composite function S must be twice continuously differentiable
#' from all boundaries or subintervals (also see \link{na_interpolation}).
#'
#' Without regression, regularization or interpolation, the user may opt to sample
#' missing values within time series data by randomly choosing a value from the
#' greatest or lowest quartile readings from each patient distribution. The R
#' function then loops over each NA element in the time series data distribution
#' of a patient for the specified parameter and randomly samples a value for the
#' chosen quartile until the data frame is complete.
#'
#' @references Jerome Friedman, Trevor Hastie, Rob Tibshirani, Balasubramanian
#' Narasimhan, Ken- neth Tay, Noah Simon, and Junyang Qian. Package ‘glmnet’.
#' Journal of Statistical Software. 2010a, 33(1), 2021.
#'
#' Hui Zou and Trevor Hastie. Regularization and variable selection via the elastic net.
#' Journal of the royal statistical society: series B (statistical methodology),
#' 67(2):301– 320, 2005.
#'
#' Steffen Moritz and Thomas Bartz-Beielstein. imputets: time series missing value
#' imputation in r. R J., 9(1):207, 2017.
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#'
#' @export
patient_list <- function (path, GitHub) {
  #The GitHub functionality was included for demo so that the user can take the
  #provided csv file from the GitHub repository
  if (missing(GitHub) || GitHub == FALSE) {
  #Prepare the raw data and get all csv. files from the indicated file path
  raw.files <- tibble(filename = list.files(path))
  #List all the files
  raw.file.paths <- raw.files  %>%
    mutate(filepath = paste0(path,"/", filename))
  raw.data <- raw.file.paths %>%
    #Combine them all together into one raw data frame
    rowwise() %>%
    do(., read_csv(file=.$filepath, show_col_types = FALSE))
  #Change object type of raw.data to data.frame
  raw.data <- as.data.frame(raw.data) }
  else {
    raw.data <- read.csv(path) #This would be a GitHub raw path to one csv file
    #only feasible for a single csv file so far
  }
  #First terminal output to present the user the available columns
  cat("\n\n")
  #Show the data frame with the first 5 lines as well as column names
  print(head(raw.data[1:5,]))
  cat("\n\n")
  #Make a small df where each column name has a corresponding number assigned
  #of which the user can decide which column number denotes the Patient_ID and
  #the time
  feedback1 <- colnames(raw.data)
  feedback2 <- 1:length(feedback1)
  feedback <- matrix(nrow = 2, ncol = length(feedback1))
  feedback[1,] <- feedback1
  feedback[2,] <- feedback2
  feedback <- as.data.frame(feedback)
  print(feedback, row.names=F, col.names=F, quote=F)
  cat("\n")
  #Ask which of the assigned column names/numbers denotes the Patient_ID
  question1 <- as.numeric(readline("Which column represents the Patient_ID (please indicate number)?: "))
  #Change the assigned column name for the indicated number to "Patient_ID"
  #for data standardization for later processing
  names(raw.data)[names(raw.data) == feedback[1,question1]] <- 'Patient_ID'
  #Ask which of the assigned column names/numbers denotes the Time
  question2 <- as.numeric(readline("Which column represents the time (please indicate number)?: "))
  #Change the assigned column name for the indicated number to "Time"
  #for data standardization for later processing
  names(raw.data)[names(raw.data) == feedback[1,question2]] <- 'Time'
  #Present the user the possible lubridate format option to standardize time formats
  cat("\n\n", "Day:Month:Year", "\t", "(1)", "\n", "Month:Day:Year", "\t", "(2)", "\n",
      "Year:Month:Day", "\t", "(3)", "\n", "Year:Day:Month", "\t", "(4)", "\n",
      "Day:Month:Year:Hour:Min", "(5)", "\n", "Month:Day:Year:Hour:Min", "(6)", "\n",
      "Year:Month:Day:Hour:Min", "(7)", "\n", "Year:Day:Month:Hour:Min", "(8)", "\n\n")
  #Store indicated format with a number and later apply lubridate function to
  #format and sort the time of each data distribution by parameter
  question3 <- readline("What's the time format (please indicate number)?: ")
  #Present the user the possible sample frequency option on time series data
  cat("\n\n", "Twice daily", "\t", "(1)", "\n", "Daily", "\t\t", "(2)", "\n",
      "Twice a week", "\t", "(3)", "\n", "Weekly", "\t", "(4)", "\n",
      "Twice a month", "\t", "(5)", "\n", "Monthly", "\t", "(6)", "\n",
      "Twice a quarter", "(7)", "\n", "Quarterly", "\t", "(8)", "\n\n")
  #Store indicated format with a number and later apply a check to recognize leaky
  #time series and eventually fill missing times up with NA
  question4 <- readline("What's the sampling frequency (please indicate number)?: ")
  #Present the user the possible options fill up missing data and let him decide
  cat("\n\n", "Sample missing values from top quantile", "\t", "(1)", "\n",
      "Sample missing values with bottom quantile", "\t", "(2)", "\n\n",
      "Interpolate missing values and apply L1 Regularization / Lasso Regression (3)", "\n",
      "Interpolate missing values and apply L2 Regulariztation / Ridge Regression (4)","\n",
      "Interpolate missing values and apply Elastic Net on Regression and Regularization (5)", "\n\n",
      "Interpolate missing values by Linear Spline (6)", "\n",
      "Interpolate missing values by Cubic C2 Spline (7)", "\n\n")
  #Store the option that the user has chosen
  usecase <- readline("Select measure for filling up NA values: ")
  #Extract the individual Patient_IDs from raw.data frame to the filter the
  #raw.data by Patient_ID and process the time series distributions according
  #to the previously specified criteria
  patnames <- names(table(raw.data[,"Patient_ID"]))
  #Preparation to for-loop to go and filter through each Patient_ID
  n <- length(patnames)
  #Create an empty list where each Patient_ID specific processed time series
  #data will be stored
  datalist <- list()
  #In case the user has chosen Interpolation/Regularizatio by Elastic Net
  if (usecase == 5) {
    cat("\n\n")
    #Let user decide which Alpha value to use for Elastic Net
    Alpha_for_elastic_net <<- as.numeric(readline("Please indicate your alpha value for elastic net (between 0 and 1): "))
    #Remark: The Alpha_for_elastic_net will be added to the environment by now
    #because after a few Iterations, there have been issues with glmnet which could
    #not find the assigned alpha value anymore when alpha is not added to the
    #environment; debugging could not resolve this so far
  }
  #Initialize a progress bar
  pb <- txtProgressBar(min = 0, max = n, style = 3)
  #Go through each Patient_ID from the raw.data df
  for (i in 1:n) {
    #Make the filter from the Patient_ID's vector for current Patient_ID
    filter <- patnames[i]
    #Filter for current Patient_ID
    patdat <- dplyr::filter(raw.data, Patient_ID %in% filter)
    #In case user has chosen dmy time format
    if (question3 == 1) {
      patdat$Time <- dmy(patdat$Time)
    }
    #In case user has chosen mdy time format
    if (question3 == 2) {
      patdat$Time <- mdy(patdat$Time)
    }
    #In case user has chosen ymd time format
    if (question3 == 3) {
      patdat$Time <- ymd(patdat$Time)
    }
    #In case user has chosen ydm time format
    if (question3 == 4) {
      patdat$Time <- ydm(patdat$Time)
    }
    #In case user has chosen dmy_hm time format
    if (question3 == 5) {
      patdat$Time <- dmy_hm(patdat$Time)
    }
    #In case user has chosen mdy_hm time format
    if (question3 == 6) {
      patdat$Time <- mdy_hm(patdat$Time)
    }
    #In case user has chosen ymd_hm time format
    if (question3 == 7) {
      patdat$Time <- ymd_hm(patdat$Time)
    }
    #In case user has chosen ydm_hm time format
    if (question3 == 8) {
      patdat$Time <- ydm_hm(patdat$Time)
    }
    #Now order by time after the time column has been formated
    patdat <- arrange(patdat, Time)
    #Find the first date of the current time series
    firstdate <- patdat[1,"Time"]
    #Find the last date of the current time series
    lastdate <- patdat[max(nrow(patdat)),"Time"]
    #Match first and last date of time series with twice daily sampling frequency
    if (question4 == 1) {
      #Create a first time sequence
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "day"))
      #Add hald a day to first and last date (so twice daily is resulting)
      firstdate_2 <- firstdate + 0.5
      lastdate_2 <- lastdate + 0.5
      #Create a second time sequence
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "day"))
      #Combine both sequences
      newdates <- full_join(newdates,newdates_2, by = "Time")
    }
    if (question4 == 2) {
      #Make a daily sequence out of first and last date
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "day"))
    }
    if (question4 == 3) {
      #Make a first weekly sequence out of first and last date
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "week"))
      #Now add 3 days (approx. 1/2 week) to first and last date for biweekly sequence
      firstdate_2 <- firstdate + 3
      lastdate_2 <- lastdate + 3
      #Make a second weekly sequence out of the new dates
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "week"))
      #Combine both sequences to a complete biweekly sequence
      newdates <- full_join(newdates,newdates_2, by = "Time")
    }
    if (question4 == 4) {
      #Make a weekly sequence out of first and last date
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "week"))
    }
    if (question4 == 5) {
      #Make a first monthly sequence out of first and last date
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "month"))
      #Add 15 days (1/2 month) to first and last date
      firstdate_2 <- firstdate + 15
      lastdate_2 <- lastdate + 15
      #Make a second monthly sequence out of these new dates
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "month"))
      #Combine both sequences to a twice monthly sequence
      newdates <- full_join(newdates,newdates_2, by = "Time")
    }
    if (question4 == 6) {
      #Make a monthly sequence out of first and last date
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "month"))
    }
    if (question4 == 7) {
      #Make a quarterly sequence out of first and last date
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "quarter"))
      #Add 45 days (1/2 quarter) to first and last date for second sequence
      firstdate_2 <- firstdate + 45
      lastdate_2 <- lastdate + 45
      #Make a second sequence out of these new date
      newdates_2 <- data.frame(Time = seq(firstdate_2, lastdate_2, by = "quarter"))
      #Combine both sequences to a twice a quarter sequence
      newdates <- full_join(newdates,newdates_2, by = "Time")
    }
    if (question4 == 8) {
      #Make a quarterly seaquence out of first and last date
      newdates <- data.frame(Time = seq(firstdate, lastdate, by = "quarter"))
    }
    #Join the current patient data by newdates and Time to check if time series is
    #incomplete; in case the time series is incomplete, a new line with this
    #time will be added and all the other cells are filled up with NAs that can
    #be filled up/interpolated later
    patdat <- full_join(patdat,newdates, by = "Time")
    #Filter current patient data by complete cases (lines with no NAs) and make a new column
    #to note with FALSE that no interpolation took place
    patdat1 <- patdat[complete.cases(patdat),]
    patdat1[,"Interpolated"] <- FALSE
    #Filter current patient data by incomplete cases (lines with NAs) and make a new column
    #to note with TRUE that an interpolation took place
    patdat2 <- patdat[!complete.cases(patdat),]
    patdat2[,"Interpolated"] <- TRUE
    #Merge both again and order again by time
    patdat <- rbind(patdat1,patdat2)[order(patdat$Time),]
    #Extract the time series data parameters from the current patient data but
    #exclude Patient_ID-, Time, Interpolation-column as they represent no time
    #series dependent data
    parameters <- (names(patdat))
    parameters <- gsub("Patient_ID",NA,parameters)
    parameters <- gsub("Time",NA,parameters)
    parameters <- gsub("Interpolated",NA,parameters)
    parameters <- na.omit(parameters)
    #In case the user opted to sample by highest quantile of data distribution
    if (usecase == 1) {
      #Go through each parameter
      for (j in 1:length(parameters)) {
        #Filter data for current parameter in complete data
        lookup <- patdat1[,parameters[j]]
        #Extract the top quantile of the data
        top_quantile <- quantile(lookup)[3]
        lookup <- subset(lookup, lookup > top_quantile)
        #Go through each NA entry
        for (k in 1:nrow(patdat2)) {
          #Randomly samople from top quantile
          patdat2[k,parameters[j]] <- sample(lookup, 1)
        }
      }
      #Combine both complete and incomplete data again
      patdat <- rbind(patdat2,patdat1)
      #Add the Patient_ID to each cell in Patient_ID column so that it remains
      #assignable
      patdat[,"Patient_ID"] <- filter
      #Sort patient data by time and add to the list as an entry wit list item
      #name is Patient_ID
      datalist[[filter]] <- patdat[order(patdat$Time),]
      #Update the progress bar once a patient distribution is processed
    }
    #Analogous approach to usecase == 1 but with bottom quantile
    if (usecase == 2) {
      for (j in 1:length(parameters)) {
        lookup <- patdat1[,parameters[j]]
        #Take here bottom quantile instead
        bottom_quantile <- quantile(lookup)[2]
        #Here smaller than...
        lookup <- subset(lookup, lookup < bottom_quantile)
        for (k in 1:nrow(patdat2)) {
          patdat2[k,parameters[j]] <- sample(lookup, 1)
        }
      }
      patdat <- rbind(patdat2,patdat1)
      patdat[,"Patient_ID"] <- filter
      datalist[[filter]] <- patdat[order(patdat$Time),]
    }
    #In case the user indicated Polynomial regression/L1 regularization
    #Remark here (as in 162); Usecase 3,4,5 were created in an isolated manner
    #even though that is not the most elegant solution; Here likewise, there has
    #been the issue that after a for-loop-runs where glmnet could not recognize
    #the defined alpha value that has been defined in an object of type numeric
    #before; These issues did not occur when alpha was indicated as a fixed number
    #within the glmnet(...) function
    if (usecase == 3) {
      #Order the current patient data by time
      patdat <- patdat[order(patdat$Time),]
      #Add a new line which indicates that the time is now an equally spaced
      #sequence to train the models via glmnet
      patdat[,"Seq"] <- 1:nrow(patdat)
      #Filter for existing data = Training data
      patdat1 <- dplyr::filter(patdat, Interpolated == FALSE)
      #Polynomial degress to try out are number of data points - 1
      pol_degrees <- length(complete.cases(patdat))-1
      #Loop through each parameters individually and find out best polynomial
      #degree and lambda for regularized model
      for (j in 1:length(parameters)) {
        #Define the current parameter
        current_par <- parameters[j]
        #Make a matrix where the line number denotes the polynomial degree and
        #the first column stores the model's lambda.min and the second column
        #the correspoinding Mean Squared Erro
        regularization_mat <- matrix(0, nrow = pol_degrees, ncol = 2)
        #Make a glmnet regularized model for each possible polynomial degree
        for (k in 1:pol_degrees) {
          #Define current polynomial degree
          PolynomialDegree <- k
          #Quiet warnings here: There will always be NA entries added automatically
          #in data matrix which would lead to many warning messages when looping
          #for each possible polynomial degree
          options(warn=-1)
          X <- matrix(patdat1[,"Seq"],ncol = 1);
          for (l in (c(2:PolynomialDegree))){
            X <- cbind(X,patdat1[,"Seq"]*X[,(l-1)])
          }
          #Disable warnings
          options(warn=0)
          #Define max of iterations for cross validation
          MaxIt<-1e+07;
          #Make glmnet for regularized model and current polynomial degree
          Y <- glmnet(X, patdat1[,current_par],family = "gaussian",alpha = 1, standardize = TRUE,intercept = TRUE, maxit = MaxIt)
          #Apply corss validation to determine best fitting lamda for min. MSE
          cvfit <- cv.glmnet(X,patdat1[,current_par],family = "gaussian",alpha = 1, standardize = TRUE, type.measure = "mse", nfolds = 10);#nfolds = n corresponds to leave-one-out cross validation. A popular value is 10.
          #Add both lamda.min and correspionding MSE to matrix
          regularization_mat[k,1] <- cvfit$lambda.min
          regularization_mat[k,2] <- min(cvfit$cvm)
        }
        #Search lowest MSE from matrix
        best_pol_deg <- which.min(regularization_mat[,2])
        #As in line 368 quiet the warning messages
        options(warn=-1)
        X <- matrix(patdat1[,"Seq"],ncol = 1);
        for (m in (c(2:best_pol_deg))){
          X <- cbind(X,patdat1[,"Seq"]*X[,(m-1)])
        }
        #Do not quiet the warning messages anymore
        options(warn=0)
        #Apply glmnet to the best fitting polynomial degree
        Y <- glmnet(X,patdat1[,current_par],family = "gaussian",alpha = 1, standardize = TRUE,intercept = TRUE, maxit = MaxIt)
        #Get the coefficients from this model with lamda.min
        C <- coef(Y, s = cvfit$lambda.min,exact = TRUE,x=X, y=patdat1[,current_par], maxit = MaxIt)
        CC <- matrix(data = c(C[2:(best_pol_deg+1),1]),ncol = 1)
        YY <- (X * as.vector(CC)) + C[1,1]
        #find the indexes within the data distribution where NA values are found
        na_in_dat <- sum(is.na(patdat[,current_par]))
        na_index <- which(is.na(patdat[,current_par]))
        #Similar as in line 368
        options(warn=-1)
        for (o in 1:na_in_dat) {
          #Go to current NA values and let predict the corresponding value at
          #this position
          seq_miss = na_index[o];
          tNA <- matrix(seq_miss, ncol = 1)
          for (p in (c(2:best_pol_deg))){
            tNA<-cbind(tNA,seq_miss*tNA[,(p-1)])
          }
          patdat[,current_par][seq_miss] <- predict(Y, newx = tNA, exact = TRUE, s = cvfit$lambda.min, x=X, y=patdat1[,current_par], maxit = MaxIt)
        }
        #Quiet warning again (see line 368)
        options(warn=-1)
        #Remove the sequence column (not relevant for later analyses anymore)
        patdat[,"Seq"] <- NULL
        #Add the Patient_ID to each cell in Patient_ID column so that it remains
        #assignable
        patdat[,"Patient_ID"] <- filter
        #Add the currently processed patient data to the list
        datalist[[filter]] <- patdat
        #Update the progress bar
      }
    }
    #In case the user indicated Polynomial regression/L1 regularization
    #Similar remarks as in line 338; The only difference here is that
    #now glmnet uses a value of 1; Besides this, the approach is similar
    if (usecase == 4) {
      patdat <- patdat[order(patdat$Time),]
      patdat[,"Seq"] <- 1:nrow(patdat)
      patdat1 <- dplyr::filter(patdat, Interpolated == FALSE)
      pol_degrees <- length(complete.cases(patdat))-1
      for (j in 1:length(parameters)) {
        current_par <- parameters[j]
        regularization_mat <- matrix(0, nrow = pol_degrees, ncol = 2)
        for (k in 1:pol_degrees) {
          PolynomialDegree <- k
          options(warn=-1)
          X <- matrix(patdat1[,"Seq"],ncol = 1);
          for (l in (c(2:PolynomialDegree))){
            X <- cbind(X,patdat1[,"Seq"]*X[,(l-1)])
          }
          options(warn=0)
          MaxIt<-1e+07;
          Y <- glmnet(X, patdat1[,current_par],family = "gaussian",alpha = 0, standardize = TRUE,intercept = TRUE, maxit = MaxIt)
          cvfit <- cv.glmnet(X,patdat1[,current_par],family = "gaussian",alpha = 0, standardize = TRUE, type.measure = "mse", nfolds = 10);#nfolds = n corresponds to leave-one-out cross validation. A popular value is 10.
          regularization_mat[k,1] <- cvfit$lambda.min
          regularization_mat[k,2] <- min(cvfit$cvm)
        }
        best_pol_deg <- which.min(regularization_mat[,2])
        options(warn=-1)
        X <- matrix(patdat1[,"Seq"],ncol = 1);
        for (m in (c(2:best_pol_deg))){
          X <- cbind(X,patdat1[,"Seq"]*X[,(m-1)])
        }
        options(warn=0)
        Y <- glmnet(X,patdat1[,current_par],family = "gaussian",alpha = 0, standardize = TRUE,intercept = TRUE, maxit = MaxIt)
        C <- coef(Y, s = cvfit$lambda.min,exact = TRUE,x=X, y=patdat1[,current_par], maxit = MaxIt)
        CC <- matrix(data = c(C[2:(best_pol_deg+1),1]),ncol = 1)
        YY <- (X * as.vector(CC)) + C[1,1]
        na_in_dat <- sum(is.na(patdat[,current_par]))
        na_index <- which(is.na(patdat[,current_par]))
        options(warn=-1)
        for (o in 1:na_in_dat) {
          seq_miss = na_index[o];
          tNA <- matrix(seq_miss, ncol = 1)
          for (p in (c(2:best_pol_deg))){
            tNA<-cbind(tNA,seq_miss*tNA[,(p-1)])
          }
          patdat[,current_par][seq_miss] <- predict(Y, newx = tNA, exact = TRUE, s = cvfit$lambda.min, x=X, y=patdat1[,current_par], maxit = MaxIt)
        }
        options(warn=-1)
        patdat[,"Seq"] <- NULL
        patdat[,"Patient_ID"] <- filter
        datalist[[filter]] <- patdat
      }
    }
    #In case the user indicated Polynomial regression/L1 regularization
    #Similar remarks as in line 338; The only difference here is that
    #now glmnet uses the indicated alpha for elastic net that has been added to
    #the environment (see line 338); Besides this, the approach is similar
    if (usecase == 5) {
      patdat <- patdat[order(patdat$Time),]
      patdat[,"Seq"] <- 1:nrow(patdat)
      patdat1 <- dplyr::filter(patdat, Interpolated == FALSE)
      pol_degrees <- length(complete.cases(patdat))-1
      for (j in 1:length(parameters)) {
        current_par <- parameters[j]
        regularization_mat <- matrix(0, nrow = pol_degrees, ncol = 2)
        for (k in 1:pol_degrees) {
          PolynomialDegree <- k
          options(warn=-1)
          X <- matrix(patdat1[,"Seq"],ncol = 1);
          for (l in (c(2:PolynomialDegree))){
            X <- cbind(X,patdat1[,"Seq"]*X[,(l-1)])
          }
          options(warn=0)
          MaxIt<-1e+07;
          Y <- glmnet(X, patdat1[,current_par],family = "gaussian",alpha = Alpha_for_elastic_net, standardize = TRUE,intercept = TRUE, maxit = MaxIt)
          cvfit <- cv.glmnet(X,patdat1[,current_par],family = "gaussian",alpha = Alpha_for_elastic_net, standardize = TRUE, type.measure = "mse", nfolds = 10);#nfolds = n corresponds to leave-one-out cross validation. A "popular" value is 10.
          regularization_mat[k,1] <- cvfit$lambda.min
          regularization_mat[k,2] <- min(cvfit$cvm)
        }
        best_pol_deg <- which.min(regularization_mat[,2])
        options(warn=-1)
        X <- matrix(patdat1[,"Seq"],ncol = 1);
        for (m in (c(2:best_pol_deg))){
          X <- cbind(X,patdat1[,"Seq"]*X[,(m-1)])
        }
        options(warn=0)
        Y <- glmnet(X,patdat1[,current_par],family = "gaussian",alpha = Alpha_for_elastic_net, standardize = TRUE,intercept = TRUE, maxit = MaxIt)
        C <- coef(Y, s = cvfit$lambda.min,exact = TRUE,x=X, y=patdat1[,current_par], maxit = MaxIt)
        CC <- matrix(data = c(C[2:(best_pol_deg+1),1]),ncol = 1)
        YY <- (X * as.vector(CC)) + C[1,1]
        na_in_dat <- sum(is.na(patdat[,current_par]))
        na_index <- which(is.na(patdat[,current_par]))
        options(warn=-1)
        for (o in 1:na_in_dat) {
          seq_miss = na_index[o];
          tNA <- matrix(seq_miss, ncol = 1)
          for (p in (c(2:best_pol_deg))){
            tNA<-cbind(tNA,seq_miss*tNA[,(p-1)])
          }
          patdat[,current_par][seq_miss] <- predict(Y, newx = tNA, exact = TRUE, s = cvfit$lambda.min, x=X, y=patdat1[,current_par], maxit = MaxIt)
        }
        options(warn=-1)
        patdat[,"Seq"] <- NULL
        patdat[,"Patient_ID"] <- filter
        datalist[[filter]] <- patdat
        setTxtProgressBar(pb, i)
      }
    }
    #In case the user opted linear interpolation
    if (usecase == 6) {
      #Order current patient data by time
      patdat <- patdat[order(patdat$Time),]
      #Perform linear interpolation on NA values
      for (j in 1:length(parameters)) {
      patdat[,parameters[j]] <- na_interpolation(patdat[,parameters[j]], option = "linear")
      }
      #Add the Patient_ID to each cell in Patient_ID column so that it remains
      #assignable
      patdat[,"Patient_ID"] <- filter
      #Add to list
      datalist[[filter]] <- patdat
    }
    #In case the user opted cubic spline interpolation
    if (usecase == 7) {
      #Order current patient data by time
      patdat <- patdat[order(patdat$Time),]
      #Perform cubic c spline interpolation on NA values
      for (j in 1:length(parameters)) {
        patdat[,parameters[j]] <- na_interpolation(patdat[,parameters[j]], option = "spline")
      }
      #Add the Patient_ID to each cell in Patient_ID column so that it remains
      #assignable
      patdat[,"Patient_ID"] <- filter
      #Add to list
      datalist[[filter]] <- patdat
    }
    setTxtProgressBar(pb, i)
  }
  datalist
}

#' Visualize patient time series data from a preprocessed in a time series plot
#' for an indicated parameter, either as normalized or non-normalized.
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list})
#' @param Patient_ID Patient_ID referring to a list element (also see function: \link{patient_list})
#' @param parameter Parameter of interest in list element
#' @param normalized TRUE/FALSE if z-normalized (TRUE by default)
#'
#' @return Visualized patient time series data in a time series plot for indicated parameter
#'
#' @import graphics
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' patient_ts_plot(list,"testpat_1","PEF")
#'
#' @export
patient_ts_plot <- function(plist, Patient_ID, parameter, normalized) {

  #Take specific data frame out of list according to specified Patient_ID
  patient_df <- as.data.frame(plist[[Patient_ID]])

  #In case "normalized" is not specified or "normalized = TRUE"
  if (missing(normalized) || normalized == TRUE) {

    #Plot time on x-axis and z-normalize specified parameter from data frame on y-axis
    graphics::plot(patient_df$Time,znorm(patient_df[, which(names(patient_df) %in% parameter)]),
         xlab = "Time", ylab = parameter, main = Patient_ID, col = "blue")

    #In case "normalized = FALSE"
  } else {
    #Plot time on x-axis and non-normalized specified parameter from data frame on y-axis
    plot(patient_df$Time,patient_df[, which(names(patient_df) %in% parameter)],
         xlab = "Time", ylab = parameter, main = Patient_ID, col = "blue")
  }
}

#' Visualize patient(s) time series from a preprocessed list with in a boxplot,
#' either as normalized or non-normalized for an indicated parameter.
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list})
#' @param patients Patient_ID(s) referring to (a) list element; can be single ID or multiple IDs (also see function: \link{patient_list})
#' @param parameter Parameter of interest in list element(s)
#' @param normalized TRUE/FALSE if z-normalized (TRUE by default)
#'
#' @return Visualized patient(s) time series data in a boxplot for indicated parameter
#'
#' @import dplyr
#' @import graphics
#'
#' @examples
#' list <- patient_list(
#' "https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' patient_boxplot(list,c("ID_2","testpat_1","testpat_2","a301"), "FEV1")
#'
#' @export
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
    dat <- patient_df %>% dplyr::filter(Patient_ID %in% patients)

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

#' Visualize patient time series data from a preprocessed in a histogram
#' for indicated parameter either normalized or non-normalized.
#'
#' @param plist List storing patient time series data (also see function: \link{patient_list})
#' @param Patient_ID Patient_ID referring to a list element (also see function: \link{patient_list})
#' @param parameter Parameter of interest in list element
#' @param normalized TRUE/FALSE if z-normalized (TRUE by default)
#'
#' @return Visualized patient time series data in a histogram for indicated parameter
#'
#' @import graphics
#'
#' @examples
#' list <- patient_list("
#' https://raw.githubusercontent.com/MrMaximumMax/FBCanalysis/master/demo/phys/data.csv",
#' GitHub = TRUE)
#' #Sampling frequency is supposed to be daily
#' patient_hist(list,"testpat_1","PEF")
#'
#' @export
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
