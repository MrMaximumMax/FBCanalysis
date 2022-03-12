#General functions to read out, process and observe time series data
#This function is still not 100% finished... (Polynomial regression + regularization required)
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
