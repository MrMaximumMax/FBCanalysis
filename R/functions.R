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
      "Sample missing values with 5 LOWEST measured", "\t", "(2)", "\n\n")
  usecase <- readline("Select measure for filling up NA values: ")

  patnames <- names(table(raw.data[,"Patient_ID"]))
  n <- length(patnames)

  datalist <- list()

  for (i in 1:n) {

    filter <- patnames[i]
    patdat <- filter(raw.data, Patient_ID %in% filter)

    #if (missing(interpolate) || interpolate == TRUE) {

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

  }
  datalist
}

tibble::tibble()
