#' @title Format Chinook salmon data from JTRAP database to run through hierarchical Bayesian models.
#'
#' @description This function reconfigures Chinook Salmon JTRAP fish observation and trap operation data into the appropriate capture-mark-recapture format needed for time stratified Bayesian  models.
#'
#' Chinook "TU" dispositioned fish are marked individuals released upstream (n) of the rotary screw trap (RST) on a given date. "RE RC" dispositioned Chinook are recaptured individuals (m) at the RST on a given date.
#' The sum of "TU", "TD", "NTT", "NTD", "NTR" & "NTS" dispositioned fish are captured unmarked individuals (u) at the RST on a given date.
#'
#' Fry and NTS disposition individuals captured between the earliest trapping date to the specified smolt cut-off are broken into "YOYnraw", "YOYm", and "YOYu". Precocial Chinook, designated with precocial disposition, are excluded from the formatting. Yearling fish migrating after the specified smolt date are excluded from the formatting.
#'
#' Years are standardized in that each year begins at the earliest trapping date and ends at the latest trapping date the RST was ever in operation since RST installation.
#'
#' This package exports a .csv file to the working directory titled "Trap name, species, today's date, formatted data.csv"
#'
#' @param trap_op_data trap operation data from IDFG JTRAP database.
#' @param fish_obs_data fish observation data from IDFG JTRAP database.
#' @param strata length of desired strata in days (e.g. 1 week -> 7 days)
#' @param trap_name character string used for titles and descriptions of reports
#' @param smolt_date "MM-DD" smolt classification date needed for formatting
#' @import tidyverse
#' @export
#' @return NULL

library(here)
library(readxl)

trap_op_data = read_xlsx("data/example_input/MARTR2 OPERATIONS.xlsx")
fish_ob_data = read_xlsx("data/example_input/MARTR2 CHINOOK.xlsx")
strata = 14
smolt_date = "07-01"

Format_Chinook <- function(trap_op_data,
                           fish_ob_data,
                           strata = 7,
                           smolt_date = "07-01",
                           trap_name){

  fish_ob_data <- fish_ob_data %>%
    mutate(FishDate = as.Date(FishDate, format = "%Y-%m-%d")) %>%
    rename("StartDate" = "FishDate")

  trap_op_data<-trap_op_data %>%
    mutate(StartDate = as.Date(StartDate, format = "%Y-%m-%d"))

  JTRAP_data <- merge(fish_ob_data, trap_op_data[,-c(1:2)], by = "StartDate")

  JTRAP_data <- JTRAP_data %>%
    rename("PT2Date" = "StartDate",
           "BroodYr" = "BroodYear",
           "Disposition_Acronym" = "Acronym",
           "NFish" = "NumberOfFish",
           "FLength" = "ForkLength")

  # Convert PT2Date to a Date object
  JTRAP_data$PT2Date <- as.POSIXct(JTRAP_data$PT2Date, format = "%Y-%m-%d %H:%M:%S")

  # Format the Date object as desired (month/day/year)
  JTRAP_data$PT2Date <- format(JTRAP_data$PT2Date, format = "%m/%d/%Y")

# merged_data$SpName <- input$data_type

  JTRAP_data$PT2Date <-as.Date(JTRAP_data$PT2Date, format = "%m/%d/%Y")
  JTRAP_data$Year <- as.numeric(format(JTRAP_data$PT2Date, format = "%Y"))
  JTRAP_data$jday <- yday(JTRAP_data$PT2Date)
  JTRAP_data$Smolt_pot <- yday(as.Date(paste(JTRAP_data$Year,"-",smolt_date, sep ="")))
  JTRAP_data$BroodYr <- as.numeric(JTRAP_data$BroodYr)
  JTRAP_data$NFish <- as.numeric(JTRAP_data$NFish)

  # -----------------------------

  day.effort <- trap_op_data %>%
    select(StartDate,Operation) %>%
    rename(date = StartDate,
           effort = Operation) %>%
    mutate(effort = as.numeric(effort)) %>%
    group_by(date) %>%
    filter(!(n() > 1 & effort == 0)) %>%
    distinct()

  # -----------------------------

  ################################
  #        n,m,u formatting      #
  ################################
  # Choose  species dataset from original raw input
  #JTRAP_data <- JTRAP_data[JTRAP_data$SpName == species , ]

  JTRAP_fry <- subset(JTRAP_data, ((JTRAP_data$jday < JTRAP_data$Smolt_pot) &
                                     ((JTRAP_data$Disposition_Acronym == "NTS") |
                                        (JTRAP_data$Disposition_Acronym == "BBY") |
                                        (JTRAP_data$Disposition_Acronym == "RE BY") |
                                        (JTRAP_data$Disposition_Acronym == "YOY")
                                     )))

  JTRAP_data <- subset(JTRAP_data, !((JTRAP_data$jday < JTRAP_data$Smolt_pot) &
                                       ((JTRAP_data$Disposition_Acronym == "NTS") |
                                          (JTRAP_data$Disposition_Acronym == "BBY") |
                                          (JTRAP_data$Disposition_Acronym == "RE BY") |
                                          (JTRAP_data$Disposition_Acronym == "YOY") |
                                          ((JTRAP_data$Year-JTRAP_data$BroodYr) == 1))))

  trap_m <- JTRAP_data[JTRAP_data$Disposition_Acronym =="RE RC" ,]
  trap_n <- JTRAP_data[JTRAP_data$Disposition_Acronym =="TU" ,]
  trap_u <- JTRAP_data[(JTRAP_data$Disposition_Acronym == "TU") | (JTRAP_data$Disposition_Acronym == "TD") | (JTRAP_data$Disposition_Acronym == "NTT") | (JTRAP_data$Disposition_Acronym == "NTD") | (JTRAP_data$Disposition_Acronym == "NTR") | (JTRAP_data$Disposition_Acronym == "NTS"), ]

  # Count the number of times the "PT2Date" shows up and weight it by "NFish" to get daily counts
  # -----------------------------

  trap_m_total <- trap_m %>%
    group_by(date = PT2Date) %>%
    summarise(m = sum(NFish))
  trap_n_total <- trap_n %>%
    group_by(date = PT2Date) %>%
    summarise(n = sum(NFish))
  trap_u_total <- trap_u %>%
    group_by(date = PT2Date) %>%
    summarise(u = sum(NFish))

  # Merge the total count data sets back together

  trap_merge <- merge(merge(trap_m_total, trap_n_total, by = "date", all = TRUE),
                      trap_u_total, by = "date", all = TRUE)

  # -----------------------------

  trap_merge$date <- as.Date(trap_merge$date, format = "%m/%d/%Y")

  minyear <- as.numeric(min(format(trap_merge$date, format= "%Y"), na.rm=TRUE))
  maxyear <- as.numeric(max(format(trap_merge$date, format= "%Y"), na.rm=TRUE))

  dates <- as.data.frame(seq(from = as.Date(paste(minyear,"-01-01",sep="")),
                             to = as.Date(paste(maxyear,"-12-31",sep="")),
                             by = "days"))
  names(dates) <- c("date")

  ###############################################
  #        Calendar fill, days, and effort      #
  ###############################################

  # merge with calendar with trap data
  trap_date_filled = merge(trap_merge, dates, by="date", all=TRUE)

  ###############################################
  #    Adding Julian Date & Year Variables      #
  ###############################################

  trap_date_filled$julian = yday(trap_date_filled$date)
  trap_date_filled$year = format(trap_date_filled$date, "%Y")

  trap_full <-trap_date_filled

  #########################################
  #                YOY                    #
  #########################################

  trap_yoym <- JTRAP_fry[JTRAP_fry$Disposition_Acronym =="RE BY" ,]
  trap_yoyn <- JTRAP_fry[JTRAP_fry$Disposition_Acronym =="BBY" ,]
  trap_yoyu <- JTRAP_fry[(JTRAP_fry$Disposition_Acronym == "YOY") | (JTRAP_fry$Disposition_Acronym == "BBY") | (JTRAP_fry$Disposition_Acronym == "NTS"), ]

  # Count the number of times the "PT2Date" shows up and weight it by "NFish" to get daily counts
  # -----------------------------

  trap_yoym_total <- trap_yoym %>%
    group_by(date = PT2Date) %>%
    summarise(yoym = sum(NFish))
  trap_yoyn_total <- trap_yoyn %>%
    group_by(date = PT2Date) %>%
    summarise(yoyn = sum(NFish))
  trap_yoyu_total <- trap_yoyu %>%
    group_by(date = PT2Date) %>%
    summarise(yoyu = sum(NFish))

  # Merge the total count data sets back together

  trap_merge_yoy <- merge(merge(trap_yoym_total, trap_yoyn_total, by = "date", all = TRUE),
                          trap_yoyu_total, by = "date", all = TRUE)

  # -----------------------------

  trap_merge_yoy$date <- as.Date(trap_merge_yoy$date, format = "%m/%d/%Y")

  trap_date_yoy_filled = merge(trap_merge_yoy, dates, by="date", all=TRUE)

  # ------------------------------------

  trap_full <- merge(trap_full, trap_date_yoy_filled, by="date", all=TRUE)

  #add a column "days" and fill with "1"s
  trap_full["days"] <- 1

  #merge an "effort" column that is filled with 1 if a fish was captured or a 0 if a fish wasn't captured
  trap_full <- merge(trap_full, day.effort, by="date", all=TRUE)

  #if the trap effort is "NA", change it to "0"
  trap_full$effort[is.na (trap_full$effort)] = 0

  #########################################
  #              Stratifying              #
  #########################################

  trap_full$strata=cut(trap_full$julian, seq(0, max(trap_full$julian, na.rm = TRUE)+7, by=strata), labels=FALSE)

  # -------------------------------------------

  trap_strata <- trap_full %>% group_by(year, strata) %>%
    summarise_at(vars(m, n, u, yoym, yoyn, yoyu, days, effort), sum, na.rm = TRUE)

  # -------------------------------------------

  #Find the first and last Jweek the trap has operated through
  #strata

  # -------------------------------------------

  stra <- trap_strata %>%
    group_by(strata) %>%
    summarise(freq = sum(effort))

  # -------------------------------------------

  trap_season=stra[stra$freq > 0 ,]
  trap_start=as.numeric(min(trap_season$strata, na.rm=TRUE))
  trap_finish=as.numeric(max(trap_season$strata, na.rm=TRUE))

  #Find the first year the trap was in operation

  # -------------------------------------------

  years <- trap_strata %>%
    group_by(year) %>%
    summarise(freq = sum(effort))

  # -------------------------------------------

  trap_year=years[years$freq > 0 ,]
  trap_start_year=as.numeric(min(trap_year$year, na.rm=TRUE))

  #Chop trap data to only the 1st year of operation, first ordinal week, and last ordinal week
  trap_complete=trap_strata[trap_strata$strata >= trap_start , ]
  trap_complete=trap_complete[trap_complete$strata <= trap_finish ,]
  trap_complete=trap_complete[trap_complete$year >= trap_start_year ,]

  trap_complete$u<- ifelse((trap_complete$effort > 0) , trap_complete$u, NA)
  trap_complete$m<- ifelse((trap_complete$m > trap_complete$n), trap_complete$n, trap_complete$m)
  trap_complete$yoyu<- ifelse((trap_complete$effort > 0) , trap_complete$yoyu, NA)
  trap_complete$yoym<- ifelse((trap_complete$yoym > trap_complete$yoyn), trap_complete$yoyn, trap_complete$yoym)

  trap_complete <- trap_complete %>%
    filter(if_any(everything(), ~ !is.na(.)))

  return(trap_complete)

  #write.csv(trap_complete, file = paste(trap_name, species, today, "formatted data.csv"))
}
