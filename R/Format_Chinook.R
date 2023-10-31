# Authors:  Bryce Oldemeyer
# Purpose:  Format the merged RST observation and operations into Capture-Mark-Recapture
# data needed for hierarchical Bayesian models. Chinook salmon data.
# Created:  10/19/2023
# Last Modified: 10/19/2023
# Published here:

Format_Chinook <- function(RST_ops_obs_data,
                           strata = 7,
                           smolt.date = "07-01",
                           species,
                           trap.name){

  RST_ops_obs_data$PT2Date <-as.Date(RST_ops_obs_data$PT2Date, format = "%m/%d/%Y")
  RST_ops_obs_data$Year <-as.numeric(format(RST_ops_obs_data$PT2Date, format = "%Y"))
  RST_ops_obs_data$jday <- yday(RST_ops_obs_data$PT2Date)
  RST_ops_obs_data$Smolt_pot <- yday(as.Date(paste(RST_ops_obs_data$Year,"-",smolt.date, sep ="")))
  RST_ops_obs_data$BroodYr <- as.numeric(RST_ops_obs_data$BroodYr)
  RST_ops_obs_data$NFish <- as.numeric(RST_ops_obs_data$NFish)

  # -----------------------------

  day.effort <- RST_ops_obs_data %>%
    select(PT2Date,Operation) %>%
    distinct() %>%
    rename("date" = "PT2Date")

  # -----------------------------

  ################################
  #        n,m,u formatting      #
  ################################
  # Choose  species dataset from original raw input

  RST_ops_obs_data <- RST_ops_obs_data[RST_ops_obs_data$SpName == species , ]

  RST_ops_obs_data_fry <- subset(RST_ops_obs_data, ((RST_ops_obs_data$jday < RST_ops_obs_data$Smolt_pot) &
                                                      ((RST_ops_obs_data$Disposition_Acronym == "NTS") |
                                                         (RST_ops_obs_data$Disposition_Acronym == "BBY") |
                                                         (RST_ops_obs_data$Disposition_Acronym == "RE BY") |
                                                         (RST_ops_obs_data$Disposition_Acronym == "YOY")
                                                      )))

  RST_ops_obs_data <- subset(RST_ops_obs_data, !((RST_ops_obs_data$jday < RST_ops_obs_data$Smolt_pot) &
                                                   ((RST_ops_obs_data$Disposition_Acronym == "NTS") |
                                                      (RST_ops_obs_data$Disposition_Acronym == "BBY") |
                                                      (RST_ops_obs_data$Disposition_Acronym == "RE BY") |
                                                      (RST_ops_obs_data$Disposition_Acronym == "YOY") |
                                                      ((RST_ops_obs_data$Year-RST_ops_obs_data$BroodYr) == 1))))

  trap_m <- RST_ops_obs_data[RST_ops_obs_data$Disposition_Acronym =="RE RC" ,]
  trap_n <- RST_ops_obs_data[RST_ops_obs_data$Disposition_Acronym =="TU" ,]
  trap_u <- RST_ops_obs_data[(RST_ops_obs_data$Disposition_Acronym == "TU") | (RST_ops_obs_data$Disposition_Acronym == "TD") | (RST_ops_obs_data$Disposition_Acronym == "NTT") | (RST_ops_obs_data$Disposition_Acronym == "NTD") | (RST_ops_obs_data$Disposition_Acronym == "NTR") | (RST_ops_obs_data$Disposition_Acronym == "NTS"), ]

  # summarize NFish by date
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

  trap_yoym <- RST_ops_obs_data_fry[RST_ops_obs_data_fry$Disposition_Acronym =="RE BY" ,]
  trap_yoyn <- RST_ops_obs_data_fry[RST_ops_obs_data_fry$Disposition_Acronym =="BBY" ,]
  trap_yoyu <- RST_ops_obs_data_fry[(RST_ops_obs_data_fry$Disposition_Acronym == "YOY") | (RST_ops_obs_data_fry$Disposition_Acronym == "BBY") | (RST_ops_obs_data_fry$Disposition_Acronym == "NTS"), ]

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
  trap_full <- merge(trap_full, day.effort, by="date", all=TRUE)

  trap_full = trap_full %>%
    mutate(Operation = ifelse(is.na(Operation), 0, Operation),
           yoyn = lag(yoyn), # lag "n" by one day for yoy and taggables
           n = lag(n)) %>%
    rename("effort" = "Operation")

  #add a column "days" and fill with "1"s
  trap_full["days"] <- 1

  #########################################
  #              Stratifying              #
  #########################################

  trap_full$strata=cut(trap_full$julian, seq(0, max(trap_full$julian, na.rm = TRUE)+strata, by=strata), labels=FALSE)

  # -------------------------------------------

  trap_week <- trap_full %>%
    group_by(year, strata) %>%
    summarise_at(vars(m, n, u, yoym, yoyn, yoyu, days, effort), sum, na.rm = TRUE)

  # -------------------------------------------


  #Find the first and last Jweek the trap has operated through
  #jweek

  # -------------------------------------------

  weeks <- trap_week %>%
    group_by(strata) %>%
    summarise(freq = sum(effort))

  # -------------------------------------------

  trap_season=weeks[weeks$freq > 0 ,]
  trap_start=as.numeric(min(trap_season$strata, na.rm=TRUE))
  trap_finish=as.numeric(max(trap_season$strata, na.rm=TRUE))

  #Find the first year the trap was in operation

  # -------------------------------------------

  years <- trap_week %>%
    group_by(year) %>%
    summarise(freq = sum(effort))

  # -------------------------------------------

  trap_year=years[years$freq > 0 ,]
  trap_start_year=as.numeric(min(trap_year$year, na.rm=TRUE))

  #Chop trap data to only the 1st year of operation, first ordinal week, and last ordinal week
  trap_complete=trap_week[trap_week$strata >= trap_start , ]
  trap_complete=trap_complete[trap_complete$strata <= trap_finish ,]
  trap_complete=trap_complete[trap_complete$year >= trap_start_year ,]

  trap_complete$u<- ifelse((trap_complete$effort > 0) , trap_complete$u, NA)
  trap_complete$m<- ifelse((trap_complete$m > trap_complete$n), trap_complete$n, trap_complete$m)
  trap_complete$yoyu<- ifelse((trap_complete$effort > 0) , trap_complete$yoyu, NA)
  trap_complete$yoym<- ifelse((trap_complete$yoym > trap_complete$yoyn), trap_complete$yoyn, trap_complete$yoym)

  trap_complete = trap_complete %>%
    filter_all(any_vars(!is.na(.)))

  today <- Sys.Date()

  return(trap_complete)

  #write.csv(trap_complete, file = paste(trap.name, species, today, "formatted data.csv"))
}
