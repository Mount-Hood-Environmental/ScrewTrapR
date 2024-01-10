#' @title Format Chinook Salmon screw trap data into capture-mark-recapture data needed for Bayesian models
#'
#' @description This function reconfigures Chinook Salmon screw trap data into the appropriate capture-mark-recapture format needed for the hierarchical Bayesian models.
#'
#' Chinook "TU" dispositioned fish are marked individuals released upstream (n) of the screw trap on for the prior day. "RE RC" dispositioned Chinook are recaptured individuals (m) at the RST on a given date.
#' The sum of "TU", "TD", "NTT", "NTD", "NTR" & "NTS" dispositioned fish are captured unmarked individuals (u) at the RST on a given date.
#'
#' Fry and NTS disposition individuals captured between the earliest trapping date to the specified smolt cut-off are broken into "YOYnraw", "YOYm", and "YOYu". Precocial Chinook, designated with precocial disposition, are excluded from the formatting. Yearling fish migrating after the specified smolt date are excluded from the formatting.
#'
#' Years are standardized in that each year begins at the earliest trapping date and ends at the latest trapping date the RST was ever in operation since RST installation.
#'
#' @param RST_ops_obs_data merged trap operation and observation data for a single species
#' @param strata length of desired strata in days (e.g. 1 week -> 7 days)
#' @param species character string needed to subset data (should be identical to the species name in the RST_ops_obs_data set)
#' @param trap.name character string used for titles and descriptions of reports
#' @param smolt.date "MM-DD" smolt/parr life stage classification date needed for formatting. Date is the beginning of the parr life stage.
#' @param parr.presmolt.date "MM-DD" parr/presmolt classification date needed for formatting. Date is the beginning of the smolt life stage.
#'
#' @import tidyverse
#' @export
#' @return NULL

Format_Chinook <- function(RST_ops_obs_data,
                           strata = 7,
                           smolt.date = "07-01",
                           parr.presmolt.date = "09-01",
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
    dplyr::rename("date" = "PT2Date") %>%
    filter(!(month(date) == 2 & day(date) == 29 & leap_year(year(date))))

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

  # merge with calendar with trap data & remove Feb. 29 during leap years for consistency.
  trap_date_filled = merge(trap_merge, dates, by="date", all=TRUE) %>%
    filter(!(month(date) == 2 & day(date) == 29 & leap_year(year(date))))

  ###############################################
  #    Adding Julian Date & Year Variables      #
  ###############################################

  trap_date_filled$julian = yday(trap_date_filled$date)

  trap_date_filled = trap_date_filled %>%
    mutate(julian = ifelse(leap_year(year(date)) & month(date) > 2 & month(date) <= 12, julian - 1, julian))

  trap_date_filled$year = format(trap_date_filled$date, "%Y")

  trap_full = trap_date_filled

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

  trap_date_yoy_filled = merge(trap_merge_yoy, dates, by="date", all=TRUE) %>%
    filter(!(month(date) == 2 & day(date) == 29 & leap_year(year(date))))

  # ------------------------------------

  trap_full <- merge(trap_full, trap_date_yoy_filled, by="date", all=TRUE)
  trap_full <- merge(trap_full, day.effort, by="date", all=TRUE)

  trap_full = trap_full %>%
    mutate(Operation = ifelse(is.na(Operation), 0, Operation),
           yoyn = lag(yoyn), # lag "n" by one day for yoy and taggables
           n = lag(n)) %>%
    dplyr::rename("effort" = "Operation")

  #add a column "days" and fill with "1"s
  trap_full["days"] <- 1

  #########################################
  #              Stratifying              #
  #########################################

  trap_smolt = trap_full %>%
    mutate(monthDay = format(trap_full$date, "%m-%d")) %>%
    filter(monthDay < smolt.date)

  # Make it so strata are grouped by days starting from the identified smolt date and working backwards to the beginning
  # of the year. This makes it so summaries for the smolt life-stage end at the defined date.

  trap_smolt$strata=cut(rev(trap_smolt$julian), seq(0, max(trap_smolt$julian, na.rm = TRUE)+strata, by=strata),
                        labels=seq(ceiling(max(trap_smolt$julian, na.rm = TRUE)/strata),1), right = TRUE)

  trap_parr = trap_full %>%
    mutate(monthDay = format(trap_full$date, "%m-%d")) %>%
    filter(monthDay >= smolt.date & monthDay< parr.presmolt.date)

  parr_seq <-c(seq(min(trap_parr$julian, na.rm = TRUE), max(trap_parr$julian, na.rm = TRUE), by=strata), max(trap_parr$julian, na.rm = TRUE)+1)

  trap_parr$strata=cut(trap_parr$julian,
                       breaks = parr_seq,
                       labels=seq(max(as.numeric(trap_smolt$strata), na.rm = TRUE)+1, max(as.numeric(trap_smolt$strata), na.rm = TRUE)+length(parr_seq)-1),
                       right = FALSE)

  # IDFG biologist set hard dates for smolt/parr and parr/presmolt life-stage designations (typically 07-01 & 09-01). These hard date designations
  # often result in the last strata of the parr life-stage not having enough days to fill a full strata.
  # e.g. Aug 12-18 = 7 days, Aug 19-25 = 7 days, Aug 26-31 = 6 days.

  # Since this typically happens at the end of August during base flow conditions, combining data during this time won't likely violate
  # the assumption of homogeneous capture probabilities if the last two strata were combined. Additionally, the large pulse of presmolts typically occurs
  # in October for most systems.

  # The multi-year model doesn't have any issues dealing with variable strata sizes but the spline model may over- or under-
  # estimate depending on data characteristics during the time period. A good compromise would be to combine the last two strata
  # if the last strata was smaller than half the designated strata size.

  # Example 1: Designated strata size = 7 days.
  # Aug 12-18 = 7 days, Aug 19-25 = 7 days, Aug 26-31 = 6 days.
  # Leave as-is

  # Example 2: Designated strata size = 10 days.
  # Aug 10-19 = 10 days, Aug 20-29 = 10 days, Aug 30-31 = 2 days.
  # Merge last two strata: Aug 10-19 = 10 days, Aug 20-31 = 12 days

  input_strata = strata  # oversite in duplicate naming in columns and function input.

  trap_parr <- trap_parr %>%
    group_by(year, strata) %>%
    mutate(row_count = n()) %>%
    ungroup() %>%
    mutate(strata = if_else(row_count < input_strata/2, max(as.numeric(as.character(trap_parr$strata)), na.rm = TRUE)-1, as.numeric(as.character(strata)))) %>%
    select(-row_count)

  trap_presmolt = trap_full %>%
    mutate(monthDay = format(trap_full$date, "%m-%d")) %>%
    filter(monthDay>= parr.presmolt.date)

  presmolt_seq <- c(seq(min(trap_presmolt$julian, na.rm = TRUE), max(trap_presmolt$julian, na.rm = TRUE), by=strata), max(trap_presmolt$julian, na.rm = TRUE)+1)

  trap_presmolt$strata=cut(trap_presmolt$julian,
                           breaks = presmolt_seq,
                           labels=seq(max(as.numeric(trap_parr$strata), na.rm = TRUE)+1, max(as.numeric(trap_parr$strata), na.rm = TRUE)+length(presmolt_seq)-1),
                           right = FALSE)

  trap_smolt = trap_smolt %>%
    mutate(strata = as.numeric(as.character(strata)))
  trap_parr = trap_parr %>%
    mutate(strata = as.numeric(as.character(strata)))
  trap_presmolt = trap_presmolt %>%
    mutate(strata = as.numeric(as.character(strata)))

  trap_strata_daily = trap_smolt %>%
    bind_rows(trap_parr, trap_presmolt) %>%
    arrange(.,date) %>%
    group_by(strata) %>%
    mutate(strata_start = min(monthDay),
           strata_end = max(monthDay),
           strata = as.numeric(as.character(strata)))

  # -------------------------------------------

  trap_strata <- trap_strata_daily %>%
    group_by(year, strata) %>%
    summarise_at(vars(m, n, u, yoym, yoyn, yoyu, days, effort), sum, na.rm = TRUE) %>%
    mutate(year = as.numeric(as.character(year)),
           strata = as.numeric(as.character(strata)))

  # Find the first and last strata the trap operated
  strata_summary <- trap_strata %>%
    group_by(strata) %>%
    summarise(freq = sum(effort))

  trap_season=strata_summary[strata_summary$freq > 0 ,]
  trap_start= min(trap_season$strata, na.rm=TRUE)
  trap_finish=max(trap_season$strata, na.rm=TRUE)

  # Find the first year the trap was in operation
  years <- trap_strata %>%
    group_by(year) %>%
    summarise(freq = sum(effort))

  trap_year=years[years$freq > 0 ,]
  trap_start_year=min(trap_year$year, na.rm=TRUE)

  # trim trap data from the first year of operation, first strata, and last strata
  trap_complete=trap_strata[trap_strata$strata >= trap_start ,]
  trap_complete=trap_complete[trap_complete$strata <= trap_finish ,]
  trap_complete=trap_complete[trap_complete$year >= trap_start_year ,]

  trap_complete$u<- ifelse((trap_complete$effort > 0) , trap_complete$u, NA)
  trap_complete$m<- ifelse((trap_complete$m > trap_complete$n), trap_complete$n, trap_complete$m)
  trap_complete$yoyu<- ifelse((trap_complete$effort > 0) , trap_complete$yoyu, NA)
  trap_complete$yoym<- ifelse((trap_complete$yoym > trap_complete$yoyn), trap_complete$yoyn, trap_complete$yoym)

  trap_complete = trap_complete %>%
    filter_all(any_vars(!is.na(.)))

  trap_complete = trap_complete %>%
    left_join(trap_strata_daily %>%
                select(strata,strata_start,strata_end) %>%
                distinct(),
              by = join_by(strata)) %>%
    arrange(., by = year, strata)


  return(trap_complete)

  #write.csv(trap_complete, file = paste(trap.name, species, today, "formatted data.csv"))
}
