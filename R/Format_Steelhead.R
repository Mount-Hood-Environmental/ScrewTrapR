#' @title Format Chinook Salmon screw trap data into capture-mark-recapture data needed for Bayesian models
#'
#' @description This function reconfigures Chinook Salmon screw trap data into the appropriate capture-mark-recapture format needed for the hierarchical Bayesian models.
#'
#' Steelhead  "TU" dispositioned fish are marked individuals released upstream (n) of the rotary screw trap (RST) on a the prior date. "RE RC" dispositioned Steelhead are recaptured individuals (m) at the RST on a given date.
#' The sum of "TU", "TD", "NTT", "NTD", and "NTR" dispositioned fish are captured unmarked individuals (u) at the RST on a given date.
#'
#' Individuals smaller than the fork length cut-off (flength.cut) are excluded from the formatting. Precocial Steelhead, designated with precocial disposition, are excluded from the formatting.
#'
#' Years are standardized in that each year begins at the earliest trapping date and ends at the latest trapping date the RST was ever in operation since RST installation.
#'
#' @param RST_ops_obs_data merged trap operation and observation data for a single species
#' @param strata length of desired strata in days (e.g. 1 week -> 7 days)
#' @param species character string needed to subset data (should be identical to the species name in the RST_ops_obs_data set)
#' @param trap.name character string used for titles and descriptions of reports
#' @param smolt.date "MM-DD" smolt classification date needed for formatting
#' @param flength.cut lower fork length cut off for steelhead to be included in the analysis
#'
#' @import tidyverse
#' @export
#' @return NULL

Format_Steelhead <- function(RST_ops_obs_data,
                             strata = 7,
                             smolt.date = "06-01",
                             flength.cut = 80,
                             species,
                             trap.name){

  RST_ops_obs_data$PT2Date <-as.Date(RST_ops_obs_data$PT2Date, format = "%m/%d/%Y")
  RST_ops_obs_data$Year <-as.numeric(format(RST_ops_obs_data$PT2Date, format = "%Y"))
  RST_ops_obs_data$jday <- yday(RST_ops_obs_data$PT2Date)
  RST_ops_obs_data$Smolt_pot <- yday(as.Date(paste(RST_ops_obs_data$Year,"-",smolt.date, sep ="")))
  RST_ops_obs_data$NFish <- as.numeric(RST_ops_obs_data$NFish)
  RST_ops_obs_data$FLength <-as.numeric(RST_ops_obs_data$FLength)

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

  #Clean data of steelhead that aren't taggable length
  RST_ops_obs_data <- subset(RST_ops_obs_data, (RST_ops_obs_data$FLength >= flength.cut) | (RST_ops_obs_data$FLength == 0))
  trap_m <- RST_ops_obs_data[RST_ops_obs_data$Disposition_Acronym =="RE RC" ,]
  trap_n <- RST_ops_obs_data[RST_ops_obs_data$Disposition_Acronym =="TU" ,]
  trap_u <- RST_ops_obs_data[(RST_ops_obs_data$Disposition_Acronym == "TU") | (RST_ops_obs_data$Disposition_Acronym == "TD") | (RST_ops_obs_data$Disposition_Acronym == "NTT") | (RST_ops_obs_data$Disposition_Acronym == "NTD") | (RST_ops_obs_data$Disposition_Acronym == "NTR"), ]

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

  trap_full <- merge(trap_full, day.effort, by="date", all=TRUE)

  trap_full = trap_full %>%
    mutate(Operation = ifelse(is.na(Operation), 0, Operation),
           n = lag(n)) %>% # lag "n" by one day for yoy and taggables
    rename("effort" = "Operation")

  #add a column "days" and fill with "1"s
  trap_full["days"] <- 1

  #########################################
  #              Stratifying              #
  #########################################

  trap_smolt = trap_full %>%
    mutate(monthDay = format(trap_date_filled$date, "%m-%d")) %>%
    filter(monthDay < smolt.date)

  # Make it so strata are grouped by days starting from the identified smolt date and working backwards to the beginning
  # of the year. This makes it so summaries for the smolt life-stage end at the defined date.

  trap_smolt$strata=cut(rev(trap_smolt$julian), seq(0, max(trap_smolt$julian, na.rm = TRUE)+strata, by=strata),
                        labels=seq(ceiling(max(trap_smolt$julian, na.rm = TRUE)/strata),1), right = TRUE)

  trap_juv = trap_full %>%
    mutate(monthDay = format(trap_date_filled$date, "%m-%d")) %>%
    filter(monthDay >= smolt.date)

  juv_seq <- c(seq(min(trap_juv$julian, na.rm = TRUE), max(trap_juv$julian, na.rm = TRUE), by=strata), max(trap_juv$julian, na.rm = TRUE)+1)

  trap_juv$strata=cut(trap_juv$julian,
                      breaks = juv_seq,
                      labels= seq(max(as.numeric(trap_smolt$strata, na.rm = TRUE))+1, max(as.numeric(trap_smolt$strata, na.rm = TRUE))+length(juv_seq)-1),
                      right = FALSE)

  # -------------------------------------------
  trap_smolt = trap_smolt %>%
    mutate(strata = as.numeric(as.character(strata)))
  trap_juv = trap_juv %>%
    mutate(strata = as.numeric(as.character(strata)))

  trap_strata_daily = trap_smolt %>%
    bind_rows(trap_juv) %>%
    arrange(.,date) %>%
    group_by(strata) %>%
    mutate(strata_start = min(monthDay),
           strata_end = max(monthDay),
           strata = as.numeric(as.character(strata)))

  trap_strata <- trap_strata_daily %>%
    group_by(year, strata) %>%
    summarise_at(vars(m, n, u, days, effort), sum, na.rm = TRUE) %>%
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
  trap_complete=trap_strata[trap_strata$strata >= trap_start , ]
  trap_complete=trap_complete[trap_complete$strata <= trap_finish ,]
  trap_complete=trap_complete[trap_complete$year >= trap_start_year ,]

  trap_complete$u<- ifelse((trap_complete$effort > 0) , trap_complete$u, NA)
  trap_complete$m<- ifelse((trap_complete$m > trap_complete$n), trap_complete$n, trap_complete$m)

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
