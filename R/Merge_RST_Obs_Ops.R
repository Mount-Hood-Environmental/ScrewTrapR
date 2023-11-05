#' @title Merge screw trap fish observation and trap operation data
#'
#' @description This function merges raw IDFG screw trap fish observation and trap operation data into one file that can then be used with the
#' `Format_Chinook()` or `Format_Steelhead()` functions that can then be used to summarize the trap data into the capture-mark-recapture format
#' needed for the bayesian models.
#' @param fish_observation trap fish observation data. Needs to include columns "FishDate","Acronym","NumberOfFish", & "ForkLength"
#' @param trap_operations trap operation data. Needs to include columns "StartDate" and "Operation"
#
#'
#' @import tidyverse
#' @export
#' @return NULL

Merge_RST_Obs_Ops <- function(fish_observation,
                              trap_operations) {

  fish <- fish_observation %>%
    mutate(FishDate = as.Date(FishDate, format = "%Y-%m-%d")) %>%
    rename("StartDate" = "FishDate")

  trap<-trap_operations %>%
    mutate(StartDate = as.Date(StartDate, format = "%Y-%m-%d"),
           Operation = as.numeric(Operation)) %>%
    group_by(StartDate) %>%
    filter(Operation %in% c(1.0, 0.5, 0.0)) %>%
    slice(which.max(Operation %in% c(1.0, 0.5, 0.0)))

  merged_data <- full_join(fish, trap[,-c(1:2)], by = "StartDate")

  fish_obs_max <- max(fish$StartDate)

  merged_data <- merged_data %>%
    rename("PT2Date" = "StartDate",
           "BroodYr" = "BroodYear",
           "Disposition_Acronym" = "Acronym",
           "NFish" = "NumberOfFish",
           "FLength" = "ForkLength") %>%
    filter(PT2Date <= fish_obs_max)

  return(merged_data)
}
