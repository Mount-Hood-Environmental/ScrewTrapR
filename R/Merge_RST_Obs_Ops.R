# Authors:  Bryce Oldemeyer
# Purpose:  Merge RST fish observation and trap operations data from IDFG databases. Same function as in ScrewTrapR package.
# Created:  10/19/2023
# Last Modified: 10/19/2023
# Published here:

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
