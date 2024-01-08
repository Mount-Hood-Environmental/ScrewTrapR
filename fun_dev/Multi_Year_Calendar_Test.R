
library(plyr)
library(lubridate)
library(R2jags)
library(coda)
library(lattice)
library(superdiag)
library(mcmcplots)
library(ggmcmc)
library(data.table)
library(gridExtra)
library(here)

getwd()

data <- read.csv(here("data/example_export/Example_data_new.csv"))

source(here("R/Multi_Year_Calendar.R"))

Multi_Year_Calendar(data,
                 effort.cor = FALSE,
                 sel.years = c(2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022 ),
                 strata.op.min = 1,
                 smolt.juv.date = "07-01",
                 species = "STHD",
                 trap.name = "CALENDAR_TEST",
                 den.plot = TRUE,
                 trace.plot = TRUE,
                 strata.length = 10,
                 burnin = 2000,
                 chains = 3,
                 iterations = 4000,
                 thin = 200,
                 boot = 5000,
                 model.params = c("p", "U", "etaP1", "etaU1", "sigmaU", "sigmaP"))
