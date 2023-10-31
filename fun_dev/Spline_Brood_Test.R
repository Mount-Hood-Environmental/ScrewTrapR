library(BTSPAS)
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

data <- read.csv(here("data/example_export/Example_data.csv"))

source("shiny_app/Spline_Brood.R")

Spline_Brood(data,
                 effort.cor = FALSE,
                 sel.years = c(2018),
                 smolt.parr.date = "07-01",
                 parr.presmolt.date = "09-01",
                 species = "CHN",
                 trap.name = "MARTRP_Brood_SPLINE_Test",
                 den.plot = TRUE,
                 trace.plot = TRUE,
                 strata.length = 10,
                 burnin = 200000,
                 chains = 3,
                 iterations = 400000,
                 thin = 200,
                 boot = 5000,
                 model.params = c("p", "U", "etaP1", "etaU1", "sigmaU", "sigmaP"))
