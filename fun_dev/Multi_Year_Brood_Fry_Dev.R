
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

data <- read.csv(here("data/example_export/Example_data_new.csv"))

effort.cor = FALSE
effort.cor.fry = FALSE
sel.years = c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 )
strata.op.min = 1
smolt.parr.date = "07-01"
parr.presmolt.date = "09-01"
species = "CHN"
trap.name = "Marsh2_TestF"
den.plot = TRUE
trace.plot = TRUE
strata.length = 10
burnin = 2000
chains = 3
iterations = 3000
thin = 200
boot = 5000
model.params = c("p", "U", "etaP1", "etaU1", "sigmaU", "sigmaP")

Multi_Year_Brood_Fry <- function(data,
                                 effort.cor = FALSE,
                                 effort.cor.fry = FALSE,
                                 sel.years = currentyear,
                                 strata.op.min = 1,
                                 smolt.parr.date = "07-01",
                                 parr.presmolt.date = "09-01",
                                 species = "",
                                 trap.name = "",
                                 den.plot = TRUE,
                                 trace.plot = TRUE,
                                 strata.length = s.length,
                                 burnin = 200000,
                                 chains = 3,
                                 iterations = 400000,
                                 thin = 200,
                                 boot = 5000,
                                 model.params = c("p", "U", "etaP1", "etaU1", "sigmaU", "sigmaP")) {

  if (max(sel.years) > (max(data$year)-1)) {
    stop(paste("ERROR: You cannot calculate the brood year summary for the years you selected (sel.years) with the input dataset.
    Check that the maximum brood year summary does not exceed the maximum calendar year in your dataset by more than one year.
    e.g. The maximum calendar year in the dataset is ",(max(data$year)),". The maximum brood year summary (which would only include parr and presmolts)
    would be brood year ",(max(data$year)-1),".", sep =""))
  }

  if (min(sel.years) < (min(data$year)-2)) {
    stop(paste("ERROR: You cannot calculate the brood year summary for the years you selected (sel.years) with the input dataset.
    Check that the minimum brood year summary does not exceed the minimum calendar year in your dataset by more than two years.
    e.g. The minimum calendar year in the dataset is ", (min(data$year)),". The minimum brood year summary (which would only include smolts)
    would be brood year ", (min(data$year)-2),".", sep =""))
  }

  currentyear <- max(data$year)-2
  s.length <- max(data$days)
  smolt.date <- paste("2010-",smolt.parr.date, sep ="")
  parr.date <- paste("2010-",parr.presmolt.date, sep ="")
  strata_key = data %>%
    select(year,strata,strata_start,strata_end)

  model <- function() {
    for(i in 1:strata){
      for(j in 1:year){
        m[i,j] ~ dbin(p[i,j],n[i,j])
        logit(p[i,j]) <- etaP[i,j]
        etaP[i,j] ~ dnorm(etaP1[i],tauP)

        u[i,j] ~ dbin(p[i,j],U[i,j])
        U[i,j] <- round(exp(etaU[i,j]))
        etaU[i,j]~ dnorm(etaU1[i],tauU)

      }
      etaP1[i] ~ dnorm(np,tp)

      etaU1[i] ~ dnorm(nu,tu)
    }
    tauP ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmaP <- 1/sqrt(tauP)
    tauU~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmaU <- 1/sqrt(tauU) #variance

    np ~ dnorm(-2,.666) #variance=1.5, precision=1/1.5->.666
    tp ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmatp <- 1/sqrt(tp) #variance

    nu ~ dnorm(10,.25) #variance=4, precision=1/4=.25 , Increased mean to 10 due to initializing issues
    tu ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmatu <- 1/sqrt(tu) #variance
  }

  #Input Data######################################################################################

  #U updated for effort correction factor
  data$cor.factor <-data$effort/data$days
  data$u.cor <- round(data$u*(1/data$cor.factor))
  data$yoyu.cor <- round(data$yoyu*(1/data$cor.factor))

  if(effort.cor == TRUE) {
    datau=matrix(data$u.cor, nrow=length(unique(data$strata)), ncol=length(unique(data$year)))
  } else {
    datau=matrix(data$u,nrow=length(unique(data$strata)), ncol=length(unique(data$year)))
  }

  datan=matrix(data$n,nrow=length(unique(data$strata)), ncol=length(unique(data$year)))
  datam=matrix(data$m,nrow=length(unique(data$strata)), ncol=length(unique(data$year)))

  strata <- as.numeric(nrow(datan))# number of unique stratum
  year <- as.numeric(ncol(datan)) # number of unique years
  n    <- as.matrix(datan) # number of marked individuals available for recapture during the stratum
  m    <- as.matrix(datam) # number of marked individuals recaptured during the stratum
  u    <- as.matrix(datau) # number of unmarked individuals captured during the stratum

  options(max.print=100000)

  model.data <- list("strata", "year", "n", "m", "u")

  model.fit <- jags(data = model.data, inits = NULL,
                    parameters.to.save = model.params, n.chains = chains, n.iter = iterations,
                    n.burnin = burnin, n.thin = thin, model.file = model)

  model.fit.mcmc <- as.mcmc(model.fit) #call model.fit from rjags
  model.fit.gg <- suppressWarnings(ggs(model.fit.mcmc)) #convert model.fit.mcmc to gg output (useable dataframe) #Suppress warnings for CRAN repo requirements

  #######################################################################################################################
  model.fry <- function() {
    for(i in 1:strata){
      for(j in 1:year){
        yoym[i,j] ~ dbin(p[i,j],yoyn[i,j])
        logit(p[i,j]) <- etaP[i,j]
        etaP[i,j] ~ dnorm(etaP1[i],tauP)

        yoyu[i,j] ~ dbin(p[i,j],U[i,j])
        U[i,j] <- round(exp(etaU[i,j]))
        etaU[i,j]~ dnorm(etaU1[i],tauU)

      }
      etaP1[i] ~ dnorm(np,tp)

      etaU1[i] ~ dnorm(nu,tu)
    }
    tauP ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmaP <- 1/sqrt(tauP)
    tauU~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmaU <- 1/sqrt(tauU) #variance

    np ~ dnorm(-2,.666) #variance=1.5, precision=1/1.5->.666
    tp ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmatp <- 1/sqrt(tp) #variance

    nu ~ dnorm(10,.25) #variance=4, precision=1/4=.25 , Increased mean to 10 due to initializing issues
    tu ~ dgamma(.001,.001) #equal to 1/epsilon^2_lambda_mu
    sigmatu <- 1/sqrt(tu) #variance
  }

  # Input Data Fry ###################################################################################

  # convert to smolt cutoff date to strata. Fry are defined as subtaggable fish captured during Smolt period
  smolt.strata <- data %>%
    filter(strata_start == smolt.parr.date) %>%
    pull(strata) %>%
    unique() %>%
    as.numeric() - (min(data$strata)-1)

  fry.data <- data[data$strata < smolt.strata+(min(data$strata)-1),] #convert smolt cutoff strata to ordinal strata

  if(effort.cor.fry == TRUE) {
    datayoyu=matrix(round(fry.data$yoyu*(1/(fry.data$effort/fry.data$days))), nrow=length(unique(fry.data$strata)), ncol=length(unique(fry.data$year)))
  } else {
    datayoyu=matrix(fry.data$yoyu,nrow=length(unique(fry.data$strata)), ncol=length(unique(fry.data$year)))
  }

  datayoyn=matrix(fry.data$yoyn,nrow=length(unique(fry.data$strata)), ncol=length(unique(fry.data$year)))
  datayoym=matrix(fry.data$yoym,nrow=length(unique(fry.data$strata)), ncol=length(unique(fry.data$year)))

  strata <- as.numeric(nrow(datayoyn))# number of unique stratum
  year <- as.numeric(ncol(datayoyn)) # number of unique years
  yoyn    <- as.matrix(datayoyn) # number of marked individuals available for recapture during the stratum
  yoym    <- as.matrix(datayoym) # number of marked individuals recaptured during the stratum
  yoyu    <- as.matrix(datayoyu) # number of unmarked individuals captured during the stratum

  options(max.print=100000)

  model.data <- list("strata", "year", "yoyn", "yoym", "yoyu")

  model.fit.fry <- jags(data = model.data, inits = NULL,
                        parameters.to.save = model.params, n.chains = chains, n.iter = iterations,
                        n.burnin = burnin, n.thin = thin, model.file = model.fry)

  #########################

  model.fit.fry.mcmc <- as.mcmc(model.fit.fry) #call model.fit from rjags
  model.fit.gg.fry <- suppressWarnings(ggs(model.fit.fry.mcmc)) #convert model.fit.mcmc to gg output (useable dataframe)

  options(scipen = 999) #change scientific notation

  ############################
  # Create folder structures #
  ############################

  main_folder <- paste(species,"_",trap.name,"_",format(Sys.Date(), "%Y_%m_%d"),sep = "")
  dir.create(main_folder)
  dir.create(paste(main_folder,"/Inputs",sep = ""))
  dir.create(paste(main_folder,"/MCMC_Chains",sep = ""))
  dir.create(paste(main_folder,"/Model_Diagnostics",sep = ""))
  dir.create(paste(main_folder,"/Results",sep = ""))
  dir.create(paste(main_folder,"/Results/General",sep = ""))

  #########################
  # Input Data & Function #
  #########################

  options(width = 10000)  # Adjust the width to fit data in .txt for printing with sink

  data<-data[,c("year","strata","strata_start","strata_end","m","n","u","yoym","yoyn","yoyu","days","effort","cor.factor","u.cor", "yoyu.cor")]

  sink(paste(main_folder,"/Inputs/Input_Data.txt",sep = ""), append=FALSE, split=FALSE)
  print(data)
  sink()

  write.csv(data, file = paste(main_folder,"/Inputs/Input_Data.csv",sep = ""))

  sink(paste(main_folder,"/Inputs/Function_Inputs.txt",sep = ""), append=FALSE, split=FALSE)

  writeLines(paste(Sys.Date(),"\n",
                   "Multi_Year_Brood_Fry()", "\n",
                   "effort.cor = ",effort.cor, "\n",
                   "effort.cor.fry = ",effort.cor.fry, "\n",
                   "sel.years = ", paste(sel.years, collapse = ", "), "\n",
                   "strata.op.min = ", strata.op.min, "\n",
                   "smolt.parr.date = ", smolt.parr.date, "\n",
                   "parr.presmolt.date = ", parr.presmolt.date, "\n",
                   "species = ", species, "\n",
                   "trap.name = ", trap.name, "\n",
                   "den.plot = ", den.plot, "\n",
                   "trace.plot =", trace.plot, "\n",
                   "strata.length =", strata.length, "\n",
                   "burnin = ", burnin,"\n",
                   "chains =", chains, "\n",
                   "iterations =", iterations, "\n",
                   "thin = ",thin,"\n",
                   "boot = ", boot,"\n",
                   "model.params = ", paste(model.params, collapse = ", "),"\n"))
  sink()

  #################
  #   Export Juv  #
  #################

  for(i in 1:chains){
    chain <- model.fit.gg[model.fit.gg$Chain == i,]
    write.table(chain, file = paste(main_folder,"/MCMC_Chains/MCMC_Chains", i,".txt", sep = ""), sep="\t")
  }

  #################
  #Juv Diagnostics#
  #################
  gd <- gelman.diag(model.fit.mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                    multivariate=FALSE)

  gdata <-as.data.frame(gd[1])
  gdata <- cbind(Parameter = rownames(gdata), gdata)
  colnames(gdata) <- c("Parameter", "Point_est", "Upper_CI")
  nonconv.gd <- gdata[gdata$Point_est > 1.1,]
  nonconv.gd <- droplevels(nonconv.gd)

  nc.gd<-as.character(nonconv.gd$Parameter)

  if (length(nc.gd) == 0) {
    nc.gd <- "   *  All Gelman-Rubin Statistics <1.1  *   "
  }

  sink(paste(main_folder,"/Model_Diagnostics/GR_Statistic.txt",sep = ""), append=FALSE, split=FALSE)

  writeLines(paste(Sys.Date(),"\n",
                   "\n",
                   "Number of chains =", chains,"\n",
                   "Number of iterations per chain =", iterations, "\n",
                   "Burn in =", burnin,"\n",
                   "\n",
                   "********** The Gelman-Rubin Statistic: **********","\n",
                   "\n",
                   "   ***  Flagged parameters with G-R point est. >1.1  ***  ", "\n",
                   "\n",
                   paste(nc.gd, collapse=', ' ),
                   "\n",
                   "\n"))
  print(gd)

  sink()

  if(den.plot == TRUE) {
    ggmcmc(model.fit.gg, file=paste(main_folder,"/Model_Diagnostics/Density_Plots.pdf",sep = ""), plot="density")
  }
  if(trace.plot == TRUE) {
    ggmcmc(model.fit.gg, file=paste(main_folder,"/Model_Diagnostics/Trace_Plots.pdf",sep = ""), plot="traceplot")
  }

  ###################
  #Manual Statistics#
  ###################

  # summary statistic by parameter for all chains & iterations
  outputsummary <- model.fit.gg %>%
    group_by(Parameter) %>%
    summarise(
      mode = as.numeric(names(which.max(table(value)))),
      mean = mean(value),
      sd = sd(value),
      naiveSE = sd / sqrt(n()),
      quantile_2.5 = quantile(value, probs = 0.025),
      quantile_25 = quantile(value, probs = 0.25),
      quantile_50 = quantile(value, probs = 0.5),
      quantile_75 = quantile(value, probs = 0.75),
      quantile_97.5 = quantile(value, probs = 0.975))

  #add year variable
  outputsummary$year<-NA
  for(j in 1:length(unique(data$year))){
    outputsummary$year<-ifelse(grepl(paste(",",j,"]$" , sep = ""), outputsummary$Parameter), min(data$year-1)+j, outputsummary$year)
  }

  #add strata variable
  outputsummary$strata<-NA
  for(s in 1:length(unique(data$strata))){
    outputsummary$strata<-ifelse(grepl(paste(s,",", sep = ""), outputsummary$Parameter), min(data$strata-1)+s, outputsummary$strata)
  }

  # clean up summary stats
  outputsummary <- outputsummary %>%
    dplyr::rename("parameter" = "Parameter",
           "mig_year"="year") %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    arrange(mig_year, parameter)

  outputsummary = outputsummary %>%
    left_join(strata_key %>%
                dplyr::rename(mig_year=year), by = c("mig_year","strata"))

  # find strata that are below the minimum strata operation threshold.
  ag.strata <- aggregate(ceiling(data$cor.factor), by=list(strata=data$strata), FUN=sum)
  excl.strata <- ag.strata[ag.strata$x < strata.op.min,]
  x.strata <- as.character(excl.strata$strata)

  #reorder output
  outputsummary <- outputsummary[,c(1,12,13,14,11,6:10,2:5 )]

  sink(paste(main_folder,"/Results/General/All_",species,"_",trap.name,"_Results.txt",sep = ""), append=FALSE, split=FALSE)
  writeLines(paste( Sys.Date(),"\n",
                    "\n",
                    species, trap.name,"\n",
                    "Number of chains =", chains,"\n",
                    "Number of iterations per chain =", iterations, "\n",
                    "Burn in =", burnin,"\n",
                    "Thin rate =",thin, "\n",
                    "\n",
                    "* Excluded strata with <", strata.op.min, "year(s) of operation (selected minimum operational threshold) * ", "\n",
                    paste(x.strata, collapse=', ' ),
                    "\n",
                    "* Flagged parameters with Gelman-Rubin statistic >1.1  *  ", "\n",
                    paste(nc.gd, collapse=', ' ),
                    "\n"
  )
  )
  print(as.data.frame(outputsummary))
  sink()

  write.csv(outputsummary, file = paste(main_folder,"/Results/General/All_",species,"_",trap.name,"_Results.csv",sep = ""))

  #########################################################################################
  #############                      Fry Summary                           #########
  #########################################################################################

  ######################
  # Export MCMC chains #
  ######################

  for(i in 1:chains){
    chain <- model.fit.gg.fry[model.fit.gg.fry$Chain == i,]
    write.table(chain, file = paste(main_folder,"/MCMC_Chains/MCMC_Chains_Fry", i,".txt", sep = ""), sep="\t")
  }

  #############
  #Diagnostics#
  #############

  gd.fry <- gelman.diag(model.fit.fry.mcmc, confidence = 0.95, transform=FALSE, autoburnin=TRUE,
                        multivariate=FALSE)

  gdata.fry <-as.data.frame(gd.fry[1])
  gdata.fry <- cbind(Parameter = rownames(gdata.fry), gdata.fry)
  colnames(gdata.fry) <- c("Parameter", "Point_est", "Upper_CI")
  nonconv.gd.fry <- gdata.fry[gdata.fry$Point_est > 1.1,]
  nonconv.gd.fry <- droplevels(nonconv.gd.fry)

  nc.gd.fry<-as.character(nonconv.gd.fry$Parameter)

  if (length(nc.gd.fry) == 0) {
    nc.gd.fry <- "   *  All Gelman-Rubin Statistics <1.1  *   "
  }

  sink(paste(main_folder,"/Model_Diagnostics/GR_Statistic_Fry.txt",sep = ""), append=FALSE, split=FALSE)

  writeLines(paste(Sys.Date(),"\n",
                   "\n",
                   "Number of chains =", chains,"\n",
                   "Number of iterations per chain =", iterations, "\n",
                   "Burn in =", burnin,"\n",
                   "\n",
                   "********** The Gelman-Rubin Statistic: **********","\n",
                   "\n",
                   "   ***  Flagged parameters with G-R point est. >1.1  ***  ", "\n",
                   "\n",
                   paste(nc.gd.fry, collapse=', ' ),
                   "\n",
                   "\n"))
  print(gd.fry)

  sink()

  if(den.plot == TRUE) {
    ggmcmc(model.fit.gg.fry, file=paste(main_folder,"/Model_Diagnostics/Density_Plots_Fry.pdf",sep = ""), plot="density")
  }
  if(trace.plot == TRUE) {
    ggmcmc(model.fit.gg.fry, file=paste(main_folder,"/Model_Diagnostics/Trace_Plots_Fry.pdf",sep = ""), plot="traceplot")
  }

  ###################
  #Manual Statistics#
  ###################

  # summary statistic by parameter for all chains & iterations
  outputsummary.fry <- model.fit.gg.fry %>%
    group_by(Parameter) %>%
    summarise(
      mode = as.numeric(names(which.max(table(value)))),
      mean = mean(value),
      sd = sd(value),
      naiveSE = sd / sqrt(n()),
      quantile_2.5 = quantile(value, probs = 0.025),
      quantile_25 = quantile(value, probs = 0.25),
      quantile_50 = quantile(value, probs = 0.5),
      quantile_75 = quantile(value, probs = 0.75),
      quantile_97.5 = quantile(value, probs = 0.975))

  #add year variable
  outputsummary.fry$year<-NA
  for(j in 1:length(unique(data$year))){
    outputsummary.fry$year<-ifelse(grepl(paste(",",j,"]$" , sep = ""), outputsummary.fry$Parameter), min(data$year-1)+j, outputsummary.fry$year)
  }

  #add strata variable
  outputsummary.fry$strata<-NA
  for(s in 1:length(unique(data$strata))){
    outputsummary.fry$strata<-ifelse(grepl(paste(s,",", sep = ""), outputsummary.fry$Parameter), min(data$strata-1)+s, outputsummary.fry$strata)
  }

  # clean up summary stats
  outputsummary.fry <- outputsummary.fry %>%
    dplyr::rename("parameter" = "Parameter",
           "mig_year"="year") %>%
    mutate(across(where(is.numeric), ~ round(., 3))) %>%
    arrange(mig_year, parameter)

  outputsummary.fry = outputsummary.fry %>%
    left_join(strata_key %>%
                dplyr::rename(mig_year=year), by = c("mig_year","strata"))

  # find strata that are below the minimum strata operation threshold.
  fry.ag.strata <- aggregate(ceiling(fry.data$cor.factor), by=list(strata=fry.data$strata), FUN=sum)
  fry.excl.strata <- fry.ag.strata[fry.ag.strata$x < strata.op.min,]
  fry.x.strata <- as.character(fry.excl.strata$strata)

  #reorder output
  outputsummary.fry <- outputsummary.fry[,c(1,12,13,14,11,6:10,2:5)]

  sink(paste(main_folder,"/Results/General/All_",species,"_",trap.name,"_Results_Fry.txt",sep = ""), append=FALSE, split=FALSE)
  writeLines(paste( Sys.Date(),"\n",
                    "\n",
                    species, trap.name,"Fry", "\n",
                    "Number of chains =", chains,"\n",
                    "Number of iterations per chain =", iterations, "\n",
                    "Burn in =", burnin,"\n",
                    "Thin rate =",thin, "\n",
                    "\n",
                    "* Excluded strata with <", strata.op.min, "year(s) of operation (selected minimum operational threshold) * ", "\n",
                    paste(fry.x.strata, collapse=', ' ),
                    "\n",
                    "* Flagged parameters with Gelman-Rubin statistic >1.1  *  ", "\n",
                    paste(nc.gd.fry, collapse=', ' ),
                    "\n"
  )
  )

  print(as.data.frame(outputsummary.fry))
  sink()

  write.csv(outputsummary.fry, file = paste(main_folder,"/Results/General/All_",species,"_",trap.name,"_Results_Fry.csv",sep = ""))


  #########################################################################################
  #############                      Brood Summary Loop                          #########
  #########################################################################################

  # convert parr cutoff date to strata
  # NOTE: smolt.strata was defined earlier because it was needed for fry data subsetting
  parr.strata <- data %>%
    filter(strata_start == parr.presmolt.date) %>%
    pull(strata) %>%
    unique() %>%
    as.numeric() - (min(data$strata)-1)

  # Temporary
  #selectyr = 2019

  for(selectyr in sel.years){

    dir.create(paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood",sep = ""))

    #remove U parameters for parr and presmolt
    pyear <- (selectyr+2) - (as.numeric(min(data$year)))

    usep <- subset(model.fit.gg,grepl("^U", Parameter)) #subset all U's
    usep.p <- subset(usep,grepl(paste(",",pyear,"]$" , sep = ""), Parameter)) #subset U's for first year
    usep.p$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep.p$Parameter) #create strata variable
    usep.p$strata <- as.numeric(substr(usep.p$strata, 2, 5)) #clean strata variable

    # adding the for loop so the code does not run for 2 years before the minimum year in the data
    if(min(data$year) - selectyr < 2){

      parr<- subset(usep.p, strata >= (smolt.strata) & strata <= (parr.strata-1))

      parr<-subset(parr, !(parr$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

      parrdt<-as.data.table(parr)

      parrUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        parrUdist[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(parrdt$strata))
      }

      parrUdist<-as.data.frame(parrUdist) #change output to dataframe
      parrUdist$parrUdist<-as.numeric(parrUdist$parrUdist) #change output to numeric

      write.table(parrUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                                          selectyr,"_",species,"_",trap.name,"_Parr_U_Dist.txt",sep = ""), sep="\t")
      #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
      parrUoutputsummary <- parrUdist %>%
        summarise(
          mode = as.numeric(names(which.max(table(parrUdist)))),
          mean = mean(parrUdist),
          sd = sd(parrUdist),
          naiveSE = sd / sqrt(n()),
          quantile_2.5 = quantile(parrUdist, probs = 0.025),
          quantile_25 = quantile(parrUdist, probs = 0.25),
          quantile_50 = quantile(parrUdist, probs = 0.5),
          quantile_75 = quantile(parrUdist, probs = 0.75),
          quantile_97.5 = quantile(parrUdist, probs = 0.975))

      parrUoutputsummary$parameter<-paste("Parr_U_Brood_",selectyr) # add Parameter variable
      parrUoutputsummary$mig_year<- selectyr+1 #add year variable
      parrUoutputsummary$strata_start<- smolt.parr.date #add year variable
      parrUoutputsummary$strata_end<- format(as.Date(parr.presmolt.date, format = "%m-%d") - 1, "%m-%d") #add year variable

      #########Presmolt###############

      presmolt<- subset(usep.p, strata >= parr.strata)

      presmolt<-subset(presmolt, !(presmolt$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

      presmoltdt<-as.data.table(presmolt)

      presmoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        presmoltUdist[i] <- sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(presmoltdt$strata))
      }

      presmoltUdist<-as.data.frame(presmoltUdist) #change output to dataframe
      presmoltUdist$presmoltUdist<-as.numeric(presmoltUdist$presmoltUdist) #change output to numeric

      write.table(presmoltUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                                              selectyr,"_",species,"_",trap.name,"_Presmolt_U_Dist.txt",sep = ""), sep="\t")

      #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
      presmoltUoutputsummary <- presmoltUdist %>%
        summarise(
          mode = as.numeric(names(which.max(table(presmoltUdist)))),
          mean = mean(presmoltUdist),
          sd = sd(presmoltUdist),
          naiveSE = sd / sqrt(n()),
          quantile_2.5 = quantile(presmoltUdist, probs = 0.025),
          quantile_25 = quantile(presmoltUdist, probs = 0.25),
          quantile_50 = quantile(presmoltUdist, probs = 0.5),
          quantile_75 = quantile(presmoltUdist, probs = 0.75),
          quantile_97.5 = quantile(presmoltUdist, probs = 0.975))

      presmoltUoutputsummary$parameter<-paste("Presmolt_U_Brood_",selectyr) # add Parameter variable
      presmoltUoutputsummary$mig_year<- selectyr+1 #add year variable
      presmoltUoutputsummary$strata_start<- parr.presmolt.date #add year variable
      presmoltUoutputsummary$strata_end<- max(data$strata_end) #add year variable
    }

    ########## Smolt #############

    # adding the for loop so the code does not run for a year before the maximum year in the data
    if(max(data$year) - selectyr > 1){

      syear <- (selectyr+3) - (as.numeric(min(data$year)))

      usep.s <- subset(usep,grepl(paste(",",syear,"]$" , sep = ""), Parameter)) #subset U's for first year
      usep.s$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep.s$Parameter) #create strata variable
      usep.s$strata <- as.numeric(substr(usep.s$strata, 2, 5)) #clean strata variable

      smolt<- subset(usep.s, strata < smolt.strata)

      smolt<-subset(smolt, !(smolt$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

      smoltdt<-as.data.table(smolt)

      smoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        smoltUdist[i] <- sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(smoltdt$strata))
      }

      smoltUdist<-as.data.frame(smoltUdist) #change output to dataframe
      smoltUdist$smoltUdist<-as.numeric(smoltUdist$smoltUdist) #change output to numeric

      write.table(smoltUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                                           selectyr,"_",species,"_",trap.name,"_Smolt_U_Dist.txt",sep = ""), sep="\t")

      #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
      smoltUoutputsummary <- smoltUdist %>%
        summarise(
          mode = as.numeric(names(which.max(table(smoltUdist)))),
          mean = mean(smoltUdist),
          sd = sd(smoltUdist),
          naiveSE = sd / sqrt(n()),
          quantile_2.5 = quantile(smoltUdist, probs = 0.025),
          quantile_25 = quantile(smoltUdist, probs = 0.25),
          quantile_50 = quantile(smoltUdist, probs = 0.5),
          quantile_75 = quantile(smoltUdist, probs = 0.75),
          quantile_97.5 = quantile(smoltUdist, probs = 0.975))

      smoltUoutputsummary$parameter<-paste("Smolt_U_Brood_",selectyr) # add Parameter variable
      smoltUoutputsummary$mig_year<- selectyr+2 #add year variable
      smoltUoutputsummary$strata_start<- min(data$strata_start) #add year variable
      smoltUoutputsummary$strata_end<- format(as.Date(smolt.parr.date, format = "%m-%d") - 1, "%m-%d") #add year variable

    }

    ######### Fry ###############

    # adding the for loop so the code does not run for 2 years before the minimum year in the data
    if(min(data$year) - selectyr < 2){

      fyear <- (selectyr+2) - (as.numeric(min(data$year)))

      usef <- subset(model.fit.gg.fry,grepl("^U", Parameter)) #subset all U's
      usep.f <- subset(usef,grepl(paste(",",fyear,"]$" , sep = ""), Parameter)) #subset U's for first year
      usep.f$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep.f$Parameter) #create strata variable
      usep.f$strata <- as.numeric(substr(usep.f$strata, 2, 5)) #clean strata variable

      fry<-subset(usep.f, !(usep.f$strata %in% fry.excl.strata$strata)) # exclude strata under the strata.op.min threshold

      frydt<-as.data.table(fry)

      fryUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        fryUdist[i] <- sum(frydt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(frydt$strata))
      }

      fryUdist<-as.data.frame(fryUdist) #change output to dataframe
      fryUdist$fryUdist<-as.numeric(fryUdist$fryUdist) #change output to numeric

      write.table(smoltUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                                           selectyr,"_",species,"_",trap.name,"_Fry_U_Dist.txt",sep = ""), sep="\t")
      #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
      fryUoutputsummary <- fryUdist %>%
        summarise(
          mode = as.numeric(names(which.max(table(fryUdist)))),
          mean = mean(fryUdist),
          sd = sd(fryUdist),
          naiveSE = sd / sqrt(n()),
          quantile_2.5 = quantile(fryUdist, probs = 0.025),
          quantile_25 = quantile(fryUdist, probs = 0.25),
          quantile_50 = quantile(fryUdist, probs = 0.5),
          quantile_75 = quantile(fryUdist, probs = 0.75),
          quantile_97.5 = quantile(fryUdist, probs = 0.975))

      fryUoutputsummary$parameter<-paste("Fry_U_Brood_",selectyr) # add Parameter variable
      fryUoutputsummary$mig_year<- selectyr+1 #add year variable
      fryUoutputsummary$strata_start<- min(data$strata_start) #add year variable
      fryUoutputsummary$strata_end<- format(as.Date(smolt.parr.date, format = "%m-%d") - 1, "%m-%d") #add year variable

    }

    #setup bootstrap to randomly draw samples from each distribution in order to obtain total U statistics

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(data$year) - selectyr > 1){

      # adding the if condition so the code only runs for 2 years before the minimum year in the data
      if(min(data$year) - selectyr == 2){

        totUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
        for (i in 1:boot){
          totUdist[i] <- sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) -
            sum(unique(smoltdt$strata))
        }

        # the code runs for years before the (maximum year-1) and except (min year - 2)
      }else{

        totUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
        for (i in 1:boot){
          totUdist[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
            sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
            sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) -
            sum(unique(parrdt$strata)) -
            sum(unique(presmoltdt$strata)) -
            sum(unique(smoltdt$strata))
        }

        totUdistfry = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
        for (i in 1:boot){
          totUdistfry[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
            sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
            sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
            sum(frydt[, value[sample.int(.N, 1, TRUE)], by = strata]) -
            sum(unique(parrdt$strata)) -
            sum(unique(presmoltdt$strata)) -
            sum(unique(smoltdt$strata)) -
            sum(unique(frydt$strata))
        }

      }

      # the code does runs for a year before the maximum year in the data (no smolt info in this case)
    }else{

      totUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        totUdist[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
          sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) -
          sum(unique(parrdt$strata)) -
          sum(unique(presmoltdt$strata))
      }

      totUdistfry = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        totUdistfry[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
          sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) +
          sum(frydt[, value[sample.int(.N, 1, TRUE)], by = strata]) -
          sum(unique(parrdt$strata)) -
          sum(unique(presmoltdt$strata)) -
          sum(unique(frydt$strata))
      }

    }

    #Get summary statistics for U bootstrapped distribution
    totUdist<-as.data.frame(totUdist) #change output to dataframe
    totUdist$totUdist<-as.numeric(totUdist$totUdist) #change output to numeric

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    totUoutputsummary <- totUdist %>%
      summarise(
        mode = as.numeric(names(which.max(table(totUdist)))),
        mean = mean(totUdist),
        sd = sd(totUdist),
        naiveSE = sd / sqrt(n()),
        quantile_2.5 = quantile(totUdist, probs = 0.025),
        quantile_25 = quantile(totUdist, probs = 0.25),
        quantile_50 = quantile(totUdist, probs = 0.5),
        quantile_75 = quantile(totUdist, probs = 0.75),
        quantile_97.5 = quantile(totUdist, probs = 0.975))

    totUoutputsummary$parameter<-paste("Total_U_Brood_",selectyr) # add Parameter variable
    totUoutputsummary$mig_year<- NA #add year variable

    # adding the if condition so the code does not run for 2 years before the minimum year in the data
    if(min(data$year) - selectyr < 2){

      #Get summary statistics for U bootstrapped distribution
      totUdistfry<-as.data.frame(totUdistfry) #change output to dataframe
      totUdistfry$totUdistfry<-as.numeric(totUdistfry$totUdistfry) #change output to numeric

      #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
      totUoutputsummaryfry <- totUdistfry %>%
        summarise(
          mode = as.numeric(names(which.max(table(totUdistfry)))),
          mean = mean(totUdistfry),
          sd = sd(totUdistfry),
          naiveSE = sd / sqrt(n()),
          quantile_2.5 = quantile(totUdistfry, probs = 0.025),
          quantile_25 = quantile(totUdistfry, probs = 0.25),
          quantile_50 = quantile(totUdistfry, probs = 0.5),
          quantile_75 = quantile(totUdistfry, probs = 0.75),
          quantile_97.5 = quantile(totUdistfry, probs = 0.975))

      totUoutputsummaryfry$parameter<-paste("Total_U_Brood_Including_Fry",selectyr) # add Parameter variable
      totUoutputsummaryfry$mig_year<- NA #add year variable

    }

    ######### Brood strata used #########

    parr1<- subset(usep.p, strata >= (smolt.strata) & strata <= (parr.strata-1))
    presmolt1<- subset(usep.p, strata >= parr.strata)
    smolt1<- subset(usep.s, strata < smolt.strata)

    if(min(data$year) - selectyr < 2) {
      fry1 <- usep.f
    }

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(data$year) - selectyr > 1){

      brood.strata <- rbind(parr1,presmolt1,smolt1)

      # the code does runs for a year before the maximum year in the data (no smolt info in this case)
    }else{

      brood.strata <- rbind(parr1,presmolt1)

    }

    boutputsummary <- brood.strata %>%
      group_by(Parameter) %>%
      summarise(
        mode = as.numeric(names(which.max(table(value)))),
        mean = mean(value),
        sd = sd(value),
        naiveSE = sd / sqrt(n()),
        quantile_2.5 = quantile(value, probs = 0.025),
        quantile_25 = quantile(value, probs = 0.25),
        quantile_50 = quantile(value, probs = 0.5),
        quantile_75 = quantile(value, probs = 0.75),
        quantile_97.5 = quantile(value, probs = 0.975))

    #add year variable
    boutputsummary$year<-NA
    for(j in 1:length(unique(data$year))){
      boutputsummary$year<-ifelse(grepl(paste(",",j,"]$" , sep = ""), boutputsummary$Parameter), min(data$year-1)+j, boutputsummary$year)
    }

    boutputsummary$strata<-NA
    for(s in 1:length(unique(data$strata))){
      boutputsummary$strata<-ifelse(grepl(paste(s,",", sep = ""), boutputsummary$Parameter), min(data$strata-1)+s, boutputsummary$strata)
    }

    # clean up summary stats
    boutputsummary <- boutputsummary %>%
      dplyr::rename("parameter" = "Parameter",
             "mig_year" = "year") %>%
      mutate(across(where(is.numeric), ~ round(., 3))) %>%
      arrange(mig_year, parameter)

    boutputsummary = boutputsummary %>%
      left_join(strata_key %>%
                  dplyr::rename(mig_year=year), by = c("mig_year","strata"))

    # adding the if condition so the code does not run for 2 years before the minimum year in the data
    if(min(data$year) - selectyr < 2){

      foutputsummary <- fry1 %>%
        group_by(Parameter) %>%
        summarise(
          mode = as.numeric(names(which.max(table(value)))),
          mean = mean(value),
          sd = sd(value),
          naiveSE = sd / sqrt(n()),
          quantile_2.5 = quantile(value, probs = 0.025),
          quantile_25 = quantile(value, probs = 0.25),
          quantile_50 = quantile(value, probs = 0.5),
          quantile_75 = quantile(value, probs = 0.75),
          quantile_97.5 = quantile(value, probs = 0.975))

      #add year variable
      foutputsummary$year<-NA
      for(j in 1:length(unique(data$year))){
        foutputsummary$year<-ifelse(grepl(paste(",",j,"]$" , sep = ""), foutputsummary$Parameter), min(data$year-1)+j, foutputsummary$year)
      }

      foutputsummary$strata<-NA
      for(s in 1:length(unique(data$strata))){
        foutputsummary$strata<-ifelse(grepl(paste(s,",", sep = ""), foutputsummary$Parameter), min(data$strata-1)+s, foutputsummary$strata)
      }

      # clean up summary stats
      foutputsummary <- foutputsummary %>%
        dplyr::rename("parameter" = "Parameter",
               "mig_year" = "year") %>%
        mutate(across(where(is.numeric), ~ round(., 3))) %>%
        arrange(mig_year, parameter) %>%
        mutate(parameter = paste("F", parameter, sep=""))

      foutputsummary = foutputsummary %>%
        left_join(strata_key %>%
                    dplyr::rename(mig_year=year), by = c("mig_year","strata"))

    }

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(data$year) - selectyr > 1){

      # adding the if condition so the code only runs for 2 years before the minimum year in the data
      if(min(data$year) - selectyr == 2){

        broodsummary<-rbind.fill(totUoutputsummary, smoltUoutputsummary, boutputsummary) #merge outputs

        # the code runs for years before the (maximum year-1) and except (min year - 2)
      }else{

        broodsummary<-rbind.fill(totUoutputsummary,totUoutputsummaryfry,
                                 fryUoutputsummary,parrUoutputsummary,
                                 presmoltUoutputsummary, smoltUoutputsummary,
                                 foutputsummary, boutputsummary) #merge outputs

      }

      # the code does runs for a year before the maximum year in the data (no smolt info in this case)
    }else{

      broodsummary<-rbind.fill(totUoutputsummary,totUoutputsummaryfry,
                               fryUoutputsummary,parrUoutputsummary,
                               presmoltUoutputsummary,foutputsummary,
                               boutputsummary) #merge outputs

    }

    broodsummary <- broodsummary[,c(10:14,5:9,1:4)]

    data_by_year <- broodsummary

    u_by_year <- broodsummary %>%
      filter(grepl("^U\\[\\d+,\\d+\\]$", parameter))

    # Define the pattern to match
    pattern <- "\\[(\\d+,\\d+)\\]"  # This pattern matches Parameter column format

    # Use regmatches to extract matching parts from the U[x,y] parameter in cohort summary to the "All" output summary
    extracted_parts1 <- regmatches(u_by_year$parameter, regexpr(pattern, u_by_year$parameter))

    extracted_parts2 <- regmatches(outputsummary$parameter, regexpr(pattern, outputsummary$parameter))

    #Extracting the rows where matches are found
    rows <- match(extracted_parts1, extracted_parts2)

    #Selecting the matched rows from the "All" output summary
    p_by_year <- outputsummary[rows,c(3,1:2,4:14)]

    ##############################################################

    # adding the if condition so the code only runs for 2 years before the minimum year in the data
    if(min(data$year) - selectyr == 2){

      broodsummary <- rbind(broodsummary, p_by_year)

      # adding the if condition so the code runs for all years except (min year - 2)
    }else{

      #Extracting the FU[x,y] parameter
      fu_by_year <- broodsummary %>%
        filter(grepl("^FU\\[\\d+,\\d+\\]$", parameter))

      # Define the pattern to match
      pattern <- "\\[(\\d+,\\d+)\\]" # This pattern matches Parameter column format

      # Use regmatches to extract matching parts from the FU[x,y] parameter in cohort summary to the "All" FRY output summary
      extracted_parts3 <- regmatches(fu_by_year$parameter, regexpr(pattern, fu_by_year$parameter))

      extracted_parts4 <- regmatches(outputsummary.fry$parameter, regexpr(pattern, outputsummary.fry$parameter))

      #Extracting the rows where matches are found
      rows <- match(extracted_parts3, extracted_parts4)

      #Selecting the matched rows from the "All" FRY output summary
      fp_by_year <- outputsummary.fry[rows,c(3,1:2,4:14)]

      #Changing p to Fp
      fp_by_year$parameter <- gsub("p\\[(\\d+),(\\d+)\\]", "Fp[\\1,\\2]", fp_by_year$parameter)

      #Merging the capture probabilities for both adult and fry to the cohort summary
      broodsummary <- rbind(broodsummary, fp_by_year, p_by_year)

    }

    broodsummary <- broodsummary %>%
      left_join(data %>%
                  select(-c("strata_start", "strata_end")) %>%
                  dplyr::rename("mig_year" ="year"), by = c("mig_year","strata"))

    #read out files
    sink(paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
               selectyr,"_",species,"_",trap.name,"_Brood_Results.txt",sep = ""), append=FALSE, split=FALSE)
    writeLines(paste( Sys.Date(),"\n",
                      "\n",
                      selectyr, species, trap.name,"\n",
                      "Boot strap iterations =", boot,"\n",
                      "Specified smolt date:", smolt.parr.date, "\n",
                      "Specified parr date:", parr.presmolt.date, "\n",
                      "\n",
                      "Strata EXCLUDED from results due to <", strata.op.min, "year(s) of operation (e.g. below defined minimum threshold) ","\n",
                      "Parr/smolt/presmolt;", paste(x.strata, collapse=', '), "\n",
                      "Fry;", paste(fry.x.strata, collapse=', '),"\n",
                      "\n"

    )
    )
    print(broodsummary)
    sink()

    write.csv(broodsummary, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                                         selectyr,"_",species,"_",trap.name,"_Brood_Results.csv",sep = ""), sep="\t")

    ##### Figs

    # adding the if condition so the code does not run for 2 years before the minimum year in the data
    if(min(data$year) - selectyr < 2){

      UandFUdata <- rbind(fu_by_year, u_by_year)

      UandFUdata$date <- format(as.Date(paste(UandFUdata$mig_year,"-",UandFUdata$strata_start, sep = "")), format = "%b %d %Y")

      # Create a new column with the first letter of Parameter
      UandFUdata <- UandFUdata %>%
        mutate(ParameterType = sub("^([FU]+).*", "\\1", parameter))

      # Define colors for 'F' and 'p' types
      parameter_colors <- c("FU" = "#E69F00", "U" = "#56B4E9")

      p1 <- ggplot(UandFUdata, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50, color = ParameterType)) +
        geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, position = position_dodge(width = 0.5)) +
        geom_point(position = position_dodge(width = 0.5)) +
        scale_color_manual(values = parameter_colors) +
        labs(y = "Abundance (U)", x = NULL, title = paste(selectyr,species,trap.name,"Brood")) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
          axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
          axis.title = element_text(size = 14),  # Increase the font size for axis titles
          plot.title = element_text(size = 16)  # Increase the font size for the plot title
        )

      PandFPdata <- rbind(fp_by_year, p_by_year)

      PandFPdata$date <- format(as.Date(paste(PandFPdata$mig_year,"-",PandFPdata$strata_start, sep = "")), format = "%b %d %Y")

      # Create a new column with the first letter of Parameter
      PandFPdata <- PandFPdata %>%
        mutate(ParameterType = sub("^([Fp]+).*", "\\1", parameter))

      # Define colors for 'F' and 'p' types
      parameter_colors <- c("Fp" = "#E69F00", "p" = "#56B4E9")

      p2 <- ggplot(PandFPdata, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50, color = ParameterType)) +
        geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, position = position_dodge(width = 0.5)) +
        geom_point(position = position_dodge(width = 0.5)) +
        scale_color_manual(values = parameter_colors) +
        labs(y = "Capture probability (p)", x = "Strata start date") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
          axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
          axis.title = element_text(size = 14),  # Increase the font size for axis titles
          plot.title = element_text(size = 16)  # Increase the font size for the plot title
        )

      ggsave(paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                   selectyr,"_",species,"_",trap.name,"_Brood_Fry_Figures.png",sep = ""),
             plot = grid.arrange(p1, p2, ncol = 1),
             width = 10,  # Adjust the width as needed
             height = 6   # Adjust the height as needed)
      )

      u_by_year$date <- format(as.Date(paste(u_by_year$mig_year,"-",u_by_year$strata_start, sep = "")), format = "%b %d %Y")

      p3 <- ggplot(u_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
        geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#56B4E9") +
        geom_point() +
        labs(y = "Abundance (U)", x = NULL, title = paste(selectyr,species,trap.name,"Brood")) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
          axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
          axis.title = element_text(size = 14),  # Increase the font size for axis titles
          plot.title = element_text(size = 16)  # Increase the font size for the plot title
        )

      p_by_year$date <- format(as.Date(paste(p_by_year$mig_year,"-",p_by_year$strata_start, sep = "")), format = "%b %d %Y")

      p4 <- ggplot(p_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
        geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#E69F00") +
        geom_point() +
        labs(y = "Capture probability (p)", x = "Strata start date") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
          axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
          axis.title = element_text(size = 14),  # Increase the font size for axis titles
          plot.title = element_text(size = 16)  # Increase the font size for the plot title
        )

      ggsave(
        filename = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                         selectyr,"_",species,"_",trap.name,"_Brood_Figures.png",sep = ""),
        plot = grid.arrange(p3, p4, ncol = 1),
        width = 10,  # Adjust the width as needed
        height = 6   # Adjust the height as needed
      )


      # adding the if condition so the code runs for 2 years before the minimum year in the data (no fry info in this case)
    }else{

      u_by_year$date <- format(as.Date(paste(u_by_year$mig_year,"-",u_by_year$strata_start, sep = "")), format = "%b %d %Y")

      p3 <- ggplot(u_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
        geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#56B4E9") +
        geom_point() +
        labs(y = "Abundance (U)", x = NULL, title = paste(selectyr,species,trap.name,"Brood")) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
          axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
          axis.title = element_text(size = 14),  # Increase the font size for axis titles
          plot.title = element_text(size = 16)  # Increase the font size for the plot title
        )

      p_by_year$date <- format(as.Date(paste(p_by_year$mig_year,"-",p_by_year$strata_start, sep = "")), format = "%b %d %Y")

      p4 <- ggplot(p_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
        geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#E69F00") +
        geom_point() +
        labs(y = "Capture probability (p)", x = "Strata start date") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
          axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
          axis.title = element_text(size = 14),  # Increase the font size for axis titles
          plot.title = element_text(size = 16)  # Increase the font size for the plot title
        )

      ggsave(
        filename = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                         selectyr,"_",species,"_",trap.name,"_Brood_Figures.png",sep = ""),
        plot = grid.arrange(p3, p4, ncol = 1),
        width = 10,  # Adjust the width as needed
        height = 6   # Adjust the height as needed
      )

    }
  }
}


Multi_Year_Brood_Fry(data,
                 effort.cor = FALSE,
                 effort.cor.fry = FALSE,
                 sel.years = c(2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021),
                 strata.op.min = 1,
                 smolt.parr.date = "07-01",
                 parr.presmolt.date = "09-01",
                 species = "Chnk",
                 trap.name = "ScrewTrap",
                 den.plot = TRUE,
                 trace.plot = TRUE,
                 strata.length = 10,
                 burnin = 2000,
                 chains = 3,
                 iterations = 4000,
                 thin = 200,
                 boot = 5000,
                 model.params = c("p","U","etaP1","etaU1","sigmaU","sigmaP"))

