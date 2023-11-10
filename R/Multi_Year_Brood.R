
#' @title Multi-year time stratified Bayesian estimator (brood summary)

#' @description This function implements a Bayesian model with a hierarchical structure for U and p by strata between years.
#' This function was intended to be used to summarize juvenile Chinook Salmon abundances by brood year.
#' e.g. Summaries for brood year 2012 consist of parr (calendar year 2013), presmolts (calendar year 2013), smolts (calendar year 2014).
#'
#'
#' @param data capture-mark-recapture data frame
#' @param burnin number of initial MCMC chain iterations to be discarded
#' @param chains number of MCMC chains (>1)
#' @param iterations number of MCMC iterations per chain
#' @param thin thin rate for MCMC chains
#' @param sel.years selected year(s) to compute abundance estimates
#' @param boot number of boot strap iterations to calculate yearly and life stage abundance estimates
#' @param species character string used for titles and descriptions of reports
#' @param model.params parameters to be returned from MCMC simulation
#' @param trap.name character string used for  titles and descriptions of reports
#' @param effort.cor expands the number of unmarked fish captured by a sampling effort
#' @param strata.length number of days in strata
#' @param smolt.parr.date "MM-DD" date to partition smolt life stage
#' @param parr.presmolt.date "MM-DD" date to partition parr life stage
#' @param strata.op.min minimum number of years data need to have been collected in a stratum to be included in summary
#' @param den.plot return density plots of MCMC chains
#' @param trace.plot return trace plots of MCMC chains#'
#' @import plyr
#' @import lubridate
#' @import R2jags
#' @import coda
#' @import lattice
#' @import superdiag
#' @import mcmcplots
#' @import ggmcmc
#' @import gridExtra
#' @importFrom data.table as.data.table
#' @export
#' @return NULL
#'
#'

Multi_Year_Brood <- function(data,
                             effort.cor = FALSE,
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

  # U updated for effort correction factor
  data$cor.factor <-data$effort/data$days
  data$u.cor <- round(data$u*(1/(data$effort/data$days)))

  data$strata_date<-format(as.Date(data$strata*strata.length,
                                   origin = paste(data$year, "-01-01", sep = "")))

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
                    n.burnin = burnin, model.file = model, n.thin = thin)

  #######################################################################################################################

  model.fit.mcmc <- as.mcmc(model.fit) #call model.fit from rjags
  model.fit.gg <- suppressWarnings(ggs(model.fit.mcmc)) #convert model.fit.mcmc to gg output (useable dataframe)

  options(scipen = 999) #change scientific notation

  main_folder <- paste(species,"_",trap.name,"_MYB_",format(Sys.Date(), "%Y_%m_%d"),sep = "")
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

  data<-data[,c(2,3,14,4:13)]

  sink(paste(main_folder,"/Inputs/Input_Data.txt",sep = ""), append=FALSE, split=FALSE)
  print(data)
  sink()

  write.csv(data, file = paste(main_folder,"/Inputs/Input_Data.csv",sep = ""))

  sink(paste(main_folder,"/Inputs/Function_Inputs.txt",sep = ""), append=FALSE, split=FALSE)
  writeLines(paste(Sys.Date(),"\n",
                   "Multi_Year_Brood()", "\n",
                   "effort.cor = ",effort.cor, "\n",
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

  ######################
  # Export MCMC chains #
  ######################

  for(i in 1:chains){
    chain <- model.fit.gg[model.fit.gg$Chain == i,]
    write.table(chain, file = paste(main_folder,"/MCMC_Chains/MCMC_Chains", i,".txt", sep = ""), sep="\t")
  }

  ###############
  # Diagnostics #
  ###############

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

  #####################
  # Manual Statistics #
  #####################

  # summary statistic by parameter for all chains & iterations
  outputsummary <- model.fit.gg %>%
    group_by(Parameter) %>%
    dplyr::summarise(
      mode = as.numeric(names(which.max(table(value)))),
      mean = mean(value),
      sd = sd(value),
      naiveSE = sd / sqrt(length(value)),
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

  outputsummary$strata_date<-format(as.Date(outputsummary$strata*strata.length,
                                            origin = paste(outputsummary$mig_year, "-01-01", sep = "")))

  # find strata that are below the minimum strata operation threshold.
  ag.strata <- aggregate(ceiling(data$cor.factor), by=list(strata=data$strata), FUN=sum)
  excl.strata <- ag.strata[ag.strata$x < strata.op.min,]
  x.strata <- as.character(excl.strata$strata)

  #reorder output
  outputsummary <- outputsummary[,c(1,12,13,11,6:10,2:5 )]

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
  #############                      Brood Summary Loop                          #########
  #########################################################################################

  smolt.strata <- round(yday(as.Date(smolt.date))/strata.length) - (min(data$strata)-1)# convert to smolt cutoff date to strata
  parr.strata <- round(yday(as.Date(parr.date))/strata.length) - (min(data$strata)-1) # convert parr cutoff date to strata

  for(selectyr in sel.years){

    dir.create(paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood",sep = ""))

    # subset U parameters for parr and presmolt
    pyear <- (selectyr+2) - (as.numeric(min(data$year)))

    usep <- subset(model.fit.gg,grepl("^U", Parameter)) #subset all U's
    usep.p <- subset(usep,grepl(paste(",",pyear,"]$" , sep = ""), Parameter)) #subset U's for first year
    usep.p$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep.p$Parameter) #create strata variable
    usep.p$strata <- as.numeric(substr(usep.p$strata, 2, 5)) #clean strata variable

    # adding the if condition so the code does not run for 2 years before the minimum year in the data
    if(min(data$year) - selectyr < 2){

      parr<- subset(usep.p, strata >= (smolt.strata) & strata <= (parr.strata-1))

      parr<-subset(parr, !(parr$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

      parrdt<-as.data.table(parr)

      parrUdist = numeric(boot); # set aside an empty vector for bootstrap
      for (i in 1:boot){
        parrUdist[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(parrdt$strata))

      }
      parrUdist<-as.data.frame(parrUdist) #change output to dataframe
      parrUdist$parrUdist<-as.numeric(parrUdist$parrUdist) #change output to numeric

      write.table(parrUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                                          selectyr,"_",species,"_",trap.name,"_Parr_U_Dist.txt",sep = ""), sep="\t")

      #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
      parrUoutputsummary <- parrUdist %>%
        dplyr::summarise(
          mode = as.numeric(names(which.max(table(parrUdist)))),
          mean = mean(parrUdist),
          sd = sd(parrUdist),
          naiveSE = sd / sqrt(length(parrUdist)),
          quantile_2.5 = quantile(parrUdist, probs = 0.025),
          quantile_25 = quantile(parrUdist, probs = 0.25),
          quantile_50 = quantile(parrUdist, probs = 0.5),
          quantile_75 = quantile(parrUdist, probs = 0.75),
          quantile_97.5 = quantile(parrUdist, probs = 0.975))

      parrUoutputsummary$parameter<-paste("Parr_U_Brood_",selectyr) # add Parameter variable
      parrUoutputsummary$mig_year<- selectyr+1 #add year variable

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
        dplyr::summarise(
          mode = as.numeric(names(which.max(table(presmoltUdist)))),
          mean = mean(presmoltUdist),
          sd = sd(presmoltUdist),
          naiveSE = sd / sqrt(length(presmoltUdist)),
          quantile_2.5 = quantile(presmoltUdist, probs = 0.025),
          quantile_25 = quantile(presmoltUdist, probs = 0.25),
          quantile_50 = quantile(presmoltUdist, probs = 0.5),
          quantile_75 = quantile(presmoltUdist, probs = 0.75),
          quantile_97.5 = quantile(presmoltUdist, probs = 0.975))

      presmoltUoutputsummary$parameter<-paste("Presmolt_U_Brood_",selectyr) # add Parameter variable
      presmoltUoutputsummary$mig_year<- selectyr+1 #add year variable
    }

    ########## Smolt #############

    # adding the if condition so the code does not run for a year before the maximum year in the data
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
        dplyr::summarise(
          mode = as.numeric(names(which.max(table(smoltUdist)))),
          mean = mean(smoltUdist),
          sd = sd(smoltUdist),
          naiveSE = sd / sqrt(length(smoltUdist)),
          quantile_2.5 = quantile(smoltUdist, probs = 0.025),
          quantile_25 = quantile(smoltUdist, probs = 0.25),
          quantile_50 = quantile(smoltUdist, probs = 0.5),
          quantile_75 = quantile(smoltUdist, probs = 0.75),
          quantile_97.5 = quantile(smoltUdist, probs = 0.975))

      smoltUoutputsummary$parameter<-paste("Smolt_U_Brood_",selectyr) # add Parameter variable
      smoltUoutputsummary$mig_year<- selectyr+2 #add year variable
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

    }

    #Get summary statistics for U bootstrapped distribution
    totUdist<-as.data.frame(totUdist) #change output to dataframe
    totUdist$totUdist<-as.numeric(totUdist$totUdist) #change output to numeric

    totUoutputsummary <- totUdist %>%
      dplyr::summarise(
        mode = as.numeric(names(which.max(table(totUdist)))),
        mean = mean(totUdist),
        sd = sd(totUdist),
        naiveSE = sd / sqrt(length(totUdist)),
        quantile_2.5 = quantile(totUdist, probs = 0.025),
        quantile_25 = quantile(totUdist, probs = 0.25),
        quantile_50 = quantile(totUdist, probs = 0.5),
        quantile_75 = quantile(totUdist, probs = 0.75),
        quantile_97.5 = quantile(totUdist, probs = 0.975))

    totUoutputsummary$parameter<-paste("Total_U_Brood_",selectyr) # add Parameter variable
    totUoutputsummary$mig_year<- NA #add year variable

    ######### Brood strata used #########
    parr1<- subset(usep.p, strata >= (smolt.strata) & strata <= (parr.strata-1))
    presmolt1<- subset(usep.p, strata >= parr.strata)
    smolt1<- subset(usep.s, strata < smolt.strata)

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(data$year) - selectyr > 1){

      brood.strata <- rbind(parr1,presmolt1,smolt1)

      # the code does runs for a year before the maximum year in the data (no smolt info in this case)
    }else{

      brood.strata <- rbind(parr1,presmolt1)

    }

    boutputsummary <- brood.strata %>%
      group_by(Parameter) %>%
      dplyr::summarise(
        mode = as.numeric(names(which.max(table(value)))),
        mean = mean(value),
        sd = sd(value),
        naiveSE = sd / sqrt(length(value)),
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

    boutputsummary$strata_date<-format(as.Date(boutputsummary$strata*strata.length,
                                               origin = paste(boutputsummary$mig_year, "-01-01", sep = "")))

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(data$year) - selectyr > 1){
      # adding the if condition so the code only runs for 2 years before the minimum year in the data
      if(min(data$year) - selectyr == 2){
        broodsummary<-rbind.fill(totUoutputsummary, smoltUoutputsummary, boutputsummary) #merge outputs
        # the code runs for years before the (maximum year-1) and except (min year - 2)
      }else{
        broodsummary<-rbind.fill(totUoutputsummary,parrUoutputsummary,presmoltUoutputsummary, smoltUoutputsummary, boutputsummary) #merge outputs
      }
      # the code does runs for a year before the maximum year in the data (no smolt info in this case)
    }else{
      broodsummary<-rbind.fill(totUoutputsummary,parrUoutputsummary,presmoltUoutputsummary, boutputsummary) #merge outputs
    }

    broodsummary <- broodsummary[,c(10:13,5:9,1:4)]

    data_by_year <- broodsummary

    u_by_year <- broodsummary %>%
      filter(grepl("^U\\[\\d+,\\d+\\]$", parameter))

    # Define the pattern to match
    pattern <- "\\[(\\d+,\\d+)\\]"  # This pattern matches Parameter column format

    # Use regmatches to extract matching parts
    extracted_parts1 <- regmatches(u_by_year$parameter, regexpr(pattern, u_by_year$parameter))

    extracted_parts2 <- regmatches(outputsummary$parameter, regexpr(pattern, outputsummary$parameter))

    rows <- match(extracted_parts1, extracted_parts2)

    p_by_year <- outputsummary[rows,c(3,1:2,4:13)]

    broodsummary <- rbind(broodsummary, p_by_year)

    broodsummary <- broodsummary %>%
      left_join(data %>%
                  select(-"strata_date") %>%
                  dplyr::rename("mig_year" ="year"), by = c("mig_year","strata"))

    #read out files
    options(width = 10000)  # Adjust the width to fit data in .txt for printing with sink

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
                      paste(x.strata, collapse=', '), "\n",
                      "\n"
    )
    )
    print(broodsummary)
    sink()

    write.csv(broodsummary, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Brood/",
                                         selectyr,"_",species,"_",trap.name,"_Brood_Results.csv",sep = ""), sep="\t")

    ##### Figs

    u_by_year$date <- format(as.Date(u_by_year$strata*strata.length, origin = paste(u_by_year$mig_year, "-01-01", sep = "")), format = "%b %d %Y")

    p1 <- ggplot(u_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
      geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#56B4E9") +
      geom_line() +
      geom_point() +
      labs(y = "Abundance (U)", x = NULL, title = paste(selectyr,species,trap.name,"Brood")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
        axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
        axis.title = element_text(size = 14),  # Increase the font size for axis titles
        plot.title = element_text(size = 16)  # Increase the font size for the plot title
      )

    p_by_year$date <- format(as.Date(p_by_year$strata*strata.length, origin = paste(p_by_year$mig_year, "-01-01", sep = "")), format = "%b %d %Y")

    p2 <- ggplot(p_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
      geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#E69F00") +
      geom_line() + geom_point() +
      labs(y = "Capture probability (p)", x = "Strata") +
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
      plot = grid.arrange(p1, p2, ncol = 1),
      width = 10,  # Adjust the width as needed
      height = 6   # Adjust the height as needed
    )

    options(scipen = 0)

  }
}
