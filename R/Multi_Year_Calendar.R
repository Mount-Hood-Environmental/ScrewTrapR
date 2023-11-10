#' @title Multi-year time stratified Bayesian estimator (calendar year summary)

#' @description This function implements a Bayesian model with a hierarchical structure for U and p by strata between years.
#' This function was intended to be used to summarize juvenile steelhead abundances by calendar year.
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
#' @param trap.name character string used for titles and descriptions of reports
#' @param effort.cor expands the number of unmarked fish captured by a sampling effort
#' @param strata.length number of days in strata
#' @param smolt.juv.date "MM-DD" date to partition smolt life stage
#' @param strata.op.min minimum number of years data need to have been collected in a stratum to be included in summary
#' @param den.plot return density plots of MCMC chains
#' @param trace.plot return trace plots of MCMC chains
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
Multi_Year_Calendar <- function(data,
                                effort.cor = FALSE,
                                sel.years = currentyear,
                                strata.op.min = 1,
                                smolt.juv.date = "06-01",
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

  if (max(sel.years) > max(data$year)) {
    stop(paste("ERROR: A selected calendar year is outside of the range of your dataset. Please modify your sel.years().", sep =""))
  }
  if (min(sel.years) < min(data$year)) {
    stop(paste("ERROR: A selected calendar year is outside of the range of your dataset. Please modify your sel.years().", sep =""))
  }

  currentyear <- max(data$year)-2
  s.length <- max(data$days)
  smolt.date <- paste("2010-",smolt.juv.date, sep ="")

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

  main_folder <- paste(species,"_",trap.name,"_MYC_",format(Sys.Date(), "%Y_%m_%d"),sep = "")
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
                   "Multi_Year_Calendar()", "\n",
                   "effort.cor = ",effort.cor, "\n",
                   "sel.years = ", paste(sel.years, collapse = ", "), "\n",
                   "strata.op.min = ", strata.op.min, "\n",
                   "smolt.juv.date = ", smolt.juv.date, "\n",
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
    rename("parameter" = "Parameter",
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

  ###########################################################
  #                      Life Stage                         #
  ###########################################################

  smolt.strata <- round(yday(as.Date(smolt.date))/strata.length) - (min(data$strata)-1) # convert to smolt cutoff date to strata

  selectyr = 2019

  for(selectyr in sel.years){

    dir.create(paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Calendar",sep = ""))

    syear <- (selectyr+1) - (as.numeric(min(data$year)))

    usep <- subset(model.fit.gg,grepl("^U", Parameter)) #subset all U's
    usep <- subset(usep,grepl(paste(",",syear,"]$" , sep = ""), Parameter)) #subset U's for first year
    usep$strata <- sub(".*?U(.*?)(,.*|$)", "\\1", usep$Parameter)
    usep$strata <- as.numeric(substr(usep$strata, 2, 5))

    juv<- subset(usep, strata >= (smolt.strata+1))

    juv<-subset(juv, !(juv$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

    juvdt<-as.data.table(juv)

    juvUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {juvUdist[i] <- sum(juvdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(juvdt$strata))
    }

    juvUdist<-as.data.frame(juvUdist) #change output to dataframe
    juvUdist$juvUdist<-as.numeric(juvUdist$juvUdist) #change output to numeric

    write.table(juvUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Calendar/",
                                       selectyr,"_",species,"_",trap.name,"_Juv_U_Dist.txt",sep = ""), sep="\t")

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    juvUoutputsummary <- juvUdist %>%
      summarise(
        mode = as.numeric(names(which.max(table(juvUdist)))),
        mean = mean(juvUdist),
        sd = sd(juvUdist),
        naiveSE = sd / sqrt(n()),
        quantile_2.5 = quantile(juvUdist, probs = 0.025),
        quantile_25 = quantile(juvUdist, probs = 0.25),
        quantile_50 = quantile(juvUdist, probs = 0.5),
        quantile_75 = quantile(juvUdist, probs = 0.75),
        quantile_97.5 = quantile(juvUdist, probs = 0.975))

    juvUoutputsummary$parameter<-paste("Juv_U_",selectyr) # add Parameter variable
    juvUoutputsummary$mig_year<- selectyr #add year variable

    ################## SMOLT ##########################

    smolt<- subset(usep, strata <= smolt.strata)

    smolt<-subset(smolt, !(smolt$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

    smoltdt<-as.data.table(smolt)

    smoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {smoltUdist[i] <- sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(smoltdt$strata))
    }

    smoltUdist<-as.data.frame(smoltUdist) #change output to dataframe
    smoltUdist$smoltUdist<-as.numeric(smoltUdist$smoltUdist) #change output to numeric

    write.table(juvUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Calendar/",
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

    smoltUoutputsummary$parameter<-paste("Smolt_U_",selectyr) # add Parameter variable
    smoltUoutputsummary$mig_year<- selectyr #add year variable

    # setup bootstrap to randomly draw samples from each distribution in order to obtain total U statistics
    usepdt<-as.data.table(usep) # turn usep into a datable for bootstrap

    usepdt<-subset(usepdt, !(usepdt$strata %in% excl.strata$strata)) # exclude strata under the strata.op.min threshold

    totUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {totUdist[i] <- sum(usepdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(usepdt$strata))
    }

    write.table(totUdist, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Calendar/",
                                       selectyr,"_",species,"_",trap.name,"_Total_U_Dist.txt",sep = ""), sep="\t")

    #Get summary statistics for U bootstrapped distribution
    totUdist<-as.data.frame(totUdist) #change output to dataframe
    totUdist$totUdist<-as.numeric(totUdist$totUdist) #change output to numeric

    #Get summary statistics for U bootstrapped distribution
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

    totUoutputsummary$parameter<-paste("Total_U_",selectyr) # add Parameter variable
    totUoutputsummary$mig_year<- selectyr #add year variable

    year_int <- outputsummary[outputsummary$mig_year == selectyr ,]
    year_int <- year_int[complete.cases(year_int),]
    year_int <- subset(year_int,grepl("^U", parameter)) #subset all U's

    cal.summary<-rbind.fill(totUoutputsummary, smoltUoutputsummary, juvUoutputsummary, year_int) #merge outputs

    cal.summary <- cal.summary[,c(10:13,5:9,1:4)]

    juv1<- subset(usep, strata >= (smolt.strata+1))
    smolt1<- subset(usep, strata <= smolt.strata)

    data_by_year <- outputsummary[outputsummary$mig_year == selectyr,]
    u_by_year <- data_by_year %>%
      filter(grepl("^U\\[\\d+,\\d+\\]$", parameter))
    p_by_year <- data_by_year %>%
      filter(grepl("^p\\[\\d+,\\d+\\]$", parameter))

    cal.summary <- rbind(cal.summary, p_by_year)

    options(width = 10000)  # Adjust the width to fit data in .txt for printing with sink

    cal.summary <- cal.summary %>%
      left_join(data %>%
                  select(-"strata_date") %>%
                  rename("mig_year" ="year"), by = c("mig_year","strata"))

    #read out files
    sink(paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Calendar/",
               selectyr,"_",species,"_",trap.name,"_Calendar_Results.txt",sep = ""), append=FALSE, split=FALSE)

    writeLines(paste( Sys.Date(),"\n",
                      "\n",
                      selectyr, species, trap.name,"\n",
                      "Boot strap iterations =", boot,"\n",
                      "Specified smolt date:", smolt.juv.date, "\n",
                      "\n",
                      "Strata EXCLUDED from results due to <", strata.op.min, "year(s) of operation (e.g. below defined minimum threshold) ","\n",
                      paste(x.strata, collapse=', '), "\n",
                      "\n"
    )
    )
    print(cal.summary)
    sink()

    write.csv(cal.summary, file = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Calendar/",
                                        selectyr,"_",species,"_",trap.name,"_Calendar_Results.csv",sep = ""))
    ##### Figs

    u_by_year$date <- format(as.Date(u_by_year$strata*strata.length, origin = paste(u_by_year$mig_year, "-01-01", sep = "")), format = "%b %d %Y")

    p1 <- ggplot(u_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
      geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#56B4E9") +
      geom_line() +
      geom_point() +
      labs(y = "Abundance (U)", x = NULL, title = paste(selectyr,species,trap.name,"Calendar")) +
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
      filename = paste(main_folder,"/Results/",selectyr,"_",species,"_",trap.name,"_Calendar/",
                       selectyr,"_",species,"_",trap.name,"_Calendar_Results.png",sep = ""),
      plot = grid.arrange(p1, p2, ncol = 1),
      width = 10,  # Adjust the width as needed
      height = 6   # Adjust the height as needed
    )

    options(scipen = 0)
  }
}
