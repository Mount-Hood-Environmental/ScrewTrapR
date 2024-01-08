#' @title Penalized-spline time stratified Bayesian estimator (calendar summary)

#' @description This function implements Bayesian methods to fit a Penalized-spline (p-spline) through U and a hierarchical structure for p between strata within a year.
#' The p-spline model was developed by Bonner & Schwarz (2011) and the R script to implement the model utilizes the TimeStratPetersenDiagError_fit()function
#' found in the BTSPAS packages.
#' This function was intended to be used to summarize juvenile steelhead abundances by calendar year.
#' At a minimum, input data should have columns "year","strata","m","n","u","days","effort","strata_start","strata_end".
#' See example data for formatting.
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
#' @import BTSPAS

#' @importFrom data.table as.data.table
#' @export
#' @return NULL
#'
#'
Spline_Calendar <- function(data,
                            effort.cor = FALSE,
                            sel.years = currentyear,
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

  #currentyear <- max(data$year)-2
  #s.length <- max(data$days)
  smolt.date <- paste("2010-",smolt.juv.date, sep ="")
  strata_key = data %>%
    select(year,strata,strata_start,strata_end)

  #U updated for effort correction factor
  data$cor.factor <-data$effort/data$days
  data$u.cor <- round(data$u*(1/(data$effort/data$days)))

  main_folder <- paste(species,"_",trap.name,"_",format(Sys.Date(), "%Y_%m_%d"),sep = "")
  dir.create(main_folder)

  #selectyr = 2017

  for(selectyr in sel.years){

    # Filter year
    year_data = data %>%
      filter(year == selectyr)

    # Identify strata the trap wasn't running
    not_op <- year_data %>%
      filter(effort <= 0.5) %>%
      pull(strata)

    # Inputs into function for strata that wasn't operating.
    not_op <- as.numeric(not_op)
    not_op_fixed_p_values <- rep(-10, length(not_op))

    # Whether to use u with effort correction or not
    if(effort.cor == TRUE) {
      u.2=year_data$u.cor
    } else {
      u.2=year_data$u
    }

    # Make the call to fit the model and generate the output files
    model.fit <- TimeStratPetersenDiagError_fit(
      title=species,
      prefix=trap.name,
      time=year_data$strata,
      n1=year_data$n,
      m2=year_data$m,
      u2=u.2,
      bad.n1=not_op,
      bad.m2=not_op,
      bad.u2=not_op,
      n.chains=chains,
      n.iter=iterations,
      n.burnin=burnin,
      n.sims = (iterations-burnin)/thin, # n.sim = "Number of simulated values to keeps for posterior distribution."
      # keep consistent with multi-year thin rate of 100.
      logitP.fixed=not_op,
      logitP.fixed.values=not_op_fixed_p_values,
      save.output.to.files=FALSE)

    auto_files_rm<-c("model.txt","inits3.txt","inits2.txt","inits1.txt", "data.txt","CODAchain1.txt", "CODAchain2.txt",
                     "CODAchain3.txt","codaIndex.txt")

    unlink(auto_files_rm, recursive = FALSE) # Delete automated files from BTSPAS functions

    model.fit.mcmc <- as.mcmc(model.fit) #call model.fit from rjags
    model.fit.gg <- suppressWarnings(ggs(model.fit.mcmc)) #convert model.fit.mcmc to gg output (useable dataframe)

    options(scipen = 999) #change scientific notation

    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/Inputs",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/MCMC_Chains",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/Model_Diagnostics",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/Results",sep = ""))

    #########################
    # Input Data & Function #
    #########################

    options(width = 10000)  # Adjust the width to fit data in .txt for printing with sink

    year_data<-year_data[,c("year","strata","strata_start","strata_end","m","n","u","yoym","yoyn","yoyu","days","effort","cor.factor","u.cor")]

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/Inputs/Input_Data.txt",sep = ""), append=FALSE, split=FALSE)
    print(year_data)
    sink()

    write.csv(year_data, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/Inputs/Input_Data.csv",sep = ""))

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/Inputs/Function_Inputs.txt",sep = ""), append=FALSE, split=FALSE)
    writeLines(paste(Sys.Date(),"\n",
                     "Spline_Calendar()", "\n",
                     "effort.cor = ",effort.cor, "\n",
                     "sel.years = ", paste(sel.years, collapse = ", "), "\n",
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
      write.table(chain, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                      "_Calendar/MCMC_Chains/MCMC_Chains", i,".txt", sep = ""), sep="\t")
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

    nonconv.gd <- nonconv.gd %>%
      filter(!Parameter %in% c("Ntot", "Utot")) %>%
      filter(!is.na(Parameter))

    nc.gd<-as.character(nonconv.gd$Parameter)

    if (length(nc.gd) == 0) {
      nc.gd <- "   *  All Gelman-Rubin Statistics <1.1  *   "
    }

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Calendar/Model_Diagnostics/GR_Statistic.txt",sep = ""), append=FALSE, split=FALSE)

    writeLines(paste(Sys.Date(),"\n",
                     "\n",
                     "Number of chains =", chains,"\n",
                     "Number of iterations per chain =", iterations, "\n",
                     "Burn in =", burnin,"\n",
                     "\n",
                     "********** The Gelman-Rubin Statistic: **********","\n",
                     "\n",
                     "   ***  Flagged parameters with point est. >1.1  ***  ", "\n",
                     "\n",
                     paste(nc.gd, collapse=', ' ),
                     "\n",
                     "\n"))
    print(gd)

    sink()

    if(den.plot == TRUE) {
      ggmcmc(model.fit.gg %>%
               filter(!Parameter %in% c("Ntot", "Utot")), file=paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                                                     "_Calendar/Model_Diagnostics/Density_Plots.pdf",sep = ""), plot="density")
    }
    if(trace.plot == TRUE) {
      ggmcmc(model.fit.gg %>%
               filter(!Parameter %in% c("Ntot", "Utot")), file=paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                                                     "_Calendar/Model_Diagnostics/Trace_Plots.pdf",sep = ""), plot="traceplot")
    }

    #####################
    # Manual Statistics #
    #####################

    usep <- subset(model.fit.gg,grepl("^U\\[", Parameter)) #subset all U's
    psep <- subset(model.fit.gg,grepl("^logitP\\[", Parameter)) #subset all U's
    psep = psep %>%
      mutate(value = 1 / (1 + exp(-value))) %>%
      mutate(Parameter = gsub("logitP\\[(\\d+)\\]", "p[\\1]", Parameter))
    model.fit.gg.renamed = rbind(usep,psep)

    #Get summary statistics for U bootstrapped distribution
    outputsummary <- model.fit.gg.renamed %>%
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
    outputsummary$mig_year = selectyr

    #add strata variable
    outputsummary$strata<-NA
    for(s in 1:length(unique(year_data$strata))){
      outputsummary$strata<-ifelse(grepl(paste(s, sep = ""), outputsummary$Parameter), min(year_data$strata-1)+s, outputsummary$strata)
    }

    # clean up summary stats
    outputsummary <- outputsummary %>%
      dplyr::rename("parameter" = "Parameter") %>%
      mutate(across(where(is.numeric), ~ round(., 3))) %>%
      arrange(mig_year, parameter)

    outputsummary = outputsummary %>%
      left_join(strata_key %>%
                  dplyr::rename(mig_year=year), by = c("mig_year","strata"))

    ###########################################################
    #                      Life Stage                         #
    ###########################################################

    smolt.strata <- data %>%
      filter(strata_start == smolt.juv.date) %>%
      pull(strata) %>%
      unique() %>%
      as.numeric() - (min(data$strata)-1)

    syear <- (selectyr+1) - (as.numeric(min(year_data$year)))

    usep <- usep %>%
      mutate(strata = as.integer(sub("U\\[(\\d+)\\]", "\\1", Parameter)))

    juv<- subset(usep, strata >= (smolt.strata+1))

    juvdt<-as.data.table(juv)

    juvUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {juvUdist[i] <- sum(juvdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(juvdt$strata))
    }

    juvUdist<-as.data.frame(juvUdist) #change output to dataframe
    juvUdist$juvUdist<-as.numeric(juvUdist$juvUdist) #change output to numeric

    write.table(juvUdist, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                       "_Calendar/Results/",selectyr,"_",species,"_",trap.name,"_Juv_U_Dist.txt", sep =""), sep="\t")

    #Get descriptive statistics mode, mean, sd, niaveSE or U bootstrap distribution
    juvUoutputsummary <- juvUdist %>%
      dplyr::summarise(
        mode = as.numeric(names(which.max(table(juvUdist)))),
        mean = mean(juvUdist),
        sd = sd(juvUdist),
        naiveSE = sd / sqrt(length(juvUdist)),
        quantile_2.5 = quantile(juvUdist, probs = 0.025),
        quantile_25 = quantile(juvUdist, probs = 0.25),
        quantile_50 = quantile(juvUdist, probs = 0.5),
        quantile_75 = quantile(juvUdist, probs = 0.75),
        quantile_97.5 = quantile(juvUdist, probs = 0.975))

    juvUoutputsummary$parameter<-paste("Juv_U_",selectyr, sep ="") # add Parameter variable
    juvUoutputsummary$mig_year<- selectyr #add year variable
    juvUoutputsummary$strata_start<- smolt.juv.date #add year variable
    juvUoutputsummary$strata_end<- max(data$strata_end) #add year variable


    ################## SMOLT ##########################

    smolt<- subset(usep, strata <= smolt.strata)

    smoltdt<-as.data.table(smolt)

    smoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {smoltUdist[i] <- sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(smoltdt$strata))
    }

    smoltUdist<-as.data.frame(smoltUdist) #change output to dataframe
    smoltUdist$smoltUdist<-as.numeric(smoltUdist$smoltUdist) #change output to numeric

    write.table(smoltUdist, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                         "_Calendar/Results/",selectyr,"_",species,"_",trap.name,"_Smolt_U_Dist.txt", sep =""), sep="\t")

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

    smoltUoutputsummary$parameter<-paste("Smolt_U_",selectyr, sep ="") # add Parameter variable
    smoltUoutputsummary$mig_year<- selectyr #add year variable
    smoltUoutputsummary$strata_start<- min(data$strata_start) #add year variable
    smoltUoutputsummary$strata_end<- format(as.Date(smolt.juv.date, format = "%m-%d") - 1, "%m-%d") #add year variable

    # setup bootstrap to randomly draw samples from each distribution in order to obtain total U statistics
    usepdt<-as.data.table(usep) # turn usep into a datable for bootstrap

    totUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
    for (i in 1:boot)
    {totUdist[i] <- sum(usepdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(usepdt$strata))
    }

    write.table(totUdist, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                       "_Calendar/Results/",selectyr,"_",species,"_",trap.name,"_Total_U_Dist.txt", sep =""), sep="\t")

    #Get summary statistics for U bootstrapped distribution
    totUdist<-as.data.frame(totUdist) #change output to dataframe
    totUdist$totUdist<-as.numeric(totUdist$totUdist) #change output to numeric

    #Get summary statistics for U bootstrapped distribution
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

    totUoutputsummary$parameter<-paste("Total_U_",selectyr, sep ="") # add Parameter variable
    totUoutputsummary$mig_year<- selectyr #add year variable

    cal.summary<-rbind.fill(totUoutputsummary, smoltUoutputsummary, juvUoutputsummary, outputsummary) #merge outputs

    cal.summary <- cal.summary[,c(10:14,5:9,1:4)]

    cal.summary <- cal.summary %>%
      left_join(year_data %>%
                  select(-c("strata_start", "strata_end")) %>%
                  dplyr::rename("mig_year" ="year"), by = c("mig_year","strata"))

    options(width = 10000)  # Adjust the width to fit data in .txt for printing with sink

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
               "_Calendar/Results/",selectyr,"_",species,"_",trap.name,"_Calendar_Results.txt", sep =""), append=FALSE, split=FALSE)
    writeLines(paste( Sys.Date(),"\n",
                      "\n",
                      selectyr, species, trap.name,"\n",
                      "Specified smolt/juvenile date:", smolt.juv.date, "\n",
                      "Number of chains =", chains,"\n",
                      "Number of iterations per chain =", iterations, "\n",
                      "Burn in =", burnin,"\n",
                      "\n",
                      "* Flagged parameters with Gelman-Rubin statistic >1.1  *  ", "\n",
                      paste(nc.gd, collapse=', ' ),
                      "\n"
    )
    )

    print(as.data.frame(cal.summary))
    sink()

    write.csv(cal.summary, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                        "_Calendar/Results/",selectyr,"_",species,"_",trap.name,"_Calendar_Results.csv", sep =""))

    ##### Figs
    u_by_year <- cal.summary %>%
      filter(grepl("^U\\[", parameter))

    u_by_year$date <- format(as.Date(paste(u_by_year$mig_year,"-",u_by_year$strata_start, sep = "")), format = "%b %d %Y")

    p1 <- ggplot(u_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
      geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#56B4E9") +
      geom_point() +
      labs(y = "Abundance (U)", x = NULL, title = paste(selectyr,species,trap.name,"Calendar")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
        axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
        axis.title = element_text(size = 14),  # Increase the font size for axis titles
        plot.title = element_text(size = 16)  # Increase the font size for the plot title
      )

    p_by_year <- cal.summary %>%
      filter(grepl("^p\\[", parameter))

    p_by_year$date <- format(as.Date(paste(p_by_year$mig_year,"-",p_by_year$strata_start, sep = "")), format = "%b %d %Y")

    p2 <- ggplot(p_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
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
      filename = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                       "_Calendar/Results/",selectyr,"_",species,"_",trap.name,"_Calendar_Figures.png", sep =""),
      plot = grid.arrange(p1, p2, ncol = 1),
      width = 10,  # Adjust the width as needed
      height = 6   # Adjust the height as needed
    )
    options(scipen = 0)
  }
}
