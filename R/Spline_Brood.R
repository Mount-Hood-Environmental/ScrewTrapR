#' @title Penalized-spline time stratified Bayesian estimator (brood year summary)

#' @description This function implements Bayesian methods to fit a Penalized-spline (p-spline) through U and a hierarchical structure for p between strata within a year.
#' The p-spline model was developed by Bonner & Schwarz (2011) and the R script to implement the model utilizes the TimeStratPetersenDiagError_fit()function
#' found in the BTSPAS packages.
#' The Spline_Brood() function was intended to be used to summarize juvenile Chinook Salmon abundances by brood year.
#' e.g. Summaries for brood year 2012 consist of parr (calendar year 2013), presmolts (calendar year 2013), smolts (calendar year 2014).
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
Spline_Brood <- function(data,
                         effort.cor = FALSE,
                         sel.years = currentyear,
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
  smolt.strata.cut <- round(yday(as.Date(smolt.date))/strata.length)# convert to smolt cutoff date to strata
  parr.strata.cut <- round(yday(as.Date(parr.date))/strata.length) # convert parr cutoff date to strata

  #U updated for effort correction factor
  data$cor.factor <-data$effort/data$days
  data$u.cor <- round(data$u*(1/(data$effort/data$days)))

  data$strata_date<-format(as.Date(data$strata*strata.length,
                                   origin = paste(data$year, "-01-01", sep = "")))

  main_folder <- paste(species,"_",trap.name,"_SB_",format(Sys.Date(), "%Y_%m_%d"),sep = "")
  dir.create(main_folder)

  selectyr = 2014

  for(selectyr in sel.years){

    # Filter out just the brood year data. Going to run the spline between years with a break in the spline
    # using jump.after() function
    brood_data_smolt = data %>%
      filter(year == (selectyr+2) &
               strata < smolt.strata.cut)

    brood_data_pp = data %>%
      filter(year == (selectyr+1) &
               strata >= smolt.strata.cut)

    brood_data = rbind(brood_data_pp,brood_data_smolt) %>%
      mutate(model_strata = c(1:nrow(.)))

    # Identify strata the trap wasn't running
    not_op <- brood_data %>%
      filter(effort <= 0.5) %>%
      pull(model_strata)

    # Inputs into function for strata that wasn't operating.
    not_op <- as.numeric(not_op)
    not_op_fixed_p_values <- rep(-10, length(not_op))

    # Whether to use u with effort correction or not
    if(effort.cor == TRUE) {
      u.2=brood_data$u.cor
    } else {
      u.2=brood_data$u
    }

    if (selectyr == min(data$year) - 2) {
      jump_after_presmolt <- NULL
    } else if (selectyr == max(data$year) - 1) {
      jump_after_presmolt <- NULL
    } else {
      jump_after_presmolt <- nrow(brood_data_pp)
    }

    # Make the call to fit the model and generate the output files
    model.fit <- TimeStratPetersenDiagError_fit(
      title=species,
      prefix=trap.name,
      time=brood_data$model_strata,
      n1=brood_data$n,
      m2=brood_data$m,
      u2=u.2,
      bad.n1=not_op,
      bad.m2=not_op,
      bad.u2=not_op,
      jump.after = jump_after_presmolt,
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

    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/Inputs",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/MCMC_Chains",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/Model_Diagnostics",sep = ""))
    dir.create(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/Results",sep = ""))

    #########################
    # Input Data & Function #
    #########################

    options(width = 10000)  # Adjust the width to fit data in .txt for printing with sink

    brood_data<-brood_data[,c(2,3,14,15,4:13)]

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/Inputs/Input_Data.txt",sep = ""), append=FALSE, split=FALSE)
    print(brood_data)
    sink()

    write.csv(brood_data, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/Inputs/Input_Data.csv",sep = ""))

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/Inputs/Function_Inputs.txt",sep = ""), append=FALSE, split=FALSE)
    writeLines(paste(Sys.Date(),"\n",
                     "Spline_Brood()", "\n",
                     "effort.cor = ",effort.cor, "\n",
                     "sel.years = ", paste(sel.years, collapse = ", "), "\n",
                     "smolt.parr.date = ", smolt.parr.date, "\n",
                     "parr.presmolt.date = ", parr.presmolt.date,"\n",
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
      write.table(chain, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/MCMC_Chains/MCMC_Chains", i,".txt", sep = ""), sep="\t")
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

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,"_Brood/Model_Diagnostics/GR_Statistic.txt",sep = ""), append=FALSE, split=FALSE)

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
      ggmcmc(model.fit.gg %>%
               filter(!Parameter %in% c("Ntot", "Utot")), file=paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                                                     "_Brood/Model_Diagnostics/Density_Plots.pdf",sep = ""), plot="density")
    }
    if(trace.plot == TRUE) {
      ggmcmc(model.fit.gg %>%
               filter(!Parameter %in% c("Ntot", "Utot")), file=paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                                                     "_Brood/Model_Diagnostics/Trace_Plots.pdf",sep = ""), plot="traceplot")
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

    #add strata variable
    outputsummary$model_strata<-NA
    for(s in 1:length(unique(brood_data$model_strata))){
      outputsummary$model_strata<-ifelse(grepl(paste(s, sep = ""), outputsummary$Parameter), min(brood_data$model_strata-1)+s, outputsummary$model_strata)
    }

    outputsummary = outputsummary %>%
      left_join(brood_data %>%
                  select(strata,model_strata,year), by = "model_strata")

    # clean up summary stats
    outputsummary <- outputsummary %>%
      rename("parameter" = "Parameter",
             "mig_year" = "year") %>%
      mutate(across(where(is.numeric), ~ round(., 3)))

    outputsummary$strata_date<-format(as.Date(outputsummary$strata*strata.length,
                                              origin = paste(outputsummary$mig_year, "-01-01", sep = "")))

    ###########################################################
    #                      Life Stage                         #
    ###########################################################

    smolt.strata <- round(yday(as.Date(smolt.date))/strata.length) - (min(data$strata)-1)# convert to smolt cutoff date to strata
    parr.strata <- round(yday(as.Date(parr.date))/strata.length) - (min(data$strata)-1) # convert parr cutoff date to strata

    usep <- usep %>%
      mutate(model_strata = as.integer(sub("U\\[(\\d+)\\]", "\\1", Parameter))) %>%
      left_join(brood_data %>%
                  select(strata,model_strata), by = "model_strata")

    if(min(brood_data$year) - selectyr < 2){

      parr <- subset(usep, strata >= (smolt.strata.cut) & strata <= (parr.strata.cut-1))

      parrdt<-as.data.table(parr)

      parrUdist = numeric(boot); # set aside an empty vector for bootstrap
      for (i in 1:boot){
        parrUdist[i] <- sum(parrdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(parrdt$strata))

      }
      parrUdist<-as.data.frame(parrUdist) #change output to dataframe
      parrUdist$parrUdist<-as.numeric(parrUdist$parrUdist) #change output to numeric

      write.table(parrUdist, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                          "_Brood/Results/",selectyr,"_",species,"_",trap.name,"_Parr_U_Dist.txt", sep =""), sep="\t")

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

      parrUoutputsummary$parameter<-paste("Parr_U_Brood_",selectyr, sep = "") # add Parameter variable
      parrUoutputsummary$mig_year<- selectyr+1 #add year variable

      #########Presmolt###############

      presmolt<- subset(usep, strata >= parr.strata.cut)

      presmoltdt<-as.data.table(presmolt)

      presmoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        presmoltUdist[i] <- sum(presmoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(presmoltdt$strata))
      }

      presmoltUdist<-as.data.frame(presmoltUdist) #change output to dataframe
      presmoltUdist$presmoltUdist<-as.numeric(presmoltUdist$presmoltUdist) #change output to numeric

      write.table(presmoltUdist, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                              "_Brood/Results/",selectyr,"_",species,"_",trap.name,"_Presmolt_U_Dist.txt", sep =""), sep="\t")

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

      presmoltUoutputsummary$parameter<-paste("Presmolt_U_Brood_",selectyr, sep = "") # add Parameter variable
      presmoltUoutputsummary$mig_year<- selectyr+1 #add year variable
    }

    ########## Smolt #############

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(brood_data$year) - selectyr > 1){

      smolt<- subset(usep, strata < smolt.strata.cut)

      smoltdt<-as.data.table(smolt)

      smoltUdist = numeric(boot); # set aside an empty vector for bootstrap disribution of the mean
      for (i in 1:boot){
        smoltUdist[i] <- sum(smoltdt[, value[sample.int(.N, 1, TRUE)], by = strata]) - sum(unique(smoltdt$strata))
      }

      smoltUdist<-as.data.frame(smoltUdist) #change output to dataframe
      smoltUdist$smoltUdist<-as.numeric(smoltUdist$smoltUdist) #change output to numeric

      write.table(smoltUdist, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                           "_Brood/Results/",selectyr,"_",species,"_",trap.name,"_Smolt_U_Dist.txt", sep =""), sep="\t")

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

      smoltUoutputsummary$parameter<-paste("Smolt_U_Brood_",selectyr, sep = "") # add Parameter variable
      smoltUoutputsummary$mig_year<- selectyr+2 #add year variable
    }

    #setup bootstrap to randomly draw samples from each distribution in order to obtain total U statistics

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(brood_data$year) - selectyr > 1){

      # adding the if condition so the code only runs for 2 years before the minimum year in the data
      if(min(brood_data$year) - selectyr == 2){

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

    totUoutputsummary$parameter<-paste("Total_U_Brood_",selectyr, sep = "") # add Parameter variable
    totUoutputsummary$mig_year<- NA #add year variable

    # adding the if condition so the code does not run for a year before the maximum year in the data
    if(max(data$year) - selectyr > 1){
      # adding the if condition so the code only runs for 2 years before the minimum year in the data
      if(min(data$year) - selectyr == 2){
        broodsummary<-rbind.fill(totUoutputsummary, smoltUoutputsummary, outputsummary) #merge outputs
        # the code runs for years before the (maximum year-1) and except (min year - 2)
      }else{
        broodsummary<-rbind.fill(totUoutputsummary,parrUoutputsummary,presmoltUoutputsummary, smoltUoutputsummary, outputsummary) #merge outputs
      }
      # the code does runs for a year before the maximum year in the data (no smolt info in this case)
    }else{
      broodsummary<-rbind.fill(totUoutputsummary,parrUoutputsummary,presmoltUoutputsummary, outputsummary) #merge outputs
    }

    broodsummary <- broodsummary[,c(10:14,5:9,1:4)]

    broodsummary <- broodsummary %>%
      left_join(brood_data %>%
                  select(-c("strata_date","model_strata")) %>%
                  rename("mig_year" ="year"), by = c("mig_year","strata"))

    #read out files
    options(width = 10000)  # Adjust the width to fit data in .txt for printing with sink

    sink(paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
               "_Brood/Results/",selectyr,"_",species,"_",trap.name,"_Brood_Results.txt", sep =""), append=FALSE, split=FALSE)
    writeLines(paste( Sys.Date(),"\n",
                      "\n",
                      selectyr, species, trap.name,"\n",
                      "Boot strap iterations =", boot,"\n",
                      "Specified smolt date:", smolt.parr.date, "\n",
                      "Specified parr date:", parr.presmolt.date, "\n",
                      "Burn in =", burnin,"\n",
                      "\n",
                      "* Flagged parameters with Gelman-Rubin statistic >1.1  *  ", "\n",
                      paste(nc.gd, collapse=', ' ),
                      "\n"
    )
    )
    print(broodsummary)
    sink()

    write.csv(broodsummary, file = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                                         "_Brood/Results/",selectyr,"_",species,"_",trap.name,"_Brood_Results.csv", sep =""))

    ##### Figs

    u_by_year <- broodsummary %>%
      filter(grepl("^U\\[", parameter))

    u_by_year$date <- format(as.Date(u_by_year$strata*strata.length, origin = paste(u_by_year$mig_year, "-01-01", sep = "")), format = "%b %d %Y")

    p1 <- ggplot(u_by_year, aes(x = reorder(reorder(date, strata), mig_year), y = quantile_50)) +
      geom_errorbar(aes(ymin = quantile_2.5, ymax = quantile_97.5), width = 0.1, col = "#56B4E9") +
      geom_line() +
      geom_point() +
      labs(y = "Abundance (U)", x = NULL,  title = paste(selectyr,species,trap.name,"Brood")) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increase the font size for x-axis text
        axis.text.y = element_text(size = 12),  # Increase the font size for y-axis text
        axis.title = element_text(size = 14),  # Increase the font size for axis titles
        plot.title = element_text(size = 16)  # Increase the font size for the plot title
      )

    p1

    p_by_year <- broodsummary %>%
      filter(grepl("^p\\[", parameter))

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

    p2

    ggsave(
      filename = paste(main_folder,"/",selectyr,"_",species,"_",trap.name,
                       "_Brood/Results/",selectyr,"_",species,"_",trap.name,"_Brood_Figures.png", sep =""),
      plot = grid.arrange(p1, p2, ncol = 1),
      width = 10,  # Adjust the width as needed
      height = 6   # Adjust the height as needed
    )

    options(scipen = 0)

  }
}
