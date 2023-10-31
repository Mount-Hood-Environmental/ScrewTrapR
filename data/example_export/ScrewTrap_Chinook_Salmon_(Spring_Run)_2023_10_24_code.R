#To install the packages that are required to run the model function 
 #If the packages are already installed, it will skill over the installation and just load them 
 package_load = function(package_list = c()){
                  options(repos = structure(c(cran= 'https://ftp.osuosl.org/pub/cran/' )))

                  if (length(setdiff(package_list, rownames(installed.packages()))) > 0) {
                        install.packages(setdiff(package_list, rownames(installed.packages())))
                  }
                 lapply(package_list, require, character.only=TRUE)

                 } 
 
 #List of the packages need to be installed and/or loaded 
 package_load(package_list = c( 'plyr' , 'lubridate' , 'R2jags' , 'coda' , 'lattice' , 'superdiag' , 'mcmcplots' , 'ggmcmc' , 'data.table' , 'ggplot2' , 'gridExtra' )) 

#Select the formatted data included in the zipped file from the Shiny export 
 data <- read.csv(file.choose()) 

#Model function to get the estimated abundance results
Multi_Year_Cohort(data, 
                     effort.cor =  FALSE , 
                     sel.years = c( 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021 ), 
                     strata.op.min =  1 , 
                     smolt.parr.date =  '07-01' , 
                     parr.presmolt.date =  '09-01' , 
                     species =  'Chinook Salmon (Spring Run)' , 
                     trap.name =  'ScrewTrap' , 
                     den.plot =  TRUE , 
                     trace.plot =  FALSE , 
                     strata.length =  10 , 
                     burnin = 100000, 
                     chains = 3, 
                     iterations = 400000, 
                     thin = 100, 
                     boot = 5000, 
                     model.params = c( 'p' , 'U' , 'etaP1' , 'etaU1' , 'sigmaU' , 'sigmaP' ) )
