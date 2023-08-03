
#                                     INITIAL COMMENTS:

# If you haven't run this code before or haven't installed the packages below in your RStudio,
# remove the "#" from in front of the "install.package(...)" and install them by running the lines.

# install.packages("plyr")
# install.packages("lubridate")
# install.packages("R2jags")
# install.packages("coda")
# install.packages("lattice")
# install.packages("superdiag")
# install.packages("mcmcplots")
# install.packages("ggmcmc")
# install.package("data.table")

# If you haven't installed MYTSBE package before, run the code below and select the "MYTSBE.zip" within the "MYTSBE Tutorial" folder
# If you have installed the MYTSBE, skip the install.package(file.choose()...) line below
install.packages(file.choose(), repos=NULL)

library(MYTSBE) #Load MYTSBE package by running this line

# Specify your working directory (where you want your output files to be exported) by running this line and selecting the file
setwd(choose.dir())
# or by specifying the locatione.g.: setwd("C:/Users/Boldemeyer/Desktop/Screw_trap_output")

rm(list = ls()) # Clear your "R environment"

# Choose your ISS.csv file by running this line and selecting the file
ISS_data <- read.table(file.choose(), header=TRUE, sep=",", fill=TRUE)
# or by specifying the location, e.g. ISS_data <- read.table("C:/Users/Boldemeyer/Desktop/MYTSBE/ISS_Example_Data.csv", header=TRUE, sep=",", fill=TRUE)

# Specify your PTAGIS.csv file by running this line and selecting the file
PTAG_data <- read.table(file.choose(), header=TRUE, sep=",", fill=TRUE)
# or by specifying the location, e.g. PTAG_data <- read.table("C:/Users/Boldemeyer/Desktop/MYTSBE/PTAGIS_Example_Data.csv", header=TRUE, sep=",", fill=TRUE)

# Call the ISSFormat function in the MYTSBE package to merge and reformat ISS_data and PTAG_data
ISSFormat(ISS_data = ISS_data, # Don't change
          PTAG_data = PTAG_data, # Don't change
          strata = 7, # specify the strata length
          species = "Chinook Salmon (Spring Run)", # specify identical name for the species in ISS_data e.g. "Chinook Salmon (Spring Run)"
          trap_name = "Marsh" # specify the name for the trap (could be any name)
          )

# Look in the working directory to see that a file titled "trap_name, species, today's date, formatted data.csv" was created


#                                            TROUBLE SHOOTING data formatting:

  #          1) Check your file directory for caps, spaces, forward slashes, and file format (needs to be .csv)
  #          3) Make sure you didn't alter the "ISS_input_file", "Strata_length", etc. titles (columns in the .csv files that need to have identical headers as example data)
  #          2) Re-run the code from the "rm(list = ls())" line until this point









#                           MYTSBE abundance estimates



# Import the output file into the R session
Input_data <- read.table(file.choose(), header=TRUE, sep=",", fill=TRUE)
# or by specifying the location e.g. Input_data <- read.table("C:/Users/Boldemeyer/Desktop/Screw_trap/Marsh Steelhead 2016-4-20 formatted data.csv", header=TRUE, sep=",", fill=TRUE)

# For steelhead run the MYSTBE_Calendar function

MYTSBE_Calendar(data = Input_data,

                effort.cor = TRUE,
                # This expands the strata "u" by an effort correction.
                # e.g. if the trap was run for 4 out of 7 days and caught 40 unmarked fish, applying effort.correction would make "u" (7/4)*40 = 70

                sel.years = 2008:2010,
                # Select the years you want output. Years must be complete.
                # Input needs to be one number (e.g. sel.years = 2013), a range of numbers (e.g. sel.years = 2010:2014),
                # or a list of numbers (e.g. sel.years = c(2010, 2012, 2014))

                strata.length = 7,
                # Specify the length of the strata in days

                smolt.date = "2012-07-01",
                # Specify the life-stage date to partition smolts & non-smolts for steelheads
                # Input needs to be formatted "YYYY-MM-DD" (e.g. smolt.date = "2012-07-01")

                species = "Steelhead",
                # Specify the species name used for titles of output (can be anything)

                trap.name = "Marsh",
                # Specify the trap name used for titles of output (can be anything)

                rm.bad.GR = TRUE
                # Remove U posterior distributions for strata with Gelman-Brooks test staistics >1.1.
                # Typically only an issue for strata that have only 1 year that has very sparse data
                # Recommended to be "TRUE"
                
)

# NOTE: This model takes upwards of several hours to run. Just hang tight and come back in a while


# For all chinook run the MYTSBE_Cohort function

MYTSBE_Cohort(data = Input_data,

                effort.cor = FALSE,
                # This expands the strata "u" by an effort correction.
                # e.g. if the trap was run for 4 out of 7 days and caught 40 unmarked fish, applying effort.correction would make "u" (7/4)*40 = 70

                sel.years = 2000:2004,
                # Select yearly output. ONLY ABLE TO CALCULATE ESTIMATE FOR COHORTS THAT HAVE LEFT FULL YEARS
                # Input needs to be one number (e.g. sel.years = 2013), a range of numbers (e.g. sel.years = 2010:2014),
                # or a list of numbers (e.g. sel.years = c(2010, 2012, 2014))

                strata.length = 7,
                # Specify the length of the strata in days

                smolt.date = "2012-07-01",
                # Specify the life-stage date to partition smolts
                # Input needs to be formatted "YYYY-MM-DD" (e.g. smolt.date = "2012-07-01")

                parr.date = "2012-09-01",
                # Specify the life-stage date to partition parr
                # Input needs to be formatted "YYYY-MM-DD" (e.g. parr.date = "2012-09-01")

                species = "Chinook",
                # Specify the species name used for titles of output (can be anything)

                trap.name = "Marsh",
                # Specify the trap name used for titles of output (can be anything)

                rm.bad.GR = FALSE
                # Remove U posterior distributions for strata with Gelman-Brooks test staistics >1.1.
                # Typically only an issue for strata that have only 1 year that has very sparse data
                #

            )

# NOTE: This model can take several hours to run. Just hang tight and come back in a while

