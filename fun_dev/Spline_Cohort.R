#######################################################################################################
#
# Authors:  *ADD AUTHORS*
# Primary contact: *PRIMARY AUTHOR EMAIL*
# Purpose:  *ADD A BRIEF DESCRIPTION OF THE PURPOSE OF THIS SHINY APP*
# Created:  *DATE CREATED*
# Last Modified: *DATE LAST MODIFIED (NOT SUPER NECESSARRY SINCE THIS IS VERSION CONTROLLED VIA GIT)*
#
#######################################################################################################
#
# Copyright © 2023 Mount Hood Environmental
# This R script is protected under copyright law. All rights reserved. No part of this script may be reproduced,
# distributed, or transmitted in any form or by any means, including photocopying, recording, or other electronic
# or mechanical methods, without the prior written permission of the copyright owner, except in the case of brief 
# quotations embodied in critical reviews and certain other noncommercial uses permitted by copyright law.
# For permissions requests or inquiries, please contact:
# *AUTHOR NAME*
# *AUTHOR EMAIL*
#
##################################################################################################################

##################################
#  Load libraries
##################################

# Some suggested libraries

library(sf)
library(here)
library(tidyverse)
library(magrittr)
library(ggmap)
library(nngeo)
library(janitor)
library(data.table)

##################################
#  Import data or functions
##################################

rm(list=ls()) #Clear global environment first

# setwd(here("FolderName")) # Set working directory using here()

# load data and/or functions

##################################
#  Start coding
##################################

# ~ DELETE COMMENTS BELOW ONCE YOU START CODING BUT JUST AS A REMINDER! ~
#
# - Whenever possible, please try to use pipe operators ( %>% , e.g. tidy/dplyr) in your code instead of 
#   nested functions (base R style).
#
# - Tidy functions > other functions
#
# - ggplot > base R plots (unless for internal data exploration)
#
# - Naming conventions:
#   + Files: underscore_separated, all lower case: e.g. numeric_version
#   + Functions: period.seperated, all lower case: e.g. my.function
#   + Variables: lowerCamelCase: e.g. juvSteelhead
#
# - Add comments and documentation throughout. Comments should try to address the "why" not the "what" for code


