#######################################################################################################
#
# Authors:  Bryce Oldemeyer
# Primary contact: Bryce.Oldemeyer@mthoodenvironmental.com
# Purpose:  Development of Bayesian P-spline model for ScrewTrapR
# Created:  8/2/2023
# Last Modified:
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


