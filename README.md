# Screw_TrapR
This packages contains the code and functions to implement three primary hierarchical Bayesian models (simple, p-spline, and multi-year models) for juvenile salmonid abundance estimation using capture-mark-recapture data collected at rotary screw traps. 

`Screw_TrapR` is an R package for estimating juvenile salmonid abundances using capture-mark-recapture data collected at rotary screw traps using three primary heirarchical Bayesian models; a simple-pooled  model (lincon-peterson model; cite), penalized-spline model (p-spline model; Bonner & Schwarz, 2011), and multi-year model (Oldemeyer et al. 2018). The `Screw_TrapR` contains two functions (`Chinook_Format()` and `Steelhead_Format()`) for formatting juvenile salmonid data pulled from the Idaho Department of Fish and Game's JTrap database into the appropriate capture-mark-recapture format needed for the hierarchical Bayesian model functions. Following, there are several functions to run the various hierarchical Bayesian models and produce different life-stage and cohort summaries dependent on the species (update later all the options later). 

Additionally, this package is complimentary to the Screw_Trap_ExploreR shiny web application (insert hyperlink later) which allows biologists and practitioners to upload raw screw trap data (Jtrap Data), visualize data, and provide recommendations based on data characteristics which model, stratification, life-stage date designation, etc. are most appropriate for the screw trap data set. The files exported from the Screw_Trap_ExploreR web app allow for seamless implementation of the model functions contained within the Screw_TrapR package.

## Getting Started

To install `Screw_TrapR` you can use Hadley Wickham's `devtools` package. To install and load the `devtools` package use:
```
install.packages("devtools")
library(devtools)
```
NOTE: To use `devtools`, you may also have to download and install Rtools (although you shouldn't). The latest version on Rtools can be found at
https://cran.r-project.org/bin/windows/Rtools/

Once `devtools` is successfully installed, use the following to install MYTSBE:
```
devtools::install_github("Mount-Hood-Environmental/Screw_TrapR")
```
If you are interested in making contributions to `Screw_TrapR`, consider getting a GitHub account, fork this repository, clone to a local directory, modify, and send me a pull request. I can then review any changes and merge.

For further information email me or, see:

Oldemeyer, B.N., Copeland, T.S., and B.P. Kennedy. (in review. 2017). A multi-year hierarchical Bayesian mark-recapture model using recurring salmonid behavior to account for sparse or missing data. Canadian Journal of Fisheries and Aquatic Sciences. 

* Will update upon publication
