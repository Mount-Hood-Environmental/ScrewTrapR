# ScrewTrapR
This packages contains the code and functions to implement three primary hierarchical Bayesian models (simple, p-spline, and multi-year models) for juvenile salmonid abundance estimation using capture-mark-recapture data collected at rotary screw traps. 

`ScrewTrapR` is an R package for estimating juvenile salmonid abundances using capture-mark-recapture data collected at rotary screw traps using three primary heirarchical Bayesian models; a simple-pooled  model (Lincoln-Peterson model; cite), penalized-spline model (p-spline model; Bonner & Schwarz, 2011), and multi-year model (Oldemeyer et al. 2018). The `ScrewTrapR` contains two functions (`Chinook_Format()` and `Steelhead_Format()`) for formatting juvenile salmonid data pulled from the Idaho Department of Fish and Game's JTrap database into the appropriate capture-mark-recapture format needed for the hierarchical Bayesian model functions. Following, there are several functions to run the various hierarchical Bayesian models and produce different life-stage and cohort summaries dependent on the species (update later once all the options later). 

Additionally, this package is complimentary to the Screw_Trap_ExploreR shiny web application (insert hyperlink later) which allows biologists and practitioners to upload raw screw trap data (Jtrap Data), visualize data, and provide recommendations based on data characteristics which model, stratification, life-stage date designation, etc. are most appropriate for the screw trap data set. The files exported from the Screw_Trap_ExploreR web app allow for seamless implementation of the model functions contained within the `ScrewTrapR` package.

## Getting Started

To install `ScrewTrapR` you can use Hadley Wickham's `devtools` package. To install and load the `devtools` package use:
```
install.packages("devtools")
library(devtools)
```
NOTE: To use `devtools`, you may also have to download and install Rtools (although you shouldn't). The latest version on Rtools can be found at
https://cran.r-project.org/bin/windows/Rtools/

Once `devtools` is successfully installed, use the following to install `ScrewTrapR`:
```
devtools::install_github("Mount-Hood-Environmental/ScrewTrapR")
```
If you are interested in making contributions to `ScrewTrapR`, consider getting a GitHub account, fork this repository, clone to a local directory, modify, and send me a pull request. I can then review any changes and merge.

