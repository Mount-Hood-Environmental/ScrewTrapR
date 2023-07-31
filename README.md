# ScrewTrapR

`ScrewTrapR` is a package for estimating juvenile salmonid abundances using capture-mark-recapture data collected at rotary screw traps using three primary hierarchical Bayesian models; a simple-pooled model, penalized-spline model (Bonner & Schwarz, 2011), and multi-year model (Oldemeyer et al. 2018). `ScrewTrapR` contains two functions (`Chinook_Format()` and `Steelhead_Format()`) for formatting raw juvenile salmonid data from the Idaho Department of Fish and Game's JTRAP database into the appropriate capture-mark-recapture format needed for the hierarchical Bayesian model functions. Following, there are several functions that implement the various hierarchical Bayesian models and produce different life-stage and cohort summaries dependent on the species (update later once all the options later).

Additionally, this package is complimentary to the Screw_Trap_ExploreR shiny web application (insert hyperlink later) which allows biologists and practitioners to upload raw rotary screw trap data (JTRAP Data), visualize data, and provide recommendations which model, stratification, life-stage date designation, etc. are most appropriate based on data characteristics. The files exported from the Screw_Trap_ExploreR web app allow for seamless implementation of the model functions contained within the `ScrewTrapR` package.

## Getting Started

To install ScrewTrapR, you first need to install the latest version of [Just Another Gibbs Sampler (JAGS)](https://mcmc-jags.sourceforge.io/). Next, make sure your R version is compatible with the version of JAGS. The README file associated with JAGS when this package was initially developed (July, 2023), states:

> If you are using R to interface with JAGS then you must ensure that you download the correct binary:
>
> If you are using R 4.2.0 or later then install JAGS-4.3.2.exe
>
> If you are using R 4.1.3 or earlier then install JAGS-4.3.0.exe

You can check your R version by simply running:

```         
R.version.string
```

Last, if you are using an Integrated Development Environment (IDE) for your coding in R (e.g. RStudio), be sure to check the IDE is using the intended version of R. You can check (and change) the selected R version in RStudio by going to the "Tools-\>Global Options-\>General-\>R Version" box.

To install `ScrewTrapR` you can use Hadley Wickham's `devtools` package. To install and load the `devtools` package use:

```         
install.packages("devtools")
library(devtools)
```

NOTE: To use `devtools`, you may also have to download and install Rtools (although you shouldn't). The latest version on Rtools can be found at <https://cran.r-project.org/bin/windows/Rtools/>

Once `devtools` is successfully installed, use the following to install `ScrewTrapR`:

```         
devtools::install_github("Mount-Hood-Environmental/ScrewTrapR")
```

If you are interested in making contributions to `ScrewTrapR`, consider getting a GitHub account, fork this repository, clone to a local directory, modify, and send me a pull request. I can then review any changes and merge.
