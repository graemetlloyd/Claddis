# Claddis

Claddis is an R package designed to import cladistic-type data sets (#NEXUS format) into R and perform disparity analysis and rate tests.

# Package status

### Travis release version [![Build Status](https://travis-ci.org/graemetlloyd/Claddis.svg?branch=CRAN)](https://travis-ci.org/graemetlloyd/Claddis)

### Travis development version [![Build Status](https://travis-ci.org/graemetlloyd/Claddis.svg?branch=master)](https://travis-ci.org/graemetlloyd/Claddis)

### CRAN [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/) [![cran version](https://www.r-pkg.org/badges/version/Claddis)](https://cran.r-project.org/package=Claddis) [![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/Claddis)](https://github.com/r-hub/cranlogs.app) ![](http://cranlogs.r-pkg.org/badges/Claddis)

# Version

Claddis is on [CRAN](https://cran.r-project.org/package=Claddis) (version 0.6.1) but is also being developed on GitHub. To get the absolute latest version you can use:

```r
if(!require(devtools)) install.packages("devtools", dependencies = TRUE)
devtools::install_github("graemetlloyd/Claddis", ref = "master")
```

However, installing a development version of a package is only recommended for expert users.

Please also consult the CHANGELOG file for all updates (including new functions, features and bug fixes) to Claddis.

# Installation

You can install Claddis in R via CRAN with:

```r
install.packages("Claddis", dependencies = TRUE)
```

Or from GitHub with:

```r
if(!require(devtools)) install.packages("devtools", dependencies = TRUE)
devtools::install_github("graemetlloyd/Claddis", ref = "CRAN")
```

And load it into memory using:

```r
library(Claddis)
```

# Help

Basic help can be found with:

```r
?Claddis
```

And clicking on the Index link at the base of the help file will reveal links to every available function.

# Tutorials

Note that I have previously linked to tutorials for the package here, but substantial reworking of the core code means these will no longer work and so currently the example code in each function's help file is the best substitute for this. New tutorials will eventually be produced and shared here.

Users should also be aware of the [dispRity](https://cran.r-project.org/package=dispRity) R package, that can form the end of a Claddis disparity pipeline.

# Citation

The first formal paper describing Claddis was published as Lloyd (2016):

Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. *Biological Journal of the Linnean Society*, **118**, 131-151.

The effects of ancestral state estimation choices on phylomorphospaces was discussed in Lloyd (2018):

Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. *Palaeontology*, **61**, 637-645.
