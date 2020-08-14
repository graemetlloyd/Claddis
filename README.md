Claddis
=======

[![Build
Status](https://travis-ci.org/graemetlloyd/Claddis.svg?branch=master)](https://travis-ci.org/graemetlloyd/Claddis)
[![cran
version](https://www.r-pkg.org/badges/version/Claddis)](https://cran.r-project.org/package=Claddis)

Claddis is an R package designed to import cladistic-type data sets (#NEXUS format) into R and perform disparity and rate tests.

# Version

Claddis is on [CRAN](https://cran.r-project.org/package=Claddis) (version 0.3.4) but is currently (13/08/20) far behind the latest GitHub version (version 0.6.0). It is recommended that you install from GitHub presently.

Please also consult the CHANGELOG file for all updates (including new functions, features and bug fixes) to Claddis.

# Installation

You can install Claddis in R via CRAN with:

```r
install.packages("Claddis", dependencies = TRUE)
```

Or from GitHub with:

```r
devtools::install_github("graemetlloyd/Claddis")
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

# Citation

The first formal paper describing Claddis was published as Lloyd (2016):

Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. *Biological Journal of the Linnean Society*, **118**, 131-151.

The effects of ancestral state estimation choices on phylomorphospaces was discussed in Lloyd (2018):

Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. *Palaeontology*, **61**, 637-645.
