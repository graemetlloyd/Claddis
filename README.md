Claddis is an R package designed to import cladistic-type data sets (#NEXUS format) into R and perform disparity and rate tests.

You can install Claddis in R using:

```r
devtools::install_github("graemetlloyd/Claddis")
```

And load it into memory using:

```r
library(Claddis)
```

Basic help can be found with:

```r
?Claddis
```
And clicking on the Index link at the base of the help file will reveal links to every available function.

Note that I have previously produced tutorials for the package [here](http://www.graemetlloyd.com/teaching/RE2014/disparity_and_rates.r), but that current substantial reworking of the core code means these will no longer function and so the example code in each function's help file is the current best substitute for this.

Citation
========

Claddis is not currently available on [CRAN](https://cran.r-project.org/package=Claddis), but version (0.3) can be installed directly from github using the code above. The first formal paper describing Claddis is now published (Lloyd 2016), with further discussion in Lloyd (2018):

Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. *Biological Journal of the Linnean Society*, **118**, 131-151.

Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. *Palaeontology*, **61**, 637-645.