Claddis is an in-development R package designed to import cladistic-type data sets (#NEXUS format) into R and perform disparity and rate tests.

Below are links to slides from three conference presentations on the package. (Note: the package has changed a lot since the first two were given!)

[Talk 1](http://www.slideshare.net/graemelloyd/a-new-r-package-for-automating-cladistic-ordination-and-the-effects-of-phylogenetic-signal-on-disparity-measures)

[Talk 2](http://www.slideshare.net/graemelloyd/claddis-a-new-r-package-for-automating-disparity-analyses-based-on-cladistic-datasets)

[Talk 3](http://www.slideshare.net/graemelloyd/new-methodologies-for-the-use-of-cladistictype-matrices-to-measure-morphological-disparity-and-evolutionary-rate)

You can install and load Claddis into R using the following:

```r
# Install the devtools package from CRAN:
install.packages("devtools")

# Load the devtools package into R:
library(devtools)

# Install the Claddis package from github:
install_github("graemetlloyd/Claddis")

# Load the Claddis package into R:
library(Claddis)
```

A brief tutorial for the package is also available [here](http://www.graemetlloyd.com/teaching/RE2014/disparity_and_rates.r).

More will be added here in future to show the basic features of the package.


Citation
========

Claddis will eventually be uploaded to [CRAN](http://cran.r-project.org/) and have a publication associated with it. Until then, please cite:

Lloyd, Graeme T. (2015). "Claddis: an R package for performing disparity and rate analysis on cladistic-type data sets." Online at GitHub, [https://github.com/graemetlloyd/Claddis](https://github.com/graemetlloyd/Claddis). Accessed (access date).
