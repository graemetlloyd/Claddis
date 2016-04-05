Claddis is an in-development R package designed to import cladistic-type data sets (#NEXUS format) into R and perform disparity and rate tests.

Below are links to slides from three conference presentations on the package. (Note: the package has changed a lot since the first two were given!)

[Talk 1](http://www.slideshare.net/graemelloyd/a-new-r-package-for-automating-cladistic-ordination-and-the-effects-of-phylogenetic-signal-on-disparity-measures)

[Talk 2](http://www.slideshare.net/graemelloyd/claddis-a-new-r-package-for-automating-disparity-analyses-based-on-cladistic-datasets)

[Talk 3](http://www.slideshare.net/graemelloyd/new-methodologies-for-the-use-of-cladistictype-matrices-to-measure-morphological-disparity-and-evolutionary-rate)

You can install and load Claddis into R using the following:

```r
# Install the devtools package from CRAN:
install.packages("devtools")

# Install the phytools and Claddis packages from github:
devtools::install_github("liamrevell/phytools")
devtools::install_github("graemetlloyd/Claddis")

# Load the Claddis package into R:
library(Claddis)
```

A brief tutorial for the package is also available [here](http://www.graemetlloyd.com/teaching/RE2014/disparity_and_rates.r).

More will be added here in future to show the basic features of the package.


Citation
========

Claddis is now on [CRAN](https://cran.r-project.org/package=Claddis), but the current version (0.1) is broken due to changes in a dependent package (phytools). Users are advised to install directly from github (version 0.2) using the code above. The first formal paper describing Claddis is now published in *Biological Journal of the Linnean Society* and can be cited as:

Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. *Biological Journal of the Linnean Society*, **118**, 131-151.
