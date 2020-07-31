#' Visualize a rate test time series
#'
#' @description
#'
#' Given the results from a rates test produces a time series visualization for a specific model.
#'
#' @param RateOutput Rate output from \link{test_rates}.
#' @param ModelNumber The number of the model you wish to visualise from the rate output.
#' @param ... Other options to be passed to \link{geoscalePlot}.
#'
#' @details
#'
#' The raw output from \link{test_rates} can be difficult to interpret without visualization and this function provides a means for doing that when the desired output is a time series (other functions will be added for other types of rate test).
#'
#' The function will only work for a single model, but in practice the user may wish to produce multiple plots in which case they simply need to rn the function multiple times or setup a multipanel window first with \link{layout}, or similar.
#'
#' Plots use the \link{geoscale} package to add geologic time to the x-axis and interested users should consult the documentation tere for a full ist of options (passed via ...) in the function (see example below).
#'
#' Calculated rates (changes per lineage million years) are plotted as filled circles and models are plotted as horizontal lines labelled by rate parameters (lambda_i).
#'
#' @return
#'
#' Nothing is returned, but a plot is produced.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' \dontrun{
#' Ages <- read.table("~/Documents/Packages/Claddis/LungfishTest/ages.txt",
#'   sep =",")
#' Matrix <-
#'   Claddis::read_nexus_matrix("~/Documents/Packages/Claddis/LungfishTest/Lloyd_etal_2012a.nex")
#' Tree <-
#'   ape::read.tree("~/Documents/Packages/Claddis/LungfishTest/Lloyd_etal_2012a.tre")
#' Tree <- Tree[sample(1:100000, 100)]
#' Tree <- lapply(Tree, function(x) strap::DatePhylo(x, Ages, rlen = 2, method = "equal"))
#' class(Tree) <- "multiPhylo"
#' TimeBins <- c(443.8, 358.9, 298.9, 251.9, 201.3, 145.0, 66.0, 0.0)
#' LungfishResults <- Claddis::test_rates(Tree[[1]], Matrix, TimeBins,
#'   TimeBinPartitionsToTest = partition_time_bins(7),
#'   CharacterPartitionsToTest = list(list(1:91), list(Cranial = 1:81,
#'   Postcranial = 82:91)))
#' plot_rates_time(LungfishResults,
#'   which(geiger::aicw(unlist(lapply(LungfishResults$TimeBinResults,
#'   function(x) x$AIC)))[, "delta"] == 0), units = "Period", cex.ts = 1,
#'   cex.age = 1, abbrev = "Period")
#' }
#'
#' @export plot_rates_time
plot_rates_time <- function(RateOutput, ModelNumber, ...) {
  
  # TO DO:
  #
  # - Make points size of evidence (i.e., amount of information) with legend.
  # - Add better example that runs.
  
  # Build vector of time bin midpoints for plotting:
  TimeBinMidpoints <- (RateOutput$TimeBinsUsed[2:length(RateOutput$TimeBinsUsed)] + RateOutput$TimeBinsUsed[1:(length(RateOutput$TimeBinsUsed) - 1)]) / 2
  
  # Get partitions used from results output:
  TimeBinPartitions <- lapply(RateOutput$TimeBinResults, function(x) lapply(strsplit(x$Partition, " \\| ")[[1]], function(y) {if (length(grep("-", y)) > 0) {z <- strsplit(y, split = "-")[[1]]; y <- paste0(z[1]:z[2])}; as.numeric(y)} ))
  
  # Get sampled rates for model:
  TimeRates <- cbind(lapply(TimeBinPartitions[ModelNumber], function(x) do.call(rbind, lapply(x, function(y) {xs <- c(RateOutput$TimeBinsUsed[y[1]], RateOutput$TimeBinsUsed[(y[length(y)] + 1)])})))[[1]], RateOutput$TimeBinResults[[ModelNumber]]$Rates, RateOutput$TimeBinResults[[ModelNumber]]$Rates)
  
  # Create base plot of rates in each time bin with any other requested options paseed as ...:
  geoscale::geoscalePlot(ages = TimeBinMidpoints, data = RateOutput$TimeRates[, "Rate"], age.lim = c(max(RateOutput$TimeBinsUsed), min(RateOutput$TimeBinsUsed)), data.lim = c(0, max(RateOutput$TimeRates[, "Rate"]) * 1.1), pch = 20, cex.pt = 2, label = "Character changes per lineage million years", ...)
  
  # Add lines representing clustering of requested model to plot:
  for(i in 1:nrow(TimeRates)) lines(x = TimeRates[i, 1:2], y = TimeRates[i, 3:4])
  
  # Add model parameters (lambda values) to plot:
  for(i in 1:nrow(TimeRates)) text(x = mean(as.numeric(TimeRates[i, 1:2])), y = as.numeric(TimeRates[i, 3]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3)
  
}

