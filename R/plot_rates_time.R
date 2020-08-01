#' Visualize a rate test time series
#'
#' @description
#'
#' Given the results from a rates test produces a time series visualization for a specific model.
#'
#' @param test_rates_output Rate output from \link{test_rates}.
#' @param model_number The number of the model you wish to visualise from the rate output.
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
#' time_tree <-
#'   ape::read.tree("~/Documents/Packages/Claddis/LungfishTest/Lloyd_etal_2012a.tre")
#' time_tree <- time_tree[sample(1:100000, 100)]
#' time_tree <- lapply(time_tree, function(x) strap::DatePhylo(x, Ages, rlen = 2, method = "equal"))
#' class(time_tree) <- "multiPhylo"
#' time_bins <- c(443.8, 358.9, 298.9, 251.9, 201.3, 145.0, 66.0, 0.0)
#' LungfishResults <- Claddis::test_rates(time_tree[[1]], Matrix, time_bins,
#'   time_partitions = partition_time_bins(7),
#'   character_partitions = list(list(1:91), list(Cranial = 1:81,
#'   Postcranial = 82:91)))
#' plot_rates_time(LungfishResults,
#'   which(geiger::aicw(unlist(lapply(LungfishResults$TimeBinResults,
#'   function(x) x$AIC)))[, "delta"] == 0), units = "Period", cex.ts = 1,
#'   cex.age = 1, abbrev = "Period")
#' }
#'
#' @export plot_rates_time
plot_rates_time <- function(test_rates_output, model_number, ...) {
  
  # TO DO:
  #
  # - Make points size of evidence (i.e., amount of information) with legend.
  # - Add better example that runs.
  
  # Build vector of time bin midpoints for plotting:
  TimeBinMidpoints <- (test_rates_output$time_binsUsed[2:length(test_rates_output$time_binsUsed)] + test_rates_output$time_binsUsed[1:(length(test_rates_output$time_binsUsed) - 1)]) / 2
  
  # Get partitions used from results output:
  TimeBinPartitions <- lapply(test_rates_output$TimeBinResults, function(x) lapply(strsplit(x$Partition, " \\| ")[[1]], function(y) {if (length(grep("-", y)) > 0) {z <- strsplit(y, split = "-")[[1]]; y <- paste0(z[1]:z[2])}; as.numeric(y)} ))
  
  # Get sampled rates for model:
  TimeRates <- cbind(lapply(TimeBinPartitions[model_number], function(x) do.call(rbind, lapply(x, function(y) {xs <- c(test_rates_output$time_binsUsed[y[1]], test_rates_output$time_binsUsed[(y[length(y)] + 1)])})))[[1]], test_rates_output$TimeBinResults[[model_number]]$Rates, test_rates_output$TimeBinResults[[model_number]]$Rates)
  
  # Create base plot of rates in each time bin with any other requested options paseed as ...:
  geoscale::geoscalePlot(ages = TimeBinMidpoints, data = test_rates_output$TimeRates[, "Rate"], age.lim = c(max(test_rates_output$time_binsUsed), min(test_rates_output$time_binsUsed)), data.lim = c(0, max(test_rates_output$TimeRates[, "Rate"]) * 1.1), pch = 20, cex.pt = 2, label = "Character changes per lineage million years", ...)
  
  # Add lines representing clustering of requested model to plot:
  for(i in 1:nrow(TimeRates)) lines(x = TimeRates[i, 1:2], y = TimeRates[i, 3:4])
  
  # Add model parameters (lambda values) to plot:
  for(i in 1:nrow(TimeRates)) text(x = mean(as.numeric(TimeRates[i, 1:2])), y = as.numeric(TimeRates[i, 3]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3)
  
}

