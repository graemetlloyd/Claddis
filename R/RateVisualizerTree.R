#' Visualize a rate test time series
#'
#' @description
#'
#' Given the results from a rates test produces a time series visualization for a specific model.
#'
#' @param RateOutput Rate output from \link{DiscreteCharacterRate}.
#' @param ModelType The type of model to plot. Must be one of "branch" or "clade".
#' @param ModelNumber The number of the model you wish to visualise from the rate output.
#' @param ... Other options to be passed to \link{geoscalePlot}.
#'
#' @details
#'
#' The raw output from \link{DiscreteCharacterRate} can be difficult to interpret without visualization and this function provides a means for doing that when the desired output is a time series (other functions will be added for other types of rate test).
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
#' # Need example(s)
#'
#' @export RateVisualizerTime
RateVisualizerTree <- function(RateOutput, ModelType, ModelNumber, ...) {
  
  # TO DO:
  #
  # - Check ModelType is a valid choice
  # - Maybe work out how to plot time trees this way too
  # - Add legend for partitions somehow?
  # - Work out how to drop/edit subtitle from legend
  # - Work out how to not round rates to 1dp (0.005 -> 0 which is useless)
  
  # Check model type is a valid option:
  if(!ModelType %in% c("branch", "clade")) stop("ModelType must be one of \"branch\" or \"clade\".")
  
  # Set resolution fro plotting (discretisation of continuous rates):
  Resolution <- 100
  
  # If requesting branch partitions extract these from rate output:
  if(ModelType == "branch") EdgePartitions <- lapply(RateOutput$BranchPartitionResults, function(x) lapply(strsplit(x$Partition, " \\| ")[[1]], function(y) unlist(lapply(y, function(z) {z <- as.list(strsplit(z, split = " ")[[1]]); unlist(lapply(z, function(p) if(length(grep("-", p)) > 0) {p <- strsplit(p, split = "-")[[1]]; as.numeric(p[1]:as.numeric(p[2]))} else {as.numeric(p)}))}))))
  
  # If requesting clade partitions extract these from rate output:
  if(ModelType == "clade") EdgePartitions <- lapply(RateOutput$CladePartitionResults, function(x) lapply(strsplit(x$Partition, " \\| ")[[1]], function(y) unlist(lapply(y, function(z) {z <- as.list(strsplit(z, split = " ")[[1]]); unlist(lapply(z, function(p) if(length(grep("-", p)) > 0) {p <- strsplit(p, split = "-")[[1]]; as.numeric(p[1]:as.numeric(p[2]))} else {as.numeric(p)}))}))))
  
  # If requesting branch rates extract these from output:
  if(ModelType == "branch") EdgeRates <- lapply(RateOutput$BranchPartitionResults, function(x) x$Rates)
  
  # If requesting clade rates extract these from output:
  if(ModelType == "clade") EdgeRates <- lapply(RateOutput$CladePartitionResults, function(x) x$Rates)
  
  # Get discretized vector of edge rates (needed for choosing plot colours):
  DiscretizedRateValues <- seq(from = 0, to = max(EdgeRates[[ModelNumber]]), length.out = Resolution)
  
  # Discretize edge rates:
  DiscretizedEdgeRates <- lapply(EdgeRates, function(x) unlist(lapply(as.list(x), function(y) DiscretizedRateValues[max(which(y >= DiscretizedRateValues))])))
  
  # Create vector of edge rate values to use in plotting:
  EdgeRateValues <- rep(0, nrow(RateOutput$Tree$edge))
  
  # Fill vector of edge rate values to use in plotting:
  for(i in 1:length(EdgePartitions[[ModelNumber]])) EdgeRateValues[EdgePartitions[[ModelNumber]][[i]]] <- DiscretizedRateValues[DiscretizedRateValues == DiscretizedEdgeRates[[ModelNumber]][i]]
  
  # Plot tree with branches colour coded by rate:
  phytools::plotBranchbyTrait(tree = RateOutput$Tree, x = EdgeRateValues, mode = "edge", xlims = c(0, max(EdgeRateValues)), title = "Changes per lineage myr", leg = max(nodeHeights(RateOutput$Tree)), ...)
  
  # Dead code attempting to basically do what phytools::plotBranchbyTrait does without calling phytools:
  #names(DiscretizedRateValues) <- hcl.colors(n = Resolution, palette = "viridis")
  #EdgeColours <- rep("white", nrow(RateOutput$Tree$edge))
  #for(i in 1:length(EdgePartitions[[ModelNumber]])) EdgeColours[EdgePartitions[[ModelNumber]][[i]]] <- names(DiscretizedRateValues[DiscretizedRateValues == DiscretizedEdgeRates[[ModelNumber]][i]])
  #ape::plot.phylo(x = RateOutput$Tree, edge.color = EdgeColours, show.tip.label = FALSE, edge.width = 3)
  #cols = names(DiscretizedRateValues)
  #tree = RateOutput$Tree
  #lwd = 4
  #lims = c(0, max(DiscretizedRateValues))
  #Rounder <- (-1 * (min(c(1, ceiling(log(lims[2], base = 10)))) - 1) + 1)
  #leg <- round(0.5 * max(nodeHeights(tree)), 2)
  #x <- max(nodeHeights(tree)) / 2
  #y <- 10
  #fsize <- 1.0
  #X <- x + cbind(0:(length(cols) - 1) / length(cols), 1:length(cols) / length(cols)) * (leg)
  #Y <- cbind(rep(y, length(cols)), rep(y, length(cols)))
  #lines(c(X[1, 1], X[nrow(X), 2]), c(Y[1, 1], Y[nrow(Y), 2]), lwd = lwd + 2, lend = 2)
  #for(i in 1:length(cols)) lines(X[i, ], Y[i, ], col = cols[i], lwd = lwd, lend = 2)
  #text(x = x, y = y, "0", pos = 3, cex = fsize)
  #text(x= x + leg, y = y, round(lims[2], Rounder), pos = 3, cex = fsize)
  #text(x = (2 * x + leg) / 2, y = y, "Changes", pos = 3, cex = fsize)
  #text(x = (2 * x + leg) / 2, y = y, "per lineage myr", pos = 1, cex = fsize)

}
