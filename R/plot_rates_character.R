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
#' # Nothing yet
#' @export plot_rates_character
plot_rates_character <- function(test_rates_output, model_number, ...) {

  # TO DO:
  #
  # - Add more options (plot colours, whetehr to show character numbers etc.)
  # - Add checks that data are present (characters were tested for)

  # Reconstruct character partitions list:
  CharacterPartitions <- lapply(test_rates_output$CharacterPartitionResults, function(x) {
    lapply(strsplit(x$Partition, " \\| ")[[1]], function(y) {
      unlist(lapply(y, function(z) {
        z <- as.list(strsplit(z, split = " ")[[1]])
        unlist(lapply(z, function(p) {
          if (length(grep("-", p)) > 0) {
            p <- strsplit(p, split = "-")[[1]]
            as.numeric(p[1]:as.numeric(p[2]))
          } else {
            as.numeric(p)
          }
        }))
      }))
    })
  })

  # Get x values for plotting partitions of model:
  ModelXValues <- lapply(apply(matrix(c(1, cumsum(unlist(lapply(CharacterPartitions[[model_number]], function(x) range(1:length(x))))[-1])), ncol = 2), 2, list), unlist)

  # Get y values for plotting partitions of model:
  ModelYValues <- lapply(as.list(test_rates_output$CharacterPartitionResults[[model_number]]$Rates), rep, 2)

  # Make vector of partition colours ready for plotting:
  PartitionColours <- unlist(unname(mapply(function(x, y) {
    rep(x, length(y))
  }, x = hcl.colors(n = length(CharacterPartitions[[model_number]]), alpha = 0.5, palette = "viridis"), y = CharacterPartitions[[model_number]])))

  # Plot character rates:
  plot(x = 1:max(unlist(CharacterPartitions[[model_number]])), y = test_rates_output$CharacterRates[unlist(CharacterPartitions[[model_number]]), "Rate"], pch = 21, bg = PartitionColours, cex = 1.5, xlab = "Character", ylab = "Changes per lineage myr", ylim = c(0, 1.1 * max(test_rates_output$CharacterRates[unlist(CharacterPartitions[[model_number]]), "Rate"])), xaxt = "n", lwd = 0.5, col = "black")

  # Add character numbrs to plot:
  text(x = 1:max(unlist(CharacterPartitions[[model_number]])), y = test_rates_output$CharacterRates[unlist(CharacterPartitions[[model_number]]), "Rate"], label = unlist(CharacterPartitions[[model_number]]), pos = 1, col = PartitionColours, cex = 0.5)

  # Add lines representing clustering of requested model to plot:
  for (i in 1:length(ModelXValues)) lines(x = c(ModelXValues[[i]][1] - 0.5, ModelXValues[[i]][2] + 0.5), y = ModelYValues[[i]])

  # Add model parameters (lambda values) to plot:
  for (i in 1:length(ModelXValues)) text(x = mean(ModelXValues[[i]]), y = mean(ModelYValues[[i]]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3, cex = 1.5)

  # Add legend to plot:
  legend(x = 0, y = 1.1 * max(test_rates_output$CharacterRates[unlist(CharacterPartitions[[model_number]]), "Rate"]), legend = paste0("Partition ", 1:length(CharacterPartitions[[model_number]])), pch = rep(21, length(CharacterPartitions[[model_number]])), pt.bg = unique(PartitionColours), col = rep("black", length(CharacterPartitions[[model_number]])), pt.lwd = 0.5, pt.cex = 1.5)
}
