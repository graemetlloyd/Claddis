#' Visualize a rate test time series
#'
#' @description
#'
#' Given the results from a rates test produces a time series visualization for a specific model.
#'
#' @param test_rates_output Rate output from \link{test_rates}.
#' @param model_number The number of the model you wish to visualise from the rate output.
#' @param ... Other options to be passed to \link{plot}.
#'
#' @details
#'
#' The raw output from \link{test_rates} can be difficult to interpret without visualization and this function provides a means for doing that when the desired output is a time series (other functions will be added for other types of rate test).
#'
#' The function will only work for a single model, but in practice the user may wish to produce multiple plots in which case they simply need to rn the function multiple times or setup a multipanel window first with \link{layout}, or similar.
#'
#' Plots use the \link{geoscale} package to add a geologic time to the x-axis and interested users should consult the documentation there for a full list of options (passed via ...) in the function (see example below).
#'
#' Calculated rates (changes per lineage million years) are plotted as filled circles and models are plotted as horizontal lines labelled by rate parameters (lambda 1, lmabda 2 etc.).
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
#' \donttest{
#' # Make time-scaled first MPT for Day 2016 data set:
#' time_tree <- ape::read.tree(text = paste0("(Biarmosuchus_tener:0.5,",
#'   "(((Hipposaurus_boonstrai:3.5,(Bullacephalus_jacksoni:0.75,",
#'   "Pachydectes_elsi:0.75):0.75):0.75,(Lemurosaurus_pricei:7.166666667,",
#'   "(Lobalopex_mordax:4.333333333,((Lophorhinus_willodenensis:3.666666667,",
#'   "(Proburnetia_viatkensis:0.8333333333,(Lende_chiweta:2,",
#'   "(Paraburnetia_sneeubergensis:1,Burnetia_mirabilis:2):1):1.833333333)",
#'   ":0.8333333333):0.8333333333,(BP_1_7098:2.25,Niuksenitia_sukhonensis:",
#'   "1.25):1.25):0.8333333333):0.8333333333):3.083333333):1.95,",
#'   "(Ictidorhinus_martinsi:15.9,(RC_20:11.6,(Herpetoskylax_hopsoni:11.3,",
#'   "Lycaenodon_longiceps:0.3):0.3):0.3):0.3):0.3);"))
#'
#' # Add root age to tree:
#' time_tree$root.time <- 269.5
#'
#' # Prune continuous block from day 2016:
#' cladistic_matrix <- prune_cladistic_matrix(
#'   cladistic_matrix = day_2016,
#'   blocks2prune = 1
#' )
#'
#' # Run test rates function for two character partitions:
#' test_rates_output <- test_rates(
#'   time_tree = time_tree,
#'   cladistic_matrix = cladistic_matrix,
#'   character_partition = list(list(1:34), list(1:17, 18:34)),
#'   time_bins = seq(from = 270, to = 252, length.out = 10)
#' )
#'
#' # Plot 2nd (arbitrary two-partition) character partition model:
#' plot_rates_character(
#'   test_rates_output = test_rates_output,
#'   model_number = 2
#' )
#' }
#' @export plot_rates_character
plot_rates_character <- function(test_rates_output, model_number, ...) {

  # TO DO:
  #
  # - Add more options (plot colours, whetehr to show character numbers etc.)
  # - Add checks that data are present (characters were tested for)

  # Reconstruct character partitions list:
  character_partitions <- lapply(X = test_rates_output$character_test_results, function(x) {
    lapply(X = strsplit(x$partition, " \\| ")[[1]], function(y) {
      unlist(x = lapply(X = y, function(z) {
        z <- as.list(x = strsplit(z, split = " ")[[1]])
        unlist(x = lapply(X = z, function(p) {
          if (length(x = grep("-", p)) > 0) {
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
  model_x_values <- lapply(X = apply(matrix(c(1, cumsum(unlist(x = lapply(X = character_partitions[[model_number]], function(x) range(1:length(x = x))))[-1])), ncol = 2), 2, list), unlist)

  # Get y values for plotting partitions of model:
  model_y_values <- lapply(X = as.list(x = test_rates_output$character_test_results[[model_number]]$rates), rep, 2)

  # Make vector of partition colours ready for plotting:
  partition_colours <- unlist(x = unname(mapply(function(x, y) {
    rep(x, length(x = y))
  }, x = hcl.colors(n = length(x = character_partitions[[model_number]]), alpha = 0.5, palette = "viridis"), y = character_partitions[[model_number]])))

  # Plot character rates:
  graphics::plot(x = 1:max(unlist(x = character_partitions[[model_number]])), y = test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"], pch = 21, bg = partition_colours, cex = 1.5, xlab = "Character", ylab = "Changes per lineage myr", ylim = c(0, 1.1 * max(test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"])), xaxt = "n", lwd = 0.5, col = "black")

  # Add character numbrs to plot:
  graphics::text(x = 1:max(unlist(x = character_partitions[[model_number]])), y = test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"], label = unlist(x = character_partitions[[model_number]]), pos = 1, col = partition_colours, cex = 0.5)

  # Add lines representing clustering of requested model to plot:
  for (i in 1:length(x = model_x_values)) graphics::lines(x = c(model_x_values[[i]][1] - 0.5, model_x_values[[i]][2] + 0.5), y = model_y_values[[i]])

  # Add model parameters (lambda values) to plot:
  for (i in 1:length(x = model_x_values)) graphics::text(x = mean(model_x_values[[i]]), y = mean(model_y_values[[i]]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3, cex = 1.5)

  # Add legend to plot:
  graphics::legend(x = 0, y = 1.1 * max(test_rates_output$character_rates[unlist(x = character_partitions[[model_number]]), "rate"]), legend = paste0("Partition ", 1:length(x = character_partitions[[model_number]])), pch = rep(21, length(x = character_partitions[[model_number]])), pt.bg = unique(x = partition_colours), col = rep("black", length(x = character_partitions[[model_number]])), pt.lwd = 0.5, pt.cex = 1.5)
}
