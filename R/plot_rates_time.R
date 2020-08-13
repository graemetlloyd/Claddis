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
#' # Run test rates function for each time bin partition:
#' test_rates_output <- test_rates(
#'   time_tree = time_tree,
#'   cladistic_matrix = cladistic_matrix,
#'   time_partitions = partition_time_bins(n_time_bins = 9),
#'   time_bins = seq(from = 270, to = 252, length.out = 10)
#' )
#'
#' # Plot 97th time bin partition model:
#' plot_rates_time(
#'   test_rates_output = test_rates_output,
#'   model_number = 97, units = "Stage", cex.ts = 1, cex.age = 1,
#'   abbrev = "Stage"
#' )
#' }
#' @export plot_rates_time
plot_rates_time <- function(test_rates_output, model_number, ...) {

  # TO DO:
  #
  # - Make points size of evidence (i.e., amount of information) with legend.
  # - Add better example that runs.

  # Build vector of time bin midpoints for plotting:
  time_bin_midpoints <- (test_rates_output$time_bins_used[2:length(x = test_rates_output$time_bins_used)] + test_rates_output$time_bins_used[1:(length(x = test_rates_output$time_bins_used) - 1)]) / 2

  # Get partitions used from results output:
  time_bin_partitions <- lapply(X = test_rates_output$time_test_results, function(x) {
    lapply(X = strsplit(x$partition, " \\| ")[[1]], function(y) {
      if (length(x = grep("-", y)) > 0) {
        z <- strsplit(y, split = "-")[[1]]
        y <- paste0(z[1]:z[2])
      }
      as.numeric(y)
    })
  })

  # Get sampled rates for model:
  time_rates <- cbind(lapply(X = time_bin_partitions[model_number], function(x) {
    do.call(what = rbind, args = lapply(X = x, function(y) {
      xs <- c(test_rates_output$time_bins_used[y[1]], test_rates_output$time_bins_used[(y[length(x = y)] + 1)])
    }))
  })[[1]], test_rates_output$time_test_results[[model_number]]$rates, test_rates_output$time_test_results[[model_number]]$rates)

  # Create base plot of rates in each time bin with any other requested options paseed as ...:
  geoscale::geoscalePlot(ages = time_bin_midpoints, data = test_rates_output$time_rates[, "rate"], age.lim = c(max(test_rates_output$time_bins_used), min(test_rates_output$time_bins_used)), data.lim = c(0, max(test_rates_output$time_rates[, "rate"]) * 1.1), pch = 20, cex.pt = 2, label = "Character changes per lineage million years", ...)

  # Add lines representing clustering of requested model to plot:
  for (i in 1:nrow(time_rates)) lines(x = time_rates[i, 1:2], y = time_rates[i, 3:4])

  # Add model parameters (lambda values) to plot:
  for (i in 1:nrow(time_rates)) text(x = mean(as.numeric(time_rates[i, 1:2])), y = as.numeric(time_rates[i, 3]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3)
}
