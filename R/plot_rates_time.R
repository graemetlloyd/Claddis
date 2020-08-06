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
#' # Get first of three MPTs for day 2016 data set:
#' tree <- ape::read.tree(text = "(Biarmosuchus_tener,
#'   (((Hipposaurus_boonstrai,(Bullacephalus_jacksoni,
#'   Pachydectes_elsi)),(Lemurosaurus_pricei,(Lobalopex_mordax,
#'   ((Lophorhinus_willodenensis,(Proburnetia_viatkensis,
#'   (Lende_chiweta,(Paraburnetia_sneeubergensis,
#'   Burnetia_mirabilis)))),(BP_1_7098,
#'   Niuksenitia_sukhonensis))))),(Ictidorhinus_martinsi,
#'   (RC_20,(Herpetoskylax_hopsoni,Lycaenodon_longiceps)))));")
#'
#' # Remove line breaks from tip names:
#' tree$tip.label <- gsub(pattern = "\n", replacement = "", x = tree$tip.label)
#'
#' # Ages for day 2016 taxa:
#' ages <- matrix(c(
#'   269, 269, 263, 263, 265, 265, 265, 265, 257, 257, 259,
#'   259, 258, 258, 260, 260, 257, 257, 257, 257, 256, 256, 259, 259, 260,
#'   260, 253, 253, 257, 257, 257, 257, 268, 268
#' ),
#' ncol = 2, byrow = TRUE,
#' dimnames = list(c(
#'   "Biarmosuchus_tener", "Hipposaurus_boonstrai",
#'   "Bullacephalus_jacksoni", "Pachydectes_elsi", "Lemurosaurus_pricei",
#'   "Lobalopex_mordax", "Lophorhinus_willodenensis",
#'   "Proburnetia_viatkensis", "Lende_chiweta",
#'   "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098",
#'   "Niuksenitia_sukhonensis", "Ictidorhinus_martinsi", "RC_20",
#'   "Herpetoskylax_hopsoni", "Lycaenodon_longiceps"
#' ), c("FAD", "LAD"))
#' )
#'
#' # Time-scale tree using equal method:
#' time_tree <- strap::DatePhylo(
#'   tree = tree, ages = ages, method = "equal",
#'   rlen = 0.5
#' )
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
    lapply(X = strsplit(x$Partition, " \\| ")[[1]], function(y) {
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
  })[[1]], test_rates_output$time_test_results[[model_number]]$Rates, test_rates_output$time_test_results[[model_number]]$Rates)

  # Create base plot of rates in each time bin with any other requested options paseed as ...:
  geoscale::geoscalePlot(ages = time_bin_midpoints, data = test_rates_output$time_rates[, "Rate"], age.lim = c(max(test_rates_output$time_bins_used), min(test_rates_output$time_bins_used)), data.lim = c(0, max(test_rates_output$time_rates[, "Rate"]) * 1.1), pch = 20, cex.pt = 2, label = "Character changes per lineage million years", ...)

  # Add lines representing clustering of requested model to plot:
  for (i in 1:nrow(time_rates)) lines(x = time_rates[i, 1:2], y = time_rates[i, 3:4])

  # Add model parameters (lambda values) to plot:
  for (i in 1:nrow(time_rates)) text(x = mean(as.numeric(time_rates[i, 1:2])), y = as.numeric(time_rates[i, 3]), labels = eval(parse(text = paste0("expression(lambda[", i, "])"))), pos = 3)
}
