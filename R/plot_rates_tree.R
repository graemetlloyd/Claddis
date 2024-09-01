#' Visualize a rate test time series
#'
#' @description
#'
#' Given the results from a rates test produces a time series visualization for a specific model.
#'
#' @param test_rates_output Rate output from \link{test_rates}.
#' @param model_type The type of model to plot. Must be one of "branch" or "clade".
#' @param model_number The number of the model you wish to visualise from the rate output.
#' @param ... Other options to be passed to \link{plot}.
#'
#' @details
#'
#' The raw output from \link{test_rates} can be difficult to interpret without visualization and this function provides a means for doing that when the desired output is a time series (other functions will be added for other types of rate test).
#'
#' The function will only work for a single model, but in practice the user may wish to produce multiple plots in which case they simply need to rn the function multiple times or setup a multipanel window first with \link{layout}, or similar.
#'
#' Plots use the \link[geoscale]{geoscale} package to add geologic time to the x-axis and interested users should consult the documentation tere for a full ist of options (passed via ...) in the function (see example below).
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
#' # Generate nine two million year time bins:
#' time_bins <- matrix(data = c(seq(from = 270, to = 252, length.out = 10)[1:9],
#'   seq(from = 270, to = 252, length.out = 10)[2:10]), ncol = 2,
#'   dimnames = list(LETTERS[1:9], c("fad", "lad")))
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Run test rates function for each clade partition:
#' test_rates_output <- test_rates(
#'   time_tree = time_tree,
#'   cladistic_matrix = cladistic_matrix,
#'   clade_partitions = as.list(x = seq(
#'     from = ape::Ntip(phy = time_tree) + 1,
#'     to = ape::Ntip(phy = time_tree) + ape::Nnode(time_tree), by = 1
#'   )),
#'   branch_partitions = lapply(X = as.list(x = seq(
#'     from = 1,
#'     to = length(x = time_tree$edge.length), by = 1
#'   )), as.list),
#'   time_bins = time_bins
#' )
#'
#' # Plot ninth branch partition model (lowest AIC value):
#' plot_rates_tree(
#'   test_rates_output = test_rates_output,
#'   model_type = "branch", model_number = 9
#' )
#'
#' # Plot third clade partition model (lowest AIC value):
#' plot_rates_tree(
#'   test_rates_output = test_rates_output,
#'   model_type = "clade", model_number = 3
#' )
#' }
#' @export plot_rates_tree
plot_rates_tree <- function(test_rates_output, model_type, model_number, ...) {

  # TO DO:
  #
  # - Check model_type is a valid choice
  # - Maybe work out how to plot time trees this way too
  # - Add legend for partitions somehow?
  # - Work out how to drop/edit subtitle from legend
  # - Work out how to not round rates to 1dp (0.005 -> 0 which is useless)
  # - Make example run faster

  # Check model type is a valid option:
  if (!model_type %in% c("branch", "clade")) stop("model_type must be one of \"branch\" or \"clade\".")

  # Set resolution for plotting (discretisation of continuous rates):
  resolution <- 100

  # If requesting branch partitions extract these from rate output:
  if (model_type == "branch") {
    edge_partitions <- lapply(X = test_rates_output$branch_test_results, function(x) {
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
  }

  # If requesting clade partitions extract these from rate output:
  if (model_type == "clade") {
    edge_partitions <- lapply(X = test_rates_output$clade_test_results, function(x) {
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
  }

  # If requesting branch rates extract these from output:
  if (model_type == "branch") edge_rates <- lapply(X = test_rates_output$branch_test_results, function(x) x$rates)

  # If requesting clade rates extract these from output:
  if (model_type == "clade") edge_rates <- lapply(X = test_rates_output$clade_test_results, function(x) x$rates)

  # Get discretized vector of edge rates (needed for choosing plot colours):
  discretized_rate_values <- seq(from = 0, to = max(edge_rates[[model_number]]), length.out = resolution)

  # Discretize edge rates:
  discretized_edge_rates <- lapply(X = edge_rates, function(x) unlist(x = lapply(X = as.list(x = x), function(y) discretized_rate_values[max(which(x = y >= discretized_rate_values))])))

  # Create vector of edge rate values to use in plotting:
  edge_rate_values <- rep(0, nrow(test_rates_output$time_tree$edge))

  # Fill vector of edge rate values to use in plotting:
  for (i in 1:length(x = edge_partitions[[model_number]])) edge_rate_values[edge_partitions[[model_number]][[i]]] <- discretized_rate_values[discretized_rate_values == discretized_edge_rates[[model_number]][i]]

  # Plot tree with branches colour coded by rate:
  phytools::plotBranchbyTrait(tree = test_rates_output$time_tree, x = edge_rate_values, mode = "edge", xlims = c(0, max(edge_rate_values)), title = "Changes per lineage myr", leg = max(nodeHeights(test_rates_output$time_tree)), ...)

  # Dead code attempting to basically do what phytools::plotBranchbyTrait does without calling phytools:
  # names(discretized_rate_values) <- hcl.colors(n = resolution, palette = "viridis")
  # EdgeColours <- rep("white", nrow(test_rates_output$time_tree$edge))
  # for(i in 1:length(x = edge_partitions[[model_number]])) EdgeColours[edge_partitions[[model_number]][[i]]] <- names(discretized_rate_values[discretized_rate_values == discretized_edge_rates[[model_number]][i]])
  # ape::plot.phylo(x = test_rates_output$time_tree, edge.color = EdgeColours, show.tip.label = FALSE, edge.width = 3)
  # cols = names(discretized_rate_values)
  # tree = test_rates_output$time_tree
  # lwd = 4
  # lims = c(0, max(discretized_rate_values))
  # Rounder <- (-1 * (min(c(1, ceiling(log(lims[2], base = 10)))) - 1) + 1)
  # leg <- round(0.5 * max(nodeHeights(tree)), 2)
  # x <- max(nodeHeights(tree)) / 2
  # y <- 10
  # fsize <- 1.0
  # X <- x + cbind(0:(length(x = cols) - 1) / length(x = cols), 1:length(x = cols) / length(x = cols)) * (leg)
  # Y <- cbind(rep(y, length(x = cols)), rep(y, length(x = cols)))
  # lines(c(X[1, 1], X[nrow(X), 2]), c(Y[1, 1], Y[nrow(Y), 2]), lwd = lwd + 2, lend = 2)
  # for(i in 1:length(x = cols)) lines(X[i, ], Y[i, ], col = cols[i], lwd = lwd, lend = 2)
  # text(x = x, y = y, "0", pos = 3, cex = fsize)
  # text(x= x + leg, y = y, round(lims[2], Rounder), pos = 3, cex = fsize)
  # text(x = (2 * x + leg) / 2, y = y, "Changes", pos = 3, cex = fsize)
  # text(x = (2 * x + leg) / 2, y = y, "per lineage myr", pos = 1, cex = fsize)
}
