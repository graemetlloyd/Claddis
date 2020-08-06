#' Visualize a rate test time series
#'
#' @description
#'
#' Given the results from a rates test produces a time series visualization for a specific model.
#'
#' @param test_rates_output Rate output from \link{test_rates}.
#' @param model_type The type of model to plot. Must be one of "branch" or "clade".
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
#' # Run test rates function for each clade partition:
#' test_rates_output <- test_rates(
#'   time_tree = time_tree,
#'   cladistic_matrix = cladistic_matrix,
#'   clade_partitions = as.list(x = seq(
#'     from = ape::Ntip(phy = tree) + 1,
#'     to = ape::Ntip(phy = tree) + ape::Nnode(tree), by = 1
#'   )),
#'   branch_partitions = lapply(X = as.list(x = seq(
#'     from = 1,
#'     to = length(x = time_tree$edge.length), by = 1
#'   )), as.list),
#'   time_bins = seq(from = 270, to = 252, length.out = 10)
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
      lapply(X = strsplit(x$Partition, " \\| ")[[1]], function(y) {
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
      lapply(X = strsplit(x$Partition, " \\| ")[[1]], function(y) {
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
  if (model_type == "branch") edge_rates <- lapply(X = test_rates_output$branch_test_results, function(x) x$Rates)

  # If requesting clade rates extract these from output:
  if (model_type == "clade") edge_rates <- lapply(X = test_rates_output$clade_test_results, function(x) x$Rates)

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
