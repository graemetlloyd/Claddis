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
#' tree$tip.label <- gsub("\n", "", tree$tip.label)
#'
#' # Ages for day 2016 taxa:
#' ages <- matrix(c(269, 269, 263, 263, 265, 265, 265, 265, 257, 257, 259,
#'   259, 258, 258, 260, 260, 257, 257, 257, 257, 256, 256, 259, 259, 260,
#'   260, 253, 253, 257, 257, 257, 257, 268, 268), ncol = 2, byrow = TRUE,
#'   dimnames = list(c("Biarmosuchus_tener", "Hipposaurus_boonstrai",
#'   "Bullacephalus_jacksoni", "Pachydectes_elsi", "Lemurosaurus_pricei",
#'   "Lobalopex_mordax", "Lophorhinus_willodenensis",
#'   "Proburnetia_viatkensis", "Lende_chiweta",
#'   "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098",
#'   "Niuksenitia_sukhonensis", "Ictidorhinus_martinsi", "RC_20",
#'   "Herpetoskylax_hopsoni", "Lycaenodon_longiceps"), c("FAD", "LAD")))
#'
#' # Time-scale tree using equal method:
#' time_tree <- strap::DatePhylo(tree = tree, ages = ages, method = "equal",
#'   rlen = 0.5)
#'
#' # Prune continuous block from day 2016:
#' cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix = day_2016,
#'   blocks2prune = 1)
#'
#' # Run test rates function for two character partitions:
#' test_rates_output <- test_rates(time_tree = time_tree,
#'   cladistic_matrix = cladistic_matrix,
#'   character_partition = list(list(1:34), list(1:17, 18:34)),
#'   time_bins = seq(from = 270, to = 252, length.out = 10))
#'
#' # Plot 2nd (arbitrary two-partition) character partition model:
#' plot_rates_character(test_rates_output = test_rates_output,
#'   model_number = 2)
#'
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
