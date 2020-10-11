#' Plot Multiple Morphopaces
#'
#' @description
#'
#' Plots multiple morphospaces up to a given number of ordination axes.
#'
#' @param pcoa_input The main input in the format outputted from \link{ordinate_cladistic_matrix}.
#' @param n_axes An integer indicating the total number of axes to plot (should minimally be three).
#' @param taxon_groups See \link{plot_morphospace}.
#' @param plot_taxon_names See \link{plot_morphospace}.
#' @param plot_convex_hulls See \link{plot_morphospace}.
#' @param plot_internal_nodes See \link{plot_morphospace}.
#' @param plot_edges See \link{plot_morphospace}.
#' @param plot_root See \link{plot_morphospace}.
#' @param root_colour See \link{plot_morphospace}.
#' @param palette See \link{plot_morphospace}.
#' @param plot_group_legend See \link{plot_morphospace}.
#'
#' @details
#'
#' Takes the output from \link{ordinate_cladistic_matrix} and uses \link{plot_morphospace} to plot the first N ordination axes.
#'
#' This allows the user a better appreciation of how variance is distributed across multiple axes and all plots are scaled the saem way to further aid visualisation. Data will seem to "shrink" towards the centre of the space on higher axes as variance decreases.
#'
#' Most of the options are simply passed to \link{plot_morphospace}, but the full range is not available as many will be inappropriate here (e.g., adding a z-axis).
#'
#' @author Emma Sherratt \email{emma.sherratt@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{assign_taxa_to_bins}, \link{plot_chronophylomorphospace}, \link{plot_morphospace_stack}, \link{plot_morphospace}, \link{ordinate_cladistic_matrix}
#'
#' @examples
#'
#' # Make PCoA for Day 2016 data set:
#' pcoa_input <- ordinate_cladistic_matrix(cladistic_matrix = day_2016)
#'
#' # Define some simple taxon groups for the data as a named list:
#' taxon_groups <- list(nonBurnetiamorpha = c("Biarmosuchus_tener",
#'    "Hipposaurus_boonstrai", "Bullacephalus_jacksoni", "Pachydectes_elsi",
#'    "Niuksenitia_sukhonensis", "Ictidorhinus_martinsi", "RC_20",
#'    "Herpetoskylax_hopsoni", "Lycaenodon_longiceps"),
#'    Burnetiamorpha = c("Lemurosaurus_pricei", "Lobalopex_mordax",
#'    "Lophorhinus_willodenensis", "Proburnetia_viatkensis", "Lende_chiweta",
#'    "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098"))
#'
#' # Set class as taxonGroups:
#' class(taxon_groups) <- "taxonGroups"
#'
#' # Plot taxon groups including convex hulls:
#' plot_multi_morphospace(pcoa_input, n_axes = 5, taxon_groups = taxon_groups,
#'   plot_convex_hulls = TRUE)
#'
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
#' # Make same plot as before but with a phylogeny:
#' plot_multi_morphospace(
#'   pcoa_input = pcoa_input,
#'   n_axes = 5,
#'   taxon_groups = taxon_groups,
#'   plot_convex_hulls = TRUE
#' )
#'
#' @export plot_multi_morphospace
plot_multi_morphospace <- function(pcoa_input, n_axes = 4, taxon_groups = NULL, plot_taxon_names = FALSE, plot_convex_hulls = FALSE, plot_internal_nodes = FALSE, plot_edges = TRUE, plot_root = TRUE, root_colour = "red", palette = "viridis", plot_group_legend = TRUE) {

  # Add zero lines or grids?
  # Make margins zero and add axes at edges somehow?
  # Conditional that changes N axes to max axes if set higher and warns user
  # Similarly, only 2 axes means plot is kinda pointless!
  # Never plot taxon names!
  
  # If not a valid taxonGroups object then stop and provide feedback to user on what is wrong:
  if (!is.null(x = taxon_groups) && !is.taxonGroups(x = taxon_groups)) stop(check_taxonGroups(taxon_groups = taxon_groups)[1])
  
  # Set x and y limits (to ensure all plots are scaled the same way):
  x_limits <- y_limits <- range(c(pcoa_input$vectors[, 1:n_axes]))

  # Work out the number of plots required:
  n_plots <- (n_axes^2 - n_axes) / 2

  # Set uo matrix that will represent plot layout:
  plot_matrix <- matrix(0, ncol = n_axes, nrow = n_axes)

  # Make plot numbers for lower triangle of layout:
  plot_matrix[which(x = lower.tri(plot_matrix) == TRUE)] <- 1:n_plots

  # Remove last empty column:
  plot_matrix <- plot_matrix[, -n_axes]

  # Add new column of PC labels:
  plot_matrix <- cbind(c(0, (n_plots + 1):(n_plots + (n_axes - 1))), plot_matrix)

  # Set first row of PC labels:
  plot_matrix[1, ] <- c(0, (n_plots + n_axes):(n_plots + n_axes + n_axes - 2))
  
  # Set final panel for legend (if used):
  if (!is.null(taxon_groups) && plot_group_legend) {
    
    # Work out available size for legend:
    legend_size <- sum(diag(x = t(x = apply(X = plot_matrix[2:nrow(x = plot_matrix), 2:ncol(x = plot_matrix)], MARGIN = 2, FUN = rev))) == 0)
    
    # Set legend block in plot_matrix:
    plot_matrix[2:(legend_size + 1), (ncol(x = plot_matrix) - legend_size + 1):ncol(x = plot_matrix)] <- max(x = plot_matrix) + 1
  }

  # Set up margins for morphospace plots:
  graphics::par(mar = c(2, 2, 0, 0))

  # Might want these to reflect actual PC size so that they are in correct relation to each other (NB: PDF will have to be square to retain the aspect ratio)
  graphics::layout(plot_matrix, widths = c(0.2, rep(1, (n_axes - 1))), heights = c(0.2, rep(1, (n_axes - 1))))

  # For each x-axis:
  for (i in 1:(n_axes - 1)) {

    # For each y-axis:
    for (j in (i + 1):n_axes) {
      
      # Plot individual morphospaces using options from input:
      plot_morphospace(pcoa_input = pcoa_input, x_axis = i, y_axis = j, z_axis = NULL, taxon_groups = taxon_groups, plot_taxon_names = plot_taxon_names, plot_convex_hulls = plot_convex_hulls, plot_internal_nodes = plot_internal_nodes, plot_root = plot_root, root_colour = root_colour, palette = palette, plot_group_legend = FALSE, inform = FALSE, x_limits = x_limits, y_limits = y_limits)
    }
  }
  
  # Create PC axis labels:
  labels <- c(paste("PC", 2:n_axes, sep = ""), paste("PC", 1:(n_axes - 1), sep = ""))

  # Set up margins for plotting PC labels:
  graphics::par(mar = c(0, 0, 0, 0))

  # Place PC labels along left
  for (i in 1:(length(x = labels) / 2)) {

    # Empty plot:
    graphics::plot(n_axes, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))

    # Add ordination axis labels:
    graphics::text(x = 0, y = 0, labels = labels[i], cex = 2, srt = 90)
  }

  # Place PC labels along top
  for (i in (((length(x = labels) / 2) + 1):length(x = labels))) {

    # Empty plot:
    graphics::plot(n_axes, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))

    # Empty plot:
    graphics::text(x = 0, y = 0, labels = labels[i], cex = 2)
  }
  
  # If plotting a legend:
  if (!is.null(taxon_groups) && plot_group_legend) {
    
    # Make empty plot:
    graphics::plot(x = c(0, 1), y = c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
    
    # Add legend in same form as plot_morphospace:
    graphics::legend(x = 0.5, y = 0.5, legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white", xjust = 0.5, yjust = 0.5)
  }
  
  # Reset plotting device so layout is not inherited by the next plot the user makes:
  graphics::layout(1)
}
