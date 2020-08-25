#' Chronophylomorphospace Plot
#'
#' @description
#'
#' Plots a three-dimensional chronophylomorphospace.
#'
#' @param pcoa_input Principal coordinate data in the format output by  \link{ordinate_cladistic_matrix} that includes a tree and ancestral states.
#' @param x_axis Which ordination axis to plot as the x-axis.
#' @param y_axis Which ordination axis to plot as the y-axis.
#' @param taxon_groups A named list of groups to which taxa are assigned (optional). This is used to plot points or convex hulls in different colours corresponding to each group. As the user names the groups these can represent any grouping of interest (e.g., taxonomic, ecological, temporal, spatial). \link{assign_taxa_to_bins} can automate temporal assignments.
#' @param plot_tips Whether or not to plot the tip nodes (defaults to TRUE).
#' @param plot_nodes Whether or not to plot the internal nodes (defaults to TRUE).
#' @param plot_taxon_names Whether or not to show the taxon nodes (defaults to TRUE).
#' @param plot_edges Whether or not to plot the branches (defaults to TRUE).
#' @param shadow Whether or not to plot a shadow (2D plot) on the bottom face of the 3D plot (defaults to TRUE).
#' @param plot_group_legend Whether or not to add a legend to identify the groups. Only relevant if using \code{"taxon_groups"}.
#' @param group_legend_position Position to plot the group legend. Must be one of \code{bottom_left}, \code{bottom_right}, \code{top_left}, or \code{top_right} (the default).
#' @param palette The palette to use for plotting each element of taxon_groups. See \link[grDevices]{palette}.
#'
#' @details
#'
#' Creates a manually repositionable three-dimensional (two ordination axes plus time) plot of a phylomorphospace.
#'
#' This function aims to mimic the data visualisation of Sakamoto and Ruta (2012; their Video S1).
#'
#' @author Emma Sherratt \email{emma.sherratt@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{assign_taxa_to_bins}, \link{plot_morphospace_stack}, \link{plot_morphospace}, \link{plot_multi_morphospace}, \link{ordinate_cladistic_matrix}
#'
#' @references
#'
#' Sakamoto, M. and Ruta, M. 2012. Convergence and divergence in the evolution of cat skulls: temporal and spatial patterns of morphological diversity. \emph{PLoS ONE}, \bold{7}, e39752.
#'
#' @examples
#'
#' \dontrun{
#' # Require rgl library to use:
#' require(rgl)
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
#' # Prune incomplete taxa from tree:
#' time_tree <- ape::drop.tip(phy = time_tree, tip = c("Lycaenodon_longiceps",
#'   "Niuksenitia_sukhonensis"))
#'
#' # Prune incomplete taxa from cladistic matrix:
#' cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix = day_2016,
#'   taxa2prune = c("Lycaenodon_longiceps", "Niuksenitia_sukhonensis"))
#'
#' # Perform a phylogenetic Principal Coordinates Analysis:
#' pcoa_input <- ordinate_cladistic_matrix(
#'   cladistic_matrix = cladistic_matrix,
#'   time_tree = time_tree
#' )
#'
#' # Define some simple taxon groups for the data as a named list:
#' taxon_groups <- list(nonBurnetiamorpha = c("Biarmosuchus_tener",
#'   "Hipposaurus_boonstrai", "Bullacephalus_jacksoni", "Pachydectes_elsi",
#'   "Niuksenitia_sukhonensis", "Ictidorhinus_martinsi", "RC_20",
#'   "Herpetoskylax_hopsoni"),
#'   Burnetiamorpha = c("Lemurosaurus_pricei", "Lobalopex_mordax",
#'   "Lophorhinus_willodenensis", "Proburnetia_viatkensis", "Lende_chiweta",
#'   "Paraburnetia_sneeubergensis", "Burnetia_mirabilis", "BP_1_7098"))
#'
#' # Plot a chronophylomorphospace:
#' plot_chronophylomorphospace(
#'   pcoa_input = pcoa_input,
#'   taxon_groups = taxon_groups,
#' )
#' }
#'
#' @export plot_chronophylomorphospace
plot_chronophylomorphospace <- function(pcoa_input, x_axis = 1, y_axis = 2, taxon_groups = NULL, plot_tips = TRUE, plot_nodes = TRUE, plot_taxon_names = TRUE, plot_edges = TRUE, shadow = TRUE, plot_group_legend = TRUE, group_legend_position = "top_right", palette = "viridis") {

  # Add top level conditionals to check for a tree etc.

  # Check user has rgl and stop and warn if not:
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop(paste0(
      "To plot three-dimensional chronophylomorphospaces, please install package rgl",
      "\n install.packages('lattice')"
    ))
  }
  
  # If using taxon_groups:
  if (methods::hasArg(name = "taxon_groups")) {
    
    # Find any taxa duplicated across taxon_group:
    duplicated_taxa <- sort(x = unique(x = unname(obj = unlist(x = taxon_groups))[duplicated(x = unname(obj = unlist(x = taxon_groups)))]))
    
    # If these exist stop and warn user:
    if (length(x = duplicated_taxa) > 0) paste0("The following taxa are duplicted in taxon_groups: ", paste(duplicated_taxa, collapse = ", "), ". Taxa can only be in one group.")
  }
  
  # Check group_legend_position is a valid value and stop and warn user if not:
  if (!group_legend_position %in% c("bottom_left", "bottom_right", "top_left", "top_right")) stop("group_legend_position must be one of \"bottom_left\", \"bottom_right\", \"top_left\", or \"top_right\".")
  
  # Set default solid colour for each taxon to black:
  solid_colours <- rep(x = "black", length.out = nrow(x = pcoa_input$vectors))
  names(solid_colours) <- rownames(x = pcoa_input$vectors)
  
  # If using taxon groups (need different colours for each one):
  if (methods::hasArg(name = "taxon_groups")) {
    
    # Build new solid colours with different colours for each group:
    solid_colours <- unlist(x = lapply(X = as.list(1:length(x = taxon_groups)), function(x) {
      
      # Build vector of group colour:
      group_solid_colours <- rep(x = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1)[x], length.out = length(x = taxon_groups[[x]]))
      
      # Add taxon names to colours:
      names(group_solid_colours) <- taxon_groups[[x]]
      
      # Return new group colours vector:
      group_solid_colours
    }))
  }

  # Default plotting parameters for a 2D morphospace. Need to change node colour for 3D using rgl
  plotting_parameters <- list(tip_colour = solid_colours, tip_symbol = 21, tip_size = 2, node_colour = "grey", node_symbol = 21, node_size = 1.5, branch_colour = grDevices::rgb(red = 0, blue = 0, green = 0, alpha = 1), branch_width = 1, tiplabel_adjustment = c(-.1, -.1), tiplabel_colour = "black", tiplabel_size = 1)

  # Isolate Tree:
  time_tree <- pcoa_input$time_tree

  # Record number of tips:
  n_tips <- ape::Ntip(phy = time_tree)

  # Isolate pcoa axes:
  pcoa_input <- pcoa_input$vectors

  # Little function to set limits for plotting (to make it cube-like):
  set_plot_limits <- function(x, scale_factor) {

    # Get range of x:
    x_range <- range(x)

    # Scale range values (i.e., place zero at centre):
    rescaled_x_range <- scale(x_range, scale = FALSE)

    # Return plot limits:
    mean(x = x_range) + scale_factor * rescaled_x_range
  }

  # Get node ages for z-axis in plotting:
  z_axis <- date_nodes(time_tree = time_tree)

  # Make x label for plot:
  xlab <- paste("PC", x_axis, sep = "")

  # Make y label for plot:
  ylab <- paste("PC", y_axis, sep = "")

  # Make empty plot:
  rgl::plot3d(
    x = pcoa_input[, x_axis],
    y = pcoa_input[, y_axis],
    z = z_axis,
    type = "n",
    xlim = set_plot_limits(x = pcoa_input[, x_axis], scale_factor = 1.5),
    ylim = set_plot_limits(x = pcoa_input[, y_axis], scale_factor = 1.5),
    zlim = rev(set_plot_limits(x = z_axis, scale_factor = 1.5)),
    #zlim = rev(x = range(z_axis)),
    asp = c(1, 1, 1),
    xlab = xlab,
    ylab = ylab,
    zlab = "Time (Ma)",
    rgl::view3d(phi = 90, fov = 30)
  )

  # If requested plots tips
  if (plot_tips) {
    rgl::points3d(
      x = pcoa_input[1:n_tips, x_axis],
      y = pcoa_input[1:n_tips, y_axis],
      z = z_axis[1:n_tips],
      col = plotting_parameters$tip_colour[rownames(pcoa_input)],
      size = plotting_parameters$tip_size * 4
    )
  }

  # If requested plot nodes
  if (plot_nodes) {
    rgl::points3d(
      x = pcoa_input[(n_tips + 1):nrow(x = pcoa_input), x_axis],
      y = pcoa_input[(n_tips + 1):nrow(x = pcoa_input), y_axis],
      z = z_axis[(n_tips + 1):nrow(x = pcoa_input)],
      col = plotting_parameters$node_colour,
      size = plotting_parameters$node_size * 4
    )
  }

  # If requested plot branches
  if (plot_edges) {
    for (i in 1:nrow(x = time_tree$edge)) {
      rgl::lines3d(
        x = pcoa_input[(time_tree$edge[i, ]), 1],
        y = pcoa_input[(time_tree$edge[i, ]), 2],
        z = z_axis[(time_tree$edge[i, ])],
        lwd = plotting_parameters$branch_width,
        col = plotting_parameters$branch_colour
      )
    }
  }

  # If requested plot taxa labels
  if (plot_taxon_names) {
    rgl::text3d(
      x = pcoa_input[time_tree$tip.label, x_axis],
      y = pcoa_input[time_tree$tip.label, y_axis],
      z = z_axis,
      texts = time_tree$tip.label,
      col = plotting_parameters$tiplabel_colour,
      cex = plotting_parameters$tiplabel_size,
      adj = plotting_parameters$tiplabel_adjustment
    )
  }

  # If plotting the shadow of x and y axes at the base of the plot:
  if (shadow) {

    # Plot branches:
    for (i in 1:nrow(x = time_tree$edge)) {
      rgl::lines3d(
        x = pcoa_input[(time_tree$edge[i, ]), x_axis],
        y = pcoa_input[(time_tree$edge[i, ]), y_axis],
        z = time_tree$root.time,
        lwd = 1,
        col = grDevices::rgb(red = 0.5, blue = 0.5, green = 0.5, alpha = 0.5)
      )
    }

    # Plot internal nodes:
    rgl::points3d(
      x = pcoa_input[(n_tips + 1):nrow(x = pcoa_input), x_axis],
      y = pcoa_input[(n_tips + 1):nrow(x = pcoa_input), y_axis],
      z = time_tree$root.time,
      size = plotting_parameters$node_size,
      col = grDevices::rgb(red = 0.5, blue = 0.5, green = 0.5, alpha = 0.5)
    )

    # Plot tips:
    rgl::points3d(
      x = pcoa_input[1:n_tips, x_axis],
      y = pcoa_input[1:n_tips, y_axis],
      z = time_tree$root.time,
      size = plotting_parameters$tip_size,
      col = grDevices::rgb(red = 0.5, blue = 0.5, green = 0.5, alpha = 0.5)
    )
  }
  
  # If plotting a group legend:
  if(methods::hasArg(name = "taxon_groups") && plot_group_legend) {
    
    # Add groups legend to plot in requested position:
    if(group_legend_position == "bottom_left") rgl::legend3d("bottomleft", legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white")
    if(group_legend_position == "bottom_right") rgl::legend3d("bottommright", legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white")
    if(group_legend_position == "top_left") rgl::legend3d("topleft", legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white")
    if(group_legend_position == "top_right") rgl::legend3d("topright", legend = names(taxon_groups), fill = grDevices::hcl.colors(n = length(x = taxon_groups), palette = palette, alpha = 1), bg = "white")
  }
}
