#' Chronophylomorphospace Plot
#'
#' @description
#'
#' Plots a three-dimensional chronophylomorphospace.
#'
#' @param pcoa_data Principal coordinate data in the format output by  \link{ordinate_cladistic_matrix} that includes a tree and ancestral states.
#' @param x_axis Which ordination axis to plot as the x-axis.
#' @param y_axis Which ordination axis to plot as the y-axis.
#' @param plot_tips Whether or not to plot the tip nodes (defaults to TRUE).
#' @param plot_nodes Whether or not to plot the internal nodes (defaults to TRUE).
#' @param show_names Whether or not to show the taxon nodes (defaults to TRUE).
#' @param plot_branches Whether or not to plot the branches (defaults to TRUE).
#' @param shadow Whether or not to plot a shadow (2D plot) on the bottom face of the 3D plot (defaults to TRUE).
#'
#' @details
#'
#' Creates a manually repositionable three-dimensional (two ordination axes plus time) plot of a phylomorphospace.
#'
#' This function aims to mimic the data visualisation of Sakamoto and Ruta (2012; their Video S1).
#'
#' @author Emma Sherratt \email{emma.sherratt@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Sakamoto, M. and Ruta, M. 2012. Convergence and divergence in the evolution of cat skulls: temporal and spatial patterns of morphological diversity. \emph{PLoS ONE}, \bold{7}, e39752.
#'
#' @examples
#'
#' \dontrun{
#' require(rgl)
#'
#' # Set random seed:
#' set.seed(4)
#'
#' # Generate a random tree for the Michaux 1989 data set:
#' time_tree <- ape::rtree(n = nrow(michaux_1989$matrix_1$matrix))
#'
#' # Set root time so latest tip terminates at the present:
#' time_tree$root.time <- max(diag(x = ape::vcv(phy = time_tree)))
#'
#' # Add taxon names to the tree:
#' time_tree$tip.label <- rownames(x = michaux_1989$matrix_1$matrix)
#'
#' # Perform a phylogenetic Principal Coordinates Analysis:
#' pcoa_data <- ordinate_cladistic_matrix(michaux_1989,
#'   time_tree = time_tree
#' )
#'
#' # Plot a chronophylomorphospace:
#' plot_chronophylomorphospace(pcoa_data)
#' }
#'
#' @export plot_chronophylomorphospace
plot_chronophylomorphospace <- function(pcoa_data, x_axis = 1, y_axis = 2, plot_tips = TRUE, plot_nodes = TRUE, show_names = TRUE, plot_branches = TRUE, shadow = TRUE) {

  # Add top level conditionals to check for a tree etc.

  # Check user has rgl and stop and warn if not:
  if (!requireNamespace("rgl", quietly = TRUE)) {
    stop(paste0(
      "To plot three-dimensional chrono-phylo-morphospaces, please install package rgl",
      "\n install.packages('lattice')"
    ))
  }

  # Default plotting parameters for a 2D morphospace. Need to change node colour for 3D using rgl
  plotting_parameters <- list(tip_colour = "black", tip_symbol = 21, tip_size = 2, node_colour = "white", node_symbol = 21, node_size = 1.25, branch_colour = "black", branch_width = 3, tiplabel_adjustment = c(-.1, -.1), tiplabel_colour = "black", tiplabel_size = 1)

  # Isolate Tree:
  time_tree <- pcoa_data$time_tree

  # Record number of tips:
  n_tips <- ape::Ntip(phy = time_tree)

  # Isolate pcoa axes:
  pcoa_data <- pcoa_data$vectors

  # Little function to set limits for plotting (make it cube-like):
  set_plot_limits <- function(x, scale_factor) {

    # Get range of x:
    x_range <- range(x)

    # Scale range values (i.e., place zero at centre):
    rescaled_x_range <- scale(x_range, scale = FALSE)

    # Return plot limits:
    mean(x_range) + scale_factor * rescaled_x_range
  }

  # Get node ages for z-axis in plotting:
  z_axis <- date_nodes(time_tree = time_tree)

  # Make x label for plot:
  xlab <- paste("PC", x_axis, sep = "")

  # Make y label for plot:
  ylab <- paste("PC", y_axis, sep = "")

  # Make empty plot:
  rgl::plot3d(
    pcoa_data,
    type = "n",
    xlim = set_plot_limits(pcoa_data[, x_axis], scale_factor = 1.5),
    ylim = set_plot_limits(pcoa_data[, y_axis], scale_factor = 1.5),
    zlim = set_plot_limits(z_axis, 0),
    asp = c(1, 1, 0.5),
    xlab = xlab, ylab = ylab,
    zlab = "Time (Ma)",
    rgl::view3d(phi = 90, fov = 30)
  )

  # If requested plots tips
  if (plot_tips) {
    rgl::points3d(
      pcoa_data[1:n_tips, 1],
      pcoa_data[1:n_tips, 2],
      z_axis[1:n_tips],
      col = plotting_parameters$tip_colour,
      size = plotting_parameters$tip_size * 4
    )
  }

  # If requested plot nodes
  if (plot_nodes) {
    rgl::points3d(
      pcoa_data[(n_tips + 1):nrow(pcoa_data), 1],
      pcoa_data[(n_tips + 1):nrow(pcoa_data), 2],
      z_axis[(n_tips + 1):nrow(pcoa_data)],
      col = plotting_parameters$node_colour,
      size = plotting_parameters$node_size * 4
    )
  }

  # If requested plot branches
  if (plot_branches) {
    for (i in 1:nrow(time_tree$edge)) {
      rgl::lines3d(
        pcoa_data[(time_tree$edge[i, ]), 1],
        pcoa_data[(time_tree$edge[i, ]), 2],
        z_axis[(time_tree$edge[i, ])],
        lwd = plotting_parameters$branch_width
      )
    }
  }

  # If requested plot taxa labels
  if (show_names) {
    rgl::text3d(
      pcoa_data[, x_axis],
      pcoa_data[, y_axis],
      z_axis,
      rownames(x = pcoa_data),
      col = plotting_parameters$tiplabel_colour,
      cex = plotting_parameters$tiplabel_size,
      adj = plotting_parameters$tiplabel_adjustment
    )
  }

  # If plotting the shadow of x and y axes at the base of the plot:
  if (shadow) {

    # Plot branches:
    for (i in 1:nrow(time_tree$edge)) {
      rgl::lines3d(
        pcoa_data[(time_tree$edge[i, ]), 1],
        pcoa_data[(time_tree$edge[i, ]), 2],
        time_tree$root.time,
        lwd = 2,
        alpha = 0.5
      )
    }

    # Plot internal nodes:
    rgl::points3d(
      pcoa_data[(n_tips + 1):nrow(pcoa_data), 1],
      pcoa_data[(n_tips + 1):nrow(pcoa_data), 2],
      time_tree$root.time,
      col = plotting_parameters$node_colour,
      size = plotting_parameters$node_size * 4,
      alpha = 0.5
    )

    # Plot tips:
    rgl::points3d(pcoa_data[1:n_tips, 1],
      pcoa_data[1:n_tips, 2],
      time_tree$root.time,
      col = plotting_parameters$tip_colour,
      size = plotting_parameters$tip_size * 4,
      alpha = 0.5
    )
  }
}
