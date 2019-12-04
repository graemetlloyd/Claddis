#' Plot Stacked Ordination Spaces
#'
#' @description
#'
#' Plots a stack of ordination spaces representing multiple time-slices.
#'
#' @param ordination_axes A matrix of the ordination axes supplied (rownames should be object names). First column should be values for first axis, second for second axis and so on.
#' @param ages A two-column matrix of the first and last apperance dates for the taxa in the same format supplied to \link{DatePhylo}.
#' @param groups A vector of colours (for plotting) with for each object name.
#' @param time_slices A vector of the boundaries for a series of time slices.
#' @param shear A single value (0 to 1) for the degree of shearing in the stacked ordination spaces.
#' @param x_axis The ordination axis to plot on the x-axis.
#' @param y_axis The ordination axis to plot nn the y-axis.
#' @param axis_label The text used to precede the axis number. Here "PC" (for principal components/coordinates) is the assumed default, but the user may wish to use something else like "RW" instead.
#'
#' @details
#'
#' This style of plot is taken from various papers by Michael Foote (Foote 1993; his Figures 2, 4, 6, 8, 10, 12, and 14; Foote 1994; his Figure 2; Foote 1995; his Figure 3; Foote 1999; his Figure 22), and can be seen elsehwere in the litarture (e.g., Friedman and Coates 2006; their Figure 2c).
#'
#' Here multiple ordination (or morpho-) spaces are plotted as a seriours of successive "stacks" representing specific intervals of time. Following geologic conventions the oldest time-slice is plotted at the base and the sequence gets younger towards the top.
#'
#' Note that the user needs to supply three pieces of data: 1) a matrix representing the ordination axes (NB: these can come from any source, they do not have to be from \link{Claddis} functions), 2) a set of ages (first adn last appearances) in the same format as required by the \link{DatePhylo} function in the \link{strap} library, and 3) a vector of ages marking the boundaries of the time-slices.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Emma Sherratt \email{emma.sherratt@@gmail.com}
#'
#' @references
#'
#' Foote, M., 1993. Discordance and concordance between morphological and taxonomic diversity. Paleobiology, 19, 185-204.
#'
#' Foote, M., 1994. Morphological disparity in Ordovician-Devonian crinoids and the early saturation of morphological space. Paleobiology, 20, 320-344.
#'
#' Foote, M., 1995. Morphological diversification of Paleozoic crinoids. Paleobiology, 21, 273-299.
#'
#' Foote, M., 1999. Morphological diversity in the evolutionary radiation of Paleozoic and post-Paleozoic crinoids. Paleobiology, 25, 1-115.
#'
#' Friedman, M. and Coates, M. I., 2006. A newly recognized fossil coelacanth highlights the early morphological diversification of the clade. Proceedings of the Royal Society of London B, 273, 245-250
#'
#' @examples
#'
#' # Create x-values that will form a grid:
#' x <- c(c(seq(0, 100, length.out = 101), seq(0, 100, length.out = 101),
#'   seq(0, 100, length.out = 101), seq(0, 100, length.out = 101)),
#'   c(rep(20, 101), rep(40, 101), rep(60, 101), rep(80, 101)))
#'
#' # Create y-values that will form grid:
#' y <- c(c(rep(20, 101), rep(40, 101), rep(60, 101), rep(80, 101)),
#'   c(seq(0, 100, length.out = 101), seq(0, 100, length.out = 101),
#'   seq(0, 100, length.out = 101), seq(0, 100, length.out = 101)))
#'
#' # Combine x and y values into
#' ordination_axes <- matrix(c(x, y), ncol = 2, dimnames =
#'   list(as.list(apply(matrix(sample(LETTERS, 8 * 8 * 101,
#'   replace = TRUE), nrow = 8 * 101), 1, paste, collapse = "")), NULL))
#'
#' # Assign ages as though taxa range through entire interval (100-0 Ma):
#' ages <- matrix(c(rep(100, 8 * 101), rep(0, 8 * 101)), ncol = 2,
#'   dimnames = list(as.list(rownames(ordination_axes)), as.list(c("FAD",
#'   "LAD"))))
#'
#' # Create five 20 million year time slices:
#' time_slices <- seq(0, 100, length.out = 6)
#'
#' # Plot grid lines to show "shearing" effect is working:
#' StackPlot(ordination_axes = ordination_axes, ages = ages, time_slices = time_slices)
#'
#' # Set random seed:
#' set.seed(17)
#'
#' # Create random values to represent ordination axes:
#' ordination_axes <- matrix(rnorm(10000), nrow = 100, dimnames =
#'   list(as.list(apply(matrix(sample(LETTERS, 8 * 100, replace = TRUE), nrow = 100),
#'   1, paste, collapse = "")), NULL))
#'
#' # Create random first and last appearance dates for objects:
#' ages <- matrix(as.vector(apply(matrix(runif(200, 0, 100), ncol = 2), 1, sort,
#'   decreasing = TRUE)), ncol = 2, byrow=TRUE, dimnames =
#'   list(as.list(rownames(ordination_axes)), as.list(c("FAD", "LAD"))))
#'
#' # Create five 20 million year long time slices:
#' time_slices <- seq(0, 100, length.out = 6)
#'
#' # Define groups for objects at random ("red" and "blue"):
#' groups <- sample(x = c("red", "blue"), size = nrow(ordination_axes), replace = TRUE)
#'
#' # Randomly assign objects to groups:
#' names(groups) <- rownames(ordination_axes)
#'
#' # Make stacked ordination plot with convex hulls for groups:
#' StackPlot(ordination_axes, ages, groups, time_slices)
#'
#' @export StackPlot
StackPlot <- function(ordination_axes, ages, groups = NULL, time_slices, shear = 0.2, x_axis = 1, y_axis = 2, axis_label = "PC") {

# Other options for treatment of age data than just ranges?
# Add spaces between stacks? Could even have this be negative to allow overlapping of highly sheared plots.
# Check axes asked for exist in data (top level conditional).
# Add geologic time at left with geoscale at some point?
# Am assuming PC axes, but could be Relative Warp...

  # Maybe let user set this later:
  plot_cushion <- 0.1

  # Put time slices in order of oldest to youngest:
  time_bins <- sort(time_slices, decreasing = TRUE)

  # Record the number of stackes to plot:
  N_stacks <- length(time_bins) - 1

  # Define x-axis label:
  xlab <- paste(axis_label, x_axis, " (", round((apply(ordination_axes, 2, var) / sum(apply(ordination_axes, 2, var)) * 100)[x_axis], 2), "% of total variance)", sep = "")

  # Define y-axis label:
  ylab <- paste(axis_label, y_axis, " (", round((apply(ordination_axes, 2, var) / sum(apply(ordination_axes, 2, var)) * 100)[y_axis], 2), "% of total variance)", sep = "")
  
  # Set some margins before plotting:
  par(mar = c(0, 0, 0, 0), oma = c(1, 1, 1, 1))
  
  # Create basic (empty) plot)
  plot(0:100, 0:100, type = "n", axes = FALSE, xlab = "", ylab = "")
  
  # Plot x-axis:
  text(x = 50, y = -0.5, pos = 1, srt = 0, labels = xlab)
  
  # Plot y-axis:
  text(x = 100.5, y = 50, pos = 1, srt = 90, labels = ylab)
  
  # For each stack:
  for(i in 1:N_stacks) {
  
    # Find taxa present in ith stack (using a range-through approach:
    taxa_in_bin <- intersect(which(ages[, "LAD"] < time_bins[i]), which(ages[, "FAD"] > time_bins[(i + 1)]))
  
    # Isolate just points present in ith stack:
    points_to_plot <- ordination_axes[taxa_in_bin, c(x_axis, y_axis)]
  
    # Shift x-axis points so minimum is zero:
    points_to_plot[, 1] <- points_to_plot[, 1] - min(ordination_axes[, x_axis])
    
    # Shift y-axis points so minimum is zero:
    points_to_plot[, 2] <- points_to_plot[, 2] - min(ordination_axes[, y_axis])
  
    # Place x-axis points on proportional scale:
    points_to_plot[, 1] <- points_to_plot[, 1] / (max(ordination_axes[, x_axis]) - min(ordination_axes[, x_axis]))
    
    # Place y-axis points on proportional scale:
    points_to_plot[, 2] <- points_to_plot[, 2] / (max(ordination_axes[, y_axis]) - min(ordination_axes[, y_axis]))
  
    # Add plot cushioning:
    points_to_plot <- (points_to_plot * (1 - plot_cushion)) + (plot_cushion / 2)
  
    # Values to add to x_axis values in order to shear them:
    x_shear_additions <- points_to_plot[, 2] * (shear * (100 - shear))

    # Get maximum y-value for plotting:
    max_y <- (100 / N_stacks) * i
    
    # Get minimum y-value for plotting:
    min_y <- (100 / N_stacks) * (i - 1)

    # Update x-axis values to final plotting value:
    points_to_plot[, 1] <- (points_to_plot[, 1] * (100 - (shear * 100))) + x_shear_additions

    # Update y-axis values to final plotting value:
    points_to_plot[, 2] <- (points_to_plot[, 2] * (max_y - min_y)) + min_y

    # Draw bounding lines for ith stack:
    polygon(x = c(0, shear * 100, 100, 100 - (shear * 100), 0), y = c((100 / N_stacks) * (i - 1), (100 / N_stacks) * i, (100 / N_stacks) * i, (100 / N_stacks) * (i - 1), (100 / N_stacks) * (i - 1)), col = "white")

    # 0, 0 lines?

    # Tick marks?
    
    # Case if groups are specified:
    if(!is.null(groups)) {
      
      # For each group:
      for(j in unique(groups)) {
          
        # Get edge points (used to plot convex hull):
        edge_points <- names(which(groups[rownames(points_to_plot)] == j))[chull(x = points_to_plot[which(groups[rownames(points_to_plot)] == j), 1], y = points_to_plot[which(groups[rownames(points_to_plot)] == j), 2])]
        
        # Plot convex hull as polygon:
        polygon(x = points_to_plot[edge_points, 1], y = points_to_plot[edge_points, 2], col = adjustcolor(j, alpha.f = 0.3) , border = 0)
        
      }
  
      # Plot data points for ith stack:
      points(x = points_to_plot[, 1], y = points_to_plot[, 2], pch = 21, bg = groups[rownames(points_to_plot)], col = groups[rownames(points_to_plot)], cex = 0.5)

    # Case if no groups are specified:
    } else {
        
      # Plot data points for ith stack:
      points(x = points_to_plot[, 1], y = points_to_plot[, 2], pch = 21, bg = "black", cex = 0.5)

    }

    # Plot ages of time slice at bottom left:
    text(x = 0, y = min_y - 2, pos = 4, labels = paste(round(time_bins[i], 1), "-", round(time_bins[(i + 1)], 1), " Ma", sep = ""))

  }
  
}

