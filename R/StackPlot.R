#' Plot Stacked Ordination Spaces
#' 
#' Plots a stack of ordination spaces representing multiple time-slices.
#' 
#' Default is oldest at bottom to youngest at top.
#'
#' @param ordination_axes A matrix of the ordination axes supplied (rownames should be object names). First column should be values for first axis, second for second axis and so on.
#' @param ages A two-column matrix of the first and last apperance dates for the taxa in the same format supplied to \link{DatePhylo}.
#' @param groups Left as NULL for now.
#' @param time_slices A vector of the boundaries for a series of time slices.
#' @param shear A single value (0 to 1) for the degree of shearing in the stacked ordination spaces.
#' @param x_axis The ordination axis to plot on the x-axis.
#' @param y_axis The ordination axis to plot nn the y-axis.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Emma Sherratt \email{emma.sherratt@@gmail.com}
#'
#' @references
#'
#' Friedman coelacanth paper and others?
#'
#' @keywords principal coordinates
#'
#' @examples
#'
#' # Nothing yet
#'
#' @export StackPlot
StackPlot <- function(ordination_axes, ages, groups = NULL, time_slices, shear = 0.2, x_axis = 1, y_axis = 2) {

# ages is FAD and LAD, so can assume range through for first start, but may want to change later.
# Add spaces between stacks? Could even have this be negative to allow overlapping of highly sheared plots.
# check axes asked for exist in data
# cushioning so points don't plot at very edges of plots.

  # Maybe let user set this later, but not yet:
  plot_cushion <- 0.1

  # Put time slices in order of oldest to youngest:
  time_bins <- sort(time_slices, decreasing = TRUE)

  # Record the number of stackes to plot:
  N_stacks <- length(time_bins) - 1

  # Add geologic time at left with geoscale at some point?

  # Define x-axis label:
  xlab <- paste("PC", x_axis, " (", round((apply(ordination_axes, 2, var) / sum(apply(ordination_axes, 2, var)) * 100)[x_axis], 2), "% of total variance)", sep = "")

  # Define y-axis label:
  ylab <- paste("PC", y_axis, " (", round((apply(ordination_axes, 2, var) / sum(apply(ordination_axes, 2, var)) * 100)[y_axis], 2), "% of total variance)", sep = "")
  
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
    lines(x = c(0, shear * 100, 100, 100 - (shear * 100), 0), y = c((100 / N_stacks) * (i - 1), (100 / N_stacks) * i, (100 / N_stacks) * i, (100 / N_stacks) * (i - 1), (100 / N_stacks) * (i - 1)), col = "black")

    # 0, 0 lines?

    # Tick marks?
    
    
    # Case if groups are specified:
    if(!is.null(groups)) {
      
      # Plot convex hulls
      
      # Plot data points for ith stack:
      points(x = points_to_plot[, 1], y = points_to_plot[, 2], pch = 21, bg = groups[rownames(points_to_plot)], col = groups[rownames(points_to_plot)], cex = 0.5)

    # Case if no groups are specified:
    } else {
        
      # Plot data points for ith stack:
      points(x = points_to_plot[, 1], y = points_to_plot[, 2], pch = 21, bg = "black", cex = 0.5)

    }

    # Plot ages of time slice at bottom left:
    text(x = 0, y = min_y - 2, pos = 4, labels = paste(time_bins[i], "-", time_bins[(i + 1)], " Ma", sep = ""))

  }
  
}

#ordination_axes <- matrix(rnorm(10000), nrow = 100)
#rownames(ordination_axes) <- apply(matrix(sample(LETTERS, 1000, replace = TRUE), nrow = 100), 1, paste, collapse = "")
#ages <- t(apply(matrix(runif(200, 0, 100), ncol = 2), 1, sort, decreasing = TRUE))
#colnames(ages) <- c("FAD", "LAD")
#rownames(ages) <- rownames(ordination_axes)
#time_slices <- seq(0, 100, length.out = 6)
#groups <- sample(x = c("red", "blue"), size = nrow(ordination_axes), replace = TRUE)
#names(groups) <- rownames(ordination_axes)

#StackPlot(ordination_axes, ages, groups = groups, time_slices, shear = 0.2, x_axis = 1, y_axis = 2)

#x <- c(c(seq(0, 100, length.out = 101), seq(0, 100, length.out = 101), seq(0, 100, length.out = 101), seq(0, 100, length.out = 101)), c(rep(20, 101), rep(40, 101), rep(60, 101), rep(80, 101)))
#y <- c(c(rep(20, 101), rep(40, 101), rep(60, 101), rep(80, 101)), c(seq(0, 100, length.out = 101), seq(0, 100, length.out = 101), seq(0, 100, length.out = 101), seq(0, 100, length.out = 101)))
#ordination_axes <- cbind(x, y)
#rownames(ordination_axes) <- apply(matrix(sample(LETTERS, 8 * 8 * 101, replace = TRUE), nrow = 8 * 101), 1, paste, collapse = "")
#ages <- cbind(rep(100, 8 * 101), rep(0, 8 * 101))
#colnames(ages) <- c("FAD", "LAD")
#rownames(ages) <- rownames(ordination_axes)
#time_slices <- seq(0, 100, length.out = 6)
#groups <- sample(x = c("red", "blue"), size = nrow(ordination_axes), replace = TRUE)
#names(groups) <- rownames(ordination_axes)

#StackPlot(ordination_axes, ages, groups = groups, time_slices, shear = 0.2, x_axis = 1, y_axis = 2)
