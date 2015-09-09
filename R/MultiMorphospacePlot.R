#' Plot Multiple Morphopaces
#' 
#' Plots morphospaces for any number of axes using the output from MorphMatrix2PCoA.
#' 
#' Input must come from \link{MorphMatrix2PCoA}. Is a wrapper for \link{MorphospacePlot}.
#'
#' @param pcoa_input Text.
#' @param N_axes Text.
#' @param plot_taxon_names Text.
#' @param plot_internal_nodes Text.
#' @param plot_root Text.
#' @param root_colour Text.
#'
#' @author Emma Sherratt \email{emma.sherratt@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords principal coordinates
#'
#' @examples
#'
#' # Nothing yet
#'
#' @export MultiMorphospacePlot
MultiMorphospacePlot <- function(pcoa_input, N_axes = 4, plot_taxon_names = FALSE, plot_internal_nodes = FALSE, plot_root = TRUE, root_colour = "grey") {

# Add zero lines or grids?
# Make margins zero and add axes at edges somehow?
# Conditional that changes N axes to max axes if set higher and warns user
# Similarly, only 2 axes means plot is kinda pointless!
# Never plot taxon names!

  # Work out the number of plots required:
  N_plots <- (N_axes ^ 2 - N_axes) / 2

  # Set uo matrix that will represent plot layout:
  A <- matrix(0, ncol = N_axes, nrow = N_axes)

  # Make plot numbers for lower triangle of layout:
  A[which(lower.tri(A) == TRUE)] <- 1:N_plots

  # Remove last empty column:
  A <- A[, -N_axes]

  # Add new column of PC labels:
  A <- cbind(c(0, (N_plots + 1):(N_plots + (N_axes - 1))), A)

  # Set first row of PC labels:
  A[1, ] <- c(0, (N_plots + N_axes):(N_plots + N_axes + N_axes - 2))

  # Set up margins for morphospace plots:
  par(mar = c(2, 2, 0, 0))

  # Might want these to reflect actual PC size so that they are in correct relation to each other (NB: PDF will have to be square to retain the aspect ratio)
  layout(A, widths = c(0.2, rep(1, (N_axes - 1))), heights = c(0.2, rep(1, (N_axes - 1))))

  # For each x-axis:
  for (i in 1:(N_axes - 1)) {
    
    # For each y-axis:
    for (j in (i + 1):N_axes) {
        
# CHANGE THIS TO EXISTING MORPHOSPACE FUNCTION:
#plot(pcoa_input$vectors[, i], pcoa_input$vectors[, j], pch = 21, bg = "black", xlab = "", ylab = "", asp = TRUE)
      MorphospacePlot(pcoa_input, x_axis = i, y_axis = j, z_axis = NULL, plot_taxon_names = plot_taxon_names, plot_internal_nodes = plot_internal_nodes, plot_root = plot_root, root_colour = root_colour)
    
    }
    
  }

  # Create PC axis labels:
  labels <- c(paste("PC", 2:N_axes, sep = "") , paste("PC", 1:(N_axes - 1), sep = ""))

  # Set up margins for plotting PC labels:
  par(mar = c(0, 0, 0, 0))

  # Place PC labels along left
  for(i in 1:(length(labels) / 2)) {
    
    # Empty plot:
    plot(N_axes, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
  
    # Add ordination axis labels:
    text(x = 0, y = 0, labels = labels[i], cex = 2, srt = 90)
    
  }

  # Place PC labels along top
  for(i in (((length(labels) / 2) + 1):length(labels))) {
    
    # Empty plot:
    plot(N_axes, type = "n", axes = FALSE, xlab = "", ylab = "", xlim = c(-1, 1), ylim = c(-1, 1))
  
    # Empty plot:
    text(x = 0, y = 0, labels = labels[i], cex = 2)
  
  }

  # Reset plotting device so layout is not inherited by the next plot the user makes:
  layout(1)

}
