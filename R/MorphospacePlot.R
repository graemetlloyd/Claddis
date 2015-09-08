#' Plot Morphopace
#' 
#' Plots a morphospace using the output from MorphMatrix2PCoA.
#' 
#' Input must come from \link{MorphMatrix2PCoA}.
#'
#' @param pcoa_input Text.
#' @param x_axis Text.
#' @param y_axis Text.
#' @param z_axis Text.
#' @param plot_taxon_names Text.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Wills? Brusatte?
#'
#' @keywords principal coordinates
#'
#' @examples
#'
#' # Nothing yet
#'
#' @export MorphospacePlot
MorphospacePlot <- function(pcoa_input, x_axis = 1, y_axis = 2, z_axis = NULL, plot_taxon_names = FALSE) {

# Option to plot names but not points

  scree_values <- apply(pcoa_input$vectors, 2, var) / sum(apply(pcoa_input$vectors, 2, var)) * 100

  x_lab <- paste("PC", x_axis, " (", round(scree_values[x_axis], 2), "% of total variance)", sep = "")
  y_lab <- paste("PC", y_axis, " (", round(scree_values[y_axis], 2), "% of total variance)", sep = "")

  plot(pcoa_input$vectors[, x_axis], pcoa_input$vectors[, y_axis], type="n", bg = "black", xlab = x_lab, ylab = y_lab)

  # Case if no z-axis chosen:
  if(is.null(z_axis)) {

    z_colours <- rep("black", nrow(pcoa_input$vectors))
    
    z_sizes <- rep(1, nrow(pcoa_input$vectors))

  # Case if a z-axis is specified:
  } else {
    
    z_colours <- rep("white", nrow(pcoa_input$vectors))
    
    z_colours[which(pcoa_input$vectors[, z_axis] > 0)] <- "black"
    
    z_sizes <- abs(pcoa_input$vectors[, z_axis]) / max(abs(pcoa_input$vectors[, z_axis])) * 3
    
  }

  points(pcoa_input$vectors[, x_axis], pcoa_input$vectors[, y_axis], pch = 21, bg = z_colours, cex = z_sizes)

  if(plot_taxon_names) {

    x_positions <- rep(2, nrow(pcoa_input$vectors))

    x_positions[which(pcoa_input$vectors[, x_axis] < 0)] <- 4

    text(x = pcoa_input$vectors[, x_axis], y = pcoa_input$vectors[, y_axis], labels = rownames(pcoa_input$vectors), pos = x_positions, cex = 0.7)

  }

  # Case if tree supplied:
  if(!is.null(pcoa_input$tree)) {
    

  # Case if no tree supplied:
  } else {


  }

}

#pcoa_input <- MorphMatrix2PCoA(Michaux1989)
