#' Plot Morphopace
#'
#' @description
#'
#' Plots a morphospace using the output from MorphMatrix2PCoA.
#'
#' @param pcoa_input The main input in the format output from \link{MorphMatrix2PCoA}.
#' @param x_axis Which ordination axis to plot as the x-axis (defaults to 1).
#' @param y_axis Which ordination axis to plot as the y-axis (defaults to 2).
#' @param z_axis Which ordination axis to plot as the z-axis (defaults to NULL, i.e., is not plotted).
#' @param plot_taxon_names Optional to plot the names of the taxa (defaults to FALSE).
#' @param plot_internal_nodes Optional to plot the internal nodes of the tree (if included in \code{pcoa_input}) (defaults to FALSE).
#' @param plot_root Optional to plot the root separately (defaults to FALSE).
#' @param root_colour If plotting the root separately (previous option) sets the root colour.
#'
#' @details
#'
#' Uses output from \link{MorphMatrix2PCoA} as input.
#'
#' Allows plotting of a third axis using the technique of Matthew Wills (Wills et al. 1994; their Figures 4 and 8; Wills 1998; his Figure 4), where black and white indicate positive and negative values respectovely, and the size of points there magnitudes.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Emma Sherratt \email{emma.sherratt@@gmail.com}
#'
#' @references
#'
#' Wills, M. A., 1998. Cambrian and Recent disparity: the picture from priapulids. Paleobiology, 24, 177-199.
#'
#' Wills, M. A., Briggs, D. E. G. and Fortey, R. A., 1994. Disparity as an evolutionary index: a comparison of Cambrian and Recent arthropods. Paleobiology, 20, 93-130.
#'
#' @keywords principal coordinates
#'
#' @examples
#'
#' # Set random seed:
#' set.seed(4)
#'
#' # Generate a random tree for the Michaux 1989 data set:
#' tree <- rtree(length(rownames(Michaux1989$Matrix_1$Matrix)))
#'
#' # Add taxon names to the tree:
#' tree$tip.label <- rownames(Michaux1989$Matrix_1$Matrix)
#'
#' # Perform a phylogenetic Principal Coordinates Analysis:
#' pcoa_input <- MorphMatrix2PCoA(Michaux1989, Tree = tree)
#'
#' # Plot the results:
#' MorphospacePlot(pcoa_input, plot_taxon_names = TRUE)
#'
#' @export MorphospacePlot
MorphospacePlot <- function(pcoa_input, x_axis = 1, y_axis = 2, z_axis = NULL, plot_taxon_names = FALSE, plot_internal_nodes = FALSE, plot_root = TRUE, root_colour = "grey") {

# Option to plot names but not points
# Add legend to z-axis if using it
# Group colours
# Group colours will conflict with z-axis so have conditional to turn z off and warn user if doing so.




#gp is a vector with n species and p groups
#is.factor() # to check if group names, otherwise is assumed to be colours

#col.gp <- rainbow(length(levels(gp))) # generates a set of different colors length p
#names(col.gp) <- levels(gp) # assign those colurs to the p groups
#col.gp <- col.gp[match(gp, names(col.gp))] # creates a vector length n with a group and colour for each




  # Get vector of values that correspond to scree plot:
  scree_values <- apply(pcoa_input$vectors, 2, var) / sum(apply(pcoa_input$vectors, 2, var)) * 100

  # Make x-axis label:
  x_lab <- paste("PC", x_axis, " (", round(scree_values[x_axis], 2), "% of total variance)", sep = "")
  
  # Make y-axis label:
  y_lab <- paste("PC", y_axis, " (", round(scree_values[y_axis], 2), "% of total variance)", sep = "")

  # Create the basic plot space (will be empty for now):
  plot(pcoa_input$vectors[, x_axis], pcoa_input$vectors[, y_axis], type="n", bg = "black", xlab = x_lab, ylab = y_lab, asp = TRUE)

  # Case if no z-axis chosen:
  if(is.null(z_axis)) {

    # Make all z-axis colours black:
    z_colours <- rep("black", nrow(pcoa_input$vectors))
    
    # Make all z-axis values equal in size:
    z_sizes <- rep(1, nrow(pcoa_input$vectors))

  # Case if a z-axis is specified:
  } else {
    
    # Create vector of colours for z-axis (default of white):
    z_colours <- rep("white", nrow(pcoa_input$vectors))
    
    # CUpdate z-axis colours for positive values to black:
    z_colours[which(pcoa_input$vectors[, z_axis] > 0)] <- "black"
    
    # Create z-axis vector for absolute size of values (up to a max of 3):
    z_sizes <- abs(pcoa_input$vectors[, z_axis]) / max(abs(pcoa_input$vectors[, z_axis])) * 3
    
  }

  # Case if tree supplied:
  if(!is.null(pcoa_input$tree)) {

    # Sort axes by node number in tree:
    pcoa_input$vectors <- pcoa_input$vectors[c(pcoa_input$tree$tip.label, setdiff(rownames(pcoa_input$vectors), pcoa_input$tree$tip.label)), ]

    # For each branch, plot branch:
    for(i in 1:nrow(pcoa_input$tree$edge)) lines(x = pcoa_input$vectors[pcoa_input$tree$edge[i, ], x_axis], y = pcoa_input$vectors[pcoa_input$tree$edge[i, ], y_axis], col = "black")

    # Establish tip node numbers:
    tip_numbers <- c(1:ape::Ntip(pcoa_input$tree))
    
    # Establish internal node numbers:
    node_numbers <- setdiff(1:nrow(pcoa_input$vectors), tip_numbers)
    
    # Establish root number:
    root_number <- ape::Ntip(pcoa_input$tree) + 1

    # If plotting internal nodes, plot internal nodes:
    if(plot_internal_nodes) points(pcoa_input$vectors[node_numbers, x_axis], pcoa_input$vectors[node_numbers, y_axis], pch = 21, bg = z_colours[node_numbers], cex = z_sizes[node_numbers])

    # If plotting root separately, plot root:
    if(plot_root) points(pcoa_input$vectors[root_number, x_axis], pcoa_input$vectors[root_number, y_axis], pch = 21, col = root_colour, bg = root_colour, cex = z_sizes[root_number])

    # Plot tip data:
    points(pcoa_input$vectors[tip_numbers, x_axis], pcoa_input$vectors[tip_numbers, y_axis], pch = 21, bg = z_colours[tip_numbers], cex = z_sizes[tip_numbers])

    # If plotting taxon names:
    if(plot_taxon_names) {
    
      # First establish a default position for names (to the left of the point):
      x_positions <- rep(2, nrow(pcoa_input$vectors))
    
      # Now changes negative values to plot on the right instead:
      x_positions[which(pcoa_input$vectors[, x_axis] < 0)] <- 4
    
      # Plot taxon names (for tips only):
      text(x = pcoa_input$vectors[tip_numbers, x_axis], y = pcoa_input$vectors[tip_numbers, y_axis], labels = rownames(pcoa_input$vectors)[tip_numbers], pos = x_positions[tip_numbers], cex = 0.7)
    
    }

  # Case if no tree supplied:
  } else {

    # Plot points:
    points(pcoa_input$vectors[, x_axis], pcoa_input$vectors[, y_axis], pch = 21, bg = z_colours, cex = z_sizes)

    # If plotting taxon names:
    if(plot_taxon_names) {
    
      # First establish a default position for names (to the left of the point):
      x_positions <- rep(2, nrow(pcoa_input$vectors))
    
      # Now changes negative values to plot on the right instead:
      x_positions[which(pcoa_input$vectors[, x_axis] < 0)] <- 4
    
      # Plot taxon names:
      text(x = pcoa_input$vectors[, x_axis], y = pcoa_input$vectors[, y_axis], labels = rownames(pcoa_input$vectors), pos = x_positions, cex = 0.7)
    
    }

  }

}
