#' Chronophylomorphospace Plot
#'
#' @description
#'
#' Plots a 3D chronophylomorphospace.
#'
#' @param pcoa_data Principal coordinate data in the format output by \link{MorphMatrix2PCoA} that includes a tree and ancestral states.
#' @param x_axis Which ordination axis to plot as the x-axis.
#' @param y_axis Which ordination axis to plot as the y-axis.
#' @param shadow Whether or not to plot a shadow (2D plot) on the bottom face of the 3D plot (defaults to TRUE).
#'
#' @details
#'
#' Creates a movable three-dimensional (two ordination axes plus time) plot of a phylomorphospace.
#'
#' Intended to mimic the data visualisation of Sakamoto and Ruta (2012; their Video S1).
#'
#' @author Emma Sherratt \email{emma.sherratt@@gmail.com} and Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Sakamoto, M. and Ruta, M. 2012. Convergence and divergence in the evolution of cat skulls: temporal and spatial patterns of morphological diversity. PLoS ONE, 7, e39752.
#'
#' @keywords principal coordinates
#'

#' @examples
#' \donttest{
#' 
#' # Set random seed:
#' set.seed(4)
#' 
#' # Generate a random tree for the Michaux 1989 data set:
#' Tree <- rtree(nrow(Michaux1989$Matrix_1$Matrix))
#' 
#' # Set root time so latest tip terminates at the present:
#' Tree$root.time <- max(diag(vcv(Tree)))
#' 
#' # Add taxon names to the tree:
#' Tree$tip.label <- rownames(Michaux1989$Matrix_1$Matrix)
#' 
#' # Perform a phylogenetic Principal Coordinates Analysis:
#' pcoa_data <- MorphMatrix2PCoA(Michaux1989, Tree = Tree)
#' 
#' # Plot a chronophylomorphospace:
#' ChronoPhyloMorphospacePlot(pcoa_data)
#' 
#' }

#' @export ChronoPhyloMorphospacePlot
ChronoPhyloMorphospacePlot <- function(pcoa_data, x_axis = 1, y_axis = 2, shadow = TRUE) {

   if (! requireNamespace("rgl", quietly = TRUE)) {
      stop(paste0(
		"To plot three-dimensional chrono-phylo-morphospaces, please install package rgl",
		"\n install.packages('lattice')"))
	  }
  
  # Add top level conditionals to check for a tree etc.

  # Default plotting parameters for a 2D morphospace. Need to change node colour for 3D using rgl
  p.p <- list()
  if(is.null(p.p$t.bg)) p.p$t.bg <- "black"
  if(is.null(p.p$t.pch)) p.p$t.pch <- 21
  if(is.null(p.p$t.cex)) p.p$t.cex <- 2
  if(is.null(p.p$n.bg)) p.p$n.bg <- "white"
  if(is.null(p.p$n.pch)) p.p$n.pch <- 21
  if(is.null(p.p$n.cex)) p.p$n.cex <- 1.25
  if(is.null(p.p$l.col)) p.p$l.col <- "black"
  if(is.null(p.p$lwd)) p.p$lwd <- 3
  if(is.null(p.p$txt.adj)) p.p$txt.adj <- c(-.1, -.1)
  if(is.null(p.p$txt.col)) p.p$txt.col <- "black"
  if(is.null(p.p$txt.cex)) p.p$txt.cex <- 1

  # Isolate Tree:
  Tree <- pcoa_data$Tree
  
  # Record number of tips:
  N <- ape::Ntip(Tree)
  
  # Isolate pcoa axes:
  pcoa_data <- pcoa_data$vectors

  # Little function to set limits for plotting (make it cube-like):
  limits <- function(x, s) {
      
    # Get range of x:
    r <- range(x)
    
    # Scale range values:
    rc <- scale(r, scale = FALSE)
    
    # ?????
    l <- mean(r) + s * rc
    
    # Return l:
    return(l)
  
  }

  # Get node ages for z-axis in plotting:
  z_axis <- GetNodeAges(Tree)

  # Make x label for plot:
  xlab <- paste("PC", x_axis, sep = "")
  
  # Make y label for plot:
  ylab <- paste("PC", y_axis, sep = "")

  # Make empty plot:
  rgl::plot3d(pcoa_data, type = "n", 
	xlim = limits(pcoa_data[, x_axis], 1.5), 
	ylim = limits(pcoa_data[, y_axis], 1.5), 
	zlim = limits(z_axis, 0), 
	asp = c(1, 1, 0.5), 
	xlab = xlab, ylab = ylab, 
	zlab = "Time (Ma)", 
	rgl::view3d(phi = 90, fov = 30))

  #plots tips
  rgl::points3d(pcoa_data[1:N, 1], 
	pcoa_data[1:N, 2], 
	z_axis[1:N], 
	col = p.p$t.bg, 
	size = p.p$t.cex * 4)

  #plots nodes
  rgl::points3d(pcoa_data[(N + 1):nrow(pcoa_data), 1], 
	pcoa_data[(N + 1):nrow(pcoa_data), 2], 
	z_axis[(N + 1):nrow(pcoa_data)], 
	col = p.p$n.bg, size = p.p$n.cex * 4)

  # plots branches
  for (i in 1:nrow(Tree$edge)) {
	rgl::lines3d(
		pcoa_data[(Tree$edge[i, ]), 1], 
		pcoa_data[(Tree$edge[i, ]), 2], 
		z_axis[(Tree$edge[i, ])], 
		lwd = 2)
	}

  # plot taxa labels
  rgl::text3d(
	pcoa_data[, x_axis], 
	pcoa_data[, y_axis], 
	z_axis, 
	rownames(pcoa_data), 
	col = p.p$txt.col, 
	cex = p.p$txt.cex, 
	adj = p.p$txt.adj)

  # If plotting the shadow of x and y axes at the base of the plot:
  if(shadow == TRUE){
      
    # Plot branches:
    for (i in 1:nrow(Tree$edge)) {
		rgl::lines3d(
			pcoa_data[(Tree$edge[i, ]), 1], 
			pcoa_data[(Tree$edge[i, ]), 2], 
			Tree$root.time, 
			lwd = 2, 
			alpha = 0.5)
		}

    # Plot internal nodes:
    rgl::points3d(
		pcoa_data[(N + 1):nrow(pcoa_data), 1], 
		pcoa_data[(N + 1):nrow(pcoa_data), 2], 
		Tree$root.time, 
		col = p.p$n.bg, 
		size = p.p$n.cex * 4, 
		alpha = 0.5)
      
    # Plot tips:
    rgl::points3d(pcoa_data[1:N, 1], 
		pcoa_data[1:N, 2], 
		Tree$root.time, 
		col = p.p$t.bg, 
		size = p.p$t.cex * 4, 
		alpha = 0.5)
      
  }

}
