#' Principal Coordinates on a Cladistic Matrix
#' 
#' Performs Principal Coordinates Analysis (PCoA) on a cladistic matrix.
#' 
#' Takes a cladistic matrix in the format imported by \link{ReadMorphNexus} and performs Principal Coordinates (Gower 1966) analysis on it.
#'
#' This function is effectively a wrapper for \link{pcoa} from the \link{ape} package and the user is referred there for some of the options (e.g., using the Caillez 1983 approach to avoiding negative eigenvalues).
#'
#' If providing a tree and inferring ancestral states then options to also infer missing or uncertain tips and whether to infer values for all characters at all internal nodes are provided (via \link{AncStateEstMatrix}).
#'
#' Other options within the function concern the distance metric to use and the transformation to be used if selecting a propotional distance (see \link{MorphDistMatrix}).
#'
#' IMPORTANT: The function can remove taxa (or if including a tree, nodes as well) if they lead to an incomplete distance matrix (see \link{TrimMorphDistMatrix}).
#'
#' @param morph.matrix A vector of mode character representing the tip names for which an ancestor is sought.
#' @param distance.method The distance method to use (one of "RED", "GED", "GC", or "MORD" - the default). See \link{MorphDistMatrix} for more details.
#' @param transform.proportional.distances The transformation to apply to propotional (0 to 1) distances (one of "none", "sqrt", or "arcsine_sqrt" - the default). See \link{MorphDistMatrix} for more details.
#' @param polymorphism.behaviour The distance behaviour for dealing with polymorphisms. Must be one of \code{"mean.difference"} or \code{"min.difference"} (the default.
#' @param uncertainty.behaviour The distance behaviour for dealing with uncertainties. Must be one of \code{"mean.difference"} or \code{"min.difference"} (the default.
#' @param correction The negative eigenvalue correction to use (one of "lingoes", "none", or "cailliez" - the default). See \link{pcoa} for more details.
#' @param tree If a phylmorphospace is desired then a tree with root age and branch-lengths must be included.
#' @param estimate.allchars If including a tree whether you want to estinate ancestral states for all characters (default is FALSE). See \link{AncStateEstMatrix} for more details.
#' @param estimate.tips If including a tree whether you want to estinate missing or polymorphic tip states (default is FALSE). See \link{AncStateEstMatrix} for more details.
#'
#' @return \item{tree}{The tree (if supplied). Note this may be pruned from the input tree by \link{TrimMorphDistMatrix}.}
#' @return \item{dist.matrix}{The distance matrix. Note this may be pruned by \link{TrimMorphDistMatrix} and thus not include all taxa.}
#' @return \item{removed.taxa}{A vector of taxa (or nodes) removed by \link{TrimMorphDistMatrix}. Returns NULL if none are removed.}
#' @return \item{note}{See \link{pcoa}.}
#' @return \item{values}{See \link{pcoa}.}
#' @return \item{vectors}{See \link{pcoa}.}
#' @return \item{trace}{See \link{pcoa}.}
#' @return \item{vectors.cor}{See \link{pcoa}.}
#' @return \item{trace.cor}{See \link{pcoa}.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Cailliez, F., 1983. The analytical solution of the additive constant problem. Psychometrika, 48, 305-308.
#'
#' Gower, J. C., 1966. Some distance properties of latent root and vector methods used in multivariate analysis. Biometrika, 53, 325-338.
#'
#' @keywords principal coordinates
#'
#' @examples
#'
#' # Run on Michaux (189) data set with defaults:
#' x <- MorphMatrix2PCoA(Michaux1989)
#'
#' # Show output:
#' x
#'
#' # Generate a (mde up) tree:
#' tree <- rtree(length(rownames(Michaux1989$matrix)))
#'
#' # Add taxon names to it:
#' tree$tip.label <- rownames(Michaux1989$matrix)
#'
#' # Run with tree:
#' y <- MorphMatrix2PCoA(Michaux1989, tree = tree)
#'
#' # Show new output:
#' y
#'
#' @export MorphMatrix2PCoA
MorphMatrix2PCoA <- function(morph.matrix, distance.method = "MORD", transform.proportional.distances = "arcsine_sqrt", polymorphism.behaviour = "min.difference", uncertainty.behaviour = "min.difference", correction = "cailliez", tree = NULL, estimate.allchars = FALSE, estimate.tips = FALSE) {
    
# Add some top level conditionsl here to check input is valid.
  
  # If no tree is supplied:
  if(is.null(tree)) {
    
    # Get morphological distances from the cladistic matrix:
    morph_distances <- MorphDistMatrix(morph.matrix, distance = distance.method, transform.proportional.distances = transform.proportional.distances, polymorphism.behaviour = polymorphism.behaviour, uncertainty.behaviour = uncertainty.behaviour)
    
    # Get trimmed distances:
    trimmed_distances <- TrimMorphDistMatrix(morph_distances$DistanceMatrix)

    # If trimming of matrix lead to taxa being removed warn user:
    if(!is.null(trimmed_distances$removed.taxa)) message(paste("The following taxa had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed.taxa, collapse = ", ")))
    
    # Perform Principal Coordinates Analysis on the data:
    pcoa_results <- pcoa(trimmed_distances$dist.matrix, correction = correction, rn = rownames(trimmed_distances$dist.matrix))
    
  # Case if a tree is included (and a phylomorphospace is requested):
  } else {
      
    # Get ancestral character states:
    ancestral_values <- AncStateEstMatrix(morph.matrix, tree, estimate.allchars = FALSE, estimate.tips = FALSE)

    # Case if estimating tip values:
    if(estimate.tips) {
    
      # Overwrite entire matrix with ancestral values as they include tips too.
      morph.matrix$matrix <- ancestral_values
    
    # Case if not estimating tip values:
    } else {
    
      # Add ancestral values into matrix with existing tip values:
      morph.matrix$matrix <- rbind(morph.matrix$matrix, ancestral_values)
    
    }

    # Get morphological distances from the cladistic matrix:
    morph_distances <- MorphDistMatrix(morph.matrix, distance = distance.method, transform.proportional.distances = transform.proportional.distances, polymorphism.behaviour = polymorphism.behaviour, uncertainty.behaviour = uncertainty.behaviour)
    
    # Get trimmed distances:
    trimmed_distances <- TrimMorphDistMatrix(morph_distances$DistanceMatrix, tree = tree)

    # If trimming of matrix lead to taxa or nodes being removed warn user:
    if(!is.null(trimmed_distances$removed.taxa)) message(paste("The following taxa or nodes had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed.taxa, collapse = ", ")))

    # Store (possibly trimmed) tree ready to be output:
    tree <- trimmed_distances$tree

    # Perform Principal Coordinates Analysis on the data:
    pcoa_results <- pcoa(trimmed_distances$dist.matrix, correction = correction, rn = rownames(trimmed_distances$dist.matrix))

  }
  
  # Compile output:
  output <- c(list(tree), list(trimmed_distances$dist.matrix), list(trimmed_distances$removed.taxa), pcoa_results)

  # Add variable name for tree:
  names(output)[[1]] <- "tree"

  # Add variable name for distaance matrix:
  names(output)[[2]] <- "dist.matrix"
  
  # Add variable name for removed taxa and nodes:
  names(output)[[3]] <- "removed.taxa"

  # Return output invisibly:
  invisible(output)
  
}
