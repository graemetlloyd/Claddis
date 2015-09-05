#' Principal Coordinates on a Cladistic Matrix
#' 
#' Performs Principal Coordinates Analysis (PCoA) on a cladistic matrix.
#' 
#' Gower (1966). Caillez (1983). Uses \link{pcoa} from the \link{ape} package.
#' 
#' @param morph.matrix A vector of mode character representing the tip names for which an ancestor is sought.
#' @param distance.method The distance method to use (one of "RED", "GED", "GC", or "MORD" - the default). See \link{MorphDistMatrix} for more details.
#' @param transform.proportional.distances The transformation to apply to propotional (0 to 1) distances (one of "none", "sqrt", or "arcsine_sqrt" - the default). See \link{MorphDistMatrix} for more details.
#' @param correction The negative eigenvalue correction to use (one of "lingoes", "none", or "cailliez" - the default). See \link{pcoa} for more details.
#' @param tree If a phylmorphospace is desired then a tree with root age and branch-lengths must be included.
#' @param estimate.allchars If including a tree whether you want to estinate ancestral states for all characters (default is FALSE). See \link {AncStateEstMatrix} for more details.
#' @param estimate.tips If including a tree whether you want to estinate missing or polymorphic tip states (default is FALSE). See \link {AncStateEstMatrix} for more details.

# Other options to go to pcoa.

#'
#' @return \item{output}{The output.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @keywords principal coordinates
#'
#' @examples
#'
#' # Nothing yet
#'
#' @export MorphMatrix2PCoA
MorphMatrix2PCoA <- function(morph.matrix, distance.method = "MORD", transform.proportional.distances = "arcsine_sqrt", correction = "cailliez", tree = NULL, estimate.allchars = FALSE, estimate.tips = FALSE) {
  
  # If no tree is supplied:
  if(is.null(tree)) {
    
    
    morph_distances <- MorphDistMatrix(morph.matrix, transform.proportional.distances = transform.proportional.distances)
    
    if(distance.method == "RED") trimmed_distances <- TrimMorphDistMatrix(morph_distances$raw.dist.matrix)

    if(distance.method == "GED") trimmed_distances <- TrimMorphDistMatrix(morph_distances$GED.dist.matrix)

    if(distance.method == "GC") trimmed_distances <- TrimMorphDistMatrix(morph_distances$gower.dist.matrix)

    if(distance.method == "MORD") trimmed_distances <- TrimMorphDistMatrix(morph_distances$max.dist.matrix)
    
    if(!is.null(trimmed_distances$removed.taxa)) message(paste("The following taxa had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed.taxa, collapse = ", ")))
    
    pcoa_results <- pcoa(trimmed_distances$dist.matrix, correction = correction, rn = rownames(trimmed_distances$dist.matrix))
    
  # Case if a tree is included (and a phylomorphospace is requested):
  } else {
      
    ancestral_values <- AncStateEstMatrix(morph.matrix, tree, estimate.allchars = FALSE, estimate.tips = FALSE)

    morph.matrix$matrix <- rbind(morph.matrix$matrix, ancestral_values)

    morph_distances <- MorphDistMatrix(morph.matrix, transform.proportional.distances = transform.proportional.distances)

    if(distance.method == "RED") trimmed_distances <- TrimMorphDistMatrix(morph_distances$raw.dist.matrix, tree = tree)
    
    if(distance.method == "GED") trimmed_distances <- TrimMorphDistMatrix(morph_distances$GED.dist.matrix, tree = tree)
    
    if(distance.method == "GC") trimmed_distances <- TrimMorphDistMatrix(morph_distances$gower.dist.matrix, tree = tree)
    
    if(distance.method == "MORD") trimmed_distances <- TrimMorphDistMatrix(morph_distances$max.dist.matrix, tree = tree)

    if(!is.null(trimmed_distances$removed.taxa)) message(paste("The following taxa or nodes had to be removed to produce a complete distance matrix:", paste(trimmed_distances$removed.taxa, collapse = ", ")))

    tree <- trimmed_distances$tree

    pcoa_results <- pcoa(trimmed_distances$dist.matrix, correction = correction, rn = rownames(trimmed_distances$dist.matrix))

  }
  
  output <- c(list(tree), list(trimmed_distances$removed.taxa), pcoa_results)

  names(output)[[1]] <- "tree"
  names(output)[[2]] <- "removed.taxa"

  invisible(output)
  
}

#tree <- rtree(length(rownames(Michaux1989$matrix)))
#tree$tip.label <- rownames(Michaux1989$matrix)

#x <- MorphMatrix2PCoA(Michaux1989)
#y <- MorphMatrix2PCoA(Michaux1989, tree = tree)

#x
#y

#tree <- rtree(length(rownames(Gauthier1986$matrix)))
#tree$tip.label <- rownames(Gauthier1986$matrix)

#x <- MorphMatrix2PCoA(Gauthier1986)
#y <- MorphMatrix2PCoA(Gauthier1986, tree = tree)

#x
#y
