#' Collapses matrix to unique character state distributions
#'
#' Collapses matrix to unique character state distributions.
#'
#' Important: not recommended for general use.
#'
#' This function is intended to make a matrix with redundant character state distributions smaller by collapsing these to single characters and upweighting them accordingly. It is intended purely for us with MRP matrices, but may have some very restricted uses elsewhere.
#'
#' The function also deletes any characters weighted zero from the matrix.
#'
#' @param clad.matrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{SafeTaxonomicReduction} and \link{MatrixPruner}
#'
#' @keywords NEXUS
#'
#' @examples
#'
#' # Examine the matrix pre-compactification:
#' Michaux1989$matrix
#'
#' # Examine the wieights pre-compactification:
#' Michaux1989$weights
#'
#' # Compactify the matrix:
#' Michaux1989compact <- CompactifyMatrix(Michaux1989)
#'
#' # Examine the matrix post-compactification:
#' Michaux1989compact$matrix
#'
#' # Examine the weights post-compactification:
#' Michaux1989compact$weights
#'
#' @export CompactifyMatrix
CompactifyMatrix <- function(clad.matrix) {
  
  # FUTURE COULD CHECK FOR UNORD AND ORD WHEN BINARY AND HENCE MEANINGLESS
  
  # List any zero weight characters:
  ZeroWeightCharacters <- which(unlist(lapply(clad.matrix[2:length(clad.matrix)], '[[', "Weights")) == 0)
  
  # If there are zero weight characters then prune them:
  if(length(ZeroWeightCharacters) > 0) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = ZeroWeightCharacters)
  
  # For each matrix block:
  for(i in 2:length(clad.matrix)) {
    
    # Get strings for each character distribution, including ordering:
    char.distrib.strings <- paste(apply(clad.matrix[[i]]$Matrix, 2, paste, collapse = ""), clad.matrix[[i]]$Ordering, sep = " ")
    
    # Case if matrix can be compactified:
    if(length(unique(char.distrib.strings)) < length(char.distrib.strings)) {
      
      # Get rle of character distribution strings:
      rle.char.distrib.strings <- rle(sort(char.distrib.strings, decreasing = TRUE))
      
      # Set ordering of newly collapsed characters:
      clad.matrix[[i]]$Ordering <- unlist(lapply(strsplit(rle.char.distrib.strings$values, " "), '[', 2))
      
      # Set weights of newly collapsed characters by aggregating weights of source characters:
      clad.matrix[[i]]$Weights <- unlist(lapply(lapply(lapply(lapply(as.list(rle.char.distrib.strings$values), '==', char.distrib.strings), which), function(x) clad.matrix[[i]]$Weights[x]), sum))
      
      # Build new collapsed matrix:
      clad.matrix[[i]]$Matrix <- matrix(unlist(lapply(lapply(strsplit(rle.char.distrib.strings$values, " "), '[', 1), strsplit, split = "")), nrow = nrow(clad.matrix[[i]]$Matrix), dimnames = list(rownames(clad.matrix[[i]]$Matrix), c()))
      
      # Get ranges of values for characters in new collapsed matrix:
      MinMax <- lapply(lapply(lapply(lapply(lapply(lapply(apply(clad.matrix[[i]]$Matrix, 2, strsplit, split = "/"), unlist), strsplit, split = "&"), unlist), unique), as.numeric), range)
      
      # Set new minimum values for collapsed matrix:
      clad.matrix[[i]]$MinVals <- unlist(lapply(MinMax, '[', 1))
      
      # Set new maximum values for collapsed matrix:
      clad.matrix[[i]]$MaxVals <- unlist(lapply(MinMax, '[', 2))
      
    # Case if matrix cannot be compactified:
    } else {
      
      # Print message to user:
      print("Matrix cannot be compactified. All character distributions are unique and weights are greater than zero.")
      
    }
    
  }

  # Output unaltered matrix:
  return(clad.matrix)
  
}
