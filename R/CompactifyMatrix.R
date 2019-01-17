#' Collapses matrix to unique character state distributions
#'
#' @description
#'
#' Collapses a cladistic matrix to just unique character state distributions.
#'
#' @param CladisticMatrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#'
#' @details
#'
#' Important: not recommended for general use.
#'
#' This function is intended to make a matrix with redundant character state distributions smaller by collapsing these to single characters and upweighting them accordingly. It is intended purely for use with MRP matrices, but may have some very restricted uses elsewhere.
#'
#' The function also deletes any characters weighted zero from the matrix.
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
#' Michaux1989$Matrix_1$Matrix
#'
#' # Examine the weights pre-compactification:
#' Michaux1989$Matrix_1$Weights
#'
#' # Compactify the matrix:
#' Michaux1989compact <- CompactifyMatrix(Michaux1989)
#'
#' # Examine the matrix post-compactification:
#' Michaux1989compact$Matrix_1$Matrix
#'
#' # Examine the weights post-compactification:
#' Michaux1989compact$Matrix_1$Weights
#'
#' @export CompactifyMatrix
CompactifyMatrix <- function(CladisticMatrix) {
  
  # FUTURE COULD CHECK FOR UNORD AND ORD WHEN BINARY AND HENCE MEANINGLESS
  
  # List any zero weight characters:
  ZeroWeightCharacters <- which(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Weights")) == 0)
  
  # If there are zero weight characters then prune them:
  if(length(ZeroWeightCharacters) > 0) CladisticMatrix <- MatrixPruner(CladisticMatrix, characters2prune = ZeroWeightCharacters)
  
  # For each matrix block:
  for(i in 2:length(CladisticMatrix)) {
    
    # Get strings for each character distribution, including ordering:
    char.distrib.strings <- paste(apply(CladisticMatrix[[i]]$Matrix, 2, paste, collapse = ""), CladisticMatrix[[i]]$Ordering, sep = " ")
    
    # Case if matrix can be compactified:
    if(length(unique(char.distrib.strings)) < length(char.distrib.strings)) {
      
      # Get rle of character distribution strings:
      rle.char.distrib.strings <- rle(sort(char.distrib.strings, decreasing = TRUE))
      
      # Set ordering of newly collapsed characters:
      CladisticMatrix[[i]]$Ordering <- unlist(lapply(strsplit(rle.char.distrib.strings$values, " "), '[', 2))
      
      # Set weights of newly collapsed characters by aggregating weights of source characters:
      CladisticMatrix[[i]]$Weights <- unlist(lapply(lapply(lapply(lapply(as.list(rle.char.distrib.strings$values), '==', char.distrib.strings), which), function(x) CladisticMatrix[[i]]$Weights[x]), sum))
      
      # Build new collapsed matrix:
      CladisticMatrix[[i]]$Matrix <- matrix(unlist(lapply(lapply(strsplit(rle.char.distrib.strings$values, " "), '[', 1), strsplit, split = "")), nrow = nrow(CladisticMatrix[[i]]$Matrix), dimnames = list(rownames(CladisticMatrix[[i]]$Matrix), c()))
      
      # Get ranges of values for characters in new collapsed matrix:
      MinMax <- lapply(lapply(lapply(lapply(lapply(lapply(apply(CladisticMatrix[[i]]$Matrix, 2, strsplit, split = "/"), unlist), strsplit, split = "&"), unlist), unique), as.numeric), range)
      
      # Set new minimum values for collapsed matrix:
      CladisticMatrix[[i]]$MinVals <- unlist(lapply(MinMax, '[', 1))
      
      # Set new maximum values for collapsed matrix:
      CladisticMatrix[[i]]$MaxVals <- unlist(lapply(MinMax, '[', 2))
      
    # Case if matrix cannot be compactified:
    } else {
      
      # Print message to user:
      print("Matrix cannot be compactified. All character distributions are unique and weights are greater than zero.")
      
    }
    
  }

  # Output unaltered matrix:
  return(CladisticMatrix)
  
}
