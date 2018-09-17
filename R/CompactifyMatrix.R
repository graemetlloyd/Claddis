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
      
      # For each character distribution string:
      for(j in 1:length(unique(char.distrib.strings))) {
        
        # Find matches for current unique character:
        matches <- which(char.distrib.strings == unique(char.distrib.strings)[j])
        
        # Sum weights of characters and store for first character:
        clad.matrix[[i]]$Weights[matches[1]] <- sum(clad.matrix[[i]]$Weights[matches])
        
        # Set additional matches to zero:
        if(length(matches) > 1) clad.matrix[[i]]$Weights[matches[-1]] <- 0
        
      }
      
    # Case if matrix cannot be compactified:
    } else {
      
      # Print message to user:
      print("Matrix cannot be compactified. All character distributions are unique and weights are greater than zero.")
      
    }
    
  }
  
  # List any zero weight characters:
  ZeroWeightCharacters <- which(unlist(lapply(clad.matrix[2:length(clad.matrix)], '[[', "Weights")) == 0)
  
  # If there are zero weight characters then prune them:
  if(length(ZeroWeightCharacters) > 0) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = ZeroWeightCharacters)

  # Output unaltered matrix:
  return(clad.matrix)
  
}
