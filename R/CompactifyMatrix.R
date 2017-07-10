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
#' # Examine the wieights post-compactification:
#' Michaux1989compact$weights
#'
#' @export CompactifyMatrix
CompactifyMatrix <- function(clad.matrix) {

# FUTURE COULD CHECK FOR UNORD AND ORD WHEN BINARY AND HENCE MEANINGLESS

    # Get strings for each character distribution, including ordering:
    char.distrib.strings <- paste(apply(clad.matrix$matrix, 2, paste, collapse = ""), clad.matrix$ordering, sep = " ")

    # Case if matrix can be compactified:
    if(length(unique(char.distrib.strings)) < length(char.distrib.strings) || sum(clad.matrix$weights == 0) > 0) {

        # If there are zero weight characters then prune them:
        if(sum(clad.matrix$weights == 0) > 0) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = which(clad.matrix$weights == 0))

        # For each character distribution string:
        for(i in 1:length(unique(char.distrib.strings))) {
    
            # Find matches for current unique character:
            matches <- which(char.distrib.strings == unique(char.distrib.strings)[i])
            
            # Sum weights of characters and store for first character:
            clad.matrix$weights[matches[1]] <- sum(clad.matrix$weights[matches])
            
            # Set additional matches to zero:
            if(length(matches) > 1) clad.matrix$weights[matches[-1]] <- 0
    
        }

        # Rerun deletion of zero weight characters:
        if(sum(clad.matrix$weights == 0) > 0) clad.matrix <- MatrixPruner(clad.matrix, characters2prune = which(clad.matrix$weights == 0))
        
        # Output compactified matrix:
        return(clad.matrix)

    # Case if matrix cannot be compactified:
    } else {
        
        # Print message to user:
        print("Matrix cannot be compactified. All character distributions are unique and weights are greater than zero.")
    
        # Output unaltered matrix:
        return(clad.matrix)
    
    }

}
