#' Collapses matrix to unique character state distributions
#' 
#' Collapses matrix to unique character state distributions.
#' 
#' Removing characters or taxa from a matrix imported using \link{ReadMorphNexus} is not simple due to asscoiated vectors for ordering, character weights etc. To save repetitively pruning each part this function takes the matrix as input and vectors of either taxon names, character numbers, or one of each and removes those from the matrix. Minimum and maximum values (used by \link{MorphDistMatrix}) are also updated.
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
#' # Remove the outgroup taxon and characters 11 and 53 from Gauthier1986:
#' prunedmatrix <- MatrixPruner(clad.matrix = Gauthier1986, taxa2prune = c("Outgroup"),
#'   characters2prune = c(11, 53))
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
        print("Matrix cannot be compactified. All character dsitributions are unique and weights are greater than zero.")
    
        # Output unaltered matrix:
        return(clad.matrix)
    
    }

}
