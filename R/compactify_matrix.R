#' Collapses matrix to unique character state distributions
#'
#' @description
#'
#' Collapses a cladistic matrix to just unique character state distributions and taxon names.
#'
#' @param cladistic_matrix The cladistic matrix in the format imported by \link{read_nexus_matrix}.
#' @param Message Logical indicating whether or not a message should be printed to the screen if the matrix cannot be compactified.
#'
#' @details
#'
#' Important: not recommended for general use.
#'
#' This function is intended to make a matrix with redundant character state distributions smaller by collapsing these to single characters and upweighting them accordingly. It is intended purely for use with MRP matrices, but may have some very restricted uses elsewhere.
#'
#' The function also deletes any characters weighted zero from the matrix and will merge duplicate taxon names into unique character strings.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Examine the matrix pre-compactification:
#' michaux_1989$matrix_1$matrix
#'
#' # Examine the weights pre-compactification:
#' michaux_1989$matrix_1$weights
#'
#' # Compactify the matrix:
#' michaux_1989compact <- compactify_matrix(michaux_1989)
#'
#' # Examine the matrix post-compactification:
#' michaux_1989compact$matrix_1$matrix
#'
#' # Examine the weights post-compactification:
#' michaux_1989compact$matrix_1$weights
#'
#' @export compactify_matrix
compactify_matrix <- function(cladistic_matrix, Message = TRUE) {
  
  # FUTURE COULD CHECK FOR UNORD AND ORD WHEN BINARY AND HENCE MEANINGLESS
  
  # List any zero weight characters:
  ZeroWeightcharacters <- which(unlist(lapply(cladistic_matrix[2:length(cladistic_matrix)], '[[', "weights")) == 0)
  
  # If there are zero weight characters then prune them:
  if (length(ZeroWeightcharacters) > 0) cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix, characters2prune = ZeroWeightcharacters)
  
  # For each matrix block:
  for(i in 2:length(cladistic_matrix)) {
    
    # Get strings for each character distribution, including ordering:
    char.distrib.strings <- paste(apply(cladistic_matrix[[i]]$matrix, 2, paste, collapse = ""), cladistic_matrix[[i]]$ordering, sep = " ")
    
    # Case if matrix can be compactified:
    if (length(unique(char.distrib.strings)) < length(char.distrib.strings) || any(duplicated(rownames(cladistic_matrix[[i]]$matrix)))) {
      
      # If collapsing characters because they are duplicated:
      if (length(unique(char.distrib.strings)) < length(char.distrib.strings)) {
        
        # Get rle of character distribution strings:
        rle.char.distrib.strings <- rle(sort(char.distrib.strings, decreasing = TRUE))
        
        # Set ordering of newly collapsed characters:
        cladistic_matrix[[i]]$ordering <- unlist(lapply(strsplit(rle.char.distrib.strings$values, " "), '[', 2))
        
        # Set weights of newly collapsed characters by aggregating weights of source characters:
        cladistic_matrix[[i]]$weights <- unlist(lapply(lapply(lapply(lapply(as.list(rle.char.distrib.strings$values), '==', char.distrib.strings), which), function(x) cladistic_matrix[[i]]$weights[x]), sum))
        
        # Build new collapsed matrix:
        cladistic_matrix[[i]]$matrix <- matrix(unlist(lapply(lapply(strsplit(rle.char.distrib.strings$values, " "), '[', 1), strsplit, split = "")), nrow = nrow(cladistic_matrix[[i]]$matrix), dimnames = list(rownames(cladistic_matrix[[i]]$Matrix), c()))
        
        # Get ranges of values for characters in new collapsed matrix:
        MinMax <- lapply(lapply(lapply(lapply(lapply(lapply(apply(cladistic_matrix[[i]]$matrix, 2, strsplit, split = "/"), unlist), strsplit, split = "&"), unlist), unique), as.numeric), range)
        
        # Set new minimum values for collapsed matrix:
        cladistic_matrix[[i]]$minimum_values <- unlist(lapply(MinMax, '[', 1))
        
        # Set new maximum values for collapsed matrix:
        cladistic_matrix[[i]]$maximum_values <- unlist(lapply(MinMax, '[', 2))
        
      }
      
      # If collapsing matrix because taxa are duplicated:
      if (any(duplicated(rownames(cladistic_matrix[[i]]$matrix)))) {
        
        # Find duplicated taxa:
        DuplicatedTaxa <- unique(rownames(cladistic_matrix[[i]]$matrix)[duplicated(rownames(cladistic_matrix[[i]]$matrix))])
        
        # For each duplicated taxon:
        for(j in DuplicatedTaxa) {
          
          # Get rows for taxon:
          DuplicateRows <- which(rownames(cladistic_matrix[[i]]$matrix) == j)
          
          # Only continue if the duplicated rows are actually variable:
          if (length(unique(apply(cladistic_matrix[[i]]$matrix[DuplicateRows, ], 1, paste, collapse = ""))) > 1) {
            
            # Build duplicated matrix from other taxa:
            TemporaryMatrix <- matrix(rep(cladistic_matrix[[i]]$matrix[-DuplicateRows, ], length(DuplicateRows)), ncol = ncol(cladistic_matrix[[i]]$matrix) * length(DuplicateRows), dimnames = list(rownames(cladistic_matrix[[i]]$matrix)[-DuplicateRows], c()))
            
            # Add duplicated taxon as single row:
            TemporaryMatrix <- rbind(TemporaryMatrix, as.vector(t(cladistic_matrix[[i]]$matrix[DuplicateRows, ])))
            
            # Add duplicated taxon name:
            rownames(TemporaryMatrix)[nrow(TemporaryMatrix)] <- j
            
            # Update stored MRP matrix:
            cladistic_matrix[[i]]$matrix <- TemporaryMatrix
            
            # Update other parts of matrix:
            cladistic_matrix[[i]]$ordering <- rep(cladistic_matrix[[i]]$ordering, length(DuplicateRows))
            cladistic_matrix[[i]]$weights <- rep(cladistic_matrix[[i]]$weights, length(DuplicateRows))
            cladistic_matrix[[i]]$minimum_values <- rep(cladistic_matrix[[i]]$minimum_values, length(DuplicateRows))
            cladistic_matrix[[i]]$maximum_values <- rep(cladistic_matrix[[i]]$maximum_values, length(DuplicateRows))
            
            # BELOW IS EEFCTIVELY A RECURSION OF THIS FUNCTION!
            
            # Get strings for each character distribution, including ordering:
            CharacterDistributionStrings <- paste(apply(cladistic_matrix[[i]]$matrix, 2, paste, collapse = ""), cladistic_matrix[[i]]$ordering, sep = " ")
            
            # If need to collapse characters because they are duplicated:
            if (length(unique(CharacterDistributionStrings)) < length(CharacterDistributionStrings)) {
              
              # Get rle of character distribution strings:
              RLECharacterDistributionStrings <- rle(sort(CharacterDistributionStrings, decreasing = TRUE))
              
              # Set ordering of newly collapsed characters:
              cladistic_matrix[[i]]$ordering <- unlist(lapply(strsplit(RLECharacterDistributionStrings$values, " "), '[', 2))
              
              # Set weights of newly collapsed characters by aggregating weights of source characters:
              cladistic_matrix[[i]]$weights <- unlist(lapply(lapply(lapply(lapply(as.list(RLECharacterDistributionStrings$values), '==', CharacterDistributionStrings), which), function(x) cladistic_matrix[[i]]$weights[x]), sum))
              
              # Build new collapsed matrix:
              cladistic_matrix[[i]]$matrix <- matrix(unlist(lapply(lapply(strsplit(RLECharacterDistributionStrings$values, " "), '[', 1), strsplit, split = "")), nrow = nrow(cladistic_matrix[[i]]$matrix), dimnames = list(rownames(cladistic_matrix[[i]]$matrix), c()))
              
              # Get ranges of values for characters in new collapsed matrix:
              MinMax <- lapply(lapply(lapply(lapply(lapply(lapply(apply(cladistic_matrix[[i]]$matrix, 2, strsplit, split = "/"), unlist), strsplit, split = "&"), unlist), unique), as.numeric), range)
              
              # Set new minimum values for collapsed matrix:
              cladistic_matrix[[i]]$minimum_values <- unlist(lapply(MinMax, '[', 1))
              
              # Set new maximum values for collapsed matrix:
              cladistic_matrix[[i]]$maximum_values <- unlist(lapply(MinMax, '[', 2))
              
            }
            
          # If duplicated rows are not variable:
          } else {
            
            # Remove all but one duplicated row from the matrix:
            cladistic_matrix[[i]]$matrix <- cladistic_matrix[[i]]$matrix[-DuplicateRows[2:length(DuplicateRows)], , drop = FALSE]
            
          }
          
        }
        
      }
      
    # Case if matrix cannot be compactified:
    } else {
      
      # Print message to user:
      if (Message) print("Matrix cannot be compactified. All character distributions are unique and weights are greater than zero.")
      
    }
    
  }

  # Output unaltered matrix:
  return(cladistic_matrix)
  
}
