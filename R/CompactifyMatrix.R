#' Collapses matrix to unique character state distributions
#'
#' @description
#'
#' Collapses a cladistic matrix to just unique character state distributions and taxon names.
#'
#' @param CladisticMatrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#' @param Message Logical indicating whether or not a message should be printed to the screen if the matrix cannot be compactified.
#'
#' @details
#'
#' Important: not recommended for general use.
#'
#' This function is intended to make a matrix with redundant character state distributions smaller by collapsing these to single characters and upweighting them accordingly. It is intended purely for use with MRP matrices, but may have some very restricted uses elsewhere.
#'
#' The function also deletes any characters weighted zero from the matrix and will merge duplicate taxon anmes into unique character strings.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{SafeTaxonomicReduction} and \link{MatrixPruner}
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
CompactifyMatrix <- function(CladisticMatrix, Message = TRUE) {
  
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
    if(length(unique(char.distrib.strings)) < length(char.distrib.strings) || any(duplicated(rownames(CladisticMatrix[[i]]$Matrix)))) {
      
      # If collapsing characters because they are duplicated:
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
        
      }
      
      # If collapsing matrix because taxa are duplicated:
      if(any(duplicated(rownames(CladisticMatrix[[i]]$Matrix)))) {
        
        # Find duplicated taxa:
        DuplicatedTaxa <- unique(rownames(CladisticMatrix[[i]]$Matrix)[duplicated(rownames(CladisticMatrix[[i]]$Matrix))])
        
        # For each duplicated taxon:
        for(j in DuplicatedTaxa) {
          
          # Get rows for taxon:
          DuplicateRows <- which(rownames(CladisticMatrix[[i]]$Matrix) == j)
          
          # Only continue if the duplicated rows are actually variable:
          if(length(unique(apply(CladisticMatrix[[i]]$Matrix[DuplicateRows, ], 1, paste, collapse = ""))) > 1) {
            
            # Build duplicated matrix from other taxa:
            TemporaryMatrix <- matrix(rep(CladisticMatrix[[i]]$Matrix[-DuplicateRows, ], length(DuplicateRows)), ncol = ncol(CladisticMatrix[[i]]$Matrix) * length(DuplicateRows), dimnames = list(rownames(CladisticMatrix[[i]]$Matrix)[-DuplicateRows], c()))
            
            # Add duplicated taxon as single row:
            TemporaryMatrix <- rbind(TemporaryMatrix, as.vector(t(CladisticMatrix[[i]]$Matrix[DuplicateRows, ])))
            
            # Add duplicated taxon name:
            rownames(TemporaryMatrix)[nrow(TemporaryMatrix)] <- j
            
            # Update stored MRP matrix:
            CladisticMatrix[[i]]$Matrix <- TemporaryMatrix
            
            # Update other parts of matrix:
            CladisticMatrix[[i]]$Ordering <- rep(CladisticMatrix[[i]]$Ordering, length(DuplicateRows))
            CladisticMatrix[[i]]$Weights <- rep(CladisticMatrix[[i]]$Weights, length(DuplicateRows))
            CladisticMatrix[[i]]$MinVals <- rep(CladisticMatrix[[i]]$MinVals, length(DuplicateRows))
            CladisticMatrix[[i]]$MaxVals <- rep(CladisticMatrix[[i]]$MaxVals, length(DuplicateRows))
            
            # BELOW IS EEFCTIVELY A RECURSION OF THIS FUNCTION!
            
            # Get strings for each character distribution, including ordering:
            CharacterDistributionStrings <- paste(apply(CladisticMatrix[[i]]$Matrix, 2, paste, collapse = ""), CladisticMatrix[[i]]$Ordering, sep = " ")
            
            # If need to collapse characters because they are duplicated:
            if(length(unique(CharacterDistributionStrings)) < length(CharacterDistributionStrings)) {
              
              # Get rle of character distribution strings:
              RLECharacterDistributionStrings <- rle(sort(CharacterDistributionStrings, decreasing = TRUE))
              
              # Set ordering of newly collapsed characters:
              CladisticMatrix[[i]]$Ordering <- unlist(lapply(strsplit(RLECharacterDistributionStrings$values, " "), '[', 2))
              
              # Set weights of newly collapsed characters by aggregating weights of source characters:
              CladisticMatrix[[i]]$Weights <- unlist(lapply(lapply(lapply(lapply(as.list(RLECharacterDistributionStrings$values), '==', CharacterDistributionStrings), which), function(x) CladisticMatrix[[i]]$Weights[x]), sum))
              
              # Build new collapsed matrix:
              CladisticMatrix[[i]]$Matrix <- matrix(unlist(lapply(lapply(strsplit(RLECharacterDistributionStrings$values, " "), '[', 1), strsplit, split = "")), nrow = nrow(CladisticMatrix[[i]]$Matrix), dimnames = list(rownames(CladisticMatrix[[i]]$Matrix), c()))
              
              # Get ranges of values for characters in new collapsed matrix:
              MinMax <- lapply(lapply(lapply(lapply(lapply(lapply(apply(CladisticMatrix[[i]]$Matrix, 2, strsplit, split = "/"), unlist), strsplit, split = "&"), unlist), unique), as.numeric), range)
              
              # Set new minimum values for collapsed matrix:
              CladisticMatrix[[i]]$MinVals <- unlist(lapply(MinMax, '[', 1))
              
              # Set new maximum values for collapsed matrix:
              CladisticMatrix[[i]]$MaxVals <- unlist(lapply(MinMax, '[', 2))
              
            }
            
          # If duplicated rows are not variable:
          } else {
            
            # Remove all but one duplicated row from the matrix:
            CladisticMatrix[[i]]$Matrix <- CladisticMatrix[[i]]$Matrix[-DuplicateRows[2:length(DuplicateRows)], , drop = FALSE]
            
          }
          
        }
        
      }
      
    # Case if matrix cannot be compactified:
    } else {
      
      # Print message to user:
      if(Message) print("Matrix cannot be compactified. All character distributions are unique and weights are greater than zero.")
      
    }
    
  }

  # Output unaltered matrix:
  return(CladisticMatrix)
  
}
