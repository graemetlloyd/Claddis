#' Prunes a character matrix of characters or taxa
#'
#' @description
#'
#' Prunes a character matrix of characters, taxa, or both.
#'
#' @param CladisticMatrix The cladistic matrix in the format imported by \link{ReadMorphNexus}.
#' @param blocks2prune A vector of number(s) of any blocks to prune.
#' @param characters2prune A vector of character numbers to prune.
#' @param taxa2prune A vector of taxon names to prune (these must be present in \code{rownames(CladisticMatrix$matrix}).
#' @param removeinvariant A logical for whether invariant characters should (TRUE) or should not (FALSE, default) be pruned.
#'
#' @details
#'
#' Removing characters or taxa from a matrix imported using \link{ReadMorphNexus} is not simple due to associated vectors for ordering, character weights etc. To save repetitively pruning each part this function takes the matrix as input and vector(s) of either block numbers, character numbers, taxon names, or any combination thereof and returns a matrix with these items removed. Minimum and maximum values (used by \link{MorphDistMatrix}) are also updated and the user has the option to remove constant characters this way as well (e.g, to reduce the memory required for a DNA matrix).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{ReadMorphNexus}
#'
#' @keywords NEXUS
#'
#' @examples
#'
#' # Remove the outgroup taxon and characters 11 and 53 from Gauthier1986:
#' prunedmatrix <- MatrixPruner(CladisticMatrix = Gauthier1986, characters2prune = c(11, 53),
#'   taxa2prune = c("Outgroup"))
#'
#' # Show priuned matrix:
#' prunedmatrix$Matrix_1$Matrix
#'
#' @export MatrixPruner
MatrixPruner <- function(CladisticMatrix, blocks2prune = c(), characters2prune = c(), taxa2prune = c(), removeinvariant = FALSE) {
  
  # Subfunction to find length of character types for each character (i.e., unique values excluding polymorphisms but included inapplicables):
  LengthFinder <- function(x) {
    
    # Convert each column of matrix to a list of numeric values:
    x <- lapply(lapply(lapply(lapply(lapply(lapply(apply(x, 2, as.list), unlist), strsplit, split = "&|/"), unlist), unique), sort), length)
    
    # Return(x):
    return(x)
    
  }
  
  # Check that something to prune has been specified:
  if(is.null(blocks2prune) && is.null(characters2prune) && is.null(taxa2prune) && removeinvariant == FALSE) stop("No blocks, taxa, or characters to prune specified.")
  
  # Check blocks specified exist and stop and warn :
  if(length(setdiff(blocks2prune, 1:(length(CladisticMatrix) - 1))) > 0) stop("Block numbers specified that are not present in data.")
  
  # Check characters specified exist and stop and warn user if not:
  if(length(setdiff(characters2prune, 1:sum(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), ncol))))) > 0) stop("Characters specified that are outside the scope of the matrix. Check and retry.")
  
  # Check taxa specified exist and stop and warn user if not:
  if(length(setdiff(taxa2prune, rownames(CladisticMatrix$Matrix_1$Matrix))) > 0) stop("Taxa specified that are not found in the matrix. Check and retry.")
  
  # Get number of characters:
  NChars <- sum(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), ncol)))
  
  # If there are characters to prune:
  if(!is.null(characters2prune)) {
    
    # Get character blocks for each character in descendng order (as want to work backwards so things match up properly):
    CharacterBlocks <- unlist(lapply(lapply(lapply(as.list(sort(characters2prune, decreasing = TRUE)), '>', cumsum(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), ncol)))), which), length)) + 1
    
    # Initial build of characters in list form:
    CharactersAsList <- lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), function(x) 1:ncol(x))
    
    # Actually form list of character numbers (i.e., renumber characters in second or higher blocks):
    if(length(CharactersAsList) > 1) for(i in 2:length(CharactersAsList)) CharactersAsList[[i]] <- CharactersAsList[[i]] + max(CharactersAsList[[(i - 1)]])
    
    # For each unique character block:
    for(i in unique(CharacterBlocks)) {
      
      # Find columns to delete in ith matrix:
      ColumnsToDelete <- match(sort(characters2prune, decreasing = TRUE)[CharacterBlocks == i], CharactersAsList[[i]])
      
      # Remove characters from matrix:
      CladisticMatrix[[(i + 1)]]$Matrix <- CladisticMatrix[[(i + 1)]]$Matrix[, -ColumnsToDelete, drop = FALSE]
      
      # Remove characters from ordering:
      CladisticMatrix[[(i + 1)]]$Ordering <- CladisticMatrix[[(i + 1)]]$Ordering[-ColumnsToDelete]
      
      # Remove characters from weights:
      CladisticMatrix[[(i + 1)]]$Weights <- CladisticMatrix[[(i + 1)]]$Weights[-ColumnsToDelete]
      
      # Remove characters from minimum values:
      CladisticMatrix[[(i + 1)]]$MinVals <- CladisticMatrix[[(i + 1)]]$MinVals[-ColumnsToDelete]
      
      # Remove characters from maximum values:
      CladisticMatrix[[(i + 1)]]$MaxVals <- CladisticMatrix[[(i + 1)]]$MaxVals[-ColumnsToDelete]
      
    }
    
  }
  
  # If there are taxa to prune:
  if(!is.null(taxa2prune)) {
    
    # Remove pruned taxa from each block:
    for(i in 2:length(CladisticMatrix)) CladisticMatrix[[i]]$Matrix <- CladisticMatrix[[i]]$Matrix[-match(taxa2prune, rownames(CladisticMatrix[[i]]$Matrix)), , drop = FALSE]
    
  }
  
  # If there are blocks to prune:
  if(!is.null(blocks2prune)) {
    
    # Remove blocks to be rpuned:
    CladisticMatrix <- CladisticMatrix[-(blocks2prune + 1)]
    
    # Rename (renumber) remaining matrix blocks:
    names(CladisticMatrix[2:length(CladisticMatrix)]) <- paste("Matrix_", 1:(length(CladisticMatrix) - 1), sep = "")
    
  }
  
  # If there are invariant characters:
  if(removeinvariant) {
    
    # Find any invariant characters:
    InvariantsAsList <- lapply(lapply(lapply(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), LengthFinder), unlist), '<', 1), which)
    
    # If there are invariant characters:
    if(length(unlist(InvariantsAsList)) > 0) {
      
      # For each matrix block:
      for(i in 1:length(InvariantsAsList)) {
        
        # Only if there are invraints for this block:
        if(length(InvariantsAsList[[i]]) > 0) {
          
          # Remove characters from matrix:
          CladisticMatrix[[(i + 1)]]$Matrix <- CladisticMatrix[[(i + 1)]]$Matrix[, -InvariantsAsList[[i]], drop = FALSE]
          
          # Remove characters from ordering:
          CladisticMatrix[[(i + 1)]]$Ordering <- CladisticMatrix[[(i + 1)]]$Ordering[-InvariantsAsList[[i]]]
          
          # Remove characters from weights:
          CladisticMatrix[[(i + 1)]]$Weights <- CladisticMatrix[[(i + 1)]]$Weights[-InvariantsAsList[[i]]]
          
          # Remove characters from minimum values:
          CladisticMatrix[[(i + 1)]]$MinVals <- CladisticMatrix[[(i + 1)]]$MinVals[-InvariantsAsList[[i]]]
          
          # Remove characters from maximum values:
          CladisticMatrix[[(i + 1)]]$MaxVals <- CladisticMatrix[[(i + 1)]]$MaxVals[-InvariantsAsList[[i]]]
          
        }
        
      }
      
    }
    
  }
  
  # Check for empty blocks and store them as blocks to delete if found:
  NewBlocksToDelete <- which(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), ncol)) == 0)
  
  # If there are new blocks to prune:
  if(length(NewBlocksToDelete) > 0) {
    
    # Remove blocks to be rpuned:
    CladisticMatrix <- CladisticMatrix[-(NewBlocksToDelete + 1)]
    
    # Rename (renumber) remaining matrix blocks:
    names(CladisticMatrix[2:length(CladisticMatrix)]) <- paste("Matrix_", 1:(length(CladisticMatrix) - 1), sep = "")
    
  }
  
  # Return pruned matrix:
  return(CladisticMatrix)
  
}
