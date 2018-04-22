#' Safe Taxonomic Reduction
#' 
#' Performs Safe Taxonomic Reduction (STR) on a character-taxon matrix.
#' 
#' Performs Safe Taxonomic Reduction (Wilkinson 1995).
#' 
#' If no taxa can be safely removed will print the text "No taxa can be safely removed", and the \code{str.list} and \code{removed.matrix} will have no rows.
#'
#' NB: If your data contains inapplicable characters these will be treated as missing data, but this is inappropriate. Thus the user is advised to double check that any removed taxa make sense in the light of inapplicable states. (As far as I am aware this same behaviour occurs in the TAXEQ3 software.)
#'
#' @param morph.matrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#'
#' @return
#'
#' \item{str.list}{A matrix listing the taxa that can be removed (\code{Junior}), the taxa which they are equivalent to (\code{Senior}) and the rule under which they can be safely removed (\code{Rule}).}
#' \item{reduced.matrix}{A character-taxon matrix excluding the taxa that can be safely removed.}
#' \item{removed.matrix}{A character-taxon matrix of the taxa that can be safely removed.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references Wilkinson, M., 1995. Coping with abundant missing entries in phylogenetic inference using parsimony. Systematic Biology, 44, 501-514.
#'
#' @keywords Safe Taxonomic Reduction
#'
#' @examples
#' 
#' # Performs STR on the Gauthier 1986 dataset used in Wilkinson (1995):
#' str.out <- SafeTaxonomicReduction(Gauthier1986)
#' 
#' # View deleted taxa:
#' str.out$str.list
#' 
#' # View reduced matrix:
#' str.out$reduced.matrix
#' 
#' # View removed matrix:
#' str.out$removed.matrix
#' 
#' @export SafeTaxonomicReduction
SafeTaxonomicReduction <- function(morph.matrix) {
  
  # Store unaltered version of matrix to return to later:
  Unaltered <- morph.matrix
  
  # Prune out any zero weight characters, if they exist:
  if(any(morph.matrix$weights == 0)) morph.matrix <- MatrixPruner(clad.matrix = morph.matrix, characters2prune = which(morph.matrix$weights == 0))
  
  # Order matrix from least to most complete taxon (as least is most likely to be removed):
  morph.matrix$matrix <- morph.matrix$matrix[order(apply(apply(morph.matrix$matrix, 1, is.na), 2, sum), decreasing = TRUE), ]
  
  # Subfunction to be used by lapply to check downwards (look at more complete taxa only) for safely removable taxa:
  CheckDownwardForMatches <- function(rownumber, morph.matrix) {
    
    # First find which characters are not scored as missing in current taxon:
    NonMissingCharacters <- !is.na(morph.matrix$matrix[rownumber, ])
    
    # Build isolated matrix block from current taxon to end of matrix only for characters coded for current taxon:
    MatrixBlockToCheck <- morph.matrix$matrix[rownumber:nrow(morph.matrix$matrix), NonMissingCharacters, drop = FALSE]
    
    # Find any taxa that have missing characters inside the block (can not be true parents):
    TaxaWithMissingCharacters <- names(which(apply(apply(MatrixBlockToCheck, 1, is.na), 2, any)))
    
    # If any taxa have missing characters inside the block (can not be true parents) remove them from the block:
    if(length(TaxaWithMissingCharacters) > 0) MatrixBlockToCheck <- MatrixBlockToCheck[-match(TaxaWithMissingCharacters, rownames(MatrixBlockToCheck)), , drop = FALSE]
    
    # Set start column (first character in matrix block):
    StartColumn <- 1
    
    # Set end column (last character in matrix block):
    EndColumn <- ncol(MatrixBlockToCheck)
    
    # As long as there is more potential seniors (more than one row in the matrix block) and there are still characters to check:
    while(nrow(MatrixBlockToCheck) > 1 && StartColumn != (EndColumn + 1)) {
      
      # Only look at rows where taxa are coded for the ith character:
      NonMissingRows <- which(!is.na(MatrixBlockToCheck[2:nrow(MatrixBlockToCheck), StartColumn]))
      
      # List any insafe taxa (i.e., those with a diffrent coding to the current taxon):
      UnsafeTaxa <- names(which(MatrixBlockToCheck[(NonMissingRows + 1), StartColumn] != MatrixBlockToCheck[1, StartColumn]))
      
      # Reduce the matrix block to just potential safe taxon parents (and the current taxon):
      if(length(UnsafeTaxa) > 0) MatrixBlockToCheck <- MatrixBlockToCheck[-match(UnsafeTaxa, rownames(MatrixBlockToCheck)), , drop = FALSE]
      
      # Iterate to next colum:
      StartColumn <- StartColumn + 1
      
    }
    
    # If safe taxonomic reduction is possible store junior and senior(s):
    if(nrow(MatrixBlockToCheck) > 1) return(cbind(rep(x = rownames(MatrixBlockToCheck)[1], times = nrow(MatrixBlockToCheck) - 1), rownames(MatrixBlockToCheck)[2:nrow(MatrixBlockToCheck)]))
    
  }
  
  # Get any safely removable taxa and their senior parents:
  SafeToRemove <- lapply(as.list(1:(nrow(morph.matrix$matrix) - 1)), CheckDownwardForMatches, morph.matrix = morph.matrix)
  
  # If no taxa can be safely deleted:
  if(all(unlist(lapply(SafeToRemove, is.null)))) {
    
    # Warn user:
    print("No taxa can be safely removed")
    
    # Create empty safe to remove matrix:
    SafeToRemove <- matrix(nrow = 0, ncol = 3, dimnames = list(c(), c("Junior", "Senior", "Rule")))
    
    # Set reduced matrix as complete matrix:
    ReducedMatrix <- Unaltered$matrix
    
    # Set removed matrix as empty matrix:
    RemovedMatrix <- Unaltered$matrix[-(1:nrow(Unaltered$matrix)), , drop = FALSE]
    
  # If there are taxa that can be deleted:
  } else {
    
    # Collapse down to just those with data:
    SafeToRemove <- SafeToRemove[which(!unlist(lapply(SafeToRemove, is.null)))]
    
    # Rebuild as matrix:
    SafeToRemove <- cbind(unlist(lapply(SafeToRemove, '[', ,1)), unlist(lapply(SafeToRemove, '[', ,2)))
    
    # Add column names:
    colnames(SafeToRemove) <- c("Junior", "Senior")
    
    # Create empty vector to store rule data:
    STRRule <- c()
    
    # For each STR pair:
    for(i in 1:nrow(SafeToRemove)) {
      
      # Check for rule 1 (symmetrical coding):
      if(paste(morph.matrix$matrix[SafeToRemove[i, 1], ], collapse = "") == paste(morph.matrix$matrix[SafeToRemove[i, 2], ], collapse = "")) {
        
        # Check for any missing data:
        if(any(is.na(morph.matrix$matrix[SafeToRemove[i, 1], ]))) {
          
          # Is rule 1b:
          STRRule <- c(STRRule, "Rule 1B")
          
        # If no missing data:
        } else {
          
          # Is rule 1a:
          STRRule <- c(STRRule, "Rule 1A")
          
        }
        
      # Must be rule 2 (asymmetric coding):
      } else {
        
        # Check for any missing data:
        if(any(is.na(morph.matrix$matrix[SafeToRemove[i, 1], ]))) {
          
          # Is rule 2B:
          STRRule <- c(STRRule, "Rule 2B")
          
        # If no missing data is rule 1a:
        } else {
          
          # Is rule 2A:
          STRRule <- c(STRRule, "Rule 2A")
          
        }
        
      }
      
    }
    
    # Add rule to output:
    SafeToRemove <- cbind(SafeToRemove, STRRule)
    
    # Update column heading for rule:
    colnames(SafeToRemove)[3] <- "Rule"
    
    # Set reduced matrix as complete matrix:
    ReducedMatrix <- Unaltered$matrix[setdiff(rownames(Unaltered$matrix), unique(SafeToRemove[, "Junior"])), ]
    
    # Set removed matrix as empty matrix:
    RemovedMatrix <- Unaltered$matrix[sort(unique(SafeToRemove[, "Junior"])), ]
    
  }
  
  # Compile output into single list:
  output <- list(SafeToRemove, ReducedMatrix, RemovedMatrix)
  
  # Add names to output:
  names(output) <- c("str.list", "reduced.matrix", "removed.matrix")
  
  # Return output inviisibly:
  return(invisible(output))

}
