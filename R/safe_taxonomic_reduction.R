#' Safe Taxonomic Reduction
#'
#' @description
#'
#' Performs Safe Taxonomic Reduction (STR) on a character-taxon matrix.
#'
#' @param cladistic.matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#'
#' @details
#'
#' Performs Safe Taxonomic Reduction (Wilkinson 1995).
#'
#' If no taxa can be safely removed will print the text "No taxa can be safely removed", and the \code{str.list} and \code{removed.matrix} will have no rows.
#'
#' NB: If your data contains inapplicable characters these will be treated as missing data, but this is inappropriate. Thus the user is advised to double check that any removed taxa make sense in the light of inapplicable states. (As far as I am aware this same behaviour occurs in the TAXEQ3 software.)
#'
#' @return
#'
#' \item{str.list}{A matrix listing the taxa that can be removed (\code{Junior}), the taxa which they are equivalent to (\code{Senior}) and the rule under which they can be safely removed (\code{Rule}).}
#' \item{reduced.matrix}{A character-taxon matrix excluding the taxa that can be safely removed.}
#' \item{removed.matrix}{A character-taxon matrix of the taxa that can be safely removed.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references Wilkinson, M., 1995. Coping with abundant missing entries in phylogenetic inference using parsimony. \emph{Systematic Biology}, \bold{44}, 501-514.
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{compactify_matrix}, \link{prune_cladistic_matrix}, \link{safe_taxonomic_reinsertion}, \link{read_nexus_matrix}, \link{write_nexus_matrix}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Performs STR on the Gauthier 1986 dataset used in Wilkinson (1995):
#' str.out <- safe_taxonomic_reduction(Gauthier1986)
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
#' @export safe_taxonomic_reduction
safe_taxonomic_reduction <- function(cladistic.matrix) {
  
  # Store unaltered version of matrix to return to later:
  Unaltered <- cladistic.matrix
  
  # If there is more than one matrix block (will need to temporarily combine into single matrix):
  if (length(cladistic.matrix) > 2) {
    
    # Get list of just the matrix blocks:
    MatrixBlocks <- lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "Matrix")
    
    # For each additional block add at end of first block:
    for(i in 2:length(MatrixBlocks)) MatrixBlocks[[1]] <- cbind(MatrixBlocks[[1]], MatrixBlocks[[i]][rownames(MatrixBlocks[[1]]), ])
    
    # Store single matrix block as first block of cladistic.matrix:
    cladistic.matrix[[2]]$Matrix <- MatrixBlocks[[1]]
    
    # Store ordering for all blocks in first block of cladistic.matrix:
    cladistic.matrix[[2]]$Ordering <- as.vector(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "Ordering")))
    
    # Store weights for all blocks in first block of cladistic.matrix:
    cladistic.matrix[[2]]$weights <- as.vector(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "weights")))
    
    # Store minimum values for all blocks in first block of cladistic.matrix:
    cladistic.matrix[[2]]$MinVals <- as.vector(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "MinVals")))
    
    # Store maximum values for all blocks in first block of cladistic.matrix:
    cladistic.matrix[[2]]$MaxVals <- as.vector(unlist(lapply(cladistic.matrix[2:length(cladistic.matrix)], '[[', "MaxVals")))
    
    # Isolate to just first matrix block:
    cladistic.matrix <- cladistic.matrix[1:2]

  }
  
  # Prune out any zero weight characters, if they exist:
  if (any(cladistic.matrix[[2]]$weights == 0)) cladistic.matrix <- prune_cladistic_matrix(cladistic.matrix = cladistic.matrix, characters2prune = which(cladistic.matrix[[2]]$weights == 0))
  
  # Order matrix from least to most complete taxon (as least is most likely to be removed):
  cladistic.matrix[[2]]$Matrix <- cladistic.matrix[[2]]$Matrix[order(apply(apply(cladistic.matrix[[2]]$Matrix, 1, is.na), 2, sum), decreasing = TRUE), ]
  
  # Subfunction to be used by lapply to check downwards (look at more complete taxa only) for safely removable taxa:
  CheckDownwardForMatches <- function(rownumber, cladistic.matrix) {
    
    # First find which characters are not scored as missing in current taxon:
    NonMissingCharacters <- !is.na(cladistic.matrix[[2]]$Matrix[rownumber, ])
    
    # Build isolated matrix block from current taxon to end of matrix only for characters coded for current taxon:
    MatrixBlockToCheck <- cladistic.matrix[[2]]$Matrix[rownumber:nrow(cladistic.matrix[[2]]$Matrix), NonMissingCharacters, drop = FALSE]
    
    # Find any taxa that have missing characters inside the block (can not be true parents):
    TaxaWithMissingCharacters <- names(which(apply(apply(MatrixBlockToCheck, 1, is.na), 2, any)))
    
    # If any taxa have missing characters inside the block (can not be true parents) remove them from the block:
    if (length(TaxaWithMissingCharacters) > 0) MatrixBlockToCheck <- MatrixBlockToCheck[-match(TaxaWithMissingCharacters, rownames(MatrixBlockToCheck)), , drop = FALSE]
    
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
      if (length(UnsafeTaxa) > 0) MatrixBlockToCheck <- MatrixBlockToCheck[-match(UnsafeTaxa, rownames(MatrixBlockToCheck)), , drop = FALSE]
      
      # Iterate to next colum:
      StartColumn <- StartColumn + 1
      
    }
    
    # If safe taxonomic reduction is possible store junior and senior(s):
    if (nrow(MatrixBlockToCheck) > 1) return(cbind(rep(x = rownames(MatrixBlockToCheck)[1], times = nrow(MatrixBlockToCheck) - 1), rownames(MatrixBlockToCheck)[2:nrow(MatrixBlockToCheck)]))
    
  }
  
  # Get any safely removable taxa and their senior parents:
  SafeToRemove <- lapply(as.list(1:(nrow(cladistic.matrix[[2]]$Matrix) - 1)), CheckDownwardForMatches, cladistic.matrix = cladistic.matrix)
  
  # If no taxa can be safely deleted:
  if (all(unlist(lapply(SafeToRemove, is.null)))) {
    
    # Warn user:
    print("No taxa can be safely removed")
    
    # Create empty safe to remove matrix:
    SafeToRemove <- matrix(nrow = 0, ncol = 3, dimnames = list(c(), c("Junior", "Senior", "Rule")))
    
    # Set reduced matrix as unaltered matrix:
    RemovedMatrix <- ReducedMatrix <- Unaltered
    
    # Set removed matrix as empty matrix:
    RemovedMatrix <- prune_cladistic_matrix(Unaltered, taxa2prune = rownames(Unaltered$Matrix_1$Matrix))
    
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
      if (paste(cladistic.matrix$Matrix_1$Matrix[SafeToRemove[i, 1], ], collapse = "") == paste(cladistic.matrix$Matrix_1$Matrix[SafeToRemove[i, 2], ], collapse = "")) {
        
        # Check for any missing data:
        if (any(is.na(cladistic.matrix$Matrix_1$Matrix[SafeToRemove[i, 1], ]))) {
          
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
        if (any(is.na(cladistic.matrix$Matrix_1$Matrix[SafeToRemove[i, 1], ]))) {
          
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
    ReducedMatrix <- prune_cladistic_matrix(Unaltered, taxa2prune = unique(SafeToRemove[, "Junior"]))
    
    # Set removed matrix as empty matrix:
    RemovedMatrix <- prune_cladistic_matrix(Unaltered, taxa2prune = setdiff(rownames(Unaltered$Matrix_1$Matrix), unique(SafeToRemove[, "Junior"])))
    
  }
  
  # Compile output into single list:
  output <- list(str.list = SafeToRemove, reduced.matrix = ReducedMatrix, removed.matrix = RemovedMatrix)
  
  # Return output inviisibly:
  return(invisible(output))
  
}
