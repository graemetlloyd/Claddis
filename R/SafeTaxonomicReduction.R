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
  
  # Store extra copy of matrix:
  full.matrix <- morph.matrix
  
  # Vector for storing constant characters:
  pars.unif <- vector(mode="numeric")
  
  # For each character:
  for(i in 1:length(morph.matrix$matrix[1, ])) {
    
    # Record constant, i.e. parsimony uninformative characters:
    if(length(unique(sort(morph.matrix$matrix[, i]))) <= 1) pars.unif[(length(pars.unif) + 1)] <- i

  }
  
  # Record zero weight characters:
  zero.wts <- which(morph.matrix$weights == 0)
  
  # Concatenate characters that are not usable in safe taxonomic reduction:
  deletes <- sort(unique(c(pars.unif, zero.wts)))
  
  # If there are characters that shold be ignored:
  if(length(deletes) > 0) {
    
    # Remove them from every section of the matrix:
    morph.matrix$matrix <- morph.matrix$matrix[, -deletes, drop = FALSE]
    morph.matrix$ordering <- morph.matrix$ordering[-deletes, drop = FALSE]
    morph.matrix$weights <- morph.matrix$weights[-deletes, drop = FALSE]
    morph.matrix$max.vals <- morph.matrix$max.vals[-deletes, drop = FALSE]
    morph.matrix$min.vals <- morph.matrix$min.vals[-deletes, drop = FALSE]

  }
  
  # Get distance matrix:
  dist.matrix <- MorphDistMatrixFast(morph.matrix)$raw.dist.matrix
  
  # Vector for storing zero value pairs:
  pairs <- c(NA, NA)
  
  # For each row:
  for(j in 1:length(dist.matrix[, 1])) {
    
    # For each column
    for(k in 1:length(dist.matrix[, 1])) {
      
      # Make sure we are only looking at one diagonal:
      if(j > k) {
        
        # For only distances that can be calculated:
        if(is.na(dist.matrix[j, k]) == FALSE) {
          
          # If the distance is zero:
          if(dist.matrix[j, k] == 0) {
            
            # Record the taxon pairing:
            pairs <- rbind(pairs, c(rownames(dist.matrix)[j], colnames(dist.matrix)[k]))

          }

        }

      }

    }

  }
  
  # Vectors to store results under each rule of Wilkinson
  rule.1a <- rule.1b <- rule.2a <- rule.2b <- matrix(nrow = 0, ncol = 2)
  
  # As long as there are zerovalue pairs:
  if(length(pairs) > 2) {
    
    # Delete first line which is empty:
    pairs <- pairs[-1, , drop = FALSE]
    
    # For each zero distance pair:
    for(j in 1:length(pairs[, 1])) {
      
      # Get scored characters for first taxon:
      missing.1 <- which(is.na(morph.matrix$matrix[match(pairs[j, 1], rownames(morph.matrix$matrix)), ]))
      
      # Get scored characters for second taxon:
      missing.2 <- which(is.na(morph.matrix$matrix[match(pairs[j, 2], rownames(morph.matrix$matrix)), ]))
      
      # STR Rule 1 test (equivalence; retain either taxon):
      if(length(setdiff(missing.1, missing.2)) == 0 && length(setdiff(missing.2, missing.1)) == 0) {
        
        # Meets Rule 1A (both taxa are known for all states):
        if(length(missing.1) == length(morph.matrix$matrix[1, ])) {
          
          # Add to Rule 1A list:
          rule.1a <- rbind(rule.1a, c(pairs[j, 1], pairs[j, 2]))
        
        # Meets rule 1B (both taxa incompletely known):
        } else {
          
          # Add to Rule 1B list:
          rule.1b <- rbind(rule.1b, c(pairs[j, 1],pairs[j, 2]))

        }
        
      # STR Rule 2 tests (asymmetric equivalence; retain more complete taxon):
      } else {
        
        # Taxon 2 redundant with respect to taxon 1:
        if(length(setdiff(missing.1, missing.2)) == 0) {
          
          # If taxon 1 is completely known (Rule 2A):
          if(length(missing.1) == length(morph.matrix$matrix[1, ])) {
            
            # Add to vector, junior first:
            rule.2a <- rbind(rule.2a, c(pairs[j, 2], pairs[j, 1]))
            
          # If taxon 1 is not completely known (Rule 2B):
          } else {
            
            # Add to vector, junior first:
            rule.2b <- rbind(rule.2b, c(pairs[j, 2], pairs[j, 1]))

          }

        }
        
        # Taxon 1 redundant with respect to taxon 2:
        if(length(setdiff(missing.2, missing.1)) == 0) {
          
          # If taxon 2 is completely known (Rule 2A):
          if(length(missing.2) == length(morph.matrix$matrix[1, ])) {
            
            # Add to vector, junior first:
            rule.2a <- rbind(rule.2a, c(pairs[j, 1], pairs[j, 2]))
            
          # If taxon 2 is not completely known (Rule 2B):
          } else {
            
            # Add to vector, junior first:
            rule.2b <- rbind(rule.2b, c(pairs[j, 1], pairs[j, 2]))

          }

		}

      }

    }

  }
  
  # List taxon pairs:
  pairs <- rbind(rule.1a, rule.1b, rule.2a, rule.2b)
  
  # List rules:
  rule <- c(rep("Rule 1A", length(rule.1a) / 2), rep("Rule 1B", length(rule.1b) / 2), rep("Rule 2A", length(rule.2a) / 2), rep("Rule 2B", length(rule.2b) / 2))
  
  # Combine into single table:
  str.list <- cbind(pairs, rule)
  
  # Name columns:
  colnames(str.list) <- c("Junior", "Senior", "Rule")

  # If there are taxa that can be safely removed:
  if(length(str.list) > 0) {
    
    # List junior taxa (i.e. those to remove):
    removes <- sort(unique(str.list[, "Junior"]))
    
    # New matrices for reduced and removed:
    reduced.matrix <- removed.matrix <- full.matrix
    
    # Remove STR taxa to create reduced matrix:
    reduced.matrix <- reduced.matrix$matrix[-match(removes, rownames(reduced.matrix$matrix)), , drop = FALSE]
    
    # Isolate removed taxa for removed matrix:
    removed.matrix <- removed.matrix$matrix[match(removes, rownames(removed.matrix$matrix)), , drop = FALSE]
    
  # If there are no taxa that can be safely removed:
  } else {
    
    # Create reduced matrix:
    reduced.matrix <- full.matrix$matrix
    
    # Create empty removed matrix:
    removed.matrix <- matrix(nrow = 0, ncol = ncol(full.matrix$matrix))
    
    # Print warning message:
    print("No taxa can be safely removed")

  }

  # Compile results into a list:
  result <- list(str.list, reduced.matrix, removed.matrix)

  # Add names to results list:
  names(result) <- c("str.list", "reduced.matrix", "removed.matrix")

  # Return results:
  return(invisible(result))

}
