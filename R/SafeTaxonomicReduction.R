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
  zero.wts <- grep(TRUE, morph.matrix$weights == 0)
  
  # Concatenate characters that are not usable in safe taxonomic reduction:
  deletes <- sort(unique(c(pars.unif, zero.wts)))
  
  # If there are characters that shold be ignored:
  if(length(deletes) > 0) {
    
    # Remove them from every section of the matrix:
    morph.matrix$matrix <- morph.matrix$matrix[, -deletes]
    morph.matrix$ordering <- morph.matrix$ordering[-deletes]
    morph.matrix$weights <- morph.matrix$weights[-deletes]
    morph.matrix$max.vals <- morph.matrix$max.vals[-deletes]
    morph.matrix$min.vals <- morph.matrix$min.vals[-deletes]

  }
  
  # Get distance matrix:
  dist.matrix <- MorphDistMatrix(morph.matrix)$raw.dist.matrix
  
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
  rule.1a <- rule.1b <- rule.2a <- rule.2b <- c(NA, NA)
  
  # As long as there are zerovalue pairs:
  if(length(pairs) > 2) {
    
    # Delete first line which is empty:
    pairs <- pairs[-1, ]
    
    # Case if only one pair:
    if(length(pairs) == 2) pairs <- t(as.matrix(pairs))
    
    # For each zero distance pair:
    for(j in 1:length(pairs[, 1])) {
      
      # Get scored characters for first taxon:
      missing.1 <- grep(TRUE, is.na(morph.matrix$matrix[match(pairs[j, 1], rownames(morph.matrix$matrix)), ]))
      
      # Get scored characters for second taxon:
      missing.2 <- grep(TRUE, is.na(morph.matrix$matrix[match(pairs[j, 2], rownames(morph.matrix$matrix)), ]))
      
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
  
  # Remove empty first rows
  if(length(rule.1a) == 2) rule.1a <- rule.1a[-(1:2)] else rule.1a <- rule.1a[-1,]
  if(length(rule.1b) == 2) rule.1b <- rule.1b[-(1:2)] else rule.1b <- rule.1b[-1,]
  if(length(rule.2a) == 2) rule.2a <- rule.2a[-(1:2)] else rule.2a <- rule.2a[-1,]
  if(length(rule.2b) == 2) rule.2b <- rule.2b[-(1:2)] else rule.2b <- rule.2b[-1,]
  
  # List taxon pairs:
  pairs <- rbind(rule.1a, rule.1b, rule.2a, rule.2b)
  
  # List rules:
  rule <- c(rep("Rule 1A", length(rule.1a) / 2), rep("Rule 1B", length(rule.1b) / 2), rep("Rule 2A", length(rule.2a) / 2), rep("Rule 2B", length(rule.2b) / 2))
  
  # Combine into single table:
  str.list <- cbind(pairs, rule)
  
  # If there are taxa that can be safely removed:
  if(length(str.list) > 0) {
    
    # Name columns:
    colnames(str.list) <- c("Junior", "Senior", "Rule")
    
    # List junior taxa (i.e. those to remove):
    removes <- sort(unique(str.list[, "Junior"]))
    
    # New matrices for reduced and removed:
    reduced.matrix <- removed.matrix <- full.matrix
    
    # Remove STR taxa to create reduced matrix:
    reduced.matrix <- reduced.matrix$matrix[-match(removes, rownames(reduced.matrix$matrix)), ]
    
    # Isolate removed taxa for removed matrix:
    removed.matrix <- removed.matrix$matrix[match(removes, rownames(removed.matrix$matrix)), ]
    
    # Compile results into a list:
    result <- list(str.list, reduced.matrix, removed.matrix)
    
    # Add names to results list:
    names(result) <- c("str.list", "reduced.matrix", "removed.matrix")
    
  # If there are no taxa that can be safely removed:
  } else {
    
    # Have warning message as result:
    result <- "No taxa can be safely removed"

  }
  
  # Return results:
  return(invisible(result))

}
