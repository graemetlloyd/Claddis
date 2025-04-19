#' Permute all combinations of x of size m with replacement
#'
#' @description
#'
#' Given a vector x, permutes all possible groups of size m ignoring order and allowing any item in x to appear multiple times.
#'
#' @param x A character vector.
#' @param m A positive integer indicating the size of the set desired.
#'
#' @details
#'
#' This is a simple combinatoric function used internally in Claddis where all possible combinations of \code{x} that are size \code{m} are permuted. Note that this ignores order (i.e., the sets \\{A,B\\} and \\{B,A\\} are considered identical) and replacements (or multiples) of an element of \code{x} are allowed (i.e., the sets \\{A,A\\} and \\{B,B\\} are both valid).
#'
#' @return A matrix of m columns where each row is a unique combination of x.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Permute all the ways the letters A-C can form a set of size 3:
#' permute_combinations_with_replacement(x = LETTERS[1:3], m = 3)
#'
#' @export permute_combinations_with_replacement
permute_combinations_with_replacement <- function(x, m) {
  
  # Get length of x (n):
  n <- length(x = x)
  
  # Subfunction to get indices for combinations:
  combination_indices <- function(n, m) {
    if (m == 1) return(value = matrix(data = 1:n, ncol = 1))
    if (m > 1) {
      x <- do.call(
        what = rbind,
        args = sapply(X = 1:n, FUN = function(i) cbind(i, i:n), simplify = FALSE)
      )
      if (m > 2) {
        while(ncol(x = x) < m) {
          x <- do.call(
            what = rbind,
            args = apply(
              X = x,
              MARGIN = 1,
              FUN = function(i) {
                size_i <- length(x = i)
                j <- as.list(x = i)
                j[[(length(x = j) + 1)]] <- i[size_i]:n
                do.call(what = cbind, args = j)
              }
            )
          )
        }
      }
      x <- unname(obj = x)
      return(value = x)
    }
  }
  
  # If there are multiple elements in x:
  if (n > 1) {
    
    # Now form all combinations of the elements of x:
    combinations <- apply(
      X = combination_indices(n = n, m = m),
      MARGIN = 1,
      FUN = function(i) x[i],
      simplify = FALSE
    )
    
    # Make into matrix:
    combinations <- do.call(what = rbind, args = combinations)
    
  # If there is only one elemnt in x:
  } else {
    
    # Make a simple matrix of x repeated m times:
    combinations <- matrix(data = rep(x = x, times = m), ncol = m)
  }
  
  # Return combinations to user:
  combinations
}
