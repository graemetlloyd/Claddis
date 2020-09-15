#' Trims marginal whitespace
#'
#' @description
#'
#' Trims any marginal whitespace from a vector of character string(s).
#'
#' @param x A character string
#'
#' @details
#'
#' Trims any marginal whitespace (spaces or tabs) from a vector of character string(s).
#'
#' @return
#'
#' A vector of character string(s) with any leading or trailing whitespace removed.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Example string:
#' x <- "   \td s f\t s  "
#'
#' # Trim only marginal whitespace:
#' trim_marginal_whitespace(x)
#'
#' @export trim_marginal_whitespace
trim_marginal_whitespace <- function(x) {
  
  # Make function to work on a scalar (single string):
  trim_scalar <- function(x) {
    
    # As long as the string has positive length:
    if (nchar(x = x) > 0) {
      
      # Split string into individual characters/whitespace:
      split_string <- strsplit(x = x, split = "")[[1]]
      
      # Find positions of any whitespace:
      position_is_whitespace <- apply(X = rbind(split_string == " ", split_string == "\t"), MARGIN = 2, FUN = any)
      
      # If everything is whitespace:
      if (all(position_is_whitespace)) {
        
        # Return empty string:
        return("")
        
        # If at least one character is not whitespace:
      } else {
        
        # Find first non-whitesapce position:
        first_non_whitespace_position <- min(which(position_is_whitespace == FALSE))
        
        # Find last non-whitesapce position:
        last_non_whitespace_position <- max(which(position_is_whitespace == FALSE))
        
        # Return trimmed string:
        return(value = paste(split_string[first_non_whitespace_position:last_non_whitespace_position], collapse = ""))
      }
      
      # If string has no length:
    } else {
      
      # Return string as is:
      return(value = x)
    }
  }
  
  # Return output with marginal whitesapce pruned:
  unlist(x = lapply(X = as.list(x = x), FUN = trim_scalar))
}

