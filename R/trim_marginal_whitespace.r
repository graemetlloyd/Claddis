#' Trims marginal whitespace
#'
#' @description
#'
#' Trims any marginal whitespace from a character string.
#'
#' @param x A character string
#'
#' @details
#'
#' Trims any marginal whitespace (spaces or tabs) from a character string.
#'
#' @return
#'
#' A character string with any leading or trailing whitespace removed.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' #
#' x <- "   \td s f\t s  "
#'
#' #
#' trim_marginal_whitespace(x)
#'
#' @export trim_marginal_whitespace
trim_marginal_whitespace <- function(x) { 

  # Split string into individual characters/whitespace:
  split_string <- strsplit(x = x, split = "")[[1]]

  # Find positions of any whitespace:
  position_is_whitespace <- apply(X = rbind(split_string == " ", split_string == "\t"), MARGIN = 2, FUN = any)
  
  # Find first non-whitesapce position:
  first_non_qwhitespace_position <- min(which(position_is_whitespace == FALSE))
  
  # Find last non-whitesapce position:
  last_non_whitespace_position <- max(which(position_is_whitespace == FALSE))
  
  # Return trimmed string:
  paste(split_string[first_non_qwhitespace_position:last_non_whitespace_position], collapse = "")

}
