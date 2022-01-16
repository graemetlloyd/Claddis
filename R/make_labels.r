#' Make unique text labels
#'
#' @description
#'
#' Given a requisite number, generates that many unique text labels.
#'
#' @param N The number of labels required,
#'
#' @details
#'
#' Where a list of unique text labels are required (i.e., where simple numbering will not suffice) it can be useful to have a simple function that generates the required amount.
#'
#' In practice, this is simple in R when N is 26 or less as the \code{LETTERS} object can be used for this purpose. For example, to get ten unique labels:
#'
#' \code{LETTERS[1:10]}
#'
#' This function works in a similar way but will add a second, third etc. letter where the value of N requires it.
#'
#' @return A character vector of N unique labels.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make 40 unique text labels:
#' make_labels(N = 40)
#'
#' @export make_labels
make_labels <- function(N) {

  # Make starting exponent:
  exponent <- 1
  
  # If necessary increase exponent to meet the demands of N:
  while (N > (26 ^ exponent)) exponent <- exponent + 1
  
  # Generate a list of letters using the exponent value:
  letters_list <- lapply(
    X = as.list(x = 1:exponent),
    FUN = function(i) LETTERS
  )
  
  # Return a character vector of unique labels:
  apply(
    X = expand.grid(letters_list),
    MARGIN = 1,
    FUN = paste,
    collapse = ""
  )[1:N]
}
