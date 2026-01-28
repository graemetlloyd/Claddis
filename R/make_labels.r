#' Make unique text labels
#'
#' @description
#'
#' Given a requisite number, generates that many unique text labels.
#'
#' @param n The number of labels required,
#' @param case The case the labels should be in. Should be "upper" for uppercase (the default), or "lower" for lower case.
#'
#' @details
#'
#' Where a list of unique text labels are required (i.e., where simple numbering will not suffice) it can be useful to have a simple function that generates the required amount.
#'
#' In practice, this is simple in R when N is 26 or less as the \code{letters} or \code{LETTERS} object can be used for this purpose. For example, to get ten unique uppercase labels we can simply use:
#'
#' \code{LETTERS[1:10]}
#'
#' This function works in a similar way but will add a second, third etc. letter where the value of \code{n} requires it.
#'
#' @return A character vector of \code{n} unique labels in upper- or lowercase.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Make 40 unique lowercase text labels:
#' make_labels(n = 40, case = "lower")
#'
#' @export make_labels
make_labels <- function(n, case = "upper") {
  
  # Check case is valid option:
  if (length(x = setdiff(x = case, y = c("upper", "lower"))) > 0) {
    stop("case must be one of \"upper\" or \"lower\".")
  }
  
  # Subfunction to recursively generated a single label for the kth value:
  make_label <- function(k) {
    if (case == "lower") letters_vector <- letters
    if (case == "upper") letters_vector <- LETTERS
    if (k <= 26) {
      letters_vector[k]
    } else {
      make_label((k - 1) %/% 26) |> paste0(letters_vector[(k - 1) %% 26 + 1])
    }
  }
  
  # Make full list of n variables:
  sapply(X = 1:n, FUN = make_label)
}
