#' Permute all possible polymorphisms for a given set of states
#'
#' @description
#'
#' Given a set of discrete states, will permute all possible polymorphic combinations of those states.
#'
#' @param single_states A vector of single states (e.g., 0, 1, 2 etc.).
#'
#' @details
#'
#' This function solves a simple phylogenetic combinatorics problem - what are all the possible outcomes for a character to be in given polymorphisms (of any size) are allowed?
#'
#' For example, for three states (0, 1, 2) there are four possible polymorphisms: 0&1, 0&2, 1&2 and 0&1&2.
#'
#' If the user is instead only interested in the size of this state space, this is simply given by 2^N - N - 1, where N is the number of single states. Thus, the first several outcomes are:
#'
#' \preformatted{----------------------------------
#' | N states | N possible outcomes |
#' ----------------------------------
#' |     2    |          1          |
#' |     3    |          4          |
#' |     4    |          11         |
#' |     5    |          26         |
#' |     6    |          57         |
#' |     7    |          120        |
#' |     8    |          247        |
#' |     9    |          502        |
#' |    10    |          1,013      |
#' |    11    |          2,036      |
#' |    12    |          4,083      |
#' |    13    |          8,178      |
#' |    14    |          16,369     |
#' ----------------------------------}
#'
#' Note that this function is really designed for internal use, but may have value to some users and so is available "visibly" here.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' A vector of all possible polymorphic states.
#'
#' @seealso
#'
#' \link{make_costmatrix} and \link{permute_all_uncertainties}
#'
#' @examples
#'
#' # Get all possible states for the character 0, 1, and 2:
#' permute_all_polymorphisms(single_states = 0:2)
#'
#' @export permute_all_polymorphisms
permute_all_polymorphisms <- function(single_states) {
  n_states <- length(x = single_states)
  if (n_states < 2) stop("single_states must contain at last two values or no polymorphisms are possible.")
  unlist(
    x = lapply(
      X = as.list(x = 2:n_states),
      FUN = function(x) {
        apply(
          X = combn(x = sort(x = single_states), m = x),
          MARGIN = 2,
          FUN = function(y) paste(x = y, collapse = "&")
        )
      }
    )
  )
}
