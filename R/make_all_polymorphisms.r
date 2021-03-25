#' Make all possible polymorphisms for a given set of states
#'
#' @description
#'
#' Given a set of discrete states will generate all possibly polymorphic comibinations of those states.
#'
#' @param single_states A vector of single states (e.g., 0, 1, 2 etc.).
#'
#' @details
#'
#' Thsi funtion solves a phylogenetic combinatorics problem - what are all the possible outcomes for a character to be in given polymorphisms are allowed?
#'
#' For example, for three states (0, 1, 2) there are seven possible states: 0, 1, 2, 0&1, 0&2, 1&2 and 0&1&2.
#'
#' If the user is instead only interested in the size of this state space this is given by 2^N - 1, where N is the number of single states. Thus, the first several sizes are:
#'
#' \preformatted{----------------------------------
#' | N states | N possible outcomes |
#' ----------------------------------
#' |     2    |          3          |
#' |     3    |          7          |
#' |     4    |          15         |
#' |     5    |          31         |
#' |     6    |          63         |
#' |     7    |          127        |
#' |     8    |          255        |
#' |     9    |          511        |
#' |    10    |          1023       |
#' |    11    |          2047       |
#' |    12    |          4095       |
#' |    13    |          8191       |
#' |    14    |          16383      |
#' ----------------------------------}
#'
#' Note that this function is really designed for internal use, but may have value to some users.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' A vector of all possible polymorphic (and non-polymorphic) states.
#'
#' @seealso
#'
#' \link{make_stepmatrix}
#'
#' @examples
#'
#' # Get all possible states for the character 0, 1, and 2:
#' make_all_polymorphisms(single_states = 0:2)
#'
#' @export make_all_polymorphisms
make_all_polymorphisms <- function(single_states) unlist(x = lapply(X = as.list(x = 1:length(x = single_states)), FUN = function(x) apply(X = combn(x = sort(x = single_states), m = x), MARGIN = 2, FUN = function(y) paste(x = y, collapse = "&"))))
