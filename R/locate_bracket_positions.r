#' Locates matching positions for sets of brackets in a text string
#'
#' @description
#'
#' Given a text string will return the positions of each matching pair of opening and closing brackets.
#'
#' @param input_string An input string containing matching brackets, such as a Newick string or character state tree.
#' @param bracket_type The type of bracket to use. Must be one of parentheses \code{()}, curly braces \code{{}}, or square brackets \code{[]}.
#'
#' @details
#'
#' This function is designed to deal with Newick strings and character state trees - ways of encoding information using nested parentheses. Although it is intended for internal use it seems sufficiently general to share as part of the package.
#'
#' The function works by traversing the string from left to right and noting the position of each opening parenthesis and then storing the corresponding position for its' matching closing parenthesis.
#'
#' It currently only works for a single string, but coud be built into a for loop or apply function if multiple strings are desired.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @return
#'
#' A two-column matrix indicating opening and closing bracket positions within \code{input_string}.
#'
#' @seealso
#'
#' \link{convert_state_tree_to_adjacency_matrix}
#'
#' @examples
#'
#' # Locate the positions of a set of parentheses in a character state tree:
#' locate_bracket_positions(
#'   input_string = "(((5)4)3,(2)1)0",
#'   bracket_type = "()"
#' )
#'
#' @export locate_bracket_positions
locate_bracket_positions <- function(input_string, bracket_type = "()") {
  
  # Check bracket_type is valid and stop and warn user if not:
  if(length(x = setdiff(x = bracket_type, y = c("()", "{}", "[]"))) > 0) stop("bracket_type must be one of \"()\", \"{}\", or \"[]\".")
  
  # Set bracket symbol based on bracket_type for use later in function:
  if(bracket_type == "()") {open <- "("; close = ")"}
  if(bracket_type == "{}") {open <- "{"; close = "}"}
  if(bracket_type == "[]") {open <- "["; close = "]"}
  
  # Check there is only one input string and stop and warn user if not:
  if(length(x = input_string) != 1) stop("input_string must be a single value (i.e., not a vector).")
  
  # Count numbers of opening and closing brackets (for checking):
  n_opens <- length(x = grep(pattern = open, x = input_string, fixed = TRUE))
  n_closes <- length(x = grep(pattern = close, x = input_string, fixed = TRUE))
  
  # Check nput string contains
  if(n_opens == 0 || n_closes == 0) stop("input_string must contain at least one pair of brackets.")
  
  # Create vector form of string for use later:
  split_string <- strsplit(x = input_string, split = "")[[1]]
  
  # Check brackets are in proportion and stop and warn user if not:
  if(n_opens != n_closes) stop("input_string does not contain equal numbers of opening and closing brackets.")
  
  # Check opening brackets always precede closing brackets and stop and warn user if not:
  if(!all(sort(x = which(x = split_string == open)) < sort(x = which(x = split_string == close)))) stop("opening brackets do not always precede closing brackets.")
  
  # Find all poistions containing a bracket (only bits we care about):
  all_bracket_positions <- which(x = ((split_string == open) + (split_string == close)) == 1)
  
  # Create empty matrix to store bracket positions (eventual output):
  bracket_positions <- matrix(nrow = 0, ncol = 2, dimnames = list(c(), c("open", "close")))
  
  # For each bracket position (can ignore rest of string:
  for(needle in all_bracket_positions) {
    
    # If an opneing bracket then start a new row in bracket positions:
    if(split_string[needle] == open) bracket_positions <- rbind(bracket_positions, c(needle, NA))
    
    # If a closing bracket then close current set of parentheses:
    if(split_string[needle] == close) bracket_positions[max(x = which(x = is.na(x = bracket_positions[, "close"]))), "close"] <- needle
    
  }
  
  # Return matrix of bracket positions:
  bracket_positions
  
}
