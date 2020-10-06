#' Time bins class
#'
#' @description
#'
#' Functions to deal with the time bins class.
#'
#' @param x A timeBins object.
#'
#' @details
#'
#' Claddis uses various classes to define specific types of data, here the use of time bins (to bin any temporal data) ae assigned the class "timeBins" and should look something like this:
#'
#' \preformatted{                  fad  lad
#'     Cenomanian    99.6 93.5
#'     Turonian      93.5 89.3
#'     Coniacian     89.3 85.8
#'     Santonian     85.8 83.5
#'     Campanian     83.5 70.6
#'     Maastrichtian 70.6 65.5}
#'
#' I.e., a matrix with two columns (fad = first appearance date and lad = last appearance date) with rows corresponding to named time bins and individual values ages in millions of years ago (Ma). The object should also have class \code{timeBins} (see example below for how to generate a valid object). Note also that the convention in Claddis is to have time bins be ordered from oldest to youngest.
#'
#' \code{is.timeBins} checks whether an object is or is not a valid timeBins object.
#'
#' @return \code{is.timeBins} returns either TRUE or FALSE.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a time bins object:
#' time_bins <- matrix(data = c(99.6, 93.5, 93.5, 89.3, 89.3, 85.8, 85.8, 83.5, 83.5, 70.6, 70.6, 65.5), ncol = 2, byrow = TRUE, dimnames = list(c("Cenomanian", "Turonian", "Coniacian", "Santonian", "Campanian", "Maastrichtian"), c("fad", "lad")))
#'
#' # Check that this is a valid timeBins object (will fail as class is not set):
#' is.timeBins(time_bins = time_bins)
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Check that this is a valid timeBins object (will succeed as format and
#' # class are correct):
#' is.timeBins(x = time_bins)
#'
#' @export is.timeBins
check_time_bins <- function(time_bins) {
  
  # Create empty vector to store error messages:
  messages <- vector(mode = "character")
  
  # Check time_bins has class timeBins and add error message to output if true:
  if (!inherits(x = time_bins, what = "timeBins")) messages <- c(messages, "time_bins must be an object of class \"timeBins\".")
  
  # Check time bins are in form of matrix add error message to output if false:
  if (!is.matrix(x = time_bins)) messages <- c(messages, "time_bins must be in the form of a matrix.")
  
  # Check time_bins has two columns and add error message to output if false:
  if (ncol(x = time_bins) != 2) messages <- c(messages, "time_bins must have exactly two columns.")
  
  # Check time_bins has mutliple rows (bins) and add error message to output if false:
  if (nrow(x = time_bins) < 2) messages <- c(messages, "time_bins must have at least two rows.")
  
  # Check time_bins column names are null and add error message to output if true:
  if (is.null(x = colnames(x = time_bins))) messages <- c(messages, "time_bins must have columns named \"fad\" and \"lad\" in that order.")
  
  # Check time_bins row names are null and add error message to output if true:
  if (is.null(x = rownames(x = time_bins))) messages <- c(messages, "time_bins must have unique row names corresponding to each time bin.")
  
  # Check time_bins column names are correct add error message to output if false:
  if (!all(colnames(x = time_bins) == c("fad", "lad"))) messages <- c(messages, "time_bins must have columns named \"fad\" and \"lad\" in that order.")
  
  # Check time_bins row names have at least one character each and add error message to output if false:
  if (!all(x = nchar(x = rownames(x = time_bins)) > 0)) messages <- c(messages, "time_bins must have row names formed from at least one character.")
  
  # Check time_bins row names are all unique and add error message to output if false:
  if (any(x = duplicated(x = rownames(x = time_bins)))) messages <- c(messages, "time_bins must have unique row names.")
  
  # Check time_bins contains only numeric values and add error message to output if false:
  if (!is.numeric(x = time_bins)) messages <- c(messages, "time_bins must contain numeric values representing millions of years ago (Ma).")
  
  # Check time_bins start before they end and have positive length and add error message to output if false:
  if (!all(x = time_bins[, "fad"] > time_bins[, "lad"])) messages <- c(messages, "time_bins must consist of \"fad\" values that exceed corresponding \"lad\" values. I.e., each time bin must have an older \"fad\" than \"lad\" and be of positive length.")
  
  # Check time_bins do not overlap each other and add error message to output if true:
  if (any(x = unlist(x = lapply(X = apply(X = time_bins, MARGIN = 1, FUN = list), FUN = function(x) length(x = which(x = apply(X = cbind(x[[1]]["fad"] > time_bins[, "lad"], x[[1]]["lad"] < time_bins[, "fad"]), MARGIN = 1, FUN = all))))) > 1)) messages <- c(messages, "time_bins must not overlap each other.")
  
  # Check time_bins are in correct order and add error message to output if false:
  if (!all(x = order(x = time_bins[, "fad"]) == seq(from = nrow(x = time_bins), to = 1, by = -1))) messages <- c(messages, "time_bins must be ordered from oldest (top) to youngest (bottom).")
  
  # Check time_bins have no gaps between them and add error message to output if true:
  if (!all(x = time_bins[1:(nrow(x = time_bins) - 1), "lad"] == time_bins[2:nrow(x = time_bins), "fad"])) messages <- c(messages, "time_bins must abut each other. I.e., there can be no gaps between successive time bins.")
  
  # Return messages:
  messages
}

is.timeBins <- function(x) {
  
  # Get any error messages for time_bins:
  messages <- check_time_bins(time_bins = x)
  
  # Return logical indicating whether object is a valid timeBins object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}
