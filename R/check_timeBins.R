#' Check timeBins object for errors
#'
#' @description
#'
#' Internal function to check timeBins object for errors.
#'
#' @param time_bins A timeBins object.
#'
#' @details
#'
#' Internal Claddis function. Nothing to see here. Carry on.
#'
#' @return An error message or empty vector if no errors found.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a time bins object:
#' time_bins <- matrix(
#'   data = c(99.6, 93.5, 93.5, 89.3, 89.3, 85.8, 85.8, 83.5, 83.5, 70.6, 70.6, 65.5),
#'   ncol = 2,
#'   byrow = TRUE,
#'   dimnames = list(
#'     c("Cenomanian", "Turonian", "Coniacian", "Santonian", "Campanian", "Maastrichtian"),
#'     c("fad", "lad")
#'   )
#' )
#'
#' # Check that this is a valid timeBins object (will return error message as class
#' # is not set):
#' check_timeBins(time_bins = time_bins)
#'
#' @export check_timeBins
check_timeBins <- function(time_bins) {
  
  # Check time_bins has class timeBins and add error message to output if true:
  if (!inherits(x = time_bins, what = "timeBins")) return("time_bins must be an object of class \"timeBins\".")
  
  # Check time bins are in form of matrix add error message to output if false:
  if (!is.matrix(x = time_bins)) return("time_bins must be in the form of a matrix.")
  
  # Check time_bins has two columns and add error message to output if false:
  if (ncol(x = time_bins) != 2) return("time_bins must have exactly two columns.")
  
  # Check time_bins has mutliple rows (bins) and add error message to output if false:
  if (nrow(x = time_bins) < 2) return("time_bins must have at least two rows.")
  
  # Check time_bins column names are null and add error message to output if true:
  if (is.null(x = colnames(x = time_bins))) return("time_bins must have columns named \"fad\" and \"lad\" in that order.")
  
  # Check time_bins row names are null and add error message to output if true:
  if (is.null(x = rownames(x = time_bins))) return("time_bins must have unique row names corresponding to each time bin.")
  
  # Check time_bins column names are correct add error message to output if false:
  if (!all(colnames(x = time_bins) == c("fad", "lad"))) return("time_bins must have columns named \"fad\" and \"lad\" in that order.")
  
  # Check time_bins row names have at least one character each and add error message to output if false:
  if (!all(x = nchar(x = rownames(x = time_bins)) > 0)) return("time_bins must have row names formed from at least one character.")
  
  # Check time_bins row names are all unique and add error message to output if false:
  if (any(x = duplicated(x = rownames(x = time_bins)))) return("time_bins must have unique row names.")
  
  # Check time_bins contains only numeric values and add error message to output if false:
  if (!is.numeric(x = time_bins)) return("time_bins must contain numeric values representing millions of years ago (Ma).")
  
  # Check time_bins start before they end and have positive length and add error message to output if false:
  if (!all(x = time_bins[, "fad"] > time_bins[, "lad"])) return("time_bins must consist of \"fad\" values that exceed corresponding \"lad\" values. I.e., each time bin must have an older \"fad\" than \"lad\" and be of positive length.")
  
  # Check time_bins do not overlap each other and add error message to output if true:
  if (any(x = unlist(x = lapply(X = apply(X = time_bins, MARGIN = 1, FUN = list), FUN = function(x) length(x = which(x = apply(X = cbind(x[[1]]["fad"] > time_bins[, "lad"], x[[1]]["lad"] < time_bins[, "fad"]), MARGIN = 1, FUN = all))))) > 1)) return("time_bins must not overlap each other.")
  
  # Check time_bins are in correct order and add error message to output if false:
  if (!all(x = order(x = time_bins[, "fad"]) == seq(from = nrow(x = time_bins), to = 1, by = -1))) return("time_bins must be ordered from oldest (top) to youngest (bottom).")
  
  # Check time_bins have no gaps between them and add error message to output if true:
  if (!all(x = time_bins[1:(nrow(x = time_bins) - 1), "lad"] == time_bins[2:nrow(x = time_bins), "fad"])) return("time_bins must abut each other. I.e., there can be no gaps between successive time bins.")
  
  # Return empty vector:
  vector(mode = "character")
}
