#' Compact display of time bins
#'
#' @description
#'
#' Displays a compact summary of a timeBins object.
#'
#' @param x An object of class \code{"timeBins"}.
#'
#' @details
#'
#' Displays some basic summary information on a time bins object, including number of bins and their names and timespans.
#'
#' @return
#'
#' Nothing is directly returned, instead a text summary describing a \code{"timeBins"} object is printed to the console.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Create a time bins object:
#' time_bins <- matrix(data = c(99.6, 93.5, 93.5, 89.3, 89.3, 85.8, 85.8, 83.5, 83.5, 70.6, 70.6, 65.5), ncol = 2, byrow = TRUE, dimnames = list(c("Cenomanian", "Turonian", "Coniacian", "Santonian", "Campanian", "Maastrichtian"), c("fad", "lad")))
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Show print.cladisticMatrix version of each included data sets:
#' print.timeBins(x = time_bins)
#' @export print.timeBins
print.timeBins <- function(x, ...) {
  
  # Check time_bins has class timeBins and stop and warn user if not:
  if (!inherits(x = time_bins, what = "timeBins")) stop("time_bins must be an object of class \"timeBins\".")
  
  # If not a valid timeBins object then stop and provide feedback to user on what is wrong:
  if (!is.timeBins(time_bins)) stop(check_time_bins(time_bins = time_bins)[1])
  
  # Return summary information about object:
  cat(paste0("timeBins object composed of ", nrow(x = time_bins), " bins:"), "\n ", paste0(unname(obj = unlist(x = apply(X = cbind(rownames(x = time_bins), time_bins), MARGIN = 1, FUN = function(x) paste0(x[1], " (", x[2], "-", x[3], " Ma)")))), collapse = "\n  "))
}
