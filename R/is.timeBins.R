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
#' # Check that this is a valid timeBins object (will fail as class is not set):
#' is.timeBins(x = time_bins)
#'
#' # Set class as timeBins:
#' class(time_bins) <- "timeBins"
#'
#' # Check that this is a valid timeBins object (will succeed as format and
#' # class are correct):
#' is.timeBins(x = time_bins)
#'
#' @export is.timeBins
is.timeBins <- function(x) {
  
  # Get any error messages for time_bins:
  messages <- check_timeBins(time_bins = x)
  
  # Return logical indicating whether object is a valid timeBins object or not:
  ifelse(test = length(x = messages) > 0, yes = FALSE, no = TRUE)
}
