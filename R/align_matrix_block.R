#' Aligns a phylogenetic matrix block
#'
#' @description
#'
#' Given a block of taxa and characters aligns text so each character block begins at same point.
#'
#' @param matrix_block The matrix block as raw input text.
#'
#' @details
#'
#' The function serves to help build NEXUS files by neatly aligning raw text blocks of taxa and characters. Or in simple terms it takes input that looks like this:
#'
#' \preformatted{Allosaurus  012100?1011
#' Abelisaurus  0100???0000
#' Tyrannosaurus  01012012010
#' Yi  10101?0????}
#'
#' And turns it into something that looks like this:
#'
#' \preformatted{Allosaurus     012100?1011
#' Abelisaurus    0100???0000
#' Tyrannosaurus  01012012010
#' Yi             10101?0????}
#'
#' I use this in building the NEXUS files on my site, \href{http://www.graemetlloyd.com/matr.html}{graemetlloyd.com}.
#'
#' @return
#'
#' Nothing is returned, instead the aligned block is sent to the clipboard ready for pasting into a text editor.
#'
#' @author
#'
#' Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Build example block from above:
#' x <- paste(c(
#'   "Allosaurus  012100?1011",
#'   "Abelisaurus  0100???0000",
#'   "Tyrannosaurus  01012012010",
#'   "Yi  10101?0????"
#' ), collapse = "\n")
#'
#' # Look at block pre-alignment:
#' x
#'
#' # Align block and place on clipboard:
#' \dontrun{
#' align_matrix_block(x)
#' }
#'
#' # To test the response open a text editor and paste the
#' # contents of the clipboard.
#' @export align_matrix_block
align_matrix_block <- function(matrix_block) {

  # Need to convert supplied block of text into list of taxon and character vectors:
  matrix_block <- lapply(X = as.list(x = strsplit(matrix_block, "\n")[[1]]), function(x) {

    # Split each line by whitespace:
    x <- unlist(x = strsplit(x, " "))

    # Return vector of name plus characters:
    x[c(1, length(x = x))]
  })

  # what is the most number of spaces to add:
  block_length <- max(unlist(x = lapply(X = matrix_block, function(x) nchar(x = x[1])))) + 2

  # Add spaces to names to align block:
  matrix_block <- lapply(X = matrix_block, function(x) {

    # Isolate taxon name:
    taxon_name <- strsplit(x[1], "")[[1]]

    # Paste line together with correct number of spaces separating taxon name and characters:
    x[1] <- paste(c(taxon_name, rep(" ", block_length - length(x = taxon_name))), collapse = "")

    # Return aligned text:
    x
  })

  # Write output to clipboard ready to paste in a text (NEXUS) file:
  clipr::write_clip(paste(unlist(x = lapply(X = matrix_block, paste, collapse = "")), collapse = "\n"))
}
