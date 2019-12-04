#' Align a phylogenetic matrix block
#'
#' @description
#'
#' Given a block of taxa and characters aligns start of character block.
#'
#' @param Block The raw input text.
#'
#' @details
#'
#' The function serves to help build NEXUS files by neatly aligning block of taxa and characters. Or in simple terms it takes inut that looks like this:
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
#' I use this in building the NEXUS on my site, \href{http://www.graemetlloyd.com/matr.html}{graemetlloyd.com}.
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
#' x <- paste(c("Allosaurus  012100?1011",
#'   "Abelisaurus  0100???0000",
#'   "Tyrannosaurus  01012012010",
#'   "Yi  10101?0????"), collapse = "\n")
#'
#' # Look at block pre-alignment:
#' x
#'
#' # Example block from above:
#' #BlockAligner(x)
#'
#' # To test the response open a text editor and paste the
#' # contents of the clipboard.
#'
#' @export BlockAligner
BlockAligner <- function(Block) {
  
  # Build matrix of input data:
  Block <- lapply(as.list(strsplit(Block, "\n")[[1]]), function(y) {y <- unlist(strsplit(y, " ")); y <- y[c(1, length(y))]; y})
  
  # Work out how many spaces to add:
  BlockLength <- max(unlist(lapply(Block, function(z) nchar(z[1])))) + 2
  
  # Add spaces to names to align block:
  Block <- lapply(Block, function(z) {TaxonName <- strsplit(z[1], "")[[1]]; z[1] <- paste(c(TaxonName, rep(" ", BlockLength - length(TaxonName))), collapse = ""); z})
  
  # Write output to clipboard ready to paste in a text (NEXUS) file:
  clipr::write_clip(paste(unlist(lapply(Block, paste, collapse = "")), collapse = "\n"))
  
}
