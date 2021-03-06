% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/align_matrix_block.R
\name{align_matrix_block}
\alias{align_matrix_block}
\title{Aligns a phylogenetic matrix block}
\usage{
align_matrix_block(matrix_block)
}
\arguments{
\item{matrix_block}{The matrix block as raw input text.}
}
\value{
Nothing is returned, instead the aligned block is sent to the clipboard ready for pasting into a text editor.
}
\description{
Given a block of taxa and characters aligns text so each character block begins at same point.
}
\details{
The function serves to help build NEXUS files by neatly aligning raw text blocks of taxa and characters. Or in simple terms it takes input that looks like this:

\preformatted{Allosaurus  012100?1011
Abelisaurus  0100???0000
Tyrannosaurus  01012012010
Yi  10101?0????}

And turns it into something that looks like this:

\preformatted{Allosaurus     012100?1011
Abelisaurus    0100???0000
Tyrannosaurus  01012012010
Yi             10101?0????}

I use this in building the NEXUS files on my site, \href{http://www.graemetlloyd.com/matr.html}{graemetlloyd.com}.
}
\examples{

# Build example block from above:
x <- paste(c(
  "Allosaurus  012100?1011",
  "Abelisaurus  0100???0000",
  "Tyrannosaurus  01012012010",
  "Yi  10101?0????"
), collapse = "\n")

# Look at block pre-alignment:
x

# Align block and place on clipboard:
\dontrun{
align_matrix_block(x)
}

# To test the response open a text editor and paste the
# contents of the clipboard.
}
\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
