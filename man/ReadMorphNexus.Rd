\name{ReadMorphNexus}

\alias{ReadMorphNexus}

\title{
Reads in a morphological #NEXUS data file
}

\description{
Reads in a morphological data file in #NEXUS format.
}

\usage{
ReadMorphNexus(file, equalise.weights=FALSE)
}

\arguments{
\item{file}{A file name specified by either a variable of mode character, or a double-quoted string.}
\item{equalise.weights}{Optional that overrides the weights specified in the file to make all characters truly equally weighted.}
}

\details{
Reads in a #NEXUS (Maddison et al. 1997) data file representing the distribution of discrete morphological characters in a set of taxa. Unlike \link{read.nexus.data} this function can handle polymorphisms (e.g., \code{(012)}).

Note that the function is generally intolerant to excursions from a standard format and it is recommended your data be formatted like the \code{morphmatrix.nex} example below.

All empty values (missing or inapplicable) are treated as NAs.
}

\value{
\item{header}{Any header text included in the file given between square brackets as character.}
\item{matrix}{A matrix of taxa (rows) and characters (columns). The matrix is in character format in order to deal with polymorphisms, which are separated by ampersands.}
\item{ordering}{A character vector of the same length as the number of morphological characters indicating whether they are ordered (\code{ORD}) or unordered (\code{UNORD}).}
\item{weights}{A numeric vector of the same length as the number of morphological characters indicating their weights.}
\item{max.vals}{A numeric vector of the same length as the number of morphological characters indicating the maximum state values.}
\item{min.vals}{A numeric vector of the same length as the number of morphological characters indicating the minimum state values.}
\item{step.matrices}{A list of any step matrices supplied in the input file. Is \code{NULL} if none are specified.}
\item{symbols}{The original symbols used in the input data (these are replaced here by the digits 0-9 and the letters A-V in order).}
}

\references{
Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. Systematic Biology, 46, 590-621.
}

\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}

\seealso{
\link{read.nexus.data}
}

\examples{
# Create example matrix (NB: creates a file in the current
# working directory called "morphmatrix.nex"):
cat("#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS  NTAX=5 NCHAR=5;\n\t
FORMAT SYMBOLS= \" 0 1 2\" MISSING=? GAP=- ;\nMATRIX\n\n
Taxon_1  010?0\nTaxon_2  021?0\nTaxon_3  02111\nTaxon_4  011-1
\nTaxon_5  001-1\n;\nEND;\n\nBEGIN ASSUMPTIONS;\n\t
OPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n\t
TYPESET * UNTITLED  = unord: 1 3-5, ord: 2;\n\t
WTSET * UNTITLED  = 1: 2, 2: 1 3-5;\nEND;", file = "morphmatrix.nex")

# Read in example matrix:
morph.matrix <- ReadMorphNexus("morphmatrix.nex")

# View example matrix in R:
morph.matrix
}

\keyword{ NEXUS }
