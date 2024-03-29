% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_tnt_matrix.R
\name{write_tnt_matrix}
\alias{write_tnt_matrix}
\title{Writes out a morphological TNT data file}
\usage{
write_tnt_matrix(cladistic_matrix, file_name, add_analysis_block = FALSE)
}
\arguments{
\item{cladistic_matrix}{A cladistic matrix in the format imported by \link{read_nexus_matrix}.}

\item{file_name}{The file name to write to. Should end in \code{.tnt}.}

\item{add_analysis_block}{Whether or not to add analysis block (i.e., tree search commands).}
}
\description{
Writes out a morphological data file in Hennig86/TNT format.
}
\details{
Writes out a TNT (Goloboff et al. 2008; Goloboff and Catalano 2016) data file representing the distribution of discrete morphological characters in a set of taxa. Data must be in the format created by importing data with \link{read_nexus_matrix}.

Note that the format can currently deal with continuous characters, sequence (DNA) data, and combinations of these and discrete morphology, but not yet the morphometric format introduced in Goloboff and Catalano (2016).
}
\examples{

# Write out Michaux 1989 to current working directory:
write_tnt_matrix(cladistic_matrix = michaux_1989, file_name = "michaux_1989.tnt")

# Remove file when finished:
file.remove(file1 = "michaux_1989.tnt")
}
\references{
Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics. \emph{Cladistics}, \bold{32}, 221-238.

Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
}
\seealso{
\link{write_nexus_matrix}

\link{build_cladistic_matrix}, \link{compactify_cladistic_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}
}
\author{
Graeme T. Lloyd \email{graemetlloyd@gmail.com}
}
