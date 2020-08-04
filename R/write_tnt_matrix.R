#' Writes out a morphological TNT data file
#'
#' @description
#'
#' Writes out a morphological data file in Hennig86/TNT format.
#'
#' @param cladistic_matrix A cladistic matrix in the format imported by \link{read_nexus_matrix}.
#' @param file_name The file name to write to. Should end in \code{.tnt}.
#' @param add.analysis.block Whether or not to add analysis block (i.e., tree search commands).
#'
#' @details
#'
#' Writes out a TNT (Goloboff et al. 2008; Goloboff and Catalano 2016) data file representing the distribution of discrete morphological characters in a set of taxa. Data must be in the format created by importing data with \link{read_nexus_matrix}.
#'
#' Note that the format can currently deal with continuous characters, sequence (DNA) data, and combinations of these and discrete morphology, but not yet the morphometric format introduced in Goloboff and Catalano (2016).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso
#'
#' \link{write_nexus_matrix}
#'
#' @references
#'
#' Goloboff, P. A. and Catalano, S. A., 2016. TNT version 1.5, including a full implementation of phylogenetic morphometrics/ \emph{Cladistics}, \bold{32}. 221-238
#'
#' Goloboff, P., Farris, J. and Nixon, K., 2008. TNT, a free program for phylogenetic analysis. \emph{Cladistics}, \bold{24}, 774-786.
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{compactify_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_nexus_matrix}
#'
#' @examples
#'
#' # Write out Michaux 1989 to current working directory:
#' write_tnt_matrix(cladistic_matrix = michaux_1989, file_name = "michaux_1989.tnt")
#'
#' # Remove file when finished:
#' file.remove("michaux_1989.tnt")
#' @export write_tnt_matrix
write_tnt_matrix <- function(cladistic_matrix, file_name, add.analysis.block = FALSE) {

  # Subfunction to convert matrices back to symbols, missing and gap characters:
  convert_matrix <- function(cladistic_matrix) {

    # If there are missing characters replace with missing symbol:
    if (any(is.na(cladistic_matrix$matrix))) cladistic_matrix$matrix[is.na(cladistic_matrix$matrix)] <- "?"

    # If there are gap characters replace with gap symbol:
    if (sum(as.vector(cladistic_matrix$matrix) == "") > 0) cladistic_matrix$matrix[cladistic_matrix$matrix == ""] <- "-"

    # If there are symbols (i.e., non-continuous data):
    if (length(cladistic_matrix$characters$symbols) > 0) {

      # If datatype is STANDARD set TNT symbols:
      if (cladistic_matrix$datatype == "STANDARD") tnt_symbols <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V")

      # If datatype is non-STANDARD (but still discrete):
      if (cladistic_matrix$datatype != "STANDARD") tnt_symbols <- cladistic_matrix$characters$symbols

      # In reverse order go through numbers:
      for (i in rev(tnt_symbols)) {

        # Replace current number with appropriate symbol:
        if (length(grep(as.character(which(tnt_symbols == i) - 1), cladistic_matrix$matrix)) > 0) cladistic_matrix$matrix <- gsub(as.character(which(tnt_symbols == i) - 1), i, cladistic_matrix$matrix)
      }
    }

    # If there are uncertainties:
    if (length(grep("/", cladistic_matrix$matrix)) > 0) {

      # Find cells that have uncertainties:
      Uncertainties <- grep("/", cladistic_matrix$matrix)

      # Replace with all possible values in curly braces:
      cladistic_matrix$matrix[Uncertainties] <- paste("[", unlist(lapply(strsplit(cladistic_matrix$matrix[Uncertainties], split = "/"), paste, collapse = "")), "]", sep = "")
    }

    # If there are polymorphisms:
    if (length(grep("&", cladistic_matrix$matrix)) > 0) {

      # Find cells with polymorphsims:
      Polymorphisms <- grep("&", cladistic_matrix$matrix)

      # Resplae with values inside parentheses:
      cladistic_matrix$matrix[Polymorphisms] <- paste("[", unlist(lapply(strsplit(cladistic_matrix$matrix[Polymorphisms], split = "&"), paste, collapse = "")), "]", sep = "")
    }

    # Get equal length taxon names (with added spaces):
    TaxonNamesWithTrailingSpaces <- paste(rownames(cladistic_matrix$matrix), unlist(lapply(lapply(as.list((max(nchar(rownames(cladistic_matrix$matrix))) + 2) - nchar(rownames(cladistic_matrix$matrix))), rep, x = " "), paste, collapse = "")), sep = "")

    # If block is continuous:
    if (cladistic_matrix$datatype == "CONTINUOUS") {

      # Format rows with spaces between values:
      cladistic_matrix$matrix <- paste(TaxonNamesWithTrailingSpaces, (apply(cladistic_matrix$matrix, 1, paste, collapse = " ")), sep = "")

      # If block is non-continuous:
    } else {

      # Format rows without spaces between values:
      cladistic_matrix$matrix <- paste(TaxonNamesWithTrailingSpaces, (apply(cladistic_matrix$matrix, 1, paste, collapse = "")), sep = "")
    }

    # Return just the newly formatted matrix (now a vector):
    cladistic_matrix$matrix
  }

  # Isolate just data blocks (i.e., cladistic_matrix without topper):
  DataBlocks <- cladistic_matrix[2:length(cladistic_matrix)]

  # Get block names:
  block_names <- unlist(lapply(DataBlocks, "[[", "block_name"))

  # Get datatypes:
  datatypes <- unlist(lapply(DataBlocks, "[[", "datatype"))

  # NEXUS versions of datatypes:
  NEXUSversion <- c("CONTINUOUS", "DNA", "NUCLEOTIDE", "PROTEIN", "RESTRICTION", "RNA", "STANDARD")

  # TNT versions of datatypes:
  TNTversion <- c("continuous", "dna", "dna", "proteins", "numeric", "dna", "numeric")

  # Replace NEXUS versions of datatypes with TNT version (as best as I can guess!):
  for (i in 1:7) datatypes <- gsub(NEXUSversion[i], TNTversion[i], datatypes)

  # Get number of taxa:
  NTaxa <- unname(unlist(lapply(lapply(DataBlocks, "[[", "matrix"), nrow))[1])

  # Get number of characters:
  Ncharacters <- unlist(lapply(lapply(DataBlocks, "[[", "matrix"), ncol))

  # Get symbols strings:
  symbols <- unlist(lapply(lapply(lapply(DataBlocks, "[[", "characters"), "[[", "symbols"), paste, collapse = " "))

  # Get missing value:
  missing <- unlist(lapply(lapply(DataBlocks, "[[", "characters"), "[[", "missing"))

  # Get gap symbol:
  gap <- unlist(lapply(lapply(DataBlocks, "[[", "characters"), "[[", "gap"))

  # Conver matrices to vectors of text strings:
  DataBlocksAsTextStrings <- lapply(DataBlocks, convert_matrix)

  # Set up header block (returns empty string if nothing there):
  headerBlock <- ifelse(length(cladistic_matrix$topper$header) > 0, paste("'", cladistic_matrix$topper$header, "'\n", sep = ""), "")

  # Set up character block (including MATRIX that will begin data):
  CharacterBlock <- paste("& [", datatypes, "]\n", sep = "")

  # Take character block and meld with matri(ces) into matrix block(s):
  MatrixBlock <- paste(paste(CharacterBlock, unlist(lapply(DataBlocksAsTextStrings, paste, collapse = "\n")), "\n", sep = ""), collapse = "")

  # Get ordering of all characters in sequence:
  ordering <- unname(unlist(lapply(DataBlocks, "[[", "ordering")))

  # Get weights of all characters in sequence:
  character_weights <- unname(unlist(lapply(DataBlocks, "[[", "character_weights")))

  # Make sure step matrices are a list if null:
  if (!is.list(cladistic_matrix$topper$step_matrices)) cladistic_matrix$topper$step_matrices <- list(NULL)

  # If there are step matrices:
  if (any(!unlist(lapply(cladistic_matrix$topper$step_matrices, is.null)))) {

    # Empty vector to store hits (characters assigned to a step matrix):
    global_hits <- vector(mode = "numeric")

    # Empty vector to store all step matrix lines:
    all_step_matrix_lines <- vector(mode = "character")

    # Check there are not too many step matrices:
    if (length(cladistic_matrix$topper$step_matrices) > 32) stop("Too many (>32) step matrices for TNT!")

    # For each step matrix:
    for (i in 1:length(cladistic_matrix$topper$step_matrices)) {

      # Set up costs vector:
      costs <- vector(mode = "character")

      # For each row state (from):
      for (j in rownames(cladistic_matrix$topper$step_matrices[[i]])) {

        # For each column state (to):
        for (k in colnames(cladistic_matrix$topper$step_matrices[[i]])) {

          # Add cost of j to k transition to costs vector:
          costs <- c(costs, paste(j, ">", k, " ", cladistic_matrix$topper$step_matrices[[i]][j, k], sep = ""))
        }
      }

      # Format top of step matrix code:
      step_matrix_top <- paste("smatrix = ", i - 1, " (", names(cladistic_matrix$topper$step_matrices)[i], ")", sep = "")

      # Get hits (characters assigned to ith step matrix):
      hits <- which(ordering == names(cladistic_matrix$topper$step_matrices)[i])

      # Add huts to global hits for all step matrices:
      global_hits <- c(global_hits, hits)

      # Stop if no hits!:
      if (length(hits) == 0) stop(paste("No characters assigned to step matrix: ", names(cladistic_matrix$topper$step_matrices)[i], ".", sep = ""))

      # Build step matrix lines:
      step_matrix_lines <- paste(c(step_matrix_top, costs, ";", paste("smatrix + ", names(cladistic_matrix$topper$step_matrices)[i], " ", paste(hits - 1, collapse = " "), ";", sep = "")), collapse = "\n")

      # Add step matrix lines to all step matrix lines:
      all_step_matrix_lines <- c(all_step_matrix_lines, step_matrix_lines)
    }

    # Make step matrix block:
    StepMatrixBlock <- paste(paste(c(all_step_matrix_lines, paste("ccode ( ", paste(sort(x = global_hits - 1), collapse = " "), ";", sep = "")), collapse = "\n"), "\n", sep = "")

    # If no step matri(ces):
  } else {

    # Create empty step matrix block:
    StepMatrixBlock <- ""
  }

  # Build ccode block:
  CCodeBlock <- paste("ccode ", paste(paste(ifelse(ifelse(ordering == "cont", "ord", ordering) == "ord", "+", "-"), ifelse(character_weights == 0, "]", "["), "/", ifelse(ordering == "cont", "1", ifelse(character_weights == 0, 1, character_weights)), " ", paste(1:sum(Ncharacters) - 1), sep = ""), collapse = " "), " ;\n", sep = "")

  # Build full string with all blocks together:
  FullString <- paste("taxname=;\nmxram 4096;\ntaxname +", max(nchar(rownames(cladistic_matrix$matrix_1$Matrix))), ";\nnstates num 32;\nxread\n", headerBlock, sum(Ncharacters), " ", NTaxa, "\n", MatrixBlock, ";\n", CCodeBlock, StepMatrixBlock, "proc/;\n", sep = "")

  # If adding analysis block:
  if (add.analysis.block) {

    # If there are few enough taxa for an exact solution:
    if (NTaxa <= 24) {

      # Use implicit enumeration:
      AnalysisBlock <- c("collapse [;", "ienum;")

      # If there are too many taxa (requiring heuristics):
    } else {

      # First iteration with new tech (where scratch.tre is created):
      AnalysisBlock <- c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "tsave scratch.tre;", "save;", "tsave /;")

      # Iterations 2-20 (where trees are appended to scratch.tre):
      AnalysisBlock <- c(AnalysisBlock, rep(c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "tsave scratch.tre +;", "save;", "tsave /;"), 19))

      # Read in new technology trees (scratch.tre) and finish with a heuristic (tbr) search for all mpts:
      AnalysisBlock <- c(AnalysisBlock, c("hold 100000;", "shortread scratch.tre;", "bbreak=tbr;"))
    }

    # Get stripped file name for use in export lines:
    out.file <- strsplit(strsplit(file_name, "/")[[1]][length(strsplit(file_name, "/")[[1]])], "\\.")[[1]][1]

    # Make name for strict consensus and MPTs tree:
    strict.name <- paste("export -", out.file, "tntmpts_plus_strict.nex;", sep = "")

    # Make MRP file name:
    mrp.name <- c("mrp;", paste("export ", out.file, "mrp.nex;", sep = ""))

    # Add analysis block to output file:
    FullString <- gsub("proc/;\n", paste(paste(AnalysisBlock, collapse = "\n"), "nelsen*;", strict.name, paste(mrp.name, collapse = "\n"), "\nproc/;\n"), FullString)
  }

  # Write to file:
  write(FullString, file_name)
}
