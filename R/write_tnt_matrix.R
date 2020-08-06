#' Writes out a morphological TNT data file
#'
#' @description
#'
#' Writes out a morphological data file in Hennig86/TNT format.
#'
#' @param cladistic_matrix A cladistic matrix in the format imported by \link{read_nexus_matrix}.
#' @param file_name The file name to write to. Should end in \code{.tnt}.
#' @param add_analysis_block Whether or not to add analysis block (i.e., tree search commands).
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
#' file.remove(file1 = "michaux_1989.tnt")
#' @export write_tnt_matrix
write_tnt_matrix <- function(cladistic_matrix, file_name, add_analysis_block = FALSE) {

  # Subfunction to convert matrices back to symbols, missing and gap characters:
  convert_matrix <- function(cladistic_matrix) {

    # If there are missing characters replace with missing symbol:
    if (any(is.na(cladistic_matrix$matrix))) cladistic_matrix$matrix[is.na(cladistic_matrix$matrix)] <- "?"

    # If there are gap characters replace with gap symbol:
    if (sum(as.vector(cladistic_matrix$matrix) == "") > 0) cladistic_matrix$matrix[cladistic_matrix$matrix == ""] <- "-"

    # If there are symbols (i.e., non-continuous data):
    if (length(x = cladistic_matrix$characters$symbols) > 0) {

      # If datatype is STANDARD set TNT symbols:
      if (cladistic_matrix$datatype == "STANDARD") tnt_symbols <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V")

      # If datatype is non-STANDARD (but still discrete):
      if (cladistic_matrix$datatype != "STANDARD") tnt_symbols <- cladistic_matrix$characters$symbols

      # In reverse order go through numbers:
      for (i in rev(tnt_symbols)) {

        # Replace current number with appropriate symbol:
        if (length(x = grep(as.character(which(x = tnt_symbols == i) - 1), cladistic_matrix$matrix)) > 0) cladistic_matrix$matrix <- gsub(pattern = as.character(which(x = tnt_symbols == i) - 1), replacement = i, x = cladistic_matrix$matrix)
      }
    }

    # If there are uncertainties:
    if (length(x = grep("/", cladistic_matrix$matrix)) > 0) {

      # Find cells that have uncertainties:
      uncertainties <- grep("/", cladistic_matrix$matrix)

      # Replace with all possible values in curly braces:
      cladistic_matrix$matrix[uncertainties] <- paste("[", unlist(x = lapply(X = strsplit(cladistic_matrix$matrix[uncertainties], split = "/"), FUN = paste, collapse = "")), "]", sep = "")
    }

    # If there are polymorphisms:
    if (length(x = grep("&", cladistic_matrix$matrix)) > 0) {

      # Find cells with polymorphsims:
      polymorphisms <- grep("&", cladistic_matrix$matrix)

      # Resplae with values inside parentheses:
      cladistic_matrix$matrix[polymorphisms] <- paste("[", unlist(x = lapply(X = strsplit(cladistic_matrix$matrix[polymorphisms], split = "&"), paste, collapse = "")), "]", sep = "")
    }

    # Get equal length taxon names (with added spaces):
    taxon_names_with_spaces <- paste(rownames(x = cladistic_matrix$matrix), unlist(x = lapply(X = lapply(X = as.list(x = (max(nchar(x = rownames(x = cladistic_matrix$matrix))) + 2) - nchar(x = rownames(x = cladistic_matrix$matrix))), FUN = rep, x = " "), FUN = paste, collapse = "")), sep = "")

    # If block is continuous:
    if (cladistic_matrix$datatype == "CONTINUOUS") {

      # Format rows with spaces between values:
      cladistic_matrix$matrix <- paste(taxon_names_with_spaces, (apply(cladistic_matrix$matrix, 1, paste, collapse = " ")), sep = "")

      # If block is non-continuous:
    } else {

      # Format rows without spaces between values:
      cladistic_matrix$matrix <- paste(taxon_names_with_spaces, (apply(cladistic_matrix$matrix, 1, paste, collapse = "")), sep = "")
    }

    # Return just the newly formatted matrix (now a vector):
    cladistic_matrix$matrix
  }

  # Isolate just data blocks (i.e., cladistic_matrix without topper):
  data_blocks <- cladistic_matrix[2:length(x = cladistic_matrix)]

  # Get block names:
  block_names <- unlist(x = lapply(X = data_blocks, FUN = "[[", "block_name"))

  # Get data_types:
  data_types <- unlist(x = lapply(X = data_blocks, FUN = "[[", "datatype"))

  # NEXUS versions of data_types:
  nexus_datatypes <- c("CONTINUOUS", "DNA", "NUCLEOTIDE", "PROTEIN", "RESTRICTION", "RNA", "STANDARD")

  # TNT versions of data_types:
  tnt_datatypes <- c("continuous", "dna", "dna", "proteins", "numeric", "dna", "numeric")

  # Replace NEXUS versions of data_types with TNT version (as best as I can guess!):
  for (i in 1:7) data_types <- gsub(pattern = nexus_datatypes[i], replacement = tnt_datatypes[i], x = data_types)

  # Get number of taxa:
  n_taxa <- unname(unlist(x = lapply(X = lapply(X = data_blocks, FUN = "[[", "matrix"), FUN = nrow))[1])

  # Get number of characters:
  n_characters <- unlist(x = lapply(X = lapply(X = data_blocks, FUN = "[[", "matrix"), FUN = ncol))

  # Get symbols strings:
  symbols <- unlist(x = lapply(X = lapply(X = lapply(X = data_blocks, FUN = "[[", "characters"), FUN = "[[", "symbols"), FUN = paste, collapse = " "))

  # Get missing value:
  missing <- unlist(x = lapply(X = lapply(X = data_blocks, FUN = "[[", "characters"), "[[", "missing"))

  # Get gap symbol:
  gap <- unlist(x = lapply(X = lapply(X = data_blocks, FUN = "[[", "characters"), FUN = "[[", "gap"))

  # Conver matrices to vectors of text strings:
  data_block_strings <- lapply(X = data_blocks, FUN = convert_matrix)

  # Set up header block (returns empty string if nothing there):
  header_block <- ifelse(length(x = cladistic_matrix$topper$header) > 0, paste("'", cladistic_matrix$topper$header, "'\n", sep = ""), "")

  # Set up character block (including MATRIX that will begin data):
  character_block <- paste("& [", data_types, "]\n", sep = "")

  # Take character block and meld with matri(ces) into matrix block(s):
  matrix_block <- paste(paste(character_block, unlist(x = lapply(X = data_block_strings, paste, collapse = "\n")), "\n", sep = ""), collapse = "")

  # Get ordering of all characters in sequence:
  ordering <- unname(unlist(x = lapply(X = data_blocks, "[[", "ordering")))

  # Get weights of all characters in sequence:
  character_weights <- unname(unlist(x = lapply(X = data_blocks, "[[", "character_weights")))

  # Make sure step matrices are a list if null:
  if (!is.list(cladistic_matrix$topper$step_matrices)) cladistic_matrix$topper$step_matrices <- list(NULL)

  # If there are step matrices:
  if (any(!unlist(x = lapply(X = cladistic_matrix$topper$step_matrices, is.null)))) {

    # Empty vector to store hits (characters assigned to a step matrix):
    global_hits <- vector(mode = "numeric")

    # Empty vector to store all step matrix lines:
    all_step_matrix_lines <- vector(mode = "character")

    # Check there are not too many step matrices:
    if (length(x = cladistic_matrix$topper$step_matrices) > 32) stop("Too many (>32) step matrices for TNT!")

    # For each step matrix:
    for (i in 1:length(x = cladistic_matrix$topper$step_matrices)) {

      # Set up costs vector:
      costs <- vector(mode = "character")

      # For each row state (from):
      for (j in rownames(x = cladistic_matrix$topper$step_matrices[[i]])) {

        # For each column state (to):
        for (k in colnames(x = cladistic_matrix$topper$step_matrices[[i]])) {

          # Add cost of j to k transition to costs vector:
          costs <- c(costs, paste(j, ">", k, " ", cladistic_matrix$topper$step_matrices[[i]][j, k], sep = ""))
        }
      }

      # Format top of step matrix code:
      step_matrix_top <- paste("smatrix = ", i - 1, " (", names(cladistic_matrix$topper$step_matrices)[i], ")", sep = "")

      # Get hits (characters assigned to ith step matrix):
      hits <- which(x = ordering == names(cladistic_matrix$topper$step_matrices)[i])

      # Add huts to global hits for all step matrices:
      global_hits <- c(global_hits, hits)

      # Stop if no hits!:
      if (length(x = hits) == 0) stop(paste("No characters assigned to step matrix: ", names(cladistic_matrix$topper$step_matrices)[i], ".", sep = ""))

      # Build step matrix lines:
      step_matrix_lines <- paste(c(step_matrix_top, costs, ";", paste("smatrix + ", names(cladistic_matrix$topper$step_matrices)[i], " ", paste(hits - 1, collapse = " "), ";", sep = "")), collapse = "\n")

      # Add step matrix lines to all step matrix lines:
      all_step_matrix_lines <- c(all_step_matrix_lines, step_matrix_lines)
    }

    # Make step matrix block:
    stepmatrix_block <- paste(paste(c(all_step_matrix_lines, paste("ccode ( ", paste(sort(x = global_hits - 1), collapse = " "), ";", sep = "")), collapse = "\n"), "\n", sep = "")

    # If no step matri(ces):
  } else {

    # Create empty step matrix block:
    stepmatrix_block <- ""
  }

  # Build ccode block:
  ccode_block <- paste("ccode ", paste(paste(ifelse(ifelse(ordering == "cont", "ord", ordering) == "ord", "+", "-"), ifelse(character_weights == 0, "]", "["), "/", ifelse(ordering == "cont", "1", ifelse(character_weights == 0, 1, character_weights)), " ", paste(1:sum(n_characters) - 1), sep = ""), collapse = " "), " ;\n", sep = "")

  # Build full string with all blocks together:
  full_string <- paste("taxname=;\nmxram 4096;\ntaxname +", max(nchar(x = rownames(x = cladistic_matrix$matrix_1$matrix))), ";\nnstates num 32;\nxread\n", header_block, sum(n_characters), " ", n_taxa, "\n", matrix_block, ";\n", ccode_block, stepmatrix_block, "proc/;\n", sep = "")

  # If adding analysis block:
  if (add_analysis_block) {

    # If there are few enough taxa for an exact solution:
    if (n_taxa <= 24) {

      # Use implicit enumeration:
      analysis_block <- c("collapse [;", "ienum;")

      # If there are too many taxa (requiring heuristics):
    } else {

      # First iteration with new tech (where scratch.tre is created):
      analysis_block <- c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "tsave scratch.tre;", "save;", "tsave /;")

      # Iterations 2-20 (where trees are appended to scratch.tre):
      analysis_block <- c(analysis_block, rep(c("rseed*;", "hold 999;", "xmult=rss fuse 50 drift 50 ratchet 50;", "tsave scratch.tre +;", "save;", "tsave /;"), 19))

      # Read in new technology trees (scratch.tre) and finish with a heuristic (tbr) search for all mpts:
      analysis_block <- c(analysis_block, c("hold 100000;", "shortread scratch.tre;", "bbreak=tbr;"))
    }

    # Get stripped file name for use in export lines:
    output_file_name <- strsplit(strsplit(file_name, "/")[[1]][length(x = strsplit(file_name, "/")[[1]])], "\\.")[[1]][1]

    # Make name for strict consensus and MPTs tree:
    scc_file_name <- paste("export -", output_file_name, "tntmpts_plus_strict.nex;", sep = "")

    # Make MRP file name:
    mrp_file_name <- c("mrp;", paste("export ", output_file_name, "mrp.nex;", sep = ""))

    # Add analysis block to output file:
    full_string <- gsub(pattern = "proc/;\n", replacement = paste(paste(analysis_block, collapse = "\n"), "nelsen*;", scc_file_name, paste(mrp_file_name, collapse = "\n"), "\nproc/;\n"), x = full_string)
  }

  # Write to file:
  write(full_string, file_name)
}
