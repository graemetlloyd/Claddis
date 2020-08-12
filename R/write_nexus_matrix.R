#' Writes out a morphological #NEXUS data file
#'
#' @description
#'
#' Writes out a morphological data file in #NEXUS format.
#'
#' @param cladistic_matrix The cladistic matrix in the format imported by \link{read_nexus_matrix}.
#' @param file_name The file name to write to. Should end in \code{.nex}.
#'
#' @details
#'
#' Writes out a #NEXUS (Maddison et al. 1997) data file representing the distribution of characters in a set of taxa. Data must be in the format created by importing data with \link{read_nexus_matrix}.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @seealso \link{write_tnt_matrix}
#'
#' @references
#'
#' Maddison, D. R., Swofford, D. L. and Maddison, W. P., 1997. NEXUS: an extensible file format for systematic information. \emph{Systematic Biology}, \bold{46}, 590-621.
#'
#' @seealso
#'
#' \link{build_cladistic_matrix}, \link{compactify_matrix}, \link{prune_cladistic_matrix}, \link{read_nexus_matrix}, \link{safe_taxonomic_reduction}, \link{write_tnt_matrix}
#'
#' @examples
#'
#' # Write out Michaux 1989 to current working directory:
#' write_nexus_matrix(cladistic_matrix = michaux_1989, file_name = "michaux_1989.nex")
#'
#' # Remove file when finished:
#' file.remove(file1 = "michaux_1989.nex")
#' @export write_nexus_matrix
write_nexus_matrix <- function(cladistic_matrix, file_name) {

  # Subfunction to convert matrices back to symbols, missing and gap characters:
  convert_matrix <- function(data_matrix) {

    # If there are missing characters replace with missing symbol:
    if (any(is.na(data_matrix$matrix))) data_matrix$matrix[is.na(data_matrix$matrix)] <- data_matrix$characters$missing

    # If there are gap characters replace with gap symbol:
    if (sum(as.vector(data_matrix$matrix) == "") > 0) data_matrix$matrix[data_matrix$matrix == ""] <- data_matrix$characters$gap

    # If there are symbols (i.e., non-continuous data):
    if (length(x = data_matrix$characters$symbols) > 0) {

      # In reverse order go through numbers:
      for (i in rev(data_matrix$characters$symbols)) {

        # Replace current number with appropriate symbol:
        if (length(x = grep(as.character(which(x = data_matrix$characters$symbols == i) - 1), data_matrix$matrix)) > 0) data_matrix$matrix <- gsub(pattern = as.character(which(x = data_matrix$characters$symbols == i) - 1), replacement = i, x = data_matrix$matrix)
      }
    }

    # If there are uncertainties:
    if (length(x = grep("/", data_matrix$matrix)) > 0) {

      # Find cells that have uncertainties:
      uncertainties <- grep("/", data_matrix$matrix)

      # Repalce with all possible values in curly braces:
      data_matrix$matrix[uncertainties] <- paste("{", unlist(x = lapply(X = strsplit(data_matrix$matrix[uncertainties], split = "/"), paste, collapse = "")), "}", sep = "")
    }

    # If there are polymorphisms:
    if (length(x = grep("&", data_matrix$matrix)) > 0) {

      # Find cells with polymorphsims:
      polymorphisms <- grep("&", data_matrix$matrix)

      # Repale with values inside parentheses:
      data_matrix$matrix[polymorphisms] <- paste("(", unlist(x = lapply(X = strsplit(data_matrix$matrix[polymorphisms], split = "&"), paste, collapse = "")), ")", sep = "")
    }

    # Get equal length taxon names (with added spaces):
    taxon_names_with_spaces <- paste(rownames(x = data_matrix$matrix), unlist(x = lapply(X = lapply(X = as.list(x = (max(nchar(x = rownames(x = data_matrix$matrix))) + 2) - nchar(x = rownames(x = data_matrix$matrix))), rep, x = " "), paste, collapse = "")), sep = "")

    # If block is continuous:
    if (data_matrix$datatype == "CONTINUOUS") {

      # Format rows with spaces between values:
      data_matrix$matrix <- paste(taxon_names_with_spaces, (apply(data_matrix$matrix, 1, paste, collapse = " ")), sep = "")

      # If block is non-continuous:
    } else {

      # Format rows without spaces between values:
      data_matrix$matrix <- paste(taxon_names_with_spaces, (apply(data_matrix$matrix, 1, paste, collapse = "")), sep = "")
    }

    # Return just the newly formatted matrix (now a vector):
    return(data_matrix$matrix)
  }

  # Subfunction for collapsing strings of character numbers (i.e., for weight and ordering in assumptions block):
  zip_string <- function(x) {

    # Set up empty zipped version of strin:
    zipped_string <- vector(mode = "character")

    # As long as there are zip-able parts to x remaining:
    while (any(diff(x) == 1)) {

      # Find start of first zippable part:
      zip_start <- which(x = diff(x) == 1)[1]

      # Find length of zippable part:
      zip_length <- rle(diff(x))$lengths[which(x = rle(diff(x))$values == 1)[1]]

      # Use length and start to find end of zip:
      zip_end <- zip_start + zip_length

      # Add zippable part to string (and any preceding unzippable parts):
      zipped_string <- c(zipped_string, setdiff(x = x[1:zip_end], y = x[zip_start:zip_end]), paste(x[zip_start], x[zip_end], sep = "-"))

      # Cut x down ready to find next zippable part:
      x <- x[-(1:zip_end)]
    }

    # If there are characters left (or there are no zippable parts):
    if (length(x = x) > 0) zipped_string <- c(zipped_string, x)

    # Collapse to a single string with spaces:
    zipped_string <- paste(zipped_string, collapse = " ")

    # Return zipped string:
    return(zipped_string)
  }

  # Isolate just data blocks (i.e., cladistic_matrix without topper):
  data_blocks <- cladistic_matrix[2:length(x = cladistic_matrix)]

  # Get block names:
  block_names <- unlist(x = lapply(X = data_blocks, "[[", "block_name"))

  # Get data_types:
  data_types <- unlist(x = lapply(X = data_blocks, "[[", "datatype"))

  # Get number of taxa:
  n_taxa <- unname(unlist(x = lapply(X = lapply(X = data_blocks, "[[", "matrix"), nrow))[1])

  # Get number of characters:
  n_characters <- unlist(x = lapply(X = lapply(X = data_blocks, "[[", "matrix"), ncol))

  # Get symbols strings:
  symbols <- unlist(x = lapply(X = lapply(X = lapply(X = data_blocks, "[[", "characters"), "[[", "symbols"), paste, collapse = " "))

  # Get missing value:
  missing <- unlist(x = lapply(X = lapply(X = data_blocks, "[[", "characters"), "[[", "missing"))

  # Get gap symbol:
  gap <- unlist(x = lapply(X = lapply(X = data_blocks, "[[", "characters"), "[[", "gap"))

  # Conver matrices to vectors of text strings:
  data_block_strings <- lapply(X = data_blocks, convert_matrix)

  # Set up header block (returns empty string if nothing there):
  header_block <- ifelse(nchar(x = cladistic_matrix$topper$header) > 0, paste("[", cladistic_matrix$topper$header, "]\n\n", sep = ""), "")

  # Set up taxa block (only required if multiple matrix blocks as sets number of taxa, will be empty string otherwise):
  taxa_block <- ifelse(length(x = data_blocks) > 1, paste("BEGIN TAXA;\n\tDIMENSIONS NTAX=", n_taxa, ";\n\tTAXLABELS\n\t\t", paste(rownames(x = cladistic_matrix$matrix_1$matrix), collapse = " "), "\n;\nEND;\n\n", sep = ""), "")

  # Set up data block (only required if a single matrix block):
  data_block <- ifelse(length(x = data_blocks) == 1, paste("BEGIN DATA;\n\tDIMENSIONS  NTAX=", n_taxa, " NCHAR=", n_characters, " ;\n\tFORMAT DATATYPE=", data_types, " SYMBOLS=\" ", symbols, "\" MISSING=", missing, " GAP=", gap, " ;\n", sep = ""), "")

  # Set up character block (including MATRIX that will begin data):
  character_block <- ifelse(rep(length(x = data_blocks), length(x = data_blocks)) > 1, paste("BEGIN CHARACTERS;\n\t", ifelse(nchar(x = block_names) > 0, paste("TITLE  ", block_names, ";\n", sep = ""), ""), "\tDIMENSIONS  NCHAR=", n_characters, ";\n\tFORMAT  DATATYPE=", ifelse(data_types == "CONTINUOUS", "CONTINUOUS ", paste(data_types, " SYMBOLS=\" ", symbols, "\" ", sep = "")), "MISSING=", missing, " GAP=", gap, " ;\nMATRIX\n\n", sep = ""), "MATRIX\n\n")

  # Take character block and meld with matri(ces) into matrix block(s):
  matrix_block <- paste(paste(character_block, unlist(x = lapply(X = data_block_strings, paste, collapse = "\n")), "\n;\nEND;\n\n", sep = ""), collapse = "")

  # Make sure step matrices are a list if null:
  if (!is.list(cladistic_matrix$topper$step_matrices)) cladistic_matrix$topper$step_matrices <- list(NULL)

  # Create step matrix block:
  stepmatrix_block <- paste(ifelse(!unlist(x = lapply(X = cladistic_matrix$topper$step_matrices, is.null)), paste(paste("\tUSERTYPE '", names(cladistic_matrix$topper$step_matrices), "' (STEPMATRIX) = ", unlist(x = lapply(X = cladistic_matrix$topper$step_matrices, ncol)), "\n", sep = ""), paste("\t", unlist(x = lapply(X = lapply(X = cladistic_matrix$topper$step_matrices, colnames), paste, collapse = " ")), "\n\t", sep = ""), unlist(x = lapply(X = lapply(X = lapply(X = cladistic_matrix$topper$step_matrices, function(x) {
    diag(x = x) <- "."
    return(x)
  }), apply, 1, paste, collapse = " "), paste, collapse = "\n\t")), "\n\t;\n", sep = ""), ""), collapse = "")

  # Get ordering of all characters in sequence:
  ordering <- unlist(x = lapply(X = data_blocks, "[[", "ordering"))

  # Get weights of all characters in sequence:
  character_weights <- unlist(x = lapply(X = data_blocks, "[[", "character_weights"))

  # Create options block (if no block names):
  if (all(is.na(block_names))) options_block <- paste(ifelse(all(ordering == "unord"), "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n", ifelse(all(ordering == "ord"), "\tOPTIONS  DEFTYPE=ord PolyTcount=MINSTEPS ;\n", "\tOPTIONS  DEFTYPE=unord PolyTcount=MINSTEPS ;\n")), ifelse(length(x = unique(x = ordering)) == 1 && length(x = setdiff(x = unique(x = ordering), y = c("ord", "unord"))) == 0, "", paste("\tTYPESET * UNTITLED  = ", paste(paste(sort(x = unique(x = ordering)), unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = sort(x = unique(x = ordering))), "==", ordering), which), zip_string)), sep = ": "), collapse = ", "), ";\n", sep = "")), sep = "")

  # Create options block (if there are block names):
  if (!all(is.na(unlist(x = block_names)))) options_block <- paste(paste("\tTYPESET * UNTITLED  (CHARACTERS = ", block_names, ")  =  ", unlist(x = lapply(X = lapply(X = data_blocks, "[[", "ordering"), function(x) paste(paste(paste(sort(x = unique(x = x)), unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = sort(x = unique(x = x))), "==", x), which), zip_string)), sep = ": "), collapse = ", "), sep = ""))), ";\n", sep = ""), collapse = "")

  # Replace cont with Squared if continuous characters present:
  if (length(x = grep(" cont: ", options_block)) > 0) options_block <- gsub(pattern = " cont: ", replacement = " Squared: ", x = options_block)

  # Convert continuous character weights to one before making weights block:
  character_weights[ordering == "cont"] <- 1

  # Create weights block (if no block names):
  if (all(is.na(block_names))) weights_block <- ifelse(all(character_weights == 1), "", paste("\tWTSET * UNTITLED  = ", paste(paste(sort(x = unique(x = character_weights)), unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = sort(x = unique(x = character_weights))), "==", character_weights), which), zip_string)), sep = ": "), collapse = ", "), ";\n", sep = ""))

  # Create weights block (if there are block names):
  if (!all(is.na(unlist(x = block_names)))) weights_block <- paste(paste("\tWTSET * UNTITLED  (CHARACTERS = ", block_names, ")  =  ", unlist(x = lapply(X = lapply(X = data_blocks, "[[", "character_weights"), function(x) paste(paste(paste(sort(x = unique(x = x)), unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = sort(x = unique(x = x))), "==", x), which), zip_string)), sep = ": "), collapse = ", "), sep = ""))), ";\n", sep = ""), collapse = "")

  # Build assumptions block:
  assumptions_block <- paste("BEGIN ASSUMPTIONS;\n", stepmatrix_block, options_block, weights_block, "END;\n", sep = "")

  # Build full string with all blocks together:
  full_string <- paste("#NEXUS\n\n", header_block, taxa_block, data_block, matrix_block, assumptions_block, sep = "")

  # Write to file:
  write(full_string, file_name)
}
