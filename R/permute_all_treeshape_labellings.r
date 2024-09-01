#' Label treeshapes
#'
#' @description
#'
#' Given a treeshape and set of labels, permutes all possible labelled phylogenetic trees.
#'
#' @param treeshapes A vector of treeshape(s) in the same format as \link{permute_treeshapes}.
#' @param labels A character vector of tip labels to use for labelling.
#'
#' @details
#'
#' A treeshape is an unlabelled phylogenetic tree and as such can be labelled to produce a phylogenetic tree. This function takes a treeshape and a set of labels and generates (permutes) all possible labellings, i.e., all phylogenetic trees which a treeshape represents.
#'
#' Note that the star tree always allows only a single labelling, whereas any more resolved treeshape will have multiple labellings.
#'
#' Here treeshapes are encoded in the same pseudo-Newick format as the \link{permute_treeshapes} function, e.g.:
#'
#' (((3),1),(1,(2)));
#'
#' (Where each pair of parentheses represents an internal node, each number the number of tips, each comma separates sets of tips, and the semicolon denotes the root clade.)
#'
#' @return A list of the same length as \code{treeshapes} composed of character vectors of labelled phylogenetic trees in the Newick format (Felsenstein 2004).
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Felsenstein, J., 2004. \emph{Inferring Phylogenies}. Sinauer Associates, Inc., Sunderland.
#'
#' @examples
#'
#' # Label some six-tip treeshapes with the letters A-F:
#' permute_all_treeshape_labellings(
#'   treeshapes = c(
#'     "(6);",
#'     "((3),(3));",
#'     "(1,(1,(1,(1,(2)))));"
#'   ),
#'   labels = LETTERS[1:6]
#' )
#'
#' @export permute_all_treeshape_labellings
permute_all_treeshape_labellings <- function(treeshapes, labels) {
  
  # TO DO:
  # - LIKE PERMUTE TIPSTATES AND PERMUTE TREESHAPES HAS REDUNDANCIES DUE TO REPEATING TREE MOTIFS (COULD BE SPED UP IF A WAY TO DEAL WITH THESE IS IDENTIFIED
  # - OTHER INPUT VALUE CHECKS COULD BE ADDED.
  
  # Check labels are all unique and stop and warn user if not:
  if (any(x = duplicated(x = labels))) stop("labels must not contain any duplicate values.")
  
  # Get numeric components for each treeshape:
  numeric_components <- lapply(
    X = treeshapes,
    FUN = function(i) as.numeric(
      x = strsplit(
        x = gsub(
          pattern = "\\)|\\(|;",
          replacement = "",
          x = i
        ),
        split = ","
      )[[1]]
    )
  )
  
  # Get tip counts for each treeshape:
  tip_counts <- unlist(x = lapply(X = numeric_components, FUN = sum))
  
  # Check treeshapes share a tip count and stop and warn user if this varies:
  if (length(x = unique(x = tip_counts)) > 1) stop("treeshapes are of different sizes (number of tips).")
  
  # Check labels are same size as the tip count and stop and warn user if not:
  if (tip_counts[1] != length(x = labels)) stop("labels must be the same length as the number of tips in treeshapes.")
  
  # Subfunction to permute labels across partition sizes:
  permute_partition_combinations <- function(partition_sizes, labels) {
    
    # Create starting labels matrix:
    labels_matrix <- t(x = utils::combn(x = labels, m = partition_sizes[1]))
    
    # As long as there is more than one clade (i.e., it isn't the star tree):
    if (length(x = partition_sizes) > 1) {
      
      # For each subsequent partition size in sequence:
      for(i in partition_sizes[2:length(x = partition_sizes)]) {
        
        # Increase labels matrix:
        labels_matrix <- do.call(
          what = rbind,
          args = lapply(
            X = apply(
             X = labels_matrix,
             MARGIN = 1,
             FUN = list
            ),
            FUN = function(y) {
              current_labels <- setdiff(x = labels, y = unlist(x = y))
              new_combinations <- t(x = utils::combn(x = current_labels, m = i))
              cbind(do.call(what = rbind, args = rep(x = y, times = nrow(x = new_combinations))), new_combinations)
            }
          )
        )
      }
    }
    
    # Return labels matrix:
    labels_matrix
  }
  
  # Subfunction to build Newick strings from a Newick Number string and associated matrix of label combinations):
  build_labelled_newicks <- function(treeshape, labels_matrix) {
    
    # Set up labelled newick strings by repeating Newick number string required number of times:
    labelled_newicks <- rep(x = treeshape, times = nrow(x = labels_matrix))
    
    # Subfunction to split treeshape into vector of components:
    split_treeshape <- function(treeshape) {
      
      # Isolate nonumeric parts:
      nonnumeric_parts <- strsplit(
        x = treeshape,
        split = "[:0-9:]+"
      )[[1]]
      
      # Isolate numeric parts:
      numeric_parts <- strsplit(
        x = gsub(
          pattern = "\\)|\\(|;",
          replacement = "",
          x = treeshape
        ),
        split = ","
      )[[1]]
     
      # Combine parts in matrix:
      combined_parts_matrix <- matrix(
        data = c(
          nonnumeric_parts,
          numeric_parts,
          NA
        ),
        nrow = 2,
        byrow = TRUE
      )
      
      # Split top row only into component parts (do not want to split numbers!):
      combined_parts_matrix[1, ] <- strsplit(x = combined_parts_matrix[1, ], split = "")
      
      # Return vector of isolated treeshape components:
      unlist(x = combined_parts_matrix[!is.na(x = combined_parts_matrix)])
    }
    
    # Create split strings in list format ready for mapply:
    split_strings <- lapply(X = as.list(x = labelled_newicks), FUN = split_treeshape)
    
    # Create split string labels as list format ready for mapply:
    string_labels <- apply(
      X = labels_matrix,
      MARGIN = 1,
      FUN = function(i) i,
      simplify = FALSE
    )
    
    # Subfunction to add labels to Newick number string:
    label_number_newicks <- function(split_string, string_labels) {
      
      # Find positions of numbers in split string:
      number_positions <- grep(pattern = "[:0-9:]", x = split_string)
      
      # Find number values for split string (N labels required):
      number_values <- as.numeric(x = split_string[number_positions])
      
      # For each number positions (CAN PROBABLY REPLACE THIS WITH A SINGLE LINE SOMEHOW):
      for(i in 1:length(x = number_positions)) {
        
        # Paste string labels over ith number position:
        split_string[number_positions[i]] <- paste(string_labels[1:as.numeric(number_values[i])], collapse = ",")
        
        # Remove used string labels from pool:
        string_labels <- string_labels[-c(1:as.numeric(x = number_values[i]))]
      }
      
      # Convert split string back into a single string:
      paste(split_string, collapse = "")
    }
    
    # Generate all labelled Newick strings from Newick number string and labels:
    labelled_newicks <- mapply(FUN = label_number_newicks, split_string = split_strings, string_labels = string_labels)
    
    # Due to repeated subgraphs (i.e., clades) may need to remove duplicate trees:
    if (length(x = labelled_newicks) > 1) labelled_newicks <- ape::write.tree(
      phy = find_unique_trees(
        trees = ape::read.tree(
          text = labelled_newicks
        )
      )
    )
    
    # Return labelled newicks:
    labelled_newicks
  }
  
  # Generate and return all labelled trees:
  mapply(
    FUN = build_labelled_newicks,
    treeshape = as.list(x = treeshapes),
    labels_matrix = lapply(
      X = numeric_components,
      FUN = permute_partition_combinations,
      labels = labels
    ),
    SIMPLIFY = FALSE
  )
}
