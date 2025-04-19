#' Permute all treeshapes of N tips
#'
#' @description
#'
#' Given a number of tips, permutes all rooted unlabelled multifurcating trees (i.e., treeshapes).
#'
#' @param n_tips The number of tips required. Note that it may be very slow or not run at all if this value is too large.
#' @param sort_by_resolution Whether or not to sort the output by number of internal nodes (from 1 to N - 1). Defaults to \code{TRUE}.
#'
#' @details
#'
#' A treeshape is essentially an unlabelled phylogenetic tree. Like other phylogenetic trees it has a root and tips, but (as you might expect) because it is unlabelled the tips have no specific identity. Thus the only information it contains is its' "shape" - the number of internal nodes and their descendants. This function permutes all \emph{unique} treeshapes and allows for multifurcations.
#'
#' Note that unique means it excludes alternative rotations of individual branch points. For example, the trees ((2),1); and (1,(2)); are identical in information content and this function would only permute one of them.
#'
#' The algorithm used here is based on the partitions approach from Felsenstein (2004), although to the best of my knowledge nobody else has formally created an algorithm to do this. (Felsenstein also lays out the expected number of such treeshapes for each value of N in his Table 3.4.)
#'
#' Here treeshapes are encoded and output in a pseudo-Newick style format where labels are replaced with the number of tips, e.g.:
#'
#' (((3),1),(1,(2)));
#'
#' Thus each pair of parentheses represents an internal node, each number the number of tips, each comma separates sets of tips, and the semicolon denotes the root clade.
#'
#' @return If \code{sort_by_resolution = TRUE} then returns a list of length N - 1, where each element is a character vector of treeshapes in Newick-style number format with that many internal nodes. I.e., the first value will always be the star tree. If \code{sort_by_resolution = FALSE} then will just be a character vector of treeshapes in Newick-style number format.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Felsenstein, J., 2004. \emph{Inferring Phylogenies}. Sinauer Associates, Inc., Sunderland.
#'
#' @examples
#'
#' # Permute all treeshapes of six tips sorted by resolution:
#' permute_treeshapes(n_tips = 6)
#'
#' @export permute_treeshapes
permute_treeshapes <- function(n_tips, sort_by_resolution = TRUE) {

  # Subfunction to add outer braces:
  add_outer_braces <- function(x) paste("(", x, ")", sep = "")
  
  # Subfunction to join elements with commas:
  comma_join <- function(x) paste(x, collapse = ",")
  
  # Subfunction to split elements by commas:
  comma_split <- function(x) {
    
    # Split x into individual characters:
    full_split <- strsplit(x = x, split = "")[[1]]
    
    # Find finishing point:
    finishing_point <- length(x = full_split)
    
    # Create empty starting vectors of clade and non clade components:
    clade_components <- non_clade_components <- nonclade <- c()

    # Initialise depth counter at zero:
    depth_counter <- 0
    
    # Step through x one character at a time:
    for(position in 1:finishing_point) {
      
      # If position is an opening parenthesis:
      if (full_split[position] == "(") {
        
        # Check if a new clade component and initialise start clade if so:
        if (depth_counter == 0) start_clade <- position
        
        # Increment depth counter:
        depth_counter <- depth_counter + 1
      }
      
      # If position is a closing parenthesis:
      if (full_split[position] == ")") {
        
        # Decrement depth:
        depth_counter <- depth_counter - 1
        
        # Check if a clade has been closed and add to clade components if so:
        if (depth_counter == 0) clade_components <- c(clade_components, paste(full_split[start_clade:position], collapse = ""))
      }
      
      # If position is a numeric:
      if (length(x = grep(pattern = "[0-9]{1}", x = full_split[position])) == 1) {
        if (depth_counter == 0) nonclade <- c(nonclade, position)
      }
    }
    
    # If there is a nonclade component then format this correctly:
    if (length(x = nonclade) > 0) non_clade_components <- c(non_clade_components, paste(full_split[nonclade], collapse = ""))
    
    # Return tree strng split by commas (of same lowest depth):
    c(clade_components, non_clade_components)
  }
  
  # Subfunction to get partitions (tree components) for n tips:
  get_tree_components <- function(n_tips) {
    apply(
      X = partitions::restrictedparts(n = n_tips, m = n_tips - 1)[, -1, drop = FALSE],
      MARGIN = 2,
      FUN = function(i) {
        
        # Prune out zero sized partitions:
        i <- i[i > 0]
      
        # Form clades from partitions of at least 2:
        component_parts <- add_outer_braces(x = i[i > 1])
      
        # Add any "grades" outside clades:
        if (any(i == 1)) component_parts <- c(component_parts, sum(x = i[i == 1]))
      
        # Return component parts:
        component_parts
      },
      simplify = FALSE
    )
  }

  # Get size of a clade component, e.g., (3) returns 3:
  get_component_size <- function(x) as.numeric(x = gsub(pattern = "\\(|\\)", replacement = "", x = x))
  
  # Subfunction to recursively resolve each partition of tree_string:
  resolve_partition <- function(tree_string) {
    
    # Split tree string apart:
    tree_components <- comma_split(x = tree_string)
    
    # Get just clade ("(N)") components from tree_components:
    clade_components <- tree_components[grep(pattern = "\\(", x = tree_components)]
    
    # Get unique clade components (only ones that need permuting):
    unique_clade_components <- unique(x = clade_components)
    
    # Get just clade ("(N)") components from tree_components:
    non_clade_components <- tree_components[grep(pattern = "\\(", x = tree_components, invert = TRUE)]

    # Get just clade ("(N)") components from tree_components:
    clade_component_sizes <- sapply(X = clade_components, FUN = get_component_size)

    # Only need to continue permuting if there are clade(s) of three or more tips:
    if (any(clade_component_sizes > 2)) {
      
      # Permuted unique_clade_components for usage later:
      permuted_unique_clade_components <- lapply(
        X = sapply(
          X = unique_clade_components,
          FUN = get_component_size,
          simplify = FALSE
        ),
        FUN = function(i) {
          add_outer_braces(
            x = c(
              paste(i, sep = ""),
              #unlist(x = lapply(X = get_tree_components(n_tips = i), FUN = comma_join))
              unlist(
                x = lapply(
                  X = get_tree_components(n_tips = i),
                  FUN = function(i) resolve_partition(tree_string = comma_join(x = i))
                )
              )
            )
          )
        }
      )
      
      # If there are any clade components of equal size:
      if (any(duplicated(x = clade_components))) {
        
        # Get just the duplicated clade components:
        duplicated_clade_components <- unique(x = clade_components[duplicated(x = clade_components)])
        
        # Get any non duplicated clade components (for seprately permuting later):
        non_duplicated_clade_components <- setdiff(x = unique(x = clade_components), y = duplicated_clade_components)
        
        # Separately permute just the duplicated components (to avoid issue of rotation symmetries):
        permuted_duplicate_components <- sapply(
          X = duplicated_clade_components,
          FUN = function(i) {
            n_duplicates <- sum(x = tree_components == i)
            combinations <- permute_combinations_with_replacement(x = permuted_unique_clade_components[[i]], m = n_duplicates)
            apply(X = combinations, MARGIN = 1, FUN = comma_join)
          },
          simplify = FALSE
        )
        
        # Permute all possible tree components:
        tree_components <- apply(
          X = expand.grid(
            x = c(
              permuted_duplicate_components,
              permuted_unique_clade_components[non_duplicated_clade_components],
              non_clade_components
            ),
            stringsAsFactors = FALSE
          ),
          MARGIN = 1,
          FUN = comma_join
        )
        
      # If no clade components are equal in size:
      } else {
        
        # Permute all possible tree components:
        tree_components <- apply(
          X = expand.grid(
            c(
              permuted_unique_clade_components[clade_components],
              non_clade_components
            ),
            stringsAsFactors = FALSE
          ),
          MARGIN = 1,
          FUN = comma_join
        )
      }
      
    # If there is nothing to permute:
    } else {
      
      # Format output as list:
      tree_components <- list(c(clade_components, non_clade_components))
    }
    
    # Return components comma separated:
    unlist(x = lapply(X = tree_components, FUN = comma_join))
  }
  
  # Get tree shapes for i tips:
  treeshapes <- resolve_partition(tree_string = add_outer_braces(x = n_tips))
  
  # Add closing semicolons:
  treeshapes <- paste(treeshapes, ";", sep = "")
  
  # If treeshapes should be output by resolution (i.e., n internal nodes):
  if (sort_by_resolution) {
    
    # Get count of internal nodes for each treeshape:
    n_nodes <- sapply(
      X = treeshapes,
      FUN = function(i) nchar(x = i) - nchar(x = gsub(pattern = "\\(", replacement = "", x = i))
    )
    
    # Overwrite treeshapes with list grouping treeshapes by numebr of internal nodes (from 1 to n tips  - 1):
    treeshapes <- sapply(X = 1:(n_tips - 1), FUN = function(i) treeshapes[n_nodes == i])
  }
  
  # Return treeshapes:
  treeshapes
}
