#' Generate all treeshapes of N tips
#'
#' @description
#'
#' Given a number of tips, generates all unlabelled multifurcating trees (i.e., treeshapes).
#'
#' @param n_tips The number of tips required. Note that it may be very slow or not run at all if this value is too large.
#' @param sort_by_resolution Whther or not to sort the output by number of internal nodes (from 1 to N - 1). Defaults to \code{TRUE}.
#'
#' @details
#'
#' A treeshape is essentially an unlabelled phylogenetic tree. Like other phylogenetic trees it has a root and tips, but (as you might expect) because it is unlabelled the tips have no specific identity. Thus the only information it contains is its' "shape" - the number of internal nodes and their descendants. This function generates all \emph{unique} treeshapes and allows for multifurcations.
#'
#' Note that unique means it excludes alternative rotations of individual branch points. For example, the trees ((2),1); and (1,(2)); are identical in information content and this function will only generate one of them.
#'
#' The algorithm used here is loosely based on the partitions approach suggested in Felsenstein (2004), although to the best of my knowledge nobody else has formally created an algorithm to do this. (Felsenstein also lays out the expected number of such treeshapes for each value of N.)
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
#' # Generate all treeshapes of six tips:
#' generate_all_treeshapes(n_tips = 6)
#'
#' @export generate_all_treeshapes
generate_all_treeshapes <- function(n_tips, sort_by_resolution = TRUE) {

  # If one of fewer tips stop and warn user:
  if (n_tips <= 1) stop("n_tips must be at least 2.")
  
  # Life is easy if only two tips were asked for:
  if (n_tips == 2) return("(2);")
  
  # Subfunction to make initial tree partitions of n_tips:
  make_partitions <- function(n_tips) {
  
    # Make all (valid as trees) partitions of n_tips:
    apply(
      X = partitions::restrictedparts(n = n_tips, m = n_tips - 1), # No point even generating the all ones partition as this can never be arranged into a valid tree
      MARGIN = 2,
      FUN = function(x) x[x > 0] # Remove the zeroes to generate a list output
    )
  }
  
  # Subfunction to add root to a tree string:
  add_root_to_treestring <- function(treestring) paste(treestring, ";", sep = "")
  
  # Convert starting partitions into starting treeshapes:
  convert_partitions_to_treeshapes <- function(partitions, add_root = FALSE) {
    lapply(
      X = partitions,
      FUN = function(x) {
      
        # Add parentheses around partitions of at least two (unless the star phylogeny)
        if (length(x = x) > 1) x[x > 1] <- paste("(", x[x > 1], ")", sep = "")
      
        # Gather singles into a paraphyletic "clump":
        outside <- as.character(x = sum(x == "1"))
      
        # If there are no singles build treeshape without them:
        if (outside == "0") x <- paste("(", paste(c(x[x != "1"]), collapse = ","), ")", sep = "")
      
        # If there are singles add them into treeshape:
        if (outside != "0") x <- paste("(", paste(c(x[x != "1"], outside), collapse = ","), ")", sep = "")
        
        # If adding root "();" then do so:
        if (add_root) x <- add_root_to_treestring(treestring = x)
      
        # Return Newick style treeshape:
        x
      }
    )
  }
  
  # Build first set of partitions:
  partitions <- make_partitions(n_tips = n_tips)
  
  # Build treeshapes from partitions:
  treeshapes <- convert_partitions_to_treeshapes(partitions = partitions, add_root = TRUE)
  
  # Only need to continue finding treeshapes if more than three tips:
  if (n_tips > 3) {
    
    # Define remaining elements to permute new treeshapes from (in decreasing size order):
    permutable_elements <- paste("(", (n_tips - 1):3, ")", sep = "")
    
    # For each permutable elements from "(N)" to "(3)":
    for(i in permutable_elements) {
    
      # Identify treeshapes that contain the ith element:
      i_hits <- grep(pattern = i, x = treeshapes, fixed = TRUE)
          
      # Get i as a numeric value:
      numeric_i <- as.numeric(x = gsub(pattern = "\\)|\\(", replacement = "", x = i))
      
      # Make all partitions of i, except the star tree (as must already have this treeshape by definition):
      i_partitions <- make_partitions(n_tips = numeric_i)[-1]
      
      # Convert partitions of i to treeshape (sans root):
      i_treeshapes <- convert_partitions_to_treeshapes(partitions = i_partitions)
      
      # Identify split sizes (will be used to identify repeated elements of i):
      split_sizes <- unlist(
        x = lapply(
          X = treeshapes[i_hits],
          FUN = function(x) length(x = strsplit(x = x, split = i, fixed = TRUE)[[1]])
        )
      )
      
      # Identify just the single hits:
      single_hits <- i_hits[split_sizes <= 2]
      
      # Identify just the multiple hits:
      multiple_hits <- i_hits[split_sizes > 2]
      
      # Create empty new treeshape lists:
      new_single_hit_treeshapes <- new_multiple_hit_treeshapes <- list()

      # If there is at least one single hit:
      if (length(x = single_hits) > 0) {
      
        # Make new single hit treeshapes:
        new_single_hit_treeshapes <- as.list(x = unlist(x = lapply(
          X = i_treeshapes,
          FUN = function(j) gsub(
            pattern = i,
            replacement = j,
            x = unlist(x = treeshapes[single_hits]),
            fixed = TRUE
          )
        )))
      }
      
      # If there is at least one multiple hit:
      if (length(x = multiple_hits) > 0) {
        
        # Isolate treeshape with multiple hits for a berry of size i:
        multihit_treeshapes <- unlist(x = treeshapes[multiple_hits])
        
        # Add a single dummy label for each tip to make it readable as a regular Newick string:
        for(j in n_tips:1) {
          multihit_treeshapes <- gsub(
            pattern = as.character(x = j),
            replacement = paste(rep(x = "A", times = j), collapse = ","),
            x = multihit_treeshapes
          )
        }

        # Read innto ape as trees:
        trees <- ape::read.tree(text = multihit_treeshapes)
        
        # If only one tree convert to a list so below wll work:
        if (class(trees) == "phylo") trees <- list(trees)
        
        
        lapply(
          X = trees,
          FUN = function(tree) {
            clade_nodes <- tree$edge[apply(
              X = do.call(
                what = cbind,
                args = lapply(
                  X = as.list(x = 1:n_tips),
                  FUN = function(tip) {
                    tree$edge[, 2] == tip
                  }
                )
              ),
              MARGIN = 1,
              FUN = any
            ), 1]
            node_freqs <- rle(x = sort(x = clade_nodes))
            berry_nodes <- node_freqs$values[node_freqs$lengths == numeric_i]
            berry_ancestors <- unlist(x = lapply(X = as.list(x = berry_nodes), FUN = function(berry) tree$edge[tree$edge[, 2] == berry, 1]))
            
            
            #make_labels(32)
            #tree$node.label <- LETTERS[2:(tree$Nnode + 1)] # SOMETHING LIKE THIS:
            
            #NEED TO LABEL NODES AND SPIT OUT THAT STRING TO KNOW WHAT REPLACEMENTS ARE!
            #actual things to permute: c(i, unlist(x = i_treeshapes))
          }
        )
        
        ### THIS IS THE HARD BIT!
        ### NEED TO DISTINGUISH BETWEEN "SYMMETRIC" AND NON-SYMMETRIC MULTIPLE HITS
        
        # If not "symmetric" then expand.grid
        # If "symmetric" need to do unordered permutations
        # If both then symmetric permutations must be done first? Add to treeshapes then do like above? Possibly bad if lots of hits...
        
        # Update below once there is actal output!:
        new_multiple_hit_treeshapes <- list()
      }
      
      # Add new treeshapes to treeshapes list:
      treeshapes <- c(treeshapes, new_single_hit_treeshapes, new_multiple_hit_treeshapes)
    }
  }
  
  # If sorting by resolution:
  if (sort_by_resolution) {
    
    # Get count of node numbers for each treeshape:
    n_internal_nodes <- unlist(
      x = lapply(
        X = treeshapes,
        FUN = function(x) length(x = strsplit(
          x = x,
          split = "(",
          fixed = TRUE
        )[[1]]) - 1
      )
    )
    
    # Gather treeshapes by resolution from 1 to N - 1 internal nodes:
    treeshapes <- lapply(
      X = 1:(n_tips - 1),
      FUN = function(x) unlist(x = treeshapes[n_internal_nodes == x])
    )
    
  # If not sorting by resolution:
  } else {
    
    # Convert treeshapes to a vector:
    treeshapes <- unlist(x = treeshapes)
  }
  
  # Return treeshapes:
  return(treeshapes)
}
