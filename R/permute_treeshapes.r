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
#' # Permute all treeshapes of six tips:
#' permute_treeshapes(n_tips = 6)
#'
#' @export permute_treeshapes
permute_treeshapes <- function(n_tips, sort_by_resolution = TRUE) {

  # If one of fewer tips stop and warn user:
  if (n_tips <= 1) stop("n_tips must be at least 2.")
  
  # Life is easy if only two tips were asked for:
  if (n_tips == 2) return("(2);")
  
  # As still broken above eight tips stop and warn user if requested:
  if (n_tips > 8) stop("Currently this function does not permute all treeshapes for tip counts above eight.")
  
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

        # Read into ape as trees:
        trees <- ape::read.tree(text = multihit_treeshapes)
        
        # If only one tree convert to a list so below wll work:
        if (inherits(x = trees, what = "phylo")) trees <- list(trees)
        
        # Permute all new multihit treeshapes:
        new_multihit_treeshapes <- lapply(
        
          # For each tree with multiple hits for a clade of size i:
          X = trees,
          FUN = function(tree) {
            
            # Isolate the nodes corresponding to clades:
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
            
            # Get frequencies of nodes (to help identify the berries of size i):
            node_freqs <- rle(x = sort(x = clade_nodes))
            
            # Isolate just the berry nodes of size i:
            berry_nodes <- node_freqs$values[node_freqs$lengths == numeric_i]
            
            # Find the ancestors of each size i berry (to use for identifying symmetries):
            berry_ancestors <- unlist(
              x = lapply(
                X = as.list(x = berry_nodes),
                FUN = function(berry) tree$edge[tree$edge[, 2] == berry, 1]
              )
            )
            
            # Label the nodes of the tree (will need later for gsub to work), excludes "A" to avoid conflict with tip labels:
            tree$node.label <- make_labels(N = (tree$Nnode + 1))[2:(tree$Nnode + 1)]
            
            # Make vector of permutation elements (what matters for symmetry is how these are permuted not what):
            permutation_elements <- c(i, i_treeshapes)
            
            # Convert permutation elements to "labelled" versions:
            for(j in n_tips:1) permutation_elements <- gsub(
              pattern = as.character(x = j),
              replacement = paste(rep(x = "A", times = j), collapse = ","),
              x = permutation_elements
            )
            
            # Set up a list to store items to expand.grid for full permutation set:
            grid_list <- list()
            
            # If any symmetries are found (i.e., cases where berries of size i share an ancestor):
            if (any(x = duplicated(x = berry_ancestors))) {
              
              # Get frequencies of each berry ancestor node:
              ancestor_freqs <- rle(x = sort(x = berry_ancestors))
              
              # Isolate just those ancestors with symmetric (i.e., multiple berry) descendants:
              symmetric_ancestors <- ancestor_freqs$values[ancestor_freqs$lengths > 1]
              
              # For each set of symmetric ancestors:
              for(j in symmetric_ancestors) {
                
                # Get the numbers of the descendant (berry) nodes:
                j_descendants_number <- berry_nodes[berry_ancestors == j]
                
                # Get the labels of those nodes:
                j_descendants_label <- tree$node.label[(j_descendants_number - n_tips)]
                
                # Permute the parts:
                permuted_parts <- do.call(
                  what = rbind,
                  args = apply(
                    X = partitions::restrictedparts(n = length(x = j_descendants_number), m = length(x = permutation_elements)),
                    MARGIN = 2,
                    FUN = function(x) {
                      multicool::allPerm(mcObj = multicool::initMC(x = x))
                    },
                    simplify = FALSE
                  )
                )
                
                # Add permutation_elements as column names to permuted parts to give them meaning:
                colnames(x = permuted_parts) <- permutation_elements
                
                # Now convert permuted parts to strings and store in symmetric_parts:
                grid_list[[(length(x = grid_list) + 1)]] <- apply(
                  X = permuted_parts,
                  MARGIN = 1,
                  FUN = function(k)
                    paste(
                      paste(
                        unlist(
                          x = lapply(
                            X = as.list(x = 1:ncol(x = permuted_parts)),
                            FUN = function(l) rep(
                              x = colnames(x = permuted_parts)[l],
                              times = k[l]
                            )
                          )
                        ),
                        j_descendants_label,
                        sep = ""
                      ),
                    collapse = "&"
                  )
                )
              }
            }
            
            # Find any asymmetric ancestors:
            asymmetric_ancestors <- setdiff(x = berry_ancestors, y = berry_ancestors[duplicated(x = berry_ancestors)])
            
            # As long as there are asymmetric ancestors:
            if (length(x = asymmetric_ancestors) > 0) {
              
              # For each set of symmetric ancestors:
              for(j in asymmetric_ancestors) {
                
                # Get the numbers of the descendant (berry) nodes:
                j_descendants_number <- berry_nodes[berry_ancestors == j]
                
                # Get the labels of those nodes:
                j_descendants_label <- tree$node.label[(j_descendants_number - n_tips)]
                
                # Add permutations to grid_list (with node label appended):
                grid_list[[(length(x = grid_list) + 1)]] <- paste(permutation_elements, j_descendants_label, sep = "")
              }
            }
            
            # Generate final set of multihit permutations:
            multihit_permutations <- apply(
              X = expand.grid(grid_list),
              MARGIN = 1,
              FUN = function(j) unname(
                obj = unlist(
                  x = strsplit(
                    x = j,
                    split = "&"
                  )
                )
              ),
              simplify = FALSE
            )
            
            # Find patterns to replace with permutations:
            pattern_to_replace <- paste(
              "(",
              paste(rep(x = "A", times = numeric_i), collapse = ","),
              ")",
              sort(x = tree$node.label[(berry_nodes - n_tips)]),
              sep = ""
            )
            
            # Generate Newick strings for each permutation:
            newick_strings <- unlist(
              x = lapply(
                X = multihit_permutations,
                FUN = function(j) {
                  names(x = j) <- unlist(x = lapply(X = as.list(x = j), FUN = function(k) rev(x = strsplit(x = k, split = "\\)")[[1]])[1]))
                  j <- j[sort(x = names(x = j))]
                  newick_string <- ape::write.tree(phy = tree)
                  for(k in 1:length(x = j)) newick_string <- gsub(
                    pattern = pattern_to_replace[k],
                    replacement = j[k],
                    x = newick_string,
                    fixed = TRUE
                  )
                  newick_string
                }
              )
            )
            
            # Prune out starting tree (just want new ones):
            newick_strings <- setdiff(x = newick_strings, ape::write.tree(phy = tree))
            
            # Read back in as phylo objects:
            new_trees <- ape::read.tree(text = newick_strings)
            
            # Strip node labels out:
            new_trees <- lapply(X = new_trees, FUN = function(k) {k$node.label <- NULL; k})
            
            # Reset class as multiPhylo:
            class(new_trees) <- "multiPhylo"
            
            # Output Newick strings:
            newick_strings <- ape::write.tree(phy = new_trees)
            
            # Convert back to numeric strings:
            for(j in n_tips:1) newick_strings <- gsub(
              pattern = paste(rep(x = "A", times = j), collapse = ","),
              replacement = as.character(x = j),
              x = newick_strings
            )
            
            # Output numeric newick strings:
            newick_strings
          }
        )
        
        # Restructure output as list of single treeshapes:
        new_multiple_hit_treeshapes <- as.list(x = unlist(x = new_multihit_treeshapes))
      }
      
      # Add new treeshapes to treeshapes list:
      treeshapes <- c(treeshapes, new_single_hit_treeshapes, new_multiple_hit_treeshapes)
    }
  }
  
  # Subfunction to encode treeshapes:
  encode_treeshape <- function(treeshape) {
    
    # Calculate number of tips:
    n_tips <- sum(x = as.numeric(x = strsplit(x = gsub(pattern = "\\)|\\(|;", replacement = "", x = treeshape), split = ",")[[1]]))
    
    # Enter a single dummy label for each tip to make it readable as regular Newick string:
    for(i in n_tips:1) treeshape <- gsub(pattern = as.character(x = i), replacement = paste(rep(x = "A", times = i), collapse = ","), x = treeshape)
    
    # Store as regular phylo object:
    tree <- ape::read.tree(text = treeshape)
    
    # Find all (candidate) penultimate nodes:
    penultimate_nodes <- rle(x = sort(x = unlist(x = lapply(X = 1:n_tips, FUN = function(i) tree$edge[tree$edge[, 2] == i, 1]))))

    # Reduce to just those with at least two tips ("berries")
    berries <- penultimate_nodes$values[penultimate_nodes$lengths > 1]
    
    # Confirm these are termini by removing any that also have internal node descendants (removes if not):
    berries <- berries[unlist(x = lapply(X = berries, FUN = function(i) !any(tree$edge[tree$edge[, 1] == i, 2] > n_tips)))]
    
    # Store root node number:
    root_node <- n_tips + 1
    
    # Get N descendants for each node in order from 1:N (used as a lookup later):
    n_descendants <- ape::node.depth(phy = tree)
    
    # Get list of termini strings:
    termini_strings <- lapply(
    
      # For each starter berry node:
      X = (n_tips + 1):(n_tips + tree$Nnode),
      FUN = function(i) {
        
        # Get initial descendant count:
        descendant_counts <- n_descendants[i]
        
        # Only need to continue of not a star tree:
        if (i > root_node) {
          
          # Set starting current node as i:
          current_node <- i
          
          # Find ancestor of current node:
          ancestor_node <- tree$edge[tree$edge[, 2] == current_node, 1]
          
          # Add descendant count of ancestor node to end of vector:
          descendant_counts <- c(descendant_counts, n_descendants[ancestor_node])
          
          # Update current node:
          current_node <- ancestor_node
          
          # Until the root is reached:
          while(current_node > root_node) {
            
            # Update ancestor ode with ancestor of current node:
            ancestor_node <- tree$edge[tree$edge[, 2] == current_node, 1]
            
            # Add descednant count to end of vector:
            descendant_counts <- c(descendant_counts, n_descendants[ancestor_node])
            
            # Update current node:
            current_node <- ancestor_node
          }
        }
        
        # If not a berry then reverse descendant_counts to indicate this:
        if (length(x = setdiff(x = i, y = berries)) > 0) descendant_counts <- rev(x = descendant_counts)
        
        # Now encode successive tip-to-root descendant counts as single string with descendant counts separated by a hyphen:
        paste(descendant_counts, collapse = "-")
      }
    )
    
    # Return a string that uniquely encodes treeshape (i.e., any rotation will give same encoding):
    paste(sort(x = unlist(x = termini_strings)), collapse = "&")
  }
  
  # Generate encodings of trees:
  encodings <- unlist(
    x = lapply(
      X = treeshapes,
      FUN = function(shape) encode_treeshape(shape)
    )
  )
  
  # Remove any duplictaes treeshapes (dont know why these exist :/):
  treeshapes <- treeshapes[!duplicated(x = encodings)]
  
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
