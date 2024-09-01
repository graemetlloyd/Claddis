#' Permute all tip states on a tree
#'
#' @description
#'
#' Given a phylogenetic tree and a set of states, permutes all possible tip state combinations.
#'
#' @param tree A phylogenetic tree in \link[ape]{ape} format.
#' @param states A vector of character states.
#' @param all_states_present A logical indicating whether or not all states should appear at least once. (Defaults to \code{TRUE}.)
#'
#' @details
#'
#' When calculating g_max or otherwise determining the limits for a given character type it can be useful to generate (permute) all possible combinations of tip states for a given tree. For example, let's imagine we have a binary character (states 0 and 1) and four tips (labelled A-D). We could simply assign each possible state to each possible label to get the following sixteen permutations:
#'
#' \preformatted{______________________
#' | #  | A | B | C | D |
#' ----------------------
#' | 1  | 0 | 0 | 0 | 0 |
#' | 2  | 1 | 0 | 0 | 0 |
#' | 3  | 0 | 1 | 0 | 0 |
#' | 4  | 0 | 0 | 1 | 0 |
#' | 5  | 0 | 0 | 0 | 1 |
#' | 6  | 1 | 1 | 0 | 0 |
#' | 7  | 1 | 0 | 1 | 0 |
#' | 8  | 1 | 0 | 0 | 1 |
#' | 9  | 0 | 1 | 1 | 0 |
#' | 10 | 0 | 1 | 0 | 1 |
#' | 11 | 0 | 0 | 1 | 1 |
#' | 12 | 1 | 1 | 1 | 0 |
#' | 13 | 1 | 1 | 0 | 1 |
#' | 14 | 1 | 0 | 1 | 1 |
#' | 15 | 0 | 1 | 1 | 1 |
#' | 16 | 1 | 1 | 1 | 1 |
#' ----------------------}
#'
#' We can know there are sixteen permutations before we even check as logically for any N states we simply multiply T times (for T tips). I.e., the answer is given by N^T. As here N is 2 and T is 4 this is 2^4 = 2 x 2 x 2 x 2 = 16.
#'
#' Technically this achieves our goal - we have generated all possible tip states, but there are at least two reasons this approach is suboptimal. Firstly, permutations 1 and 16 are invariant (all tip states are the same) and so are unlikely to be of interest. Or, to generalise this, we might not want to consider permutations where not all possible states are sampled - i.e., we may wish to stipulate that every state appear at least once. Secondly, this approach does not consider the structure of the phylogenetic tree of interest. It might not be immediately obvious why this matters, so let's consider another example, Sticking with our binary character and four tips let's consider the tree:
#'
#' \preformatted{A   B   C   D
#'  \   \   \ /
#'   \   \   /
#'    \   \ /
#'     \   /
#'      \ /}
#'
#' And the following permutations (from above):
#'
#' \preformatted{______________________
#' | #  | A | B | C | D |
#' ----------------------
#' | 4  | 0 | 0 | 1 | 0 |
#' | 5  | 0 | 0 | 0 | 1 |
#' | 7  | 1 | 0 | 1 | 0 |
#' | 8  | 1 | 0 | 0 | 1 |
#' | 9  | 0 | 1 | 1 | 0 |
#' | 10 | 0 | 1 | 0 | 1 |
#' | 12 | 1 | 1 | 1 | 0 |
#' | 13 | 1 | 1 | 0 | 1 |
#' ----------------------}
#'
#' Here are the eight permutations where C and D are assigned different states. However, in practical terms there is no need to generate half of these permutations as the order CD/DC doesn't matter (due to the principal of free rotation). In other words, the "identity" (i.e., the label of a tip) is not as important as the structure of the tree.
#'
#' Thus, there are two ways in which the simple ordered permutation method inefficiently generates permutations that are redundant with respect to a given tree. This is important as scaling the simple assignment of states to labels (N^T) can rapidly lead to very large numbers of permutations:
#'
#' \preformatted{           ___________________________________________________________________
#'            |                          Number of tips                         |
#' ------------------------------------------------------------------------------
#' | N States |  5    |    10     |       20      |      50      |      100     |
#' ------------------------------------------------------------------------------
#' |     2    |    32 |     1,024 |     1,048,576 | 1.13 x 10^15 | 1.27 x 10^30 |
#' |     3    |   243 |    59,049 | 3,486,784,401 | 7.18 x 10^23 | 5.15 x 10^47 |
#' |     4    | 1,024 | 1,048,576 |  1.10 x 10^12 | 1.27 x 10^30 | 1.61 x 10^60 |
#' |     5    | 3,125 | 9,765,625 |  9.54 x 10^13 | 8.88 x 10^34 | 7.89 x 10^69 |
#' ------------------------------------------------------------------------------}
#'
#' Thus for many empirically realistic examples the computational expense is too great and a more complex, less redundant, permutation algorithm is desirable that directly considers both tree structure and the desire for full variance sampling.
#'
#' As this is not a simple algorithm let's begin by consider a simple and limiting case - the star tree:
#'
#' \preformatted{A   B   C   D
#' |   |   |   |
#' |   |   |   |
#' -------------}
#'
#' Here there is no "normal" tree structure, but in practical terms this means there is no reason to consider \emph{which} label a tip state is assigned to. If we again consider our binary character from above this means by pruning invariant and redundant partitions we get just three permutations:
#'
#' \preformatted{______________________
#' | #  | A | B | C | D |
#' ----------------------
#' | 2  | 1 | 0 | 0 | 0 |
#' | 6  | 1 | 1 | 0 | 0 |
#' | 12 | 1 | 1 | 1 | 0 |
#' ----------------------}
#'
#' Put simply, these are the cases where there is one, two, or three of state 1, and three, two or one of state 0. Thus we have dramatically reduced the number of permutations, here by over 80%. The number of permutations in this case is now just T - 1, which scales much better:
#'
#' \preformatted{           __________________________
#'            |      Number of tips    |
#' -------------------------------------
#' | N States | 5 | 10 | 20 | 50 | 100 |
#' -------------------------------------
#' |     2    | 4 |  9 | 19 | 49 |  99 |
#' -------------------------------------}
#'
#' However, things become a little more complex when we consider a third state. It also makes more sense to start tabulating these permutations by the number of times a state appears instead of arbitrarily assigning states to a particular label. So our binary example would look like this:
#'
#' \preformatted{_________
#' | 0 | 1 |
#' ---------
#' | 3 | 1 |
#' | 2 | 2 |
#' | 1 | 3 |
#' ---------}
#'
#' If we now add a third possible state (2) we can permute all the possible scenarios (excluding those where a state doesn't appear at all) as:
#'
#' \preformatted{_____________
#' | 0 | 1 | 2 |
#' -------------
#' | 2 | 1 | 1 |
#' | 1 | 2 | 1 |
#' | 1 | 1 | 2 |
#' -------------}
#'
#' Thus for this example there are still only three options. This is because stipulating that every state appears at least once means there is only one "free" assignment that can be allocated as either state 0, state 1 or state 2 - accounting for the three permutations. In other words, the assignments here had to be 2, 1, and 1, the only distinction is the "order" of these assignments. Let's consider a fifth tip but stick with a star tree. Now the assignments can be 3, 1, and 1 or 2, 2, and 1:
#'
#' \preformatted{_____________
#' | 0 | 1 | 2 |
#' -------------
#' | 3 | 1 | 1 |
#' | 1 | 3 | 1 |
#' | 1 | 1 | 3 |
#' | 2 | 2 | 1 |
#' | 2 | 1 | 2 |
#' | 1 | 2 | 2 |
#' -------------}
#'
#' Thus by increasing the number of tips by one in this instance we have also doubled the permutations. Note, though, that we do not need an additional column to represent this, just additional rows. Adding a sixth tip we get:
#'
#' \preformatted{_____________
#' | 0 | 1 | 2 |
#' -------------
#' | 4 | 1 | 1 |
#' | 1 | 4 | 1 |
#' | 1 | 1 | 4 |
#' | 3 | 2 | 1 |
#' | 3 | 1 | 2 |
#' | 2 | 3 | 1 |
#' | 1 | 3 | 2 |
#' | 1 | 2 | 3 |
#' | 2 | 1 | 3 |
#' | 2 | 2 | 2 |
#' -------------}
#'
#' Now we have three possible assignments (4-1-1, 3-2-1 and 2-2-2), but find the ways to permute these are no longer equal (3, 6 and 1, respectively). This makes it harder to derive an equation that will allow us to scale the problem to any value of N or T. We can repeat our table from above though to get a general sense of scale. Remember that this is for the star tree with the additional stipulation that every state appears at least once:
#'
#' \preformatted{           __________________________________________
#'            |             Number of tips             |
#' -----------------------------------------------------
#' | N States |  5 |  10 |    20 |      50 |       100 |
#' -----------------------------------------------------
#' |     2    |  4 |   9 |    19 |      49 |        99 |
#' |     3    |  6 |  36 |   171 |   1,176 |     4,851 |
#' |     4    |  4 |  84 |   969 |  18,424 |   156,849 |
#' |     5    |  1 | 126 | 3,876 | 211,876 | 3,764,376 |
#' -----------------------------------------------------}
#'
#' We see, then, that this scales much more reasonably than our initial approach, although this is the "best case" scenario of no tree structure. We can also relax things a little by allowing states to not always appear:
#'
#' \preformatted{           ______________________________________________
#'            |               Number of tips               |
#' ---------------------------------------------------------
#' | N States |   5 |    10 |     20 |      50 |       100 |
#' ---------------------------------------------------------
#' |     2    |   6 |    11 |     21 |      51 |       101 |
#' |     3    |  21 |    66 |    231 |   1,326 |     5,151 |
#' |     4    |  56 |   286 |  1,771 |  23,426 |   176,851 |
#' |     5    | 126 | 1,001 | 10,626 | 316,251 | 4,598,126 |
#' ---------------------------------------------------------}
#'
#' These are larger, but not intolerably so.
#'
#' In practice things can become more complex when the full variety of possible treeshapes are considered and in fact the present algorithm is not without the possibility of redundancies. Specifically, because trees can contain repeating "motifs" that share a common ancestor fewer permutations would still cover the full range of possibilities. For example, consider the tree:
#'
#' \preformatted{A   B   C   D   E   F
#'  \   \ /     \   \ /
#'   \   /       \   /
#'    \ /         \ /
#'     \           /
#'      \         /
#'       \       /
#'        \     /
#'         \   /
#'          \ /}
#'
#' Here the motifs:
#'
#' \preformatted{A   B   C          D   E   F
#'  \   \ /            \   \ /
#'   \   /     and      \   /
#'    \ /                \ /}
#'
#' Are the same as such so would be the permutations:
#'
#' \preformatted{_________________________
#' | A | B | C | D | E | F |
#' -------------------------
#' | 0 | 1 | 1 | 1 | 0 | 0 |
#' | 1 | 0 | 0 | 0 | 1 | 1 |
#' -------------------------}
#'
#' Currently the function does not account for these, nor is a closed form equation known that would generate the true smallest number of unique permtations.
#'
#' Instead, it simply permutes each "berry" (a node where all descendants are terminals, an extension of the notion pf phylogenetic "cherries") as though it is the star tree, and permuting all other tips as simple labellings.
#'
#' @return A matrix where each row is a unique suite of states and each column a lebelled tip of the input tree.
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @examples
#'
#' # Permute four tip states (A, C, G, T) on a six-taxon tree:
#' permute_tipstates(
#'   tree = ape::read.tree(text = "((A,B),(C,(D,(E,F))));"),
#'   states = c("A", "C", "G", "T"),
#'   all_states_present = TRUE
#' )
#'
#' # Permute three tip states (0, 1, 2) on a ten-taxon star tree:
#' permute_tipstates(
#'   tree = ape::read.tree(text = "(A,B,C,D,E,F,G,H,I);"),
#'   states = c("0", "1", "2"),
#'   all_states_present = TRUE
#' )
#'
#' @export permute_tipstates
permute_tipstates <- function(tree, states, all_states_present = TRUE) {

  # TO DO:
  # - Add additional input checks including for tree, tips, and any additional input variables.
  # - Can we make an equation that generates number of permutations?
  
  # Check input takes the form of a vector and stop and warn user if not:
  if (!is.vector(x = states)) stop("states must be a vector.")
  
  # Check if states are fromatted as characters and if not then coerce to be characters:
  if (!is.character(x = states)) states <- as.character(x = states)
  
  # Check tip labels are unique and stop and warn user if not:
  if (any(duplicated(x = tree$tip.label))) stop("Tip labels must be unique.")

  # Store number of tips:
  n_tips <- ape::Ntip(tree)
  
  # Sanity check that there are at least two tips:
  if (n_tips < 2) stop("Not a valid tree! Must have at east two tips.")
  
  # Store number of nodes:
  n_nodes <- tree$Nnode
  
  # Store node numbers:
  node_numbers <- (n_tips + 1):(n_tips + n_nodes)
  
  # Store number of states:
  n_states <- length(x = states)
  
  # Little check that all_states_present is even possible:
  if (all_states_present && n_states > n_tips) {
    
    # Warn user that they asked for something impossible:
    cat("It is not possible to return all states as number of states is greater than number of tips!\n")
    
    # Turn all_states_present to FALSE (as impossible to comply):
    all_states_present <- FALSE
  }
  
  # Case if only one state (return tip state matrix immediately - no need for further analysis):
  if (n_states == 1) return(matrix(data = states, nrow = 1, ncol = n_tips, dimnames = list(c(), tree$tip.label)))
  
  # Case if a star tree (no real reason to consider tree structure beyond lack of order to tip labels):
  if (n_nodes == 1) {
    
    # If tip and states counts are equal and all_states_present = TRUE then only one permutation possible:
    if (n_tips == n_states && all_states_present) state_assignments <- matrix(
      data = 1,
      nrow = 1,
      ncol = n_states,
      dimnames = list(c(), states)
    )

    # If numbers of tips and states differ and all_states_present = TRUE then perform limited permutations:
    if (n_tips != n_states && all_states_present) state_assignments <- permute_restricted_compositions(
      n = n_tips - n_states,
      m_labels = states
    ) + 1

    # If numbers of tips and states differ but all_states_present = FALSE then perform all permutations:
    if (n_tips != n_states && !all_states_present) state_assignments <- permute_restricted_compositions(
      n = n_tips,
      m_labels = states
    )

    # Convert to a tip_state_assignments matrix:
    tip_state_assignments <- do.call(
      what = rbind,
      args = apply(
        X = state_assignments,
        MARGIN = 1,
        FUN = function(state_permutation) unlist(
          x = lapply(
            X = as.list(x = 1:n_states),
            FUN = function(i) rep(
              x = names(x = state_permutation)[i],
              times = state_permutation[i]
            )
          )
        ),
        simplify = FALSE
      )
    )
    
    # Add tip labels as column names:
    colnames(x = tip_state_assignments) <- tree$tip.label
    
    # Return tip state assignments:
    return(tip_state_assignments)
  }
  
  # If not a star tree:
  if (n_nodes > 1) {
    
    # First identify any berries (i.e., any node whose descendants are all terminals):
    berries <- node_numbers[unlist(
      x = lapply(
        X = as.list(node_numbers),
        FUN = function(node) {
          all(tree$edge[tree$edge[, 1] == node, 2] <= n_tips)
        }
      )
    )]
    
    # Establish sizes of each berry:
    berry_sizes <- unlist(x = lapply(X = as.list(x = berries), FUN = function(i) sum(tree$edge[, 1] == i)))
    
    # And unique sizes:
    unique_berry_sizes <- unique(x = berry_sizes)
    
    # Generate a berry permutation for each unique berry size:
    berry_permutations <- lapply(
      X = as.list(x = unique_berry_sizes),
      FUN = function(i) {
        state_assignments <- permute_restricted_compositions(n = i, m_labels = states)
        tip_state_assignments <- do.call(
          what = rbind,
          args = apply(
            X = state_assignments,
            MARGIN = 1,
            FUN = function(state_permutation) unlist(
              x = lapply(
                X = as.list(x = 1:n_states),
                FUN = function(i) rep(
                  x = names(x = state_permutation)[i],
                  times = state_permutation[i]
                )
              )
            ),
            simplify = FALSE
          )
        )
      }
    )
    
    # Start to build a grid list to later expand for full permutation set:
    grid_list <- lapply(
      X = as.list(x = 1:length(x = berries)),
      FUN = function(i) {
        i_tip_labels <- tree$tip.label[tree$edge[tree$edge[, 1] == berries[i], 2]]
        i_permutations <- berry_permutations[[which(x = unique_berry_sizes == berry_sizes[i])]]
        i_permutations <- unlist(
          x = apply(
            X = i_permutations,
            MARGIN = 1,
            FUN = paste,
            collapse = "&",
            simplify = FALSE
          )
        )
        names(x = i_permutations)[1:length(x = i_permutations)] <- paste(i_tip_labels, collapse = "&")
        i_permutations
      }
    )
    
    # Find any non berry tips:
    nonberry_tip_labels <- setdiff(
     x = tree$tip.label,
     y = unlist(
       x = lapply(
         X = grid_list,
         FUN = function(i) strsplit(
           x = names(x = i)[1],
           split = "&"
          )[[1]]
        )
      )
    )
    
    # If there are any non berry tips:
    if (length(x = nonberry_tip_labels) > 0) {
      
      # Add non berry tips to grid_list:
      grid_list <- c(
        grid_list,
        lapply(
          X = nonberry_tip_labels,
          FUN = function(tip) {
            tip_states <- states
            names(x = tip_states)[1:n_states] <- tip
            tip_states
          }
        )
      )
    }
    
    # Make into tip state assignments matrix by expanding grid list:
    tip_state_assignments <- do.call(
      what = rbind,
      args = apply(
        X = expand.grid(grid_list),
        MARGIN = 1,
        FUN = function(i) {
          
          # Make vector of tip states:
          tip_states <- unname(
            obj = unlist(
              x = strsplit(
                x = i,
                split = "&"
              )
            )
          )
          
          # If all_states_present is TRUE then prune out any where they are not all present:
          if (all_states_present && length(x = setdiff(x = states, y = tip_states)) > 0) tip_states <- c()
          
          # Output tip states:
          tip_states
        },
        simplify = FALSE
      )
    )
    
    # Add tip labels as column names:
    colnames(x = tip_state_assignments) <- unlist(x = lapply(X = grid_list, FUN = function(i) strsplit(x = names(x = i)[1], split = "&")))
    
    # Return tip state assignments:
    return(tip_state_assignments[, tree$tip.label])
  }
}
