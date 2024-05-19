#' Determine maximum parsimony ancestral state reconstruction(s)
#'
#' @description
#'
#' Given a tree, or set of trees, and a cladistic matrix returns all most parsimonious reconstruction(s) for every character.
#'
#' @param trees A tree (phylo object) or set of trees (multiPhylo object).
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}. These should be discrete with labels matching the tip labels of \code{tree}.
#' @param estimate_all_nodes Logical that allows the user to make estimates for all ancestral values. The default (\code{FALSE}) will only make estimates for nodes that link coded terminals (recommended).
#' @param estimate_tip_values Logical that allows the user to make estimates for tip values. The default (\code{FALSE}) will only makes estimates for internal nodes (recommended).
#' @param inapplicables_as_missing Argument passed to \link{calculate_tree_length}.
#' @param polymorphism_behaviour Argument passed to \link{calculate_tree_length}.
#' @param uncertainty_behaviour Argument passed to \link{calculate_tree_length}.
#' @param polymorphism_geometry Argument passed to \link{make_costmatrix}.
#' @param polymorphism_distance Argument passed to \link{make_costmatrix}.
#' @param state_ages Argument passed to \link{make_costmatrix}.
#' @param dollo_penalty Argument passed to \link{make_costmatrix}.
#'
#' @details
#'
#' Although most phylogenetic inference considers character evolution on a tree, they rarely explicitly assign states to internal nodes in the process, and hence this is normally done \emph{a posteriori} using some form of ancestral state estimation algorithm. This function uses the parsimony optimality criterion (i.e., maximum parsimony, or the process that involves the fewest possible transitions, or, more accurately, the lowest total transition cost).
#'
#' \bold{Algorithm}
#'
#' Like \link{calculate_tree_length}, this function is built around the generalised Swofford and Maddison (1992) algorithm, but adds the second pass as this is the one that explicitly generates \emph{all} most parsimonious reconstructions. As it is a generalised approach it can work with the broadest range of character types and (rooted) trees of any degree of resolution - treating polytomies as hard.
#'
#'

# estimate all nodes/estimate tips cite Lloyd 2018
# REFERENCE castor asr_max_parsimony and phangorn acctran as alternates

#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Lloyd, G. T., 2018. Journeys through discrete-character morphospace: synthesizing phylogeny, tempo, and disparity. \emph{Palaeontology}, \bold{61}, 637-645.
#'
#' Swofford, D. L. and Maddison, W. P., 1992. Parsimony, character-state reconstructions, and evolutionary inferences. \emph{In} R. L. Mayden (ed.) Systematics, Historical Ecology, and North American Freshwater Fishes. Stanford University Press, Stanford, p187-223.
#'
#' @return
#'
#' A list with multiple components, including:
#'
#' \item{ITEM}{ITEM EXPLANATION.}
#'
#' @seealso
#'
#' \link{calculate_tree_length}, \link{estimate_ancestral_states}
#'
#' @examples
#'
#' # Generate two trees:
#' trees <- ape::read.tree(text = c("(A,(B,(C,(D,E))));", "(A,((B,C),(D,E)));"))
#'
#' # Build a simple 5-by-10 binary character matrix:
#' cladistic_matrix <- build_cladistic_matrix(
#'   character_taxon_matrix = matrix(
#'     data = sample(
#'       x = c("0", "1"), # ADD MISSING, POLYMORPHISM, INAPPLICABLE ETC. HERE LATER TO TEST
#'       size = 50,
#'       replace = TRUE
#'     ),
#'     nrow = 5,
#'     ncol = 10,
#'     dimnames = list(
#'       LETTERS[1:5],
#'       c()
#'     )
#'   )
#' )
#'
#' # Reconstruct ancestral states (limiting output to all most parsimonious
#' # ancestral state reconstruction for every tree and character combination):
#' reconstruct_ancestral_states(
#'   trees = trees,
#'   cladistic_matrix = cladistic_matrix,
#'   estimate_all_nodes = FALSE,
#'   estimate_tip_values = FALSE,
#'   inapplicables_as_missing = FALSE,
#'   polymorphism_behaviour = "uncertainty",
#'   uncertainty_behaviour = "uncertainty",
#'   polymorphism_geometry = "simplex",
#'   polymorphism_distance = "euclidean",
#'   dollo_penalty = 999,
#'   state_ages = c(2, 1)
#' )$node_estimates
#'
#' @export reconstruct_ancestral_states
reconstruct_ancestral_states <- function(
  trees,
  cladistic_matrix,
  estimate_all_nodes = FALSE,
  estimate_tip_values = FALSE,
  inapplicables_as_missing = FALSE,
  polymorphism_behaviour = "uncertainty",
  uncertainty_behaviour = "uncertainty",
  polymorphism_geometry,
  polymorphism_distance,
  state_ages,
  dollo_penalty
) {
  
  # COLLAPSE OPTIONS (MAYBE ANOTHER FUNCTION?):
  # - ACCTRAN
  # - DELTRAN
  # - MINF (Swofford and Maddison 1987) - I THINK NOT! FREQUENCIES IN GENERAL THOUGH? MOST EQUAL DISTRIBUTION OF STATES?
  # - MINSTATE *shrug emoji*
  # - MAXSTATE *shrug emoji*
  # - WEIGHTED BY BRANCH DURATION? (E.G., ((0:1,1:4)); SHOULD FAVOUR A 0 ROOT STATE). NEED WAY TO "SCORE" THESE SUCH THAT CAN FIND MPR WITH BEST SCORE.
  # - ALLOW ROOT TO LEAN TOWARDS LOWEST STATE? (e.g., 0 > 1).
  
  ### PARSIMONY ALGORITHM SHOULDN'T EVEN CONSIDER COSTS FOR ANYTHING WHERE "TO" COST IS ALWAYS INFINITY - E.G. UNCERTAINTIES
  
  # ALLOW RECONSTRUCTIONS OF SOME COST LONGER SOMEHOW? PROBABLY A BACKBURNER/LONG TERM THING
  
  # TO ADD INTO NODE ESTIMATE VERSION OF FUNCTION:
  #estimate_all_nodes <- FALSE
  #estimate_tip_values <- FALSE # Should not happen for inapplicables when they are treated as inapplicables

  # WILL NEED TO MODIFY BELOW TO DEAL WITH UNCERTAINTIES AND POLYMORPHISMS AT TIPS
  
  # Perform Swofford and Maddison (1992) first pass of tree by calling calculate_tree_length:
  first_pass_output <- calculate_tree_length(
    trees = trees,
    cladistic_matrix = cladistic_matrix,
    inapplicables_as_missing = inapplicables_as_missing,
    polymorphism_behaviour = polymorphism_behaviour,
    uncertainty_behaviour = uncertainty_behaviour,
    polymorphism_geometry = polymorphism_geometry,
    polymorphism_distance = polymorphism_distance,
    state_ages = state_ages,
    dollo_penalty = dollo_penalty
  )
  
  # Subfunction to perform Swofford and Maddison (1992) second pass (generates all most parsimonious reconstructions):
  second_pass <- function(tree, node_values, tip_count, node_count, costmatrix) {
    
    # Build node estimates from single state results:
    node_estimates <- matrix(
      data = apply(
        X = node_values,
        MARGIN = 1,
        FUN = function(x) {
          x <- names(x = x[x == min(x = x)])
          ifelse(
            test = length(x = x) == 1,
            yes = x,
            no = NA
          )
        }
      ),
      ncol = 1,
      nrow = node_count
    )
    
    # Update tip states with their original values:
    #node_estimates[1:tip_count] <- tip_states
    
    ### ABOVE NEEDS TO BE INPUT TIP STATES NOT MODIFIED VERSION! BUT LATER FOR OUTPUT? NEED TO RECORD SOMEWHERE WHAT HAS HAPPENED THOUGH SO CAN PASS TO A CHARACTER MAP FUNCTION AND TREAT EVERYTHING CORRECTLY.
    ### ACTUALLY SHOULD PROBABLY BE
    
    # Only continue if there is any ambiguity at internal nodes:
    if (any(x = is.na(x = node_estimates[(tip_count + 1):node_count]))) {
      
      # Find all possible root states:
      possible_root_states <- colnames(x = node_values)[node_values[tip_count + 1, ] == min(x = node_values[tip_count + 1, ])]
      
      # Make new node estimates with every possible root state:
      node_estimates <- do.call(
        what = cbind,
        args = lapply(
          X = as.list(x = possible_root_states),
          FUN = function(x) {
            y <- node_estimates
            y[tip_count + 1, ] <- x
            y
          }
        )
      )
      
      # For each internal node above the root to the tips:
      for(needle in (tip_count + 2):(tip_count + node_count - tip_count)) {
        
        # Establish ancestor of current node:
        ancestor_node <- tree$edge[tree$edge[, 2] == needle, 1]
        
        # Reformat node_estimates as a list to enable lapply:
        node_estimates <- split(x = node_estimates, f = col(x = node_estimates))
        
        # Permute all possible values and reformat as matrix:
        node_estimates <- do.call(
          what = cbind,
          args = lapply(
            X = node_estimates,
            FUN = function(x) {
          
              # Get updated tree lengths for current node as per Swofford and Maddison 1992 second pass:
              updated_tree_lengths <- costmatrix$costmatrix[x[ancestor_node], ] + node_values[needle, ]
          
              # Store all possible most parsimonious state(s) for current node:
              possible_states <- names(x = updated_tree_lengths[updated_tree_lengths == min(x = updated_tree_lengths)])
          
              # Store total number of possible states:
              n_states <- length(x = possible_states)
          
              # Create new node estimates to store possible state(s):
              new_estimates <- matrix(data = rep(x = x, times = n_states), ncol = n_states)
          
              # Add possible state(s)
              new_estimates[needle, ] <- possible_states
          
              # Return new estimates only:
              new_estimates
            }
          )
        )
      }
    }
    
    # Return node estimates:
    node_estimates
  }
  
  # Perform second pass and add to first pass output:
  first_pass_output$node_estimates <- lapply(
    X = as.list(x = 1:length(x = first_pass_output$input_trees)),
    FUN = function(tree_n) {
      tree <- first_pass_output$input_trees[[tree_n]]
      tip_count <- ape::Ntip(tree)
      node_count <- tip_count + tree$Nnode
      lapply(
        X = as.list(1:ncol(x = first_pass_output$character_matrix$matrix)),
        FUN = function(character_n) {
          costmatrix <- first_pass_output$costmatrices[[character_n]]
          node_values <- first_pass_output$node_values[[tree_n]][[character_n]]
          second_pass(
            tree = tree,
            node_values = node_values,
            tip_count = tip_count,
            node_count = node_count,
            costmatrix = costmatrix
          )
        }
      )
    }
  )
  
  ### NEED TO DEAL WITH MISSING/INAPPLICABLE AND GENERAL THIRD PASS STUFF IN HERE
  ### ALSO WHAT TO SET TIP STATES AS (MAYBE RAW INPUT BUT THINK ABOUT STOCHASTIC CHARACTER MAP IMPLICATIONS)?
  
  # Return already compiled output:
  first_pass_output
}

# WHAT IF ALL CHARACTERS ARE UNCERTAIN!!!!!! (WHAT ARE ANCESTRAL STATES?)
# MANUAL: Swofford and Maddison (1992; p210): "Both [ACCTRAN and DELTRAN] start by either assuming that the ancestral state is known or by choosing a state from the set of optimal assignments to the root node" - i.e., always an arbitrainess problem at root!


# ALSO NEEDS TO BE RECORDED AS SUCH FOR STOCHASTIC CHARACTER MAPS AND HOMOPLASY METRICS.

# ADD SWOFFORD REF TO DESCRIPTION FILE (PROLLY NO DOI!)



# Record in manual that "weighted" discrete parsimony is folly (unlike square changed parsimony) - will lead to more changes than "unweighted" parsimony, but could use to chose between MPRs (maybe?)




# OUTPUT:
# - Ambiguities?
# - Tip states used
# - How root value was chosen (arbitrary forced, unambiguous?) - ADD ROOT OPTION A LA PAUP*
# - Algorithm (parsimony, ML, etc.) - for later integration into single ASR function with other options (likelihood, rerooting etc.)
# - Character map? (need to think carefully what this format should be and hence how stuff like ASE needs to inform it.)
# - MAKE THIS A CLASS TOO!

# NEED TO COVER SWOFFORD AND MADDISON 1992 ROOT VALUE CAVEAT SOMEWHERE, MAYBE TREE LENGTH FUNCTION? BUT RELATES HERE TOO

# MORE OUTPUT OR SUMMARY OPTIONS COULD BE N CHANGES ON EACH BRANCH ACROSS MPRS, MEAN OF SAME, MIN AND MAX OF SAME. SIMILAR FOR CHANGE TYPES I.E., COSTMATRIX BUT ACTUALLY FREQUENCY OF CHANGES FOR EACH TRANSITION (THIS IS REALLY AN STOCHASTIC CHARACTER MAP OUTPUT).
# AS RATE TENDS TOWARD INFINITY THEN RECONSTRUCTION AT A NODE TENDS TOWARDS EQUAL FREQUENCY OF EACH STATE (I.E., MAXIMAL UNCERTAINTY).
# "This illustrates the fact that ACCTRAN and DELTRAN do not always choose a single one of the most parsimonious reconstructions." MacClade 4 manual, page 99).
# "Note that MacClade, unlike PAUP*, does not choose the lowest-valued state at the root to begin these processes. Thus MacClade's ACCTRAN and DELTRAN may not fully resolve ambiguity in the ancestral state reconstruction."
# INTERMEDIATES ARE GONNA MATTER IF MOVING TO CHARACTER MAPS. E.G., IF ONLY 0 AND 2 ARE SAMPLED BUT A CHARACTER IS ORDERED THEN THERE ARE TWO CHANGES ALONG THE BRANCH NOT ONE TWO-COST CHANGE.
# TEST WEIGHTING OF STRATOCLADISTICS BY USING A STRATIGRAPHIC CHARACTER AND WEIGHTING IT MULTIPLE WAYS (AS A SLIDER) AND SHOW HOW IT EFFECTS PARSIMONY VS STRAT CONGRUENCE TREE LANDSCAPES
# MONOFURCATIONS IDEA IN CASTOR IS INTERESTING AS ALLOWS POINT ESTIMATES ALONG BRANCHES (E.G., FOR STOCHASTIC CHARACTER MAP OF A CONTINUOUS CHARACTER)
