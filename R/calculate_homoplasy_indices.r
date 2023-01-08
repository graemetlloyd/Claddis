#' Calculates common homoplasy metrics for a set of phylogenetic tree(s)
#'
#' @description
#'
#' Given a tree, or set of trees, and a cladistic matrix calculates common homoplasy metrics for each tree.
#'
#' @param trees A tree (\code{phylo} object) or set of trees (\code{multiPhylo} object).
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}. These should be discrete with rownames (taxon labels) matching the tip labels of \code{trees}.
#' @param inapplicables_as_missing Logical that decides whether or not to treat inapplicables as missing (TRUE) or not (FALSE, the default and recommended option).
#' @param polymorphism_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#' @param uncertainty_behaviour One of either "missing", "uncertainty", "polymorphism", or "random". See details.
#' @param polymorphism_geometry Argument passed to \link{make_costmatrix}.
#' @param polymorphism_distance Argument passed to \link{make_costmatrix}.
#' @param state_ages Argument passed to \link{make_costmatrix}.
#' @param dollo_penalty Argument passed to \link{make_costmatrix}.
#'
#' @details
#'
#' Several metrics have been developed to assess the degree of homoplasy - broadly the proportion of parallelisms plus reversals to first acquisitions of derived states - implied by a tree. This function calculates the common forms of these, excluding those that require computationally intensive permutations to which R is not well suited.
#'
#' \bold{Basic terms}
#'
#' Common metrics are built around the minimum number of steps - or really the minimum cost as steps refers to only the discrete case - a character can have on the tree (m), the maximum number (g), and those sampled from the empirical topolog(ies) (s). When summed across all characters these are given as the capitals M, G, and S, respectively (Kitching et al. 1998). The function also tabulates h/H which is given as s - m/S - M (for a character or entire matrix, respectively), and was proposed in Hoyal Cuthill (2015) and is identical to Fisher's (1994) "parsimony debt".
#'
#' \emph{Minimum steps (cost)}
#'
#' The minimum number of steps (cost) is given by either the length of the minimum spanning tree for a symmetric character, or the shortest arboresence for an asymmetric character, on the costmatrix of \emph{sampled} tip states. More specifically, although many taxa may be involved, we can imagine for the minimum possible number of steps (cost) any taxa that share a character state must form a subtree and that any arrangement of that subtree would result in zero inferred steps (cost). Thus each such subtree may be collapsed to a single node (vertex) representing that character state and hence the problem collapses to a minimum spanning tree on a costmatrix (i.e., a graph representation) problem. Fortunately, such graph theory problems have wide applicability and so a range of solutions are already available in the literature. For example, for a symmetric costmatrix - the case of an undirected graph - a minimum spanning tree can be found with Borůvka's (Borůvka 1926a,b; Choquet 1938; Florek et al. 1951; Nešetřil et al. 2001), DJP's (Jarnik 1930; Prim 1957; Dijkstra 1959), or Kruskal's (Kruskal 1956) algorithm. For an asymmetric costmatrix the problem is the same as finding a minimum spanning tree on a directed graph (digraph) and because such a tree must have a tree-like shape (including a root and subsequent direction of branching) this is more technically called an arboresence. Again, though, the wide applicablity of solutions means multiple algorithms for finding such a tree are already available (e.g., Chu and Liu 1965; Edmonds 1967; Bock 1971).
#'
#' \emph{Maximum steps (cost)}
#'
#' Mathematical proofs for maximum tree lengths are considerably more complex to solve (Hoyal Cuthill et al. 2010), and to the best of my knowledge there is no general costmatrix solution currently available. Instead, this function simply uses the approach common to other software and uses the minimum steps required on the star phylogeny as the maximum steps value. Note: this doesn't mean this is not the true maximum value, simply that this is not known for certain beyond the simpler cases, such as binary and unordered characters (see Hoyal Cuthill et al. 2010 and references therein).
#'





#' \bold{The Consistency Index}
#'
#' CI = M/S character version vs ensemble version (Kitching et al. 1998) [WITHOUT UNINFORMATIVE CHARACTERS? - RESCALED?]
#'
#' \bold{The Retention Index}
#'
#' Farris (1989) character version vs ensemble version (Kitching et al. 1998)
#'


# CAN GET LONGER TREE LENGTHS HERE THAN FROM SOFTWARE SUCH AS TNT DUE TO POLYTOMIES BEING TREATED AS HARD HERE
# OTHER OPTIONS AND THEIR EFFECTS ON HOMOPLASY?

#' \bold{Permutation methods}
#'
#' Cannot do HER as this requires inferring additional trees which Claddis is not setup to do. [REFS of alternates such as Hoyal Cuthill 2015]
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com}
#'
#' @references
#'
#' Bock, F., 1971. An algorithm to construct a minimum directed spanning tree in a directed Network. In, \emph{Developments in Operations Research}. Gordon and Breach, New York, pp. 29-44.
#'
#' Borůvka, O., 1926a. O jistém problému minimálním [About a certain minimal problem]. \emph{Práce Moravské Přírodovědecké Společnosti V Brně III}, \bold{3}, 37-58.
#'
#' Borůvka, O., 1926b. Příspěvek k řešení otázky ekonomické stavby elektrovodních sítí [Contribution to the solution of a problem of economical construction of electrical networks]. \emph{Elektronický Obzor}, \bold{15}, 153-154.
#'
#' Choquet, G., 1938. Étude de certains réseaux de routes. \emph{Comptes Rendus de l'Académie des Sciences}, \bold{206}, 310-313.
#'
#' Chu, Y. J. and Liu, T. H., 1965. On the shortest arborescence of a directed graph. \emph{Science Sinica}, \bold{14}, 1396-1400.
#'
#' Dijkstra, E. W., 1959. A note on two problems in connexion with graphs. \emph{Numerische Mathematik}, \bold{1}, 269–271.
#'
#' Edmonds, J., 1967. Optimum branchings. \emph{Journal of Research of the National Bureau of Standards Section B}, \bold{71B}, 233-240.
#'
#' Farris, J. S., 1989. The retention index and the rescaled consistency index.\emph{Cladistics}, \bold{5}, 417-419.
#'
#' Fisher, D. C., 1994. Stratocladistics: morphological and temporal patterns and their relation to phylogenetic process. In L. Grande and O. Rieppel (eds.), \emph{Interpreting the Hierarchy of Nature}. Academic Press, San Diego. pp133–171.
#'
#' Florek, K., Łukaszewicz, J., Perkal, J., Steinhaus, H. and Zubrzycki, S., 1951. Sur la liaison et la division des points d'un ensemble fini. \emph{Colloquium Mathematicae}, \bold{2}, 282-285.
#'
#' Hoyal Cuthill, J., 2015. The size of the character state space affects the occurrence and detection of homoplasy: modelling the probability of incompatibility for unordered phylogenetic characters. \emph{Journal ofTheoretical Biology}, \bold{366}, 24-32.
#'
#' Hoyal Cuthill, J. F., Braddy, S. J. and Donoghue, P. C. J., 2010. A formula for maximum possible steps in multistate characters: isolating matrix parameter effects on measures of evolutionary convergence. \emph{Cladistics}, \bold{26}, 98-102.
#'
#' Jarník, V., 1930, O jistém problému minimálním [About a certain minimal problem]. \emph{Práce Moravské Přírodovědecké Společnosti}, \bold{6}, 57–63.
#'
#' Kitching, I. J., Forey, P. L., Humphries, C. J. and Williams, D. M., 1998. \emph{Cladistics: The Theory and Practice of Parsimony Analysis (Second Edition)}. Oxford University Press, Oxford. 229pp.
#'
#' Kruskal, J. B., 1956. On the shortest spanning subtree of a graph and the traveling salesman problem. \emph{Proceedings of the American Mathematical Society}, \bold{7}, 48-50.
#'
#' Nešetřil, J., Milková, E. and Nešetřilová, H., 2001. Otakar Borůvka on minimum spanning tree problem: translation of both the 1926 papers, comments, history. \emph{Discrete Mathematics}, \bold{233}, 3-36.
#'
#' Prim, R. C., 1957. Shortest connection networks and some generalizations, \emph{Bell System Technical Journal}, \bold{36}, 1389-1401.
#'
#' @return
#'
#' Text.
#'
#' @seealso
#'
#' \link{calculate_tree_length}
#'
#' @examples
#'
#' # TO DO
#'
#' @export calculate_homoplasy_indices
calculate_homoplasy_indices <- function(
  trees,
  cladistic_matrix,
  inapplicables_as_missing = FALSE,
  polymorphism_behaviour,
  uncertainty_behaviour,
  polymorphism_geometry,
  polymorphism_distance,
  state_ages,
  dollo_penalty
) {
  
  ### ADD IMPLIED WEIGHTS TO OUTPUT k=0.5; plot(x =0:100, y = (k + 1) / (1:101 + k + 1 - 1), type = "l", xlab = "Homoplasy (extra steps)", ylab = "Weight", ylim = c(0, 1))
  ### IF CALLING find_stategraph_minimum_span REMEMBER TO NOT JUST TAKE THE ARC WEIGHTS BUT TO CONSIDEr CHARACTER WEIGHT TOO!
  
  # HOMOPLASY FUNCTION:
  # - Option to define outgroup (root state?) as will affect min cost value
  # - Otherwise need to do path both ways, e.g., 0->1->2 and 2->1->0 in case of asymmetric characters
  # - Keep check for variance as no changes at all if invariant.
  # - How to handle uncertainties?
  # - Continuous characters?
  # - Polymorphism option of splitting into multiple taxa (increase t) as a berry in tree - faster than expanding costmatrix!
  
  # CALCULATE G_MIN WITH PERMUTATIONS LIKE WE DO WITH G_MAX (BEST AND WORST CHARACTER STATE FREQUENCIES ON STAR TREE?)
  
  # CHECK FOR STRATIGRAPHY AS SHOULD NOT INCLUDE THESE IN HOMOPLASY MEASURE!
  
  # ANY BINARY:
  # MAKE GMAX ETC. WORK LIKE MINIMUM COST EQUATION (LOTS OF CONDITIONALS TO DROP OUT EARLY IF SIMPLE)
  
  # NEED PRUNE FUNCTION
  ### FUNCTION TO PRUNE COSTMATRIX SHOULD ALLOW PRUNING OF POLYMORPHISMS AND UNCERTAINTIES?
  # SAFELY PRUNE COSTMATRIX TO SAMPLED STATES (INCLUDING POLYMORPHISMS AND UNCERTAINTIES

#NEED TO FORMAT ALL CLADISTIC MATRICES WITH COSTMATRICES? I.E., ON IMPORT?
#MIGHT BE TRICKY AND POINTLESS?
  
  # TO INCLUDE:
  # - Add sanity check that minimum is less than maximum and sampled fits inbetween and also that nothing is greater than g_maximum
  # - RI? (Need to know max value which is just fit on star tree)
  # - N reversals (Will depend on root and "direction") Only works for an MPR - an SCM feature not a homoplasy one?)
  # - N parallelisms (Will depend on root and "direction") Only works for an MPR - an SCM feature not a homoplasy one?)
  # - Record character types too (Dollo important as really want to count single steps but will affect min and max values)
  
  # EXAMPLE DATA
  inapplicables_as_missing = TRUE
  polymorphism_behaviour = "uncertainty"
  uncertainty_behaviour = "uncertainty"
  polymorphism_geometry = "hypersphere"
  polymorphism_distance = "great_circle"
  state_ages = c()
  dollo_penalty = 100
  cladistic_matrix <- Claddis::read_nexus_matrix("~/Documents/Homepage/www.graemetlloyd.com/nexus/Gheerbrant_etal_2014a.nex")
  trees <- ape::read.tree("~/Documents/Homepage/www.graemetlloyd.com/mpts/Gheerbrant_etal_2014a.tre")
  trees <- lapply(trees, function(i) {
    i$tip.label[i$tip.label == "Perissodactyl"] <- "Perissodactyla"
    i$tip.label[i$tip.label == "Phosphatheriu"] <- "Phosphatherium"
    i
  })
  class(trees) <- "multiPhylo"
  cladistic_matrix$topper$costmatrices <- lapply(cladistic_matrix$topper$costmatrices, function(x) {if (any(x$costmatrix == 10)) {x$costmatrix[x$costmatrix == 10] <- Inf}; x})
  
  
  # Modify:
  #tips_to_remove <- trees[[1]]$tip.label[1:13]
  #trees <- lapply(trees, ape::drop.tip, tip = tips_to_remove)
  #class(trees) <- "multiPhylo"
  #cladistic_matrix <- prune_cladistic_matrix(cladistic_matrix, characters2prune = 11:175, taxa2prune = tips_to_remove)
  #cladistic_matrix$matrix_1$matrix[, 1:10] <- sample(c("0", "1"), size = 120, replace = TRUE)
  #cladistic_matrix$matrix_1$ordering <- paste("step_", LETTERS[1:10], sep = "")
  #cladistic_matrix$matrix_1$maximum_values <- rep(1, 10)
  #for(i in 1:10) {
  #  cladistic_matrix$topper$costmatrices[[i]]$size <- 2
  #  cladistic_matrix$topper$costmatrices[[i]]$type <- "custom"
  #  if(i >= 5) j <- i + 1
  #  if(i < 5) j <- i
  #  cladistic_matrix$topper$costmatrices[[i]]$costmatrix <- matrix(c(0, 5, j, 0), nrow = 2, dimnames = list(c(0, 1), c(0, 1)))
  #  cladistic_matrix$topper$costmatrices[[i]]$symmetry <- "Asymmetric"
  #}
  
  
  
  #prunechars <- grep("step", cladistic_matrix$matrix_1$ordering)
  #cladistic_matrix  <- prune_cladistic_matrix(cladistic_matrix=cladistic_matrix, characters2prune=prunechars)
  
  
  # NEED TO EDIT PRUNE FUNCTION TO REMOVE STEP MATRICES IF THESE ARE REMOVED BY CHARACTERS TO PRUNE OR MODIFIED BY REMOVING TAXA THAT ALSO ELIMINATE TIP STATES! (ALSO DOES IT EDIT CHARACTER FOR OUTPUT?
  
  
  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")
  
  # Checks trees has class phylo or multiPhylo and stop and warn user if not:
  if (!any(c(inherits(x = trees, what = "phylo"), inherits(x = trees, what = "multiPhylo")))) stop("trees must be an object of class \"phylo\" or \"multiPhylo\".")

  # If only a single tree:
  if (class(trees) == "phylo") {
    
    # Make trees a list:
    trees <- list(trees)
    
    # Set class as multiPhylo:
    class(trees) <- "multiPhylo"
  }
  
  # Perform first pass of data to get character length information:
  sampled_costs <- calculate_tree_length(
    trees = trees,
    cladistic_matrix = cladistic_matrix,
    inapplicables_as_missing = inapplicables_as_missing,
    polymorphism_behaviour = polymorphism_behaviour,
    uncertainty_behaviour = uncertainty_behaviour,
    polymorphism_geometry = polymorphism_geometry,
    polymorphism_distance = polymorphism_distance,
    state_ages = state_ages,
    dollo_penalty = dollo_penalty
  )$character_lengths
 
  # Subfunction to convert input matrix to a costmatrix format that enables min/max cost calculations:
  build_costmatrix_matrix <- function(cladistic_matrix) {
  
    # Build combined matrix:
    combined_matrix <- list(
    
      # Add matrix part:
      matrix = do.call(
        what = cbind,
        args = lapply(
          X = cladistic_matrix[names(x = cladistic_matrix)[-1]],
          FUN = function(i) i$matrix
        )
      ),
      
      # Add ordering part:
      ordering = unname(
        obj = do.call(
          what = c,
          args = lapply(
            X = cladistic_matrix[names(x = cladistic_matrix)[-1]],
            FUN = function(i) i$ordering
          )
        )
      ),
      
      # Add character weights part:
      character_weights = unname(
        obj = do.call(
          what = c,
          args = lapply(
            X = cladistic_matrix[names(x = cladistic_matrix)[-1]],
            FUN = function(i) i$character_weights
          )
        )
      ),
      
      # Add minimum values part:
      minimum_values = unname(
        obj = do.call(
          what = c,
          args = lapply(
            X = cladistic_matrix[names(x = cladistic_matrix)[-1]],
            FUN = function(i) i$minimum_values
          )
        )
      ),
      
      # Add maximum values part:
      maximum_values = unname(
        obj = do.call(
          what = c,
          args = lapply(
            X = cladistic_matrix[names(x = cladistic_matrix)[-1]],
            FUN = function(i) i$maximum_values
          )
        )
      )
    )
    
    # Get number of characters:
    n_characters <- ncol(x = combined_matrix$matrix)
    
    # Add costmatrices part (the real meat of the function):
    combined_matrix[["costmatrices"]] <- lapply(
      X = as.list(x = 1:n_characters),
      FUN = function(i) {
        
        # POSSIBLE ORDERING TYPES:
        # - continuous !!!
        # - dollo !!!!
        # - stratigraphic !!!
        # - NOT custom as these must be coded as named cost matrices!
        
        # NEED TO CHECK FOR ABOVE AND ACT ACCORDINGLY

        # If character is ordered:
        if (combined_matrix$ordering[i] == "ordered") {
          costmatrix <- make_costmatrix(
            min_state = combined_matrix$minimum_values[i],
            max_state = combined_matrix$maximum_values[i],
            character_type = "ordered",
            include_polymorphisms = FALSE,
            polymorphism_geometry = polymorphism_geometry,
            polymorphism_distance = polymorphism_distance,
            state_ages = c(),
            dollo_penalty = dollo_penalty
          )
        }
        
        # If character is unordered:
        if (combined_matrix$ordering[i] == "unordered") {
          costmatrix <- make_costmatrix(
            min_state = combined_matrix$minimum_values[i],
            max_state = combined_matrix$maximum_values[i],
            character_type = "unordered",
            include_polymorphisms = FALSE,
            polymorphism_geometry = polymorphism_geometry,
            polymorphism_distance = polymorphism_distance,
            state_ages = c(),
            dollo_penalty = dollo_penalty
          )
        }
        
        # If character is a custom step (cost) matrix:
        if (length(x = grep(pattern = "step", x = combined_matrix$ordering[i])) == 1) {
          costmatrix <- cladistic_matrix$topper$costmatrices[[combined_matrix$ordering[i]]]
        }
        
        # REMOVE UNSAMPLED STATES FROM COSTMATRIX
        # NEED SUBFUNCTION FOR PRUNING A COSTMATRIX! - MUST CHECK PRUNED MATRIX MAKES SENSE!
        
        # Output costmatrix:
        costmatrix
      }
    )
    
    # ADD SAMPLED STATES (HAVE TO THINK WHAT TO DO WITH POLYMORPHISMS AND UNCERTANTIES!)
    
    # Output combined
    combined_matrix
  }
  
  # Convert cladistic_matrix to costmatrix_matrix:
  costmatrix_matrix <- build_costmatrix_matrix(cladistic_matrix = cladistic_matrix)
  
  # Get total number of taxa:
  n_taxa <- nrow(x = costmatrix_matrix$matrix)
  
  # Now calculate (m) minimum costs for each character:
  minimum_costs <- unlist(
    x = lapply(
      X = costmatrix_matrix$costmatrices,
      FUN = function(i) find_costmatrix_minimum_span(i)
    )
  )
  
  # Reduce each cost matrix to a single character value to check for unique types:
  costmatrix_encodings <- unlist(
    x = lapply(
      X = costmatrix_matrix$costmatrices,
      function(i) paste(
        unlist(x = i),
        collapse = "&"
      )
    )
  )
  
  # Find unique encodings:
  unique_costmatrix_encodings <- unique(x = costmatrix_encodings)
  
  # Find first character number for each unique encoding:
  unique_character_numbers <- unlist(
    x = lapply(
      X = as.list(x = unique_costmatrix_encodings),
      FUN = function(i) which(x = costmatrix_encodings == i)[1]
    )
  )
  
  # Calculate g_max costs for each unique character:
  unique_g_maximum_costs <- unlist(
    x = lapply(
      X = as.list(x = unique_character_numbers),
      FUN = function(i) {
      
        # Define tip states for character i:
        tip_states <- as.character(x = costmatrix_matrix$minimum_values[i]:costmatrix_matrix$maximum_values[i])
        
        # If only one state then build single permutation of n_taxa:
        if (length(x = tip_states) == 1) permuted_tipstate_frequencies <- matrix(data = n_taxa, nrow = 1, dimnames = list(c(), c(tip_states)))
      
        # If there is at least two states then permute all ways to assign n_taxa to them:
        if (length(x = tip_states) > 1) permuted_tipstate_frequencies <- permute_restricted_compositions(
          n = n_taxa,
          m_labels = tip_states,
          allow_zero = FALSE
        )
        
        # Permute costs on star tree:
        permuted_costs <- do.call(
          what = cbind,
          args = lapply(
            X = as.list(x = 1:costmatrix_matrix$costmatrices[[i]]$size),
            FUN = function(j) {
          
              # For each row of the cost matrix multiply by character frequency:
              cost_by_frequency <- do.call(
                what = cbind,
                args = lapply(
                  X = as.list(x = 1:costmatrix_matrix$costmatrices[[i]]$size),
                  function(k) costmatrix_matrix$costmatrices[[i]]$costmatrix[j, k] * permuted_tipstate_frequencies[, k]
                )
              )
     
              # Return total cost:
              apply(
                X = cost_by_frequency,
                MARGIN = 1,
                FUN = sum
              )
            }
          )
        )
      
        # Take minimum cost (i.e., preferred ancestral state under maximum parsimony):
        minimum_permuted_costs <- apply(
          X = permuted_costs,
          MARGIN = 1,
          FUN = min
        )
     
        # Return maximum cost of these (g max):
        max(x = minimum_permuted_costs)
      }
    )
  )
  
  # Now return full g max vector by filling out duplicated values:
  g_maximum_costs <- unlist(
    x = lapply(
      X = as.list(x = costmatrix_encodings),
      FUN = function(i) unique_g_maximum_costs[unique_costmatrix_encodings == i]
    )
  )
  
  
  
  
  
  
  
  # TEST ABOVE WITH CUSTOM ASYMMETRIC BINARY COST MATRICES
  # Jen's equation for g_max on custom asymmetric binary: gmax = floor((c0,1 * c1,0 * t) / (c1,0 + c0,1))
  for(i in 1:10) {
    print(floor((cladistic_matrix$topper$costmatrices[[i]]$costmatrix["0", "1"] * cladistic_matrix$topper$costmatrices[[i]]$costmatrix["1", "0"] * 12) / (cladistic_matrix$topper$costmatrices[[i]]$costmatrix["0", "1"] + cladistic_matrix$topper$costmatrices[[i]]$costmatrix["1", "0"])))
    #print(cladistic_matrix$topper$costmatrices[[i]]$costmatrix)
  }
  g_maximum_costs

  
  plot(g_maximum_costs, type = "n")
  for(i in 1:length(g_maximum_costs)) {
    lines(x = c(i - 0.35, i + 0.35), y = c(g_maximum_costs[i], g_maximum_costs[i]))
    lines(x = c(i - 0.35, i + 0.35), y = c(minimum_costs[i], minimum_costs[i]))
    lines(x = c(i, i), y = c(minimum_costs[i], g_maximum_costs[i]))
    points(x = i, y = sampled_costs[i, 1], pch = 18, col = "red", cex = 0.7)
  }
  
  # HAVE TO MULTIPLY THROUGH BY CHARACTER WEIGHTS AT SOME POINT!
  # PROBABLY BEST TO DO THIS AT END AND MAYBE EVEN ADD AN "UNWEIGHTED" OPTION TO OUTPUT
  
  
  
  
  
  
  character_cis <- apply(
    X = character_lengths,
    MARGIN = 2,
    FUN = function(s) min_lengths / s
  )
  
  
  
  
  tree_cis <- apply(
    X = character_lengths,
    MARGIN = 2,
    FUN = function(s) sum(x = min_lengths) / sum(x = s)
  )
  
  ### NOT SURE IF NaN REPLACEMENT IN BELOW SHOULD BE 1, 0 OR SOMETHING ELSE? (NA?) RANDOM SLIDE ON INTERNET SUGGESTS 0 BUT SHOULD FIND REF AND EXPLAIN IN MANUAL! (NB: IS TRIGGERED BY PARSIMONY UNINFORMATIVE CHARACTERS, I.E., AUTAPOMORPHIES AND INVARIANTS)
  
  character_ris <- apply(
    X = character_lengths,
    MARGIN = 2,
    FUN = function(s) as.numeric(x = gsub(pattern = NaN, replacement = 0, x = (max_lengths - s) / (max_lengths - min_lengths)))
  )
  
  tree_ris <- apply(
    X = character_lengths,
    MARGIN = 2,
    FUN = function(s) as.numeric(x = gsub(pattern = NaN, replacement = 0, x = (sum(x = max_lengths) - sum(x = s)) / (sum(x =max_lengths) - sum(x = min_lengths))))
  )
  
  # Calculate H (= S - M) of Hoyal Cuthill (2015) - Fisher's "parsimony debt":
  character_parsimony_debt <- apply(
    X = character_lengths,
    MARGIN = 2,
    FUN = function(s) s - min_lengths
  )
  tree_parsimony_debt <- apply(
    X = character_lengths,
    MARGIN = 2,
    FUN = function(s) sum(x = s) - sum(x = min_lengths)
  )
  
  # "Hmax = Gmax - M, where Gmax is the sum of the individual gmax values for the characters: gmax = t–[t=n_e]" (Hoyal Cuthill 2015)
  # Add other metrics from Hoyal CUthill 2015/Hoyal Cuthill et al. 2010
  # "A character will be parsimony uninformative if it allocates all the taxa the same state (known as an invariant or constant   character), if it allocates each taxon a different state (all states are unique, or autapomorphies) (Steel and Penny, 2005) or if all states except one are unique (Naylor and Kraus, 1995)" (Hoyal Cuthill 2015) did not know this and somewehere should be labelling autapomorphic/parsimony unninformative characters (maybe add some kind of summariser to the cladisticMatrix class? (is summary a primitive function for this like print and plot?)
  
  # Need to add g, s, m output
  
  # Return compiled output:
  list(
    tree_values = rbind(
      "S (sample cost)" = tree_lengths,
      "M (minimum cost)" = rep(x = sum(x = min_lengths), length.out = length(x = tree_lengths)),
      "G (maximum cost)" = rep(x = sum(x = max_lengths), length.out = length(x = tree_lengths)),
      "H (S - M; parsimony debt)" = tree_parsimony_debt,
      "Ensemble Consistency Index (CI)" = tree_cis,
      "Ensemble Retention Index (RI)" = tree_ris,
      "Rescaled Consistency Index (RCI)" = tree_cis * tree_ris,
      "Homoplasy Index (HI)" = 1 - tree_cis
    ),
    character_lengths = character_lengths,
    character_parsimony_debt = character_parsimony_debt,
    character_cis = character_cis,
    character_ris = character_ris,
    characters_rcis = character_cis * character_ris,
    character_his = 1 - character_cis
  )

}






prune_costmatrix <- function(costmatrix, sampled_states, message = TRUE) {
  
  # DATA CHECKS
  
  # POLYMORPHISMS CAN NEVER BE BRIDGE VERTICES (I AM PRETTY SURE)
  # SAME FOR UNCERTAINTIES AS THEY CAN NEVER BE A FROM STATE
  
  # Isolate states in costmatrix:
  costmatrix_states <- rownames(x = costmatrix$costmatrix)
  
  # Find any states absent from the costmatrix:
  states_not_in_costmatrix <- setdiff(x = sampled_states, y = costmatrix_states)
  
  # If any such states found stop and warn user:
  if (length(x = states_not_in_costmatrix) > 0) stop("Some sampled state(s) are not present in the costmatrix. Either prune sampled_states or add these to the costmatrix.")
  
  # Isolate any states to prune:
  states_to_prune <- setdiff(x = costmatrix_states, y = sampled_states)
  
  # If nothing to prune:
  if (length(x = states_to_prune) == 0) {
    
    # Message user that nothing needs to be pruned:
    if (message) print("No states to prune. Returning original costmatrix.")
    
    # Return unchanged costmatrix:
    return(object = costmatrix)
  }
  
  # If reached this point can safely state costmatrix is pruned:
  costmatrix$pruned <- TRUE

#' This is the type of graph where a state (0) - if unsampled - could not be pruned as it would lead to incorrect inferences. I.e., 0 is a cut vertex.
#'  1     2
#'   \   /
#'    \ /
#'     0
#'     |
#'     |
#'     3
#'
#'
#' Costmatrix representation:
#'
#' \preformatted{    -----------------
#'     | 0 | 1 | 2 | 3 |
#' ---------------------
#' | 0 | 0 | 1 | 1 | 1 |
#' ---------------------
#' | 1 | 1 | 0 | 2 | 2 |
#' ---------------------
#' | 2 | 1 | 2 | 0 | 2 |
#' ---------------------
#' | 3 | 1 | 2 | 2 | 0 |
#' ---------------------}



  
  # SEPARATE SINGLE STATES (POSSIBLE BRIDGE VERTICES) AND POLYMORPHISMS/UNCERTAINTIES (NON BRIDGE VERTICES)
  
  # If costmatrix is of type ordered:
  if (costmatrix$type == "ordered") {}
  
  # If costmatrix is of type unordered:
  if (costmatrix$type == "unordered") {}
  
  # If costmatrix is of type dollo:
  if (costmatrix$type == "dollo") {}
  
  # If costmatrix is of type irreversible:
  if (costmatrix$type == "irreversible") {}
  
  # If costmatrix is of type stratigraphy:
  if (costmatrix$type == "stratigraphy") {}
  
  # If costmatrix is of type custom:
  if (costmatrix$type == "custom") {}
  
  # UPDATE IF NECESSARY
  #costmatrix$single_states
  #costmatrix$costmatrix
  #costmatrix$symmetry
  #costmatrix  $includes_polymorphisms
  #costmatrix$includes_uncertainties
  #costmatrix$base_age <- ?

}


#costmatrix <- make_costmatrix(
#  min_state = 0,
#  max_state = 3,
#  character_type = "unordered"
#)
#sampled_states <- c("0", "1", "3")
#prune_costmatrix(
#  costmatrix = costmatrix,
#  sampled_states = sampled_states
#)



