#' Discrete character rates across trees, time, and character types
#'
#' @description
#'
#' Given a tree and a cladistic-type matrix uses likelihood ratio tests to compare N-rate and 1-rate models across branches, clades, time bins, or character partitions.
#'
#' @param tree A tree (phylo object) with branch lengths that represents the relationships of the taxa in \code{CladisticMatrix}.
#' @param CladisticMatrix A character-taxon matrix in the format imported by \link{ReadMorphNexus}.
#' @param TimeBins A vector of ages (in millions of years) indicating the boundaries of a series of time bins in order from oldest to youngest.
#' @param BranchPartitionsToTest A list of branch(es) (edge numbers) to test for a 2-rate parameter model (i.e., one rate for the edge and another for the rest of the tree). If NULL (the default) then no partition test(s) will be made.
#' @param CharacterPartitionsToTest A list of character partition(s) (character numbers) to test for a 2-rate parameter model (i.e., one rate for the partition and another for the remaining characters). If NULL (the default) then no partition test(s) will be made.
#' @param CladePartitionsToTest A list of clade partition(s) (node numbers) to test for a 2-rate parameter model (i.e., one rate for the clade and another for the rest of the tree). If NULL (the default) then no partition test(s) will be made.
#' @param TimeBinPartitionsToTest A list of time bin partition(s) (numbered 1 to N) to test for a 2-rate parameter model (i.e., one rate for the time bin(s) and another for the remaining time bins). If NULL (the default) then no partition test(s) will be made.
#' @param ChangeTimes The time at which to record the character changes. One of \code{"midpoint"} (changes occur at the midpoint of the branch), \code{"spaced"} (changes equally spaced along branch), or \code{"random"} (change times drawn at random from a uniform distribution; the default and recommended option). Note: this is only meaningful if testing for time bin partitions.
#' @param Alpha The alpha value to be used for the significance tests. The default is 0.01.
#' @param MultipleComparisonCorrection One of \code{"BenjaminiHochberg"} (the Benjamini and Hochberg 1995 false discovery rate approach; default and recommended) or \code{"Bonferroni"} (the Bonferroni correction).
#' @param PolymorphismState One of \code{"missing"} (converts polymorphic values to NA; the default) or \code{"random"} (picks one of the possible polymorphic states at random).
#' @param UncertaintyState One of \code{"missing"} (converts uncertain values to NA; the default) or \code{"random"} (picks one of the possible uncertain states at random).
#' @param InapplicableState The only current option is \code{"missing"} (converts value to NA).
#' @param TimeBinApproach One of \code{"Close"} or \code{"Lloyd"} (the default).
#' @param EnsureAllWeightsAreIntegers Logical for whether (\code{TRUE}) to reweight non-integer weights until all weights are integers or to leave them as they are (\code{FALSE}; the default).
#' @param EstimateAllNodes Option passed to internal use of \link{AncStateEstMatrix}.
#' @param EstimateTipValues Option passed to internal use of \link{AncStateEstMatrix}.
#' @param InapplicablesAsMissing Option passed to internal use of \link{AncStateEstMatrix}.
#' @param PolymorphismBehaviour Option passed to internal use of \link{AncStateEstMatrix}.
#' @param UncertaintyBehaviour Option passed to internal use of \link{AncStateEstMatrix}.
#' @param Threshold Option passed to internal use of \link{AncStateEstMatrix}.
#'
#' @details
#'
#' \bold{Introduction}
#'
#' Morphological change can be captured by discrete characters and their evolution modelled as occurring along the branches of a phylogenetic tree. This function takes as primary input a character-taxon matrix of discrete characters (in the format imported by \link{ReadMorphNexus}) and a time-scaled phylogenetic tree (in the format of \pkg{paleotree} or \pkg{strap}) and begins by inferring ancestral states at the tree's internal nodes using the \link{AncStateEstMatrix} function. From here changes along individual branches can be estimated (only the minimum number of changes are inferred; see \link{GetAllStateChanges} for an alternative but unfinished approach) and hence rates can be calculated.
#'
#' A discrete character rate can be expressed as the mean number of changes per million years (users may wish to normalise this by the number of characters for interpretation) and can be calculated for a branch (edge) of the tree, a clade (a mean rate for the edges descended from a single node), a character partition (the mean rate for a subset of the characters across all edges), or, most complex (see Lloyd 2016), the mean rate across the edges (or parts of edges) present in a time bin (defined by two values denoting the beginning and end of the time bin). In an ideal scenario these rates could be compared at face value, but that would require a large number of characters and very minimal (or zero) missing data. I.e., at an extreme of missing data if only one character can be observed along a branch it will either change (the maximum possible rate of evolution) or it will not (the minimum possible rate of evolution). In such cases it would be unwise to consider either outcome as being a significant departure from the mean rate.
#'
#' Because of these complications Lloyd et al. (2012) devised tests by which the significance of an edge (or other paritioning of the data, i.e., a clade, time bin etc.) could be considered to be significantly high or low in comparison to the mean rate for the whole tree (i.e., whether a two-rate model could be considered more likely than a one-rate model). This is achieved through a likelihood ratio test:
#'
#' \deqn{LR = value of likehood function under the null (one-rate) hypothesis / maximum possible value of likehood function under the alternative (two-rate) hypotheses}
#'
#' Typically we might expect the two hypotheses to be well defined a priori. E.g., an expectation that a specific branch of the tree might have a higher or lower rate than background due to some evolutionary shift. However, Lloyd et al. (2012) instead provided an exploratory approach whereby every possible one edge value was compared with the rate for the rest of the tree (and the equivalent with clades and time bins). This was the default in Claddis up to version 0.2, but this has now been replaced (since version 0.3) with a more customisable set of options that allows different types of hypotheses (e.g., partitioning the data by character), as well as more complex hypotheses (e.g., a three-rate model), to be tested.
#'
#' \bold{The four types of rate hypothesis}
#'
#' Following Cloutier (1991), Lloyd (2016) extended the two main types of rate hypotheses to four:
#'
#' \enumerate{
#'   \item A branch rate (available here with the \code{BranchPartitionsToTest} option).
#'   \item A clade rate (available here with the \code{CladePartitionsToTest} option).
#'   \item A time bin rate (available here with the \code{TimeBinPartitionsToTest} option).
#'   \item A character partition rate (available here with the \code{CharacterPartitionsToTest} option).
#' }
#'
#' In Claddis (>=0.3) these partitions are defined as a list of lists of vectors where only the first N - 1 partitions need be defined. E.g., if comparing the first edge value to the rest of the tree then the user only needs to define the value "1" and the function will automatically add a second partition containing all other edges. This can be set with the option \code{BranchPartitionsToTest = list(list(1))}. Similarly, to do what Lloyd et al. (2012) did and repeat the test for every edge in the tree (and assuming this variable is already named "tree") you could use, \code{BranchPartitionsToTest = lapply(as.list(1:nrow(tree$edge)), as.list)}.
#'
#' Because of the flexibility of this function the user can define any set of edges. For example, they could test whether terminal branches have a different rate from internal branches with \code{BranchPartitionsToTest = list(list(match(1:Ntip(tree), tree$edge[, 2])))}. The \code{CladePartitionsToTest} is really just a special subset of this type of hypothesis, but with edges being defined as descending from a specific internal node in the tree. Once again, an exploratory approach like that of Lloyd et al. (2012) can be used with: \code{CladePartitionsToTest = lapply(as.list(Ntip(tree) + (2:Nnode(tree))), as.list)}. Note that this excludes the root node as this would define a single partition and hence would represent the null hypothesis (a single rate model for the whole tree). More generally clades must be defined by the node numbers they correspond to. In R an easy way to identify these is with: \code{plot(tree); nodelabels()}.
#'
#' Time bin partitions are defined in a similar way, but are numbered 1:N starting from the oldest time bin. So if wanting to do an exploratory test of single bin partitions (and only four time bins were specified) you could use: \code{TimeBinPartitionsToTest = lapply(as.list(1:4), as.list)}. Bins can be combined too, just as edges are above. For example, time bins 1 and 2 could form a single partition with: \code{TimeBinPartitionsToTest = list(list(1:2))}. Or if looking to test a model where each bin has its' own rate value you could use: \code{TimeBinPartitionsToTest = list(as.list(1:3))}. Note, as before we do not need to specify the fourth bin as this will be automatically done by the function, however, \code{TimeBinPartitionsToTest = list(as.list(1:4))} will also work. Some caution needs to be applied with N-rate models (where N is three or larger) as a result favouring such models does not necessarily endorse N-separate rates. I.e., it could simply be that one bin has such a large excursion that overall the N-rate model fits better than the 1-rate model, but some 2-rate models might be better still. It is up to the user to check this themselves by exploring smaller combinations of bins.
#'
#' Finally, character partitions allow the user to explore whether rates vary across different character types, e.g., skeletal characters versus soft tissue characters, or cranial characters versus postcranial characters. Here characters are simply numbered 1:N, but here single character partitions are less likely to be of interest. As an example of use lets say the first ten characters are what we are interested in as a partition (the second partition being the remaining characters), we could use: \code{CharacterPartitionsToTest = list(list(1:10))} to test for a two-rate model.
#'
#' Note that the list of lists structure is critical to defining partitions as it allows them to be of different sizes and numbers. For example, one partition of three and another of six, or one set of two partitions and another set of four partitions - structures not possible using vectors or matrices. However, it may not be intuitive to some users so it is recommended that the user refers to the examples above as a guide.
#'
#' Additionally, it should be noted that the user can test multiple types of hypotheses simultaneously with the function. For example, performing several branch tests whilst also performing clade tests. However, it is not necessary to perform all types simultaneously (as was the case in version 0.2) and unused partition types can be set to NULL, the default in each case.
#'
#' \bold{Other options}
#'
#' Since Claddis version 0.3 this function has allowed the user greater control with many more options than were offered previously and these should be considered carefully before running any tests.
#'
#' Firstly, the user can pick an option for \code{ChangeTimes} which sets the times character changes are inferred to occur. This is only relevant when the user is performing time bin partition test(s) as this requires some inference to be made about when changes occur on branches that may span multiple time bins. The current options are: \code{"midpoint"} (all changes are inferred to occur midway along the branch, effectively mimicking the approach of Ruta et al. 2006), \code{"spaced"} (all changes are inferred to occur equally spaced along the branch, with changes occurring in character number order), or \code{"random"} (changes are assigned a random time by drawing from a uniform distribution between the beginning and end of each branch). The first of these is likely to lead to unrealistically "clumped" changes and by extension implies completely correlated character change that would violate the assumptions of the Poisson distribution that underlies the significance tests here (Lloyd et al. 2012). At the other extreme, the equally spaced option will likely unrealistically smooth out time series and potentially make it harder to reject the single-rate null. For these reasons, the random option is recommended and is set as the default. However, because it is random this makes the function stochastic (the answer can vary each time it is run) and so the user should therefore run the function multiple times if using this option (i.e., by using a for loop) and aggregating the results at the end (e.g., as was done by previous authors; Lloyd et al. 2012; Close et al. 2015).
#'
#' Secondly, the \code{Alpha} value sets the significance threshold by which the likelihood ratio test's resulting p-value is compared. Following Lloyd et al. (2012) this is set lower (0.01) than the standard 0.05 value by default as those authors found rates to be highly heterogenous in their data set (fossil lungfish). However, this should not be adopted as a "standard" value without question. Note that the function also corrects for multiple comparisons (using the \code{MultipleComparisonCorrection} option) to avoid Type I errors (false positives). It does so (following Lloyd et al. 2012) using the Benjamini-Hochberg (Benjamini and Hochberg 1995) False Discovery Rate approach (see Lloyd et al. 2012 for a discussion of why), but the Bonferroni correction is also offered.
#'
#' Thirdly, polymorphisms and uncertainities create complications for assessing character changes along branches. These can occur at the tips (true polymorphisms or uncertainties in sampled taxa) and internal nodes (uncertainty over the estimated ancestral state). There are two options presented here, and applicable to both \code{PolymorphismState} and \code{UncertaintyState} (allowing these to be set separately). These are to convert such values to missing (NA) or to pick one of the possible states at random. Using missing values will increase overall uncertainty and potentially lead to Type II errors (false negatives), but represents a conservative solution. The random option is an attempt to avoid Type II errors, but can be considered unrealistic, especially if there are true polymorphisms. Additionally, the random option will again make the function stochastic meaning the user should run it multiple times amd aggregate the results. Note that if there are no polymorphisms or uncertainties in the character-taxon matrix the latter can still occur with ancestral state estimates, especially if the threshold value is set at a high value (see \link{AncStateEstMatrix} for details).
#'
#' Fourthly, inapplicable characters can additionally complicate matters as they are not quite the same as missing data. I.e., they can mean that change in a particular character is not even possible along a branch. However, there is no easy way to deal with such characters at present so the user is not presented with a true option here - currently all inapplicable states are simply converted to missing values by the function. In future, though, other options may be available here. For now it is simply noted that users should be careful in making inferences if there are inapplicable characters in their data and should perhaps consider removing them with \link{MatrixPruner} to gauge their effect.
#'
#' Fifthly, there are currenty two further options for assessing rates across time bins. As noted above a complication here is that character changes (the rate numerator) and character completeness (part of the rate denominator) are typically assessed on branches. However, branches will typically span time bin boundaries and hence many bins will contain only some portion of particular branches. The exact portion can be easily calculated for branch durations (the other part of the rate denominator) and the \code{ChangeTimes} option above is used to set the rate numerator, however, the completeness remains complex to deal with. The first attempt to deal with this was made by Close et al. (2015) who simply did weighted mean completeness by using the proportion of a branch in each bin as the weight and multiplying this by each branch's completeness (the \code{"Close"} option here). However, this may lead to unrealistic "smoothing" of the data and perhaps more importantly makes no account of which characters are known in a bin. Lloyd (2016) proposed an alternative "subtree" approach which assesses completeness by considering each character to be represented by a subtree where only branches that are complete are retained then branch durations in each bin are summed across subtrees such that the duration term automatically includes completeness (the \code{"Lloyd"} option here). Here the latter is strongly recommended, for example, because this will lead to the same global rate across the whole tree as the branch, clade or character partitions, whereas the Close approach will not.
#'
#' Sixthly, all character changes are weighted according to the weights provided in the input character-taxon matrix. In many cases these will simply all be one, although see the equalise weights option in \link{ReadMorphNexus}. However, when weights vary they can create some issues for the function. Specifically, changes are expected to be in the same (integer) units, but if weights vary then they have to be modelled accordingly. I.e., a character twice the weight of another may lead to a single change being counted as two changes. This is most problematic when the user has continuous characters which are automatically converted to gap-weighted (Thiele 1993) characters. However, this conversion creates drastically down-weighted characters and hence the user may wish to set the \code{EnsureAllWeightsAreIntegers} option to TRUE. Note that reweighting will affect the results and hence shifting the weights of characters up or down will necessarily lead to shifts in the relative Type I and II errors. This is an unexplored aspect of such approaches, but is something the user should be aware of. More broadly it is recommended that continuous (or gap-weighted) characters be avoided when using these approaches.
#'
#' Finally, the remaining options (\code{EstimateAllNodes}, \code{EstimateTipValues}, \code{InapplicablesAsMissing}, \code{PolymorphismBehaviour}, \code{UncertaintyBehaviour}, and \code{Threshold}) are all simply passed directly to \link{AncStateEstMatrix} for estimating the ancestral states and users should consult the help file for that function for further details.
#'
#' Note that currently the function cannot deal with step matrices and that the terminal versus internal option from Brusatte et al. (2014) is yet to be implemented.
#'
#' \bold{Output}
#'
#' The output for each test (i.e., the components of the \code{BranchPartitionResults}, \code{CharacterPartitionResults}, \code{CladePartitionResults} and \code{TimeBinResults} parts of the output) includes three main parts:
#'
#' \enumerate{
#'   \item Rates.
#'   \item PValue.
#'   \item CorrectedAlpha.
#' }
#'
#' For each rate test the \code{Rates} part of the output is a vector of the absolute rate (number of changes per million years) for each partition in the test. So, for example, a branch rate test for the sixth edge in a tree would be the rate for the sixth edge followed by the pooled rate for all other edges. The length of the vector is the length of the number of partitions.
#'
#' The PValue is a single value indicating the probability that the likelihood ratio (see above and Lloyd et al. 2012) is one, i.e., that the likelihoods of the one-rate and N-rate models are the same.
#'
#' The CorrectedAlpha is the alpha-value that should be used to determine the significance of the current partition test (i.e., The PValue, above). If the PValue exceeds the CorrectedAlpha then the null (single-rate) hypothesis should be accepted, if lower then the null should be rejected in favour of the N-rate hypthesis. Note that the CorrectedAlpha will not typically be the same for each partition and will also typically be different from the input \code{Alpha} value due to the \code{MultipleComparisonCorrection} option used.
#'
#' @return
#'
#' \item{InferredCharacterChanges}{Matrix of inferred character changes.}
#' \item{IntrinsicCharacterRate}{The intrinsic (global) character rate in changes per million years.}
#' \item{ContinuousCharactersConvertedToDiscrete}{Whether or not continuous characters were converted to discrete characters (important for handling the data in downstream analys(es)).}
#' \item{BranchPartitionResults}{List of branch partition results (corresponding to \code{BranchPartitionsToTest}. NULL if not requested.}
#' \item{CharacterPartitionResults}{List of character partition results (corresponding to \code{CharacterPartitionsToTest}. NULL if not requested.}
#' \item{CladePartitionResults}{List of clade partition results (corresponding to \code{CladePartitionsToTest}. NULL if not requested.}
#' \item{TimeBinResults}{List of time bin partition results (corresponding to \code{TimeBinPartitionsToTest}. NULL if not requested.}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Steve C. Wang \email{scwang@@swarthmore.edu}
#'
#' @references
#'
#' Benjamini, Y. and Hochberg, Y., 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society, Series B, 57, 289-300.
#'
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. Current Biology, 24, 2386-2392.
#'
#' Close, R. A., Friedman, M., Lloyd, G. T. and Benson, R. B. J., 2015. Evidence for a mid-Jurassic adaptive radiation in mammals. Current Biology, 25, 2137-2142.
#'
#' Cloutier, R., 1991. Patterns, trends, and rates of evolution within the Actinistia. Environmental Biology of Fishes, 32, 23–58.
#'
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. Biological Journal of the Linnean Society, 118, 131-151.
#'
#' Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). Evolution, 66, 330-348.
#'
#' Ruta, M., Wagner, P. J. and Coates, M. I., 2006. Evolutionary patterns in early tetrapods. I. Rapid initial diversification followed by decrease in rates of character change. Proceedinsg of the Royal Society of London B, 273, 2107–2111.
#'
#' Thiele, K.. 1993. The Holy Grail of the perfect character: the cladistic treatment of morphometric data. Cladistics, 9, 275-304.
#'
#' @examples
#' 
#' # Set random seed:
#' set.seed(17)
#' 
#' # Generate a random tree for the Michaux data set:
#' tree <- rtree(nrow(Michaux1989$Matrix_1$Matrix))
#' 
#' # Update taxon names to match those in the data matrix:
#' tree$tip.label <- rownames(Michaux1989$Matrix_1$Matrix)
#' 
#' # Set root time by making youngest taxon extant:
#' tree$root.time <- max(diag(vcv(tree)))
#' 
#' # Get discrete character rates:
#' x <- DiscreteCharacterRate(tree = tree, CladisticMatrix =
#'   Michaux1989, TimeBins = seq(from = tree$root.time,
#'   to = 0, length.out = 5), BranchPartitionsToTest =
#'   lapply(as.list(1:nrow(tree$edge)), as.list),
#'   CharacterPartitionsToTest = lapply(as.list(1:3),
#'   as.list), CladePartitionsToTest =
#'   lapply(as.list(Ntip(tree) + (2:Nnode(tree))),
#'   as.list), TimeBinPartitionsToTest =
#'   lapply(as.list(1:4), as.list), ChangeTimes =
#'   "random", Alpha = 0.01, PolymorphismState =
#'   "missing", UncertaintyState = "missing",
#'   InapplicableState = "missing", TimeBinApproach =
#'   "Lloyd")
#'
#' @export DiscreteCharacterRate
DiscreteCharacterRate <- function(tree, CladisticMatrix, TimeBins, BranchPartitionsToTest = NULL, CharacterPartitionsToTest = NULL, CladePartitionsToTest = NULL, TimeBinPartitionsToTest = NULL, ChangeTimes = "random", Alpha = 0.01, MultipleComparisonCorrection = "BenjaminiHochberg", PolymorphismState = "missing", UncertaintyState = "missing", InapplicableState = "missing", TimeBinApproach = "Lloyd", EnsureAllWeightsAreIntegers = FALSE, EstimateAllNodes = FALSE, EstimateTipValues = FALSE, InapplicablesAsMissing = FALSE, PolymorphismBehaviour = "equalp", UncertaintyBehaviour = "equalp", Threshold = 0.01) {
  
  # DESIDERATA (STUFF IT WOULD BE NICE TO ADD IN FUTURE):
  #
  # WRITE SEARCH VERSION FOR FINDING RATE SHIFTS? SHOULD THIS EVEN BE AN OPTION? DOES THIS REQUIRE MODIFYING LRT TO COMPARE E.G. 2-RATE DIRECTLY WITH 3-RATE MODEL? OR CAN USE OUTPUT P-VALUES AND CONVERT MODELS TO AICS? WOULD NEED TO PERMUTE ALL POSSIBLE COMBOS AND NOT SURE HOW LARGE THESE MIGHT GET.
  # MAYBE MAKE ANCESTRAL STATE UNCERTAINTY DIFFERENT FOR TIPS THAN NODES? I.E., HOW IT GETS RESOLVED CAN BE DIFFERENT (MORE OPTIONS TO FUNCTION)
  # THESE TWO ARE RELATED: 1. ADD TERMINAL VERSUS INTERNAL OPTION SOMEHOW/SOMEWHERE, 2. ALLOW OPTION TO IGNORE SOME PARTS OF THE TREE FOR PARTITION TESTS? MAKES CALCULATING THE MEAN RATE TRICKIER BUT MIGHT MAKE SENSE E.G. FOR INGROUP ONLY TESTS. EXCLUDE EDGES AFTER DOING ANCESTRAL STATES? OR SET THESE TO ALL NAS TO ENSURE OTHER THINGS WORK FINE?
  # EXTRA FUNCTION(S) TO VISUALISE RESULTS MOST LIKELY
  # CHECK FOR AUTAPOMORPHIES AND INFORM USER IF FOUND?
  # ADD CONTRIVED EXAMPLES TO SHOW HOW FUNCTION WORKS, E.G. RATE OF ONE CHANGES PER MILLION YEARS THEN DUPLICATED BLOCK WITH CHARCTER PARTITION TEST.
  # PROBABLY NEED MORE CAREFUL CHECKS FOR ZERO VALUES GENERALLY, E.G., CHARACTER WITH ALL MISSING DATA
  # ALLOW REWEIGHTING OF INAPPLICABLES ZERO AS AN OPTION FOR THEM?
  # HOW TO FORMAT OUTPUT? GET CIS FOR EACH PARTITION FOR VISUALISATION (E.G., BARPLOT OF PARTITION VALUES WITH DASHED LINE FOR MEAN AND ERROR BARS FOR CIS)? STICK WITH LIST OR COLLAPSE TO A MATRIX SOMEHOW?
  # TIME BINS WITH NOTHING IN WILL CAUSE ISSUES AS DIVIDE BY ZEROES WILL OCCUR - ADD CHECK FOR THIS.
  # GET CLADE NUMBERS BACK FOR OUTPUTTING?
  # WHAT IS SIGNIFICANTLY HIGH OR LOW IF THERE ARE THREE OR MORE PARTITIONS? THIS IS NOT EVEN IN OUTPUT YET.

  # Check for step matrices and stop and warn user if found:
  if(is.list(CladisticMatrix$Topper$StepMatrices)) stop("Function cannot currently deal with step matrices.")

  # Check tree has branch lengths:
  if(is.null(tree$edge.length)) stop("Tree does not have branch lengths (durations). Try timescaling the tree, e.g., with DatePhylo.")
  
  # Check tree has root age:
  if(is.null(tree$root.time)) stop("Tree is missing $root.time. Try setting this before continuing, e.g., tree$root.time <- 104.2.")
  
  # Check ChangeTimes is correctly formatted or stop and warn user:
  if(length(setdiff(ChangeTimes, c("midpoint", "spaced", "random"))) > 0) stop("ChangeTimes must be one of \"midpoint\", \"spaced\", or \"random\".")
  
  # Check MultipleComparisonCorrection is correctly formatted or stop and warn user:
  if(length(setdiff(MultipleComparisonCorrection, c("BenjaminiHochberg", "Bonferroni"))) > 0) stop("MultipleComparisonCorrection must be one of \"BenjaminiHochberg\" or \"Bonferroni\".")
  
  # Check PolymorphismState is correctly formatted or stop and warn user:
  if(length(setdiff(PolymorphismState, c("missing", "random"))) > 0) stop("PolymorphismState must be one of \"missing\" or \"random\".")
  
  # Check UncertaintyState is correctly formatted or stop and warn user:
  if(length(setdiff(UncertaintyState, c("missing", "random"))) > 0) stop("UncertaintyState must be one of \"missing\" or \"random\".")
  
  # Check InapplicableState is correctly formatted or stop and warn user:
  if(length(setdiff(InapplicableState, c("missing"))) > 0) stop("InapplicableState must be \"missing\".")
  
  # Check TimeBinApproach is correctly formatted or stop and warn user:
  if(length(setdiff(TimeBinApproach, c("Close", "Lloyd"))) > 0) stop("TimeBinApproach must be one of \"Close\" or \"Lloyd\".")
  
  # Check partitions are not all NULL values:
  if(is.null(BranchPartitionsToTest) && is.null(CharacterPartitionsToTest) && is.null(CladePartitionsToTest) && is.null(TimeBinPartitionsToTest)) stop("No partitions are requested. Set at least one of BranchPartitionsToTest, CharacterPartitionsToTest, CladePartitionsToTest, or TimeBinPartitionsToTest to a list of appropriate values. Type \"?DiscreteCharacterRate\" for help.")
  
  # Get internal node numbers:
  InternalNodeNumbers <- 1:ape::Nnode(tree) + ape::Ntip(tree)
  
  # Get edge numbers:
  EdgeNumbers <- 1:nrow(tree$edge)
  
  # Get character numbers:
  CharacterNumbers <- 1:sum(unlist(lapply(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Matrix"), ncol)))
  
  # Ensure time bins are in correct order:
  TimeBins <- sort(unique(TimeBins), decreasing = TRUE)
  
  # Find the Time bin midpoints:
  TimeBinMidpoints <- (TimeBins[2:length(TimeBins)] + TimeBins[1:(length(TimeBins) - 1)]) / 2
  
  # Get the numbers for each time bins:
  TimeBinNumbers <- 1:length(TimeBinMidpoints)
  
  # Subfunction to ensure partitions are formatted correctly:
  PartitionFormatter <- function(PartitionsToTest, ValidValues, PartitionName) {
    
    # Check partitions are in the form of a list of lists:
    if(!all(c(all(unlist(lapply(PartitionsToTest, is.list))), is.list(PartitionsToTest)))) stop(paste(PartitionName, " must be in the form of a list of lists.", sep = ""))
    
    # Get a vector of any non-present valid values:
    NonPresentValues <- setdiff(unique(unlist(PartitionsToTest)), ValidValues)
    
    # Check valid values have been used and if not stop and warn user:
    if(length(NonPresentValues) > 0) stop(paste(PartitionName, "Partitions to test must be defined using the valid range of values (", paste(range(ValidValues), collapse = " to "), ") only.", sep = ""))
    
    # Check partitions never overlap and stop and warn user if they do:
    Check <- lapply(PartitionsToTest, function(x) if(any(duplicated(sort(unlist(x))))) stop(paste("Each partition of ", PartitionName, " must not contain overlapping values (e.g., can not have 1:3 and 3:5 as both contain 3).", sep = "")))
    
    # Subfunction to ad the missing partition (if exists):
    AddMissingPartitions <- function(x, ValidValues) {
      
      # Define any missing values:
      MissingValues <- setdiff(ValidValues, unlist(x))
      
      # If there are missing values add them to list at end:
      if(length(MissingValues) > 0) x[[(length(x) + 1)]] <- MissingValues
      
      # Return x:
      return(x)
      
    }
    
    # Add in missing partitions (if any):
    PartitionsToTest <- lapply(PartitionsToTest, AddMissingPartitions, ValidValues = ValidValues)
    
    # Check partitions are all at least two in size or else no comparison can be made:
    if(any(unlist(lapply(PartitionsToTest, length)) == 1)) stop("Partitions must divide the available data into at least two parts.")

    # Return formatted partitions to test:
    return(PartitionsToTest)

  }

  # If performing branch partition test(s) check and reformat branch partitions:
  if(!is.null(BranchPartitionsToTest)) BranchPartitionsToTest <- PartitionFormatter(PartitionsToTest = BranchPartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "BranchPartitionsToTest")
  
  # If performing character partition test(s) check and reformat character partitions:
  if(!is.null(CharacterPartitionsToTest)) CharacterPartitionsToTest <- PartitionFormatter(PartitionsToTest = CharacterPartitionsToTest, ValidValues = CharacterNumbers, PartitionName = "CharacterPartitionsToTest")
  
  # If performing clade partition test(s)
  if(!is.null(CladePartitionsToTest)) {
    
    # Convert clade partitions to edge partitions:
    CladePartitionsToTest <- lapply(CladePartitionsToTest, lapply, GetDescendantEdges, tree = tree)
    
    # Check and reformat clade partitions:
    CladePartitionsToTest <- PartitionFormatter(PartitionsToTest = CladePartitionsToTest, ValidValues = EdgeNumbers, PartitionName = "CladePartitionsToTest")
    
  }

  # If performing time bin partition test(s) check and reformat time bin partitions:
  if(!is.null(TimeBinPartitionsToTest)) TimeBinPartitionsToTest <- PartitionFormatter(PartitionsToTest = TimeBinPartitionsToTest, ValidValues = TimeBinNumbers, PartitionName = "TimeBinPartitionsToTest")
  
  # Subfunction to calculate maximum likelihood p value:
  GetMaximumLikelihoodPValue <- function(MeanRate, SampledRates, SampledChanges, SampledCompleteness, SampledTime) {
    
    # Set maximum likelihood numerator:
    MaximumLikelihoodNumerator <- MeanRate
    
    # Set maximum likelihood denominator:
    MaximumLikelihoodDenominator <- SampledRates
    
    # Get log numerator:
    LogNumerator <- sum(log(dpois(round(SampledChanges), MaximumLikelihoodNumerator * SampledCompleteness * SampledTime)))
    
    # Get log denominator:
    LogDenominator <- sum(log(dpois(round(SampledChanges), MaximumLikelihoodDenominator * SampledCompleteness * SampledTime)))
    
    # Get test statistic:
    TestStatistic <- -2 * (LogNumerator - LogDenominator)
    
    # Calculate position of test statistic in chi-square distribution to get probability:
    PValue <- pchisq(TestStatistic, length(SampledRates) - 1, lower.tail = FALSE)
    
    # Output probability for later alpha comparison:
    return(PValue)
    
  }
 
  # Get ages for each (tip and internal) node:
  NodeAges <- GetNodeAges(tree)

  # Get branch ages (from and to):
  BranchAges <- unname(cbind(NodeAges[as.character(tree$edge[, 1])], NodeAges[as.character(tree$edge[, 2])]))

  # Build edge list from node numbers (from-to) for each branch:
  EdgeList <- lapply(apply(tree$edge, 1, list), function(x) {names(x) <- "NodeNumberFromTo"; return(x)})
  
  # Add node ages to edge list:
  for(i in 1:length(EdgeList)) EdgeList[[i]]$NodeAgeFromTo <- BranchAges[i, ]
  
  # Add node ages (from-to) to each edge in list:
  EdgeList <- lapply(EdgeList, function(x) {x$BranchDuration <- x$NodeAgeFromTo[1] - x$NodeAgeFromTo[2]; return(x)})
  
  # Get vector of branch types:
  BranchTypes <- gsub("0", "Internal", gsub("1", "Terminal", as.numeric(tree$edge[, 2] <= ape::Ntip(tree))))
  
  # Add branch type to edge list:
  for(i in 1:length(EdgeList)) EdgeList[[i]]$BranchType <- BranchTypes[i]
  
  # Find descendant edges for each internal node:
  DescendantEdgesForEachInternalNode <- lapply(as.list(InternalNodeNumbers), GetDescendantEdges, tree = tree)
  
  # Get ancestral character states:
  AncestralStates <- AncStateEstMatrix(CladisticMatrix = CladisticMatrix, Tree = tree, EstimateAllNodes = EstimateAllNodes, EstimateTipValues = EstimateTipValues, InapplicablesAsMissing = InapplicablesAsMissing, PolymorphismBehaviour = PolymorphismBehaviour, UncertaintyBehaviour = UncertaintyBehaviour, Threshold = Threshold)
  
  # Build single matrix of all states in tip label then node number order:
  AllStates <- do.call(cbind, lapply(lapply(AncestralStates[2:length(AncestralStates)], '[[', "Matrix"), function(x) x[c(tree$tip.label, 1:ape::Nnode(tree) + ape::Ntip(tree)), , drop = FALSE]))
  
  # Make vector of ordering of characters:
  Ordering <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Ordering")))
  
  # Make vector of weights of characters:
  Weights <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "Weights")))
  
  # Make vector of minimum values:
  MinVals <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MinVals")))
  
  # Make vector of maximum values:
  MaxVals <- unname(unlist(lapply(CladisticMatrix[2:length(CladisticMatrix)], '[[', "MaxVals")))
  
  # Find positions in matrix with polymorphisms:
  PolymorphismPositions <- grep("&", AllStates)
  
  # Find positions in matrix with uncertainties:
  UncertaintyPositions <- grep("/", AllStates)
  
  # Find positions in matrix with inapplicables:
  InapplicablePositions <- which(AllStates == "")

  # If polymorphisms were found:
  if(length(PolymorphismPositions) > 0) {
    
    # If replacing polymorphsims with missing do so:
    if(PolymorphismState == "missing") AllStates[PolymorphismPositions] <- NA
    
    # If replacing polymorphisms with random values draw and replace:
    if(PolymorphismState == "random") AllStates[PolymorphismPositions] <- unlist(lapply(strsplit(AllStates[PolymorphismPositions], "&"), sample, size = 1))
    
  }
  
  # If uncertainties were found:
  if(length(UncertaintyPositions) > 0) {
    
    # If replacing uncertainties with missing do so:
    if(UncertaintyState == "missing") AllStates[UncertaintyPositions] <- NA
    
    # If replacing uncertainties with random values draw and replace:
    if(UncertaintyState == "random") AllStates[UncertaintyPositions] <- unlist(lapply(strsplit(AllStates[UncertaintyPositions], "/"), sample, size = 1))
    
  }
  
  # If inapplicable states were found:
  if(length(InapplicablePositions) > 0) {
    
    # If replacing inapplicables with missing do so:
    if(InapplicableState == "missing") AllStates[InapplicablePositions] <- NA
    
  }
  
  # Set default converted continuous characters to FALSE:
  ContinuousCharactersConverted <- FALSE
  
  # Check for continuous characters as these will need to be modified for modelling to discrete characters:
  if(any(Ordering == "cont")) {
    
    # Tell user this is happening:
    cat("Continuous characters found. Converting to gap-weighted discrete characters.\n")
    
    # Set default converted continuous characters to TRUE:
    ContinuousCharactersConverted <- TRUE
    
    # Find out which characters are continuous:
    ContinuousCharactersFound <- which(Ordering == "cont")
    
    # Rescale continous characters as zero to one values:
    ListOfContinuousValuesRescaledZeroToOne <- lapply(lapply(lapply(apply(AllStates[, ContinuousCharactersFound, drop = FALSE], 2, list), unlist), as.numeric), function(x) {x <- x - min(sort(x)); x <- x / max(sort(x)); return(x)})
    
    # Now discretize and store these characters (0 to 31 scale):
    AllStates[, ContinuousCharactersFound] <- do.call(cbind, lapply(lapply(lapply(ListOfContinuousValuesRescaledZeroToOne, function(x) as.list(x)), lapply, function(x) ifelse(is.na(x), NA, max(which(x >= (0:31) / 31)) - 1)), unlist))
    
    # Convert character type to ordered:
    Ordering[ContinuousCharactersFound] <- "ord"
    
    # Convert weights to 1/31:
    Weights[ContinuousCharactersFound] <- 1 / 31
    
    # Set minimum value to zero:
    MinVals[ContinuousCharactersFound] <- 0
    
    # Set maximum value to 31:
    MaxVals[ContinuousCharactersFound] <- 31
    
  }
  
  # If EnsureAllWeightsAreIntegers is TRUE rescale weights until they are all integers so can model appropriately with Poisson later:
  if(EnsureAllWeightsAreIntegers) while(is.character(all.equal(sum(Weights %% 1), 0))) Weights <- (1 / (Weights %% 1)[(Weights %% 1) > 0])[1] * Weights

  # Add from-to node states for each character to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$CharacterStatesFromTo <- matrix(AllStates[x$NodeNumberFromTo, , drop = FALSE], nrow = 2, dimnames = list(c("From", "To"))); return(x)})
  
  # Subfunction to define character changes:
  BuildChangesMatrix <- function(x) {
    
    # Find only comparable characters (those scored for both from and to states):
    ComparableCharacters <- which(apply(!apply(x$CharacterStatesFromTo, 2, is.na), 2, all))
    
    # Isolate comparable ordering:
    ComparableOrdering <- Ordering[ComparableCharacters]
    
    # Isolate comparable weights:
    ComparableWeights <- Weights[ComparableCharacters]
    
    # Isolate only characters that actually differ (change):
    CharacterDifferences <- which(x$CharacterStatesFromTo["From", ComparableCharacters] != x$CharacterStatesFromTo["To", ComparableCharacters])
    
    # Build character change matrix:
    CharacterChanges <- matrix(nrow = 0, ncol = 5, dimnames = list(c(), c("Character", "From", "To", "Steps", "Weight")))
    
    # If characters change then make a matrix from them:
    if(length(CharacterDifferences) > 0) CharacterChanges <- rbind(CharacterChanges, cbind(as.numeric(ComparableCharacters[CharacterDifferences]), as.numeric(x$CharacterStatesFromTo["From", ComparableCharacters[CharacterDifferences]]), as.numeric(x$CharacterStatesFromTo["To", ComparableCharacters[CharacterDifferences]]), ifelse(ComparableOrdering[CharacterDifferences] == "unord", 1, abs(as.numeric(x$CharacterStatesFromTo["To", ComparableCharacters[CharacterDifferences]]) - as.numeric(x$CharacterStatesFromTo["From", ComparableCharacters[CharacterDifferences]]))), ComparableWeights[CharacterDifferences]))
    
    # Store character changes as new sublist for x:
    x$CharacterChanges <- CharacterChanges
    
    # Store comparable characters as new sublist of x:
    x$ComparableCharacters <- ComparableCharacters
    
    # Return x:
    return(x)
    
  }
  
  # Get character changes and comparable characters and add to edge list:
  EdgeList <- lapply(EdgeList, BuildChangesMatrix)
  
  # Check whether time bins are being compared (otherwise no need to assign character changes):
  if(!is.null(TimeBinPartitionsToTest)) {
    
    # Subfunction to add change times to character changes:
    AddChangeTimes <- function(x, ChangeTimes) {
      
      # Isolate character changes:
      CharacterChanges <- x$CharacterChanges
      
      # If any changes involve two or more steps (requiring replacement with multiple changes):
      if(any(CharacterChanges[, "Steps"] > 1)) {
        
        # Get multistep character changes:
        MultiStepCharacters <- which(CharacterChanges[, "Steps"] > 1)
        
        # For each multistep character change:
        for(i in rev(MultiStepCharacters)) {
          
          # Isolate other rows:
          OtherRowNumbers <- setdiff(1:nrow(CharacterChanges), i)
          
          # Get unpacked changes (X:Y, e.g., 0:2 would become 0 1 2):
          UnpackedChanges <- CharacterChanges[i, "From"]:CharacterChanges[i, "To"]
          
          # Update character changes with multistep changes unpacked:
          CharacterChanges <- rbind(CharacterChanges[OtherRowNumbers, ], unname(cbind(rep(CharacterChanges[i, "Character"], length.out = length(UnpackedChanges) - 1), UnpackedChanges[1:(length(UnpackedChanges) - 1)], UnpackedChanges[2:length(UnpackedChanges)], rep(1, length.out = length(UnpackedChanges) - 1), rep(CharacterChanges[i, "Weight"], length.out = length(UnpackedChanges) - 1))))
          
        }
        
        # Resort character changes by character number:
        CharacterChanges <- CharacterChanges[order(CharacterChanges[, "Character"]), ]
        
      }
      
      # If using midpoint option set character change times as midpoint of branch:
      if(ChangeTimes == "midpoint") CharacterChanges <- cbind(CharacterChanges, rep(x$NodeAgeFromTo[1] - (x$BranchDuration / 2), length.out = nrow(CharacterChanges)))
      
      # If using spaced then set character change times as equally spaced along branch:
      if(ChangeTimes == "spaced") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - (seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1)[1:nrow(CharacterChanges)] + (diff(seq(from = 0, to = x$BranchDuration, length.out = nrow(CharacterChanges) + 1))[1] / 2)))
      
      # If using random then set character change times as random draws from a uniform distribution:
      if(ChangeTimes == "random") CharacterChanges <- cbind(CharacterChanges, x$NodeAgeFromTo[1] - runif(n = nrow(CharacterChanges), min = 0, max = x$BranchDuration))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Time"
      
      # Subfunction to re-sort character change times so they occur in correct order:
      SortChangeTimes <- function(CharacterChanges) {
        
        # Sort change time for each character from oldest (first) to youngest (last) and store it:
        CharacterChanges[, "Time"] <- unname(unlist(lapply(as.list(unique(CharacterChanges[, "Character"])), function(x) sort(CharacterChanges[which(CharacterChanges[, "Character"] == x), "Time"], decreasing = TRUE))))
        
        # Return sorted character changes:
        return(CharacterChanges)
        
      }
      
      # Re-sort character change times so they occur in correct order:
      CharacterChanges <- SortChangeTimes(CharacterChanges)
      
      # Add bin for character change as last column:
      CharacterChanges <- cbind(CharacterChanges, unlist(lapply(as.list(CharacterChanges[, "Time"]), function(x) max(which(x <= TimeBins)))))
      
      # Add column name to change time column:
      colnames(CharacterChanges)[ncol(CharacterChanges)] <- "Bin"
      
      # Overwrite character changes with new version with changes added:
      x$CharacterChanges <- CharacterChanges
      
      # Return x:
      return(x)
    
    }
    
    # Add character change times to edge list:
    EdgeList <- lapply(EdgeList, AddChangeTimes, ChangeTimes = ChangeTimes)
    
  }
  
  # Subfunction to get edge sections in time bins:
  EdgeSectionsInBins <- function(x, TimeBins = TimeBins) {
    
    # Set first appearance datum of edge:
    FAD <- x$NodeAgeFromTo[1]
    
    # Set last appearance datum of edge:
    LAD <- x$NodeAgeFromTo[2]
    
    # Get any time bin boundaries crossed (an be empty if none are):
    BoundariesCrossed <- TimeBins[2:(length(TimeBins) - 1)][intersect(which(TimeBins[2:(length(TimeBins) - 1)] > LAD), which(TimeBins[2:(length(TimeBins) - 1)] < FAD))]
    
    # If boundaries are crossed:
    if(length(BoundariesCrossed) > 0) {
      
      # Break up branch into binned sections as vector of FADs:
      FAD <- c(FAD, BoundariesCrossed)
      
      # Break up branch into binned sections as vector of LADs:
      LAD <- c(BoundariesCrossed, LAD)
      
    }
    
    # Build matrix of branch sections with FADs and LADs:
    BranchSections <- rbind(FAD, LAD)
    
    # Add bin number present in to column names:
    colnames(BranchSections) <- unlist(lapply(lapply(lapply(as.list(BranchSections["FAD", ]), '<=', TimeBins), which), max))
    
    # Add new list section for branch (edge) sections binned by time:
    x$BinnedEdgeSections <- BranchSections
    
    # Return output:
    return(x)
    
  }
  
  # Get edge sections in time bins:
  EdgeList <- lapply(EdgeList, EdgeSectionsInBins, TimeBins = TimeBins)
  
  # Add binned branch durations to edge list:
  EdgeList <- lapply(EdgeList, function(x) {BranchDurations <- rep(0, length(TimeBins) - 1); BranchDurations[as.numeric(colnames(x$BinnedEdgeSections))] <- abs(apply(x$BinnedEdgeSections, 2, diff)); x$BinnedBranchDurations <- BranchDurations; return(x)})
  
  # Add proportional binned branch lengths to edge list:
  EdgeList <- lapply(EdgeList, function(x) {x$ProportionalBinnedEdgeDurations <- x$BinnedBranchDurations / sum(x$BinnedBranchDurations); return(x)})
  
  # Start to build matrix of all changes with list of character changes:
  AllChanges <- lapply(EdgeList, function(x) x$CharacterChanges)
  
  # Add edge number to each matrix of character changes:
  for(i in 1:length(AllChanges)) AllChanges[[i]] <- cbind(rep(i, times = nrow(AllChanges[[i]])), AllChanges[[i]])
  
  # Combine all changes into a single matrix:
  AllChanges <- do.call(rbind, lapply(AllChanges, function(x) {colnames(x)[1] <- "Edge"; x}))
  
  # Remove silly rownames from all changes:
  rownames(AllChanges) <- NULL

  # If doing some kind of edge test (branch or clade):
  if(!is.null(BranchPartitionsToTest) || !is.null(CladePartitionsToTest)) {
    
    # Get (weighted) number of changes on each edge:
    EdgeChanges <- unlist(lapply(EdgeList, function(x) sum(x$CharacterChanges[, "Steps"] * x$CharacterChanges[, "Weight"])))

    # Get completeness for each edge:
    EdgeCompleteness <- unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters]) / sum(Weights)))
    
    # Get duration of each edge:
    EdgeDurations <- unlist(lapply(EdgeList, function(x) x$BranchDuration))
    
    # Set global rate:
    GlobalRate <- sum(EdgeChanges) / sum(EdgeCompleteness * EdgeDurations)
    
    # If performing branch partition tests:
    if(!is.null(BranchPartitionsToTest)) {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
      PartitionedData <- lapply(BranchPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(EdgeChanges[y]), sum(EdgeCompleteness[y] * EdgeDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
      
      # Add sampled rate to paritioned data matrices:
      PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
      
      # Get P-Values and combine output as edge test results:
      BranchPartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
    
    # If not performing branch partition tests:
    } else {
      
      # Make empty branch partition result output:
      BranchPartitionTestResults <- NULL
      
    }
    
    # If performing clade partition tests:
    if(!is.null(CladePartitionsToTest)) {
      
      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
      PartitionedData <- lapply(CladePartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(EdgeChanges[y]), sum(EdgeCompleteness[y] * EdgeDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
      
      # Add sampled rate to paritioned data matrices:
      PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
      
      # Get P-Values and combine output as edge test results:
      CladePartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
      
    # If not performing clade partition tests:
    } else {
      
      # Make empty clade partition result output:
      CladePartitionTestResults <- NULL
      
    }
    
  # If not doing clade OR branch tests:
  } else {
    
    # Make empty branch partition result output:
    BranchPartitionTestResults <- NULL

    # Make empty clade partition result output:
    CladePartitionTestResults <- NULL

  }
  
  # If performing branch partition tests:
  if(!is.null(CharacterPartitionsToTest)) {
    
    # Get vector of (weighted) changes for each character:
    CharacterChanges <- unlist(lapply(as.list(CharacterNumbers), function(x) {CharacterRows <- which(AllChanges[, "Character"] == x); sum(AllChanges[CharacterRows, "Steps"] * AllChanges[CharacterRows, "Weight"])}))
    
    # Get vector of weighted durations for each character:
    CharacterDurations <- (Weights / sum(Weights)) * sum(tree$edge.length)
    
    # Get vector of completness (opportunity to observe changes) for each character:
    CharacterCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) {CharacterPresence <- rep(0, times = length(CharacterNumbers)); CharacterPresence[x$ComparableCharacters] <- 1; CharacterPresence * x$BranchDuration})), 2, sum) / sum(tree$edge.length)
    
    # Set global rate:
    GlobalRate <- sum(CharacterChanges) / sum(CharacterCompleteness * CharacterDurations)
    
    # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
    PartitionedData <- lapply(CharacterPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(CharacterChanges[y]), sum(CharacterCompleteness[y] * CharacterDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
    
    # Add sampled rate to paritioned data matrices:
    PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
    
    # Get P-Values and combine output as edge test results:
    CharacterPartitionTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})

  # If performing branch partition tests:
  } else {
    
    # Make empty character partition result output:
    CharacterPartitionTestResults <- NULL
    
  }
  
  # If performing time bin partition tests:
  if(!is.null(TimeBinPartitionsToTest)) {
    
    # Get weighted number of changes from each time bin:
    TimeBinChanges <- unlist(lapply(as.list(1:(length(TimeBins) - 1)), function(x) {ChangeRows <- AllChanges[, "Bin"] == x; sum(AllChanges[ChangeRows, "Steps"] * AllChanges[ChangeRows, "Weight"])}))
    
    # If using the Close time bin completeness approach get completeness value for each time bin:
    if(TimeBinApproach == "Close") TimeBinCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) x$ProportionalBinnedEdgeDurations * (sum(Weights[x$ComparableCharacters]) / sum(Weights)))), 2, sum) / apply(do.call(rbind, lapply(EdgeList, function(x) x$ProportionalBinnedEdgeDurations)), 2, sum)
    
    # If using the Lloyd time bin completeness approach get completeness value for each time bin::
    if(TimeBinApproach == "Lloyd") TimeBinCompleteness <- apply(do.call(rbind, lapply(EdgeList, function(x) apply(matrix(Weights[x$ComparableCharacters], ncol = 1) %*% x$BinnedBranchDurations, 2, sum))), 2, sum) / apply(do.call(rbind, lapply(EdgeList, function(x) apply(matrix(Weights, ncol = 1) %*% x$BinnedBranchDurations, 2, sum))), 2, sum)
    
    # Get durations of edges in each time bin:
    TimeBinDurations <- apply(do.call(rbind, lapply(EdgeList, function(x) x$BinnedBranchDurations)), 2, sum)
    
    # Set global rate (NB: will differ between Close and Lloyd approaches, but Lloyd approach will match edge or character global rate):
    GlobalRate <- sum(TimeBinChanges) / sum(TimeBinCompleteness * TimeBinDurations)
    
    # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
    PartitionedData <- lapply(TimeBinPartitionsToTest, function(x) matrix(unlist(lapply(x, function(y) c(sum(TimeBinChanges[y]), sum(TimeBinCompleteness[y] * TimeBinDurations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("Changes", "Completeness", "Duration"))))
    
    # Add sampled rate to paritioned data matrices:
    PartitionedData <- lapply(PartitionedData, function(x) {x <- cbind(as.numeric(gsub(NaN, 0, c(x[, "Changes"] / (x[, "Completeness"] * x[, "Duration"])))), x); colnames(x)[1] <- "Rate"; x})
    
    # Get P-Values and combine output as edge test results:
    TimeBinTestResults <- lapply(PartitionedData, function(x) {x <- list(x[, "Rate"], GetMaximumLikelihoodPValue(MeanRate = GlobalRate, SampledRates = x[, "Rate"], SampledChanges = x[, "Changes"], SampledCompleteness = x[, "Completeness"], SampledTime = x[, "Duration"])); names(x) <- c("Rates", "PValue"); x})
    
  # If not performing time bin partition tests:
  } else {
    
    # Make empty time bin partition result output:
    TimeBinTestResults <- NULL
    
  }
  
  # Set global rate for output:
  GlobalRate <- sum(unlist(lapply(EdgeList, function(x) sum(x$CharacterChanges[, "Steps"] * x$CharacterChanges[, "Weight"])))) / sum(unlist(lapply(EdgeList, function(x) sum(Weights[x$ComparableCharacters]) / sum(Weights))) * unlist(lapply(EdgeList, function(x) x$BranchDuration)))

  # Subfunction to calculate adjusted alphas for multiple comparison corrections:
  AddMultipleComparisonCorrectionCutoffs <- function(TestResults, Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection) {
    
    # Get number of comparisons performed:
    NComparisons <- length(TestResults)
    
    # If using the Benjamini-Hochberg false discovery rate approach:
    if(MultipleComparisonCorrection == "BenjaminiHochberg") {
      
      # Set cutoff values:
      CutoffValues <- ((1:NComparisons) / NComparisons) * Alpha
      
      # Get actual p-values found:
      PValues <- unlist(lapply(TestResults, '[[', "PValue"))
      
      # Order cutoffs by p-value rank:
      CutoffValues <- CutoffValues[rank(PValues, ties.method = "random")]
      
    }
    
    # If using the Bonferroni correction set cutoff values as alpha over N:
    if(MultipleComparisonCorrection == "Bonferroni") CutoffValues <- Alpha / NComparisons
    
    # Add cutoffs to output:
    for(i in 1:length(TestResults)) TestResults[[i]]$CorrectedAlpha <- CutoffValues[i]
    
    # Return modified test results:
    return(TestResults)
    
  }
  
  # If doing branch partition tests then add multiple comparison alpha cutoffs:
  if(!is.null(BranchPartitionsToTest)) BranchPartitionTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = BranchPartitionTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
  
  # If doing character partition tests then add multiple comparison alpha cutoffs:
  if(!is.null(CharacterPartitionsToTest)) CharacterPartitionTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = CharacterPartitionTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
  
  # If doing clade partition tests then add multiple comparison alpha cutoffs:
  if(!is.null(CladePartitionsToTest)) CladePartitionTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = CladePartitionTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
  
  # If doing time bin partition tests then add multiple comparison alpha cutoffs:
  if(!is.null(TimeBinPartitionsToTest)) TimeBinTestResults <- AddMultipleComparisonCorrectionCutoffs(TestResults = TimeBinTestResults, Alpha = Alpha, MultipleComparisonCorrection = MultipleComparisonCorrection)
  
  # Compile output:
  Output <- list(AllChanges, GlobalRate, ContinuousCharactersConverted, BranchPartitionTestResults, CharacterPartitionTestResults, CladePartitionTestResults, TimeBinTestResults)
  
  # Add naems to output:
  names(Output) <- c("InferredCharacterChanges", "IntrinsicCharacterRate", "ContinuousCharactersConvertedToDiscrete", "BranchPartitionResults", "CharacterPartitionResults", "CladePartitionResults", "TimeBinResults")
  
  # Return output:
  return(Output)

}
