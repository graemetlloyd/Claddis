#' Discrete character rates across trees, time, and character types
#'
#' @description
#'
#' Given a tree and a cladistic-type matrix uses either likelihood ratio tests or the Akaike Information Criterion to compare rate models across branches, clades, time bins, or character partitions.
#'
#' @param time_tree A tree (phylo object) with branch durations that represents the relationships of the taxa in \code{cladistic_matrix}.
#' @param cladistic_matrix A character-taxon matrix in the format imported by \link{read_nexus_matrix}.
#' @param time_bins A vector of ages (in millions of years) indicating the boundaries of a series of time bins in order from oldest to youngest.
#' @param branch_partitions A list of branch(es) (edge number) partitions to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param character_partitions A list of character partition(s) (character numbers) to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param clade_partitions A list of clade partition(s) (node numbers) to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param time_partitions A list of time bin partition(s) (numbered 1 to N) to test as N-rate parameter model (where N is the total number of partitions). If NULL (the default) then no partition test(s) will be made.
#' @param change_times The time at which to record the character changes. One of \code{"midpoint"} (changes occur at the midpoint of the branch), \code{"spaced"} (changes equally spaced along branch), or \code{"random"} (change times drawn at random from a uniform distribution; the default and recommended option). Note: this is only meaningful if testing for time bin partitions.
#' @param test_type Whether to apply an Akaike Information Criterion (\code{"aic"}; the default) or likelihood ratio test (\code{"lrt"}).
#' @param alpha The alpha value to be used for the significance tests. The default is 0.01. Thsi is only relevant if using likelihood ratio tests.
#' @param multiple_comparison_correction Current options are: 1. \code{"benjaminihochberg"} (the Benjamini and Hochberg 1995 false discovery rate approach; default and recommended), or 2. \code{"bonferroni"} (the Bonferroni correction). This is only relevant if using likelihood ratio tests.
#' @param polymorphism_state Current options are: 1. \code{"missing"} (converts polymorphic values to NA; the default), or 2. \code{"random"} (picks one of the possible polymorphic states at random).
#' @param uncertainty_state Current options are: 1. \code{"missing"} (converts uncertain values to NA; the default), or 2. \code{"random"} (picks one of the possible uncertain states at random).
#' @param inapplicable_state The only current option is \code{"missing"} (converts value to NA).
#' @param time_binning_approach One of \code{"close"} or \code{"lloyd"} (the latter is the default and recommended option).
#' @param all_weights_integers Logical for whether (\code{TRUE}) to reweight non-integer weights until all weights are integers or to leave them as they are (\code{FALSE}; the default).
#' @param estimate_all_nodes Option passed to \link{estimate_ancestral_states}.
#' @param estimate_tip_values Option passed to \link{estimate_ancestral_states}.
#' @param inapplicables_as_missing Option passed to \link{estimate_ancestral_states}.
#' @param polymorphism_behaviour Option passed to \link{estimate_ancestral_states}.
#' @param uncertainty_behaviour Option passed to \link{estimate_ancestral_states}.
#' @param threshold Option passed to \link{estimate_ancestral_states}.
#' @param all_missing_allowed Option passed to \link{estimate_ancestral_states}.
#'
#' @details
#'
#' \bold{Introduction}
#'
#' Morphological change can be captured by discrete characters and their evolution modelled as occurring along the branches of a phylogenetic tree. This function takes as primary input a character-taxon matrix of discrete characters (in the format imported by \link{read_nexus_matrix}) and a time-scaled phylogenetic tree (in the format of \pkg{paleotree} or \pkg{strap}) and begins by inferring ancestral states at the tree's internal nodes using the \link{estimate_ancestral_states} function. From here changes along individual branches can be estimated (only the minimum number of changes are inferred; see \link{map_stochastic_changes} for an alternative but unfinished approach) and hence rates can be calculated.
#'
#' A discrete character rate can be expressed as the mean number of changes per million years (users may wish to normalise this by the number of characters for interpretation) and can be calculated for a branch (edge) of the tree, a clade (a mean rate for the edges descended from a single node), a character partition (the mean rate for a subset of the characters across all edges), or, most complex (see Lloyd 2016), the mean rate across the edges (or parts of edges) present in a time bin (defined by two values denoting the beginning and end of the time bin). In an ideal scenario these rates could be compared at face value, but that would require a large number of characters and very minimal (or zero) missing data. I.e., at an extreme of missing data if only one character can be observed along a branch it will either change (the maximum possible inferrable rate of evolution) or it will not (the minimum possible inferrable rate of evolution). In such cases it would be unwise to consider either outcome as being a significant departure from the mean rate.
#'
#' Because of these complications Lloyd et al. (2012) introduced tests by which the significance of an edge (or other partitioning of the data, i.e., a clade, time bin etc.) could be considered to be significantly high or low in comparison to the mean rate for the whole tree (i.e., whether a two-rate model could be considered more likely than a one-rate model). This is achieved through a likelihood ratio test:
#'
#' \deqn{LR = value of likehood function under the null (one-rate) hypothesis / maximum possible value of likehood function under the alternative (two-rate) hypotheses}
#'
#' Typically we might expect the two hypotheses to be well defined a priori. E.g., an expectation that a specific branch of the tree might have a higher or lower rate than background due to some evolutionary shift. However, Lloyd et al. (2012) instead provided an exploratory approach whereby every possible one edge value was compared with the rate for the rest of the tree (and the equivalent with clades and time bins). This was the default in Claddis up to version 0.2, but this has now been replaced (since version 0.3) with a more customisable set of options that allows different types of hypotheses (e.g., partitioning the data by character), as well as more complex hypotheses (e.g., a three-rate model), to be tested. Since version 0.4 the option to replace likelihood ratio tests with the Akaike Information Criterion has also been added.
#'
#' \bold{The four types of rate hypothesis}
#'
#' Following Cloutier (1991), Lloyd (2016) extended the two main types of rate hypotheses to four:
#'
#' \enumerate{
#'   \item A branch rate (available here with the \code{branch_partitions} option).
#'   \item A clade rate (available here with the \code{clade_partitions} option).
#'   \item A time bin rate (available here with the \code{time_partitions} option).
#'   \item A character partition rate (available here with the \code{character_partitions} option).
#' }
#'
#' In Claddis (>=0.3) these partitions are defined as a list of lists of vectors where only the first N - 1 partitions need be defined. E.g., if comparing the first edge value (based on \pkg{ape} numbering, i.e., \code{plot(tree); edgelabels()}) to the rest of the tree then the user only needs to define the value "1" and the function will automatically add a second partition containing all other edges. This can be set with the option \code{branch_partitions = list(list(1))}. Similarly, to do what Lloyd et al. (2012) did and repeat the test for every edge in the tree (and assuming this variable is already named "tree") you could use, \code{branch_partitions = lapply(X = as.list(x = 1:nrow(tree$edge)), as.list)}.
#'
#' Because of the flexibility of this function the user can define any set of edges. For example, they could test whether terminal branches have a different rate from internal branches with \code{branch_partitions = list(list(match(1:ape::Ntip(phy = tree), tree$edge[, 2])))}. The \code{clade_partitions} is really just a special subset of this type of hypothesis, but with edges being defined as descending from a specific internal node in the tree. Once again, an exploratory approach like that of Lloyd et al. (2012) can be used with: \code{clade_partitions = lapply(X = as.list(x = ape::Ntip(phy = tree) + (2:Nnode(tree))), as.list)}. Note that this excludes the root node as this would define a single partition and hence would represent the null hypothesis (a single rate model for the whole tree). (If using \code{test_type = "aic"} then the user typically \emph{will} want a value for a single partition.) More generally clades must be defined by the node numbers they correspond to. In R an easy way to identify these is with: \code{plot(tree); nodelabels()}.
#'
#' Time bin partitions are defined in a similar way, but are numbered 1:N starting from the oldest time bin. So if wanting to do an exploratory test of single bin partitions (and only four time bins were specified) you could use: \code{time_partitions = lapply(X = as.list(x = 1:4), as.list)}. Bins can be combined too, just as edges are above. For example, time bins 1 and 2 could form a single partition with: \code{time_partitions = list(list(1:2))}. Or if looking to test a model where each bin has its' own rate value you could use: \code{time_partitions = list(as.list(x = 1:3))}. Note, as before we do not need to specify the fourth bin as this will be automatically done by the function, however, \code{time_partitions = list(as.list(x = 1:4))} will also work. Some caution needs to be applied with N-rate models (where N is three or larger) and \code{test_type = "lrt"} as a result favouring such models does not necessarily endorse N-separate rates. I.e., it could simply be that one bin has such a large excursion that overall the N-rate model fits better than the 1-rate model, but some 2-rate models might be better still. It is up to the user to check this themselves by exploring smaller combinations of bins and more genrally if exploring partitions of three or more use of the Akaike Information Criterion (\code{test_type = "aic"}) is recommended.
#'
#' Finally, character partitions allow the user to explore whether rates vary across different character types (numbers), e.g., skeletal characters versus soft tissue characters, or cranial characters versus postcranial characters. Here characters are simply numbered 1:N (across all blocks of a matrix), but single character partitions are less likely to be of interest. As an example of use lets say the first ten characters are what we are interested in as a partition (the second partition being the remaining characters), we could use: \code{character_partitions = list(list(1:10))} to test for a two-rate model with \code{test_type = "lrt"}.
#'
#' Note that the list of lists structure is critical to defining partitions as it allows them to be of different sizes and numbers. For example, one partition of three and another of six, or one set of two partitions and another set of four partitions - structures not easily realizable using vectors or matrices. However, it may not be intuitive to some users so it is recommended that the user refers to the examples above as a guide.
#'
#' Additionally, it should be noted that the user can test multiple types of hypotheses simultaneously with the function. For example, performing several branch tests whilst also performing clade tests. However, it is not necessary to perform all types simultaneously (as was the case up to version 0.2) and unused partition types can be set to NULL, the default in each case.
#'
#' \bold{AIC vs LRT}
#'
#' Since Claddis version 0.4 the option to use the Akaike Information Criterion (AIC) instead of likelihood ratio tests (LRTs) has been added (although it was not properly functional until version 0.4.2). Practically speaking the AIC uses something similar to the denominator term from the LRT (see equation above) and adds a penalty for the number of parameters (partitions). However, it also fundamentally alters the comparative framework applied and hence needs more careful attention from the user to be applied correctly. Specifically, the LRT is by its nature comparative, always comparing an N-rate partition with a one-rate partition. By contrast the AIC does not directly apply any comparison and so the user must logically supply multiple partitionings of the data in order for the results to be meaningful. It might be assumed that a user will always want to apply a single partition that pools all the data for each type of test, whether this is all the edges (branches), time bins, or characters. This will thus be the obvious comparator for any multiple partition supplied, ensuring that any more complex partitioning is minimally superior to this. Additionally, it is also logical to consider each possible way of joining partitions simpler than the most complex partition being considered. E.g., if considering a four-partition model then the user should also consider all possible three-partition and two-partition combinations of that four-partition model. This can obviously lead to some complexity in supplying partitions on the user's part and so some automating of this process is planned in future (but is not fully available yet). For an example, see \link{partition_time_bins}.
#'
#' Additionally, AIC values are not simple to compare as there is no direct equivalent of the alpha value from the LRT. Instead the user can modify the AIC values returned themselves to get delta-AIC or Akaike weight values (e.g., with \code{geiger::aicw}). (NB: I will not explain these here as there are better explanations online.) Furthermore, since version 0.4.2 sample-size corrected AIC (AICc) is also available in the output. Note that some caution should be used in applying this if the number of partitions is equal to the sample size or only one fewer. E.g., if you have ten time bins and ten partitions you can get a negative value due to the denominator term in the AICc calculation. Thus it is advised to use the raw AIC values as first approximation and be wary if the AICc flips the preferred model to a more complex model or models (i.e., those with more partitions) as this is the opposite of the intent of the AICc.
#'
#' \bold{High versus low rates}
#'
#' Prior to Claddis version 0.3, rate results were classified as significantly high or significntly low as part of the output. This was done simply on whether the estimated p-value fell above or below the corrected alpha level (the significance threshold) and whether the first part of the two-rate partition had a higher or lower rate than the second part (the pooling of all other values). This simple interpretation is no longer valid here as the function can consider more than two partitions (high versus low is meaningless) and the allowing of AIC values means a significance test need not be performed either. Although the same interepretation can still be applied manually when only using two-partition tests and the LRT, other partition sizes and use of the AIC complicate this interpretation (in the same way an ANOVA is more complex than a two-sample t-test). This will also affect visualisation of the data (i.e. the simple pie chart coluring of non-significant, significntly high or significantly low rates seeen since Lloyd et al. 2012 will no longer apply). Instead the user should isolate the best model(s) and attempt to visualise these, perhaps using something like a heat map, with the mean rate for each partition being represented by an appropriate colour.
#'
#' \bold{Other options}
#'
#' Since Claddis version 0.3 this function has allowed the user greater control with many more options than were offered previously and these should be considered carefully before running any tests.
#'
#' Firstly, the user can pick an option for \code{change_times} which sets the times character changes are inferred to occur. This is only relevant when the user is performing time bin partition test(s) as this requires some inference to be made about when changes occur on branches that may span multiple time bins. The current options are: \code{"midpoint"} (all changes are inferred to occur midway along the branch, effectively mimicking the approach of Ruta et al. 2006), \code{"spaced"} (all changes are inferred to occur equally spaced along the branch, with changes occurring in character number order), or \code{"random"} (changes are assigned a random time by drawing from a uniform distribution between the beginning and end of each branch). The first of these is likely to lead to unrealistically "clumped" changes and by extension implies completely correlated character change that would violate the assumptions of the Poisson distribution that underlies the significance tests here (Lloyd et al. 2012). At the other extreme, the equally spaced option will likely unrealistically smooth out time series and potentially make it harder to reject the single-rate null (leading to Type II errors). For these reasons, the random option is recommended and is set as the default. However, because it is random this makes the function stochastic (the answer can vary each time it is run) and so the user should therefore run the function multiple times if using this option (i.e., by using a for loop) and aggregating the results at the end (e.g., as was done by previous authors; Lloyd et al. 2012; Close et al. 2015).
#'
#' Secondly, the \code{alpha} value sets the significance threshold by which the likelihood ratio test's resulting p-value is compared (i.e., it is only reelevant when \code{test_type = "lrt"}. Following Lloyd et al. (2012) this is set lower (0.01) than the standard 0.05 value by default as those authors found rates to be highly heterogenous in their data set (fossil lungfish). However, this should not be adopted as a "standard" value without question (just as 0.05 shouldn't). Note that the function also corrects for multiple comparisons (using the \code{multiple_comparison_correction} option) to avoid Type I errors (false positives). It does so (following Lloyd et al. 2012) using the Benjamini-Hochberg (Benjamini and Hochberg 1995) False Discovery Rate approach (see Lloyd et al. 2012 for a discussion of why), but the Bonferroni correction is also offered here (albeit not recommended).
#'
#' Thirdly, polymorphisms and uncertainities create complications for assessing character changes along branches. These can occur at the tips (true polymorphisms or uncertainties in sampled taxa) and internal nodes (uncertainty over the estimated ancestral state). There are two options presented here, and applicable to both \code{polymorphism_state} and \code{uncertainty_state} (allowing these to be set separately). These are to convert such values to missing (NA) or to pick one of the possible states at random. Using missing values will increase overall uncertainty and potentially lead to Type II errors (false negatives), but represents a conservative solution. The random option is an attempt to avoid Type II errors, but can be considered unrealistic, especially if there are true polymorphisms. Additionally, the random option will again make the function stochastic meaning the user should run it multiple times amd aggregate the results. Note that if there are no polymorphisms or uncertainties in the character-taxon matrix the latter can still occur with ancestral state estimates, especially if the threshold value is set at a high value (see \link{estimate_ancestral_states} for details).
#'
#' Fourthly, inapplicable characters can additionally complicate matters as they are not quite the same as missing data. I.e., they can mean that change in a particular character is not even possible along a branch. However, there is no easy way to deal with such characters at present so the user is not presented with a true option here - currently all inapplicable states are simply converted to missing values by the function. In future, though, other options may be available here. For now it is simply noted that users should be careful in making inferences if there are inapplicable characters in their data and should perhaps consider removing them with \link{prune_cladistic_matrix} to gauge their effect.
#'
#' Fifthly, there are currenty two further options for assessing rates across time bins. As noted above a complication here is that character changes (the rate numerator) and character completeness (part of the rate denominator) are typically assessed on branches. However, branches will typically span time bin boundaries and hence many bins will contain only some portion of particular branches. The exact portion can be easily calculated for branch durations (the other part of the rate denominator) and the \code{change_times} option above is used to set the rate numerator, however, completeness remains complex to deal with. The first attempt to deal with this was made by Close et al. (2015) who simply used weighted mean completeness by calculating the proportion of a branch in each bin as the weight and multiplying this by each branch's completeness (the \code{"close"} option here). However, this may lead to unrealistic "smoothing" of the data and perhaps more importantly makes no account of which characters are known in a bin. Lloyd (2016) proposed an alternative "subtree" approach which assesses completeness by considering each character to be represented by a subtree where only branches that are complete are retained then branch durations in each bin are summed across subtrees such that the duration term automatically includes completeness (the \code{"lloyd"} option here). Here the latter is strongly recommended, for example, because this will lead to the same global rate across the whole tree as the branch, clade or character partitions, whereas the Close approach will not.
#'
#' Sixthly, all character changes are weighted according to the weights provided in the input character-taxon matrix. In many cases these will simply all be one, although see the equalise weights option in \link{read_nexus_matrix}. However, when weights vary they can create some issues for the function. Specifically, changes are expected to be in the same (integer) units, but if weights vary then they have to be modelled accordingly. I.e., a character twice the weight of another may lead to a single change being counted as two changes. This is most problematic when the user has continuous characters which are automatically converted to gap-weighted (Thiele 1993) characters. However, this conversion creates drastically down-weighted characters and hence the user may wish to set the \code{all_weights_integers} option to TRUE. Note that reweighting will affect the results and hence shifting the weights of characters up or down will necessarily lead to shifts in the relative Type I and II errors. This is an unexplored aspect of such approaches, but is something the user should be aware of. More broadly it is recommended that continuous (or gap-weighted) characters be avoided when using this function.
#'
#' Finally, the remaining options (\code{estimate_all_nodes}, \code{estimate_tip_values}, \code{inapplicables_as_missing}, \code{polymorphism_behaviour}, \code{uncertainty_behaviour}, and \code{threshold}) are all simply passed directly to \link{estimate_ancestral_states} for estimating the ancestral states and users should consult the help file for that function for further details.
#'
#' Note that currently the function cannot deal with step matrices and that the terminal versus internal option from Brusatte et al. (2014) is yet to be implemented.
#'
#' \bold{Output}
#'
#' The output for each LRT test (i.e., the components of the \code{branch_test_results}, \code{character_test_results}, \code{clade_test_results} and \code{time_test_results} parts of the output) includes three main parts:
#'
#' \enumerate{
#'   \item Rates.
#'   \item p_value.
#'   \item CorrectedAlpha.
#' }
#'
#' Or for each AIC test there are:
#'
#' \enumerate{
#'   \item Rates.
#'   \item AIC.
#'   \item AICc.
#' }
#'
#' For each rate test the \code{Rates} part of the output is a vector of the absolute rate (number of changes per million years) for each partition in the test (in the order they were supplied to the function). So, for example, a branch rate for the sixth edge in a tree would be the rate for the sixth edge followed by the pooled rate for all other edges. The length of the vector is the length of the number of partitions.
#'
#' The p_value is a single value indicating the probability that the likelihood ratio (see above and Lloyd et al. 2012) is one, i.e., that the likelihoods of the one-rate and N-rate models are the same.
#'
#' The CorrectedAlpha is the alpha-value that should be used to determine the significance of the current partition test (i.e., The p_value, above). If the p_value exceeds the CorrectedAlpha then the null (single-rate) hypothesis should be accepted, if lower then the null should be rejected in favour of the N-rate hypthesis. Note that the CorrectedAlpha will not typically be the same for each partition and will also typically be different from the input \code{alpha} value due to the \code{multiple_comparison_correction} option used.
#'
#' The AIC is the Akaike Information Criterion, and is relatively meaningless on its own and can only really be used to compare with the AIC values for other partitions of the data. The AICc is simply the sample-size corrected version of the AIC and is preferable when sample sizes are small.
#'
#' @return
#'
#' \item{time_bins_used}{The time binning used (NB: May be slightly altered from the input values).}
#' \item{inferred_character_changes}{Matrix of inferred character changes.}
#' \item{mean_rate}{The global (mean) character rate in changes per million years. I.e, the average rate across all characters, branches and time bins, effectively the null hypothesis for any test performed.}
#' \item{continuous_characters_discretized}{Whether or not continuous characters were converted to discrete characters (important for handling the data in downstream analys(es)).}
#' \item{branch_test_results}{List of branch partition results (corresponding to \code{branch_partitions}. NULL if not requested.}
#' \item{character_test_results}{List of character partition results (corresponding to \code{character_partitions}. NULL if not requested.}
#' \item{clade_test_results}{List of clade partition results (corresponding to \code{clade_partitions}. NULL if not requested.}
#' \item{time_test_results}{List of time bin partition results (corresponding to \code{time_partitions}. NULL if not requested.}
#' \item{branch_rates}{Matrix showing calculated rates for each branch. NULL if \code{branch_partitions} is not requested.}
#' \item{character_rates}{Matrix showing calculated rates for each character. NULL if \code{character_partitions} is not requested.}
#' \item{clade_rates}{Matrix showing calculated rates for each clade. NULL if \code{clade_partitions} is not requested.}
#' \item{time_rates}{Matrix showing calculated rates for each time bin. NULL if \code{time_partitions} is not requested.}
#' \item{time_tree}{The time-scaled input tree used as input (provided as output for use with visualisation functions).}
#'
#' @author Graeme T. Lloyd \email{graemetlloyd@@gmail.com} and Steve C. Wang \email{scwang@@swarthmore.edu}
#'
#' @references
#'
#' Benjamini, Y. and Hochberg, Y., 1995. Controlling the false discovery rate: a practical and powerful approach to multiple testing. \emph{Journal of the Royal Statistical Society, Series B}, \bold{57}, 289-300.
#'
#' Brusatte, S. L., Lloyd, G. T., Wang, S. C. and Norell, M. A., 2014. Gradual assembly of avian body plan culminated in rapid rates of evolution across dinosaur-bird transition. \emph{Current Biology}, \bold{24}, 2386-2392.
#'
#' Close, R. A., Friedman, M., Lloyd, G. T. and Benson, R. B. J., 2015. Evidence for a mid-Jurassic adaptive radiation in mammals. \emph{Current Biology}, \bold{25}, 2137-2142.
#'
#' Cloutier, R., 1991. Patterns, trends, and rates of evolution within the Actinistia. \emph{Environmental Biology of Fishes}, \bold{32}, 23-58.
#'
#' Lloyd, G. T., 2016. Estimating morphological diversity and tempo with discrete character-taxon matrices: implementation, challenges, progress, and future directions. \emph{Biological Journal of the Linnean Society}, \bold{118}, 131-151.
#'
#' Lloyd, G. T., Wang, S. C. and Brusatte, S. L., 2012. Identifying heterogeneity in rates of morphological evolution: discrete character change in the evolution of lungfish (Sarcopterygii; Dipnoi). \emph{Evolution}, \bold{66}, 330-348.
#'
#' Ruta, M., Wagner, P. J. and Coates, M. I., 2006. Evolutionary patterns in early tetrapods. I. Rapid initial diversification followed by decrease in rates of character change. \emph{Proceedinsg of the Royal Society of London B}, \bold{273}, 2107-2111.
#'
#' Thiele, K.. 1993. The Holy Grail of the perfect character: the cladistic treatment of morphometric data. \emph{Cladistics}, \bold{9}, 275-304.
#'
#' @examples
#'
#' # Set random seed:
#' set.seed(seed = 17)
#'
#' # Generate a random tree for the Michaux data set:
#' time_tree <- ape::rtree(n = nrow(michaux_1989$matrix_1$matrix))
#'
#' # Update taxon names to match those in the data matrix:
#' time_tree$tip.label <- rownames(x = michaux_1989$matrix_1$matrix)
#'
#' # Set root time by making youngest taxon extant:
#' time_tree$root.time <- max(diag(x = ape::vcv(phy = time_tree)))
#'
#' # Get discrete character rates:
#' x <- test_rates(
#'   time_tree = time_tree, cladistic_matrix =
#'     michaux_1989, time_bins = seq(
#'     from = time_tree$root.time,
#'     to = 0, length.out = 5
#'   ), branch_partitions =
#'     lapply(X = as.list(x = 1:nrow(time_tree$edge)), as.list),
#'   character_partitions = lapply(
#'     X =
#'       as.list(x = 1:3),
#'     as.list
#'   ), clade_partitions =
#'     lapply(
#'       X =
#'         as.list(x = ape::Ntip(phy = time_tree) + (2:ape::Nnode(phy = time_tree))),
#'       as.list
#'     ), time_partitions =
#'     lapply(X = as.list(x = 1:4), as.list), change_times =
#'     "random", alpha = 0.01, polymorphism_state =
#'     "missing", uncertainty_state = "missing",
#'   inapplicable_state = "missing", time_binning_approach =
#'     "lloyd"
#' )
#' @export test_rates
test_rates <- function(time_tree, cladistic_matrix, time_bins, branch_partitions = NULL, character_partitions = NULL, clade_partitions = NULL, time_partitions = NULL, change_times = "random", test_type = "aic", alpha = 0.01, multiple_comparison_correction = "benjaminihochberg", polymorphism_state = "missing", uncertainty_state = "missing", inapplicable_state = "missing", time_binning_approach = "lloyd", all_weights_integers = FALSE, estimate_all_nodes = FALSE, estimate_tip_values = FALSE, inapplicables_as_missing = FALSE, polymorphism_behaviour = "equalp", uncertainty_behaviour = "equalp", threshold = 0.01, all_missing_allowed = FALSE) {
  
  # MAKE RATES OUTPUT MAKE SENSE BY PROELRY SEPARATING COMPLETENESS AND DURATION

  # ADD EXAMPLES OF VISUALISED OUTPUT (MENTION IN MANUAL AN ADD PLOT FUNCTIONS TO SEALSO_

  # TO EXCLUDE OUTGROUP CAN JUST SET UP PARTITIONS THAT EXCLUDE THESE (REQUIRES REMOVING PARTITION FIXER AS STANDARD)
  # COULD ALSO DO INTERNAL AND TERMINAL BRANCHES SEPARATELY THIS WAY
  # ALLOW NAMING OF PARTITIONS AND POTENTIALLY CALLING NAME FOR PLOTS

  # AICc BREAKS IF MORE THAN N-2 PARAMETERS (INFINITY OR WRONGLY NEGATIVE OUTPUT CAN OCCUR THIS WAY)
  # GLOBAL RATE NEEDS TO GO INTO LRT OPTION IF NOT USED IN AIC
  # SOMEHOW NEED TO ALLOW MULTIPLE VERSIONS OF MAIN PIEPLINE IF DOING RANDOM ASSIGNMENTS OF CHANGE TIMES
  # NEED TO ADD PARTITIONS USED TO OUTPUT SOMEHOW...
  # STOP REQUIRING TIME BINS TO BE SET IF NOT TESTING TIME PARTITIONS.

  # MAYBE DO SOME KIND OF WEIGHTING FOR AIC AS SOME PARTITIONS WILL CONTAIN VERY LITTLE DATA AND BE EXTREME OUTLIERS?
  # NEED TO CHECK FOR SINGLE PARTITION WITH LRT (ALLOWED WITH AIC)
  # NEED TO CHECK FOR FULLY BIFURCATING IF IMPLEMENTING WANG STUDENTS APPROACH? OR IS THAT DONE ANYWAY?
  # IF USING AIC NEED TO CHECK FOR EACH TEST TYPE AT LEAST TWO PARTITIONS ARE SUPPLIED OR ELSE THE RESULTS ARE MEANINGLESS
  # ADD CHECK TO INCLUDE ALL LESS COMPLEX SUBPARTITIONS OF ANY 3 OR MORE SIZE PARTITOINING OF THE DATA

  # DESIDERATA (STUFF IT WOULD BE NICE TO ADD IN FUTURE):
  #
  # WRITE SEARCH VERSION FOR FINDING RATE SHIFTS? SHOULD THIS EVEN BE AN OPTION? DOES THIS REQUIRE MODIFYING LRT TO COMPARE E.G. 2-RATE DIRECTLY WITH 3-RATE MODEL? WOULD NEED TO PERMUTE ALL POSSIBLE COMBOS AND NOT SURE HOW LARGE THESE MIGHT GET (VERY FOR EDGES).
  # MAYBE MAKE ANCESTRAL STATE UNCERTAINTY DIFFERENT FOR TIPS THAN NODES? I.E., HOW IT GETS RESOLVED CAN BE DIFFERENT (MORE OPTIONS TO FUNCTION)
  # THESE TWO ARE RELATED: 1. ADD TERMINAL VERSUS INTERNAL OPTION SOMEHOW/SOMEWHERE (BRANCH TYPE ALREADY RECORDED ON EDGE LIST!), 2. ALLOW OPTION TO IGNORE SOME PARTS OF THE TREE FOR PARTITION TESTS? MAKES CALCULATING THE MEAN RATE TRICKIER BUT MIGHT MAKE SENSE E.G. FOR INGROUP ONLY TESTS. EXCLUDE EDGES AFTER DOING ANCESTRAL STATES? OR SET THESE TO ALL NAS TO ENSURE OTHER THINGS WORK FINE? FOR EXAMPLE, USE OUTGROUPS TO SET ANCESTOR THEN EXCLUDE THEM FROM THE PROCESS.
  # EXTRA FUNCTION(S) TO VISUALISE RESULTS MOST LIKELY. DEFO NEEDED! HEAT MAP WITH EDGE BLOCKS?
  # CHECK FOR AUTAPOMORPHIES AND INFORM USER IF FOUND?
  # ADD CONTRIVED EXAMPLES (UNIT TESTS) TO SHOW HOW FUNCTION WORKS, E.G. RATE OF ONE CHANGES PER MILLION YEARS THEN DUPLICATED BLOCK WITH CHARACTER PARTITION TEST.
  # PROBABLY NEED MORE CAREFUL CHECKS FOR ZERO VALUES GENERALLY, E.G., CHARACTER WITH ALL MISSING DATA
  # ALLOW REWEIGHTING OF INAPPLICABLES ZERO AS AN OPTION FOR THEM?
  # HOW TO FORMAT OUTPUT? GET CIS FOR EACH PARTITION FOR VISUALISATION (E.G., BARPLOT OF PARTITION VALUES WITH DASHED LINE FOR MEAN AND ERROR BARS FOR CIS)? STICK WITH LIST OR COLLAPSE TO A MATRIX SOMEHOW?
  # TIME BINS WITH NOTHING IN WILL CAUSE ISSUES AS DIVIDE BY ZEROES WILL OCCUR - ADD CHECK FOR THIS.
  # WHAT IS SIGNIFICANTLY HIGH OR LOW IF THERE ARE THREE OR MORE PARTITIONS? THIS IS NOT EVEN IN OUTPUT YET. PROLLY CANNOT DO FULL STOP NOW PARTITIONS ARE MORE COMPLEX

  # Check cladistic_matrix has class cladisticMatrix and stop and warn user if not:
  if (!inherits(x = cladistic_matrix, what = "cladisticMatrix")) stop("cladistic_matrix must be an object of class \"cladisticMatrix\".")

  # Check for step matrices and stop and warn user if found:
  if (is.list(cladistic_matrix$topper$step_matrices)) stop("Function cannot currently deal with step matrices.")

  # Check tree has branch lengths:
  if (is.null(time_tree$edge.length)) stop("time_tree does not have branch lengths (durations). Try timescaling the tree, e.g., with DatePhylo.")

  # Check tree has root age:
  if (is.null(time_tree$root.time)) stop("time_tree is missing $root.time. Try setting this before continuing, e.g., tree$root.time <- 104.2.")

  # Check change_times is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = change_times, y = c("midpoint", "spaced", "random"))) > 0) stop("change_times must be one of \"midpoint\", \"spaced\", or \"random\".")

  # Check multiple_comparison_correction is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = multiple_comparison_correction, y = c("benjaminihochberg", "bonferroni"))) > 0) stop("multiple_comparison_correction must be one of \"benjaminihochberg\" or \"bonferroni\".")

  # Check polymorphism_state is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = polymorphism_state, y = c("missing", "random"))) > 0) stop("polymorphism_state must be one of \"missing\" or \"random\".")

  # Check uncertainty_state is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = uncertainty_state, y = c("missing", "random"))) > 0) stop("uncertainty_state must be one of \"missing\" or \"random\".")

  # Check inapplicable_state is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = inapplicable_state, y = c("missing"))) > 0) stop("inapplicable_state must be \"missing\".")

  # Check time_binning_approach is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = time_binning_approach, y = c("close", "lloyd"))) > 0) stop("time_binning_approach must be one of \"close\" or \"lloyd\".")

  # Check partitions are not all NULL values:
  if (is.null(branch_partitions) && is.null(character_partitions) && is.null(clade_partitions) && is.null(time_partitions)) stop("No partitions are requested. Set at least one of branch_partitions, character_partitions, clade_partitions, or time_partitions to a list of appropriate values. Type \"?test_rates\" for help.")

  # Get tip number:
  n_tips <- ape::Ntip(phy = time_tree)

  # Get node number:
  n_nodes <- ape::Nnode(phy = time_tree)

  # Get internal node numbers:
  node_numbers <- 1:n_nodes + n_tips

  # Get edge numbers:
  edge_numbers <- 1:nrow(time_tree$edge)

  # Get character numbers:
  character_numbers <- 1:sum(unlist(x = lapply(X = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "matrix"), ncol)))

  # Ensure time bins are in correct order:
  time_bins <- sort(x = unique(x = time_bins), decreasing = TRUE)

  # Find the Time bin midpoints:
  time_bin_midpoints <- (time_bins[2:length(x = time_bins)] + time_bins[1:(length(x = time_bins) - 1)]) / 2

  # Get the numbers for each time bins:
  time_bin_numbers <- 1:length(x = time_bin_midpoints)

  # Subfunction to ensure partitions are formatted correctly:
  format_partition <- function(partitions_to_test, valid_values, partition_name) {

    # Check partitions are in the form of a list of lists:
    if (!all(c(all(unlist(x = lapply(X = partitions_to_test, is.list))), is.list(partitions_to_test)))) stop(paste(partition_name, " must be in the form of a list of lists.", sep = ""))

    # Get a vector of any absent valid values:
    absent_values <- setdiff(x = unique(x = unlist(x = partitions_to_test)), y = valid_values)

    # Check valid values have been used and if not stop and warn user:
    if (length(x = absent_values) > 0) stop(paste(partition_name, "Partitions to test must be defined using the valid range of values (", paste(range(valid_values), collapse = " to "), ") only.", sep = ""))

    # Check partitions never overlap and stop and warn user if they do:
    overlap_check <- lapply(X = partitions_to_test, function(x) if (any(duplicated(sort(x = unlist(x = x))))) stop(paste("Each partition of ", partition_name, " must not contain overlapping values (e.g., can not have 1:3 and 3:5 as both contain 3).", sep = "")))

    # Subfunction to add the missing partition (if exists):
    add_missing_partitions <- function(x, valid_values) {

      # Define any missing values:
      missing_values <- setdiff(x = valid_values, y = unlist(x = x))

      # If there are missing values add them to list at end:
      if (length(x = missing_values) > 0) x[[(length(x = x) + 1)]] <- missing_values

      # Return x:
      x
    }

    # Add in missing partitions (if any):
    partitions_to_test <- lapply(X = partitions_to_test, add_missing_partitions, valid_values = valid_values)

    # Check partitions are all at least two in size or else no comparison can be made:
    if (any(unlist(x = lapply(X = partitions_to_test, length)) == 1) && test_type == "lrt") stop("Partitions must divide the available data into at least two parts if performing likelihood ratio tests.")

    # Return formatted partitions to test:
    return(partitions_to_test)
  }

  # Subfunction to pack partitions to short format for output:
  pack_partitions <- function(formatted_partitions) {
    unlist(x = lapply(X = formatted_partitions, function(x) {
      paste(unlist(x = lapply(X = x, function(y) {

        # First make sure y is sorted:
        y <- sort(x = y)

        # Covnvert y to a list (splitting if gaps greater than 1 are found):
        y <- unname(split(y, cumsum(c(TRUE, diff(y) > 1))))

        # Collapse gaps of one with hyphens:
        paste0(unlist(x = lapply(X = y, function(z) {
          res <- as.character(z)
          if (length(x = z) > 1) {
            r <- rle(c(1, pmin(diff(z), 2)))
            res <- paste0(z[c(1, cumsum(r$lengths))], c("-", " ")[r$values], collapse = "")
            res <- substr(res, 1, nchar(x = res) - 1)
          }
          res
        })), collapse = " ")
      })), collapse = " | ")
    }))
  }

  # If performing branch partition test(s) check and reformat branch partitions:
  if (!is.null(branch_partitions)) branch_partitions <- format_partition(partitions_to_test = branch_partitions, valid_values = edge_numbers, partition_name = "branch_partitions")

  # If performing character partition test(s) check and reformat character partitions:
  if (!is.null(character_partitions)) character_partitions <- format_partition(partitions_to_test = character_partitions, valid_values = character_numbers, partition_name = "character_partitions")

  # If performing clade partition test(s)
  if (!is.null(clade_partitions)) {

    # Convert clade partitions to edge partitions:
    clade_partitions <- lapply(X = clade_partitions, lapply, find_descendant_edges, tree = time_tree)

    # Check and reformat clade partitions:
    clade_partitions <- format_partition(partitions_to_test = clade_partitions, valid_values = edge_numbers, partition_name = "clade_partitions")
  }

  # If performing time bin partition test(s) check and reformat time bin partitions:
  if (!is.null(time_partitions)) time_partitions <- format_partition(partitions_to_test = time_partitions, valid_values = time_bin_numbers, partition_name = "time_partitions")

  # Check test_type is correctly formatted or stop and warn user:
  if (length(x = setdiff(x = test_type, y = c("aic", "lrt"))) > 0) stop("test_type must be one of \"aic\" or \"lrt\".")

  # Subfunction to calculate maximum likelihood p value:
  get_likelihood_p <- function(mean_rate, sampled_rates, sampled_changes, sampled_completeness, sampled_time) {

    # Get log numerator:
    log_numerator <- sum(log(dpois(round(sampled_changes), mean_rate * sampled_completeness * sampled_time)))

    # Get log denominator:
    log_denominator <- sum(log(dpois(round(sampled_changes), sampled_rates * sampled_completeness * sampled_time)))

    # Get test statistic:
    test_statistic <- -2 * (log_numerator - log_denominator)

    # Calculate position of test statistic in chi-square distribution to get probability and return:
    pchisq(test_statistic, length(x = sampled_rates) - 1, lower.tail = FALSE)
  }

  # Subfunction to calculate AIC:
  calculate_aic <- function(sampled_rates, sampled_changes, sampled_completeness, sampled_time) {

    # Get log maximum likelihood estimate:
    log_mle <- sum(log(dpois(x = round(sampled_changes), lambda = sampled_rates * sampled_completeness * sampled_time)))

    # Calculate and return AIC:
    (2 * length(x = sampled_rates)) - (2 * log_mle)
  }

  # Subfunction to calculate AIC from partition (with columns labelled partition, rate, completeness, duration):
  calculate_partition_aic <- function(partition, aicc = FALSE) {

    # Get log maximum likelihood estimate:
    log_mle <- sum(log(dpois(round(partition[, "changes"]), partition[, "rate"] * partition[, "completeness"] * partition[, "duration"])))

    # Get k (number of parameters) term:
    k <- max(partition[, "partition"])

    # Calculate AIC:
    aic <- (2 * k) - (2 * log_mle)

    # If AICc is desired then calculate this and overwrite AIC with it:
    if (aicc) aic <- aic + (((2 * (k^2)) + (2 * k)) / (nrow(partition) - k - 1))

    # Return AIC:
    aic
  }

  # Get ages for each (tip and internal) node:
  node_dates <- date_nodes(time_tree = time_tree)

  # Get branch ages (from and to):
  branch_ages <- unname(cbind(node_dates[as.character(time_tree$edge[, 1])], node_dates[as.character(time_tree$edge[, 2])]))

  # Build edge list from node numbers (from-to) for each branch:
  edge_list <- lapply(X = apply(time_tree$edge, 1, list), function(x) {
    names(x) <- "node_number_from_to"
    x
  })

  # Add node ages to edge list:
  for (i in 1:length(x = edge_list)) edge_list[[i]]$node_age_from_to <- branch_ages[i, ]

  # Add node ages (from-to) to each edge in list:
  edge_list <- lapply(X = edge_list, function(x) {
    x$branch_duration <- x$node_age_from_to[1] - x$node_age_from_to[2]
    x
  })

  # Get vector of branch types:
  branch_types <- gsub(pattern = "0", replacement = "internal", x = gsub(pattern = "1", replacement = "terminal", x = as.numeric(time_tree$edge[, 2] <= n_tips)))

  # Add branch type to edge list:
  for (i in 1:length(x = edge_list)) edge_list[[i]]$branch_type <- branch_types[i]

  # Find descendant edges for each internal node:
  find_descendant_edges_for_each_internal_node <- lapply(X = as.list(x = node_numbers), find_descendant_edges, tree = time_tree)

  # Get ancestral character states:
  ancestral_states <- estimate_ancestral_states(cladistic_matrix = cladistic_matrix, time_tree = time_tree, estimate_all_nodes = estimate_all_nodes, estimate_tip_values = estimate_tip_values, inapplicables_as_missing = inapplicables_as_missing, polymorphism_behaviour = polymorphism_behaviour, uncertainty_behaviour = uncertainty_behaviour, threshold = threshold, all_missing_allowed = all_missing_allowed)

  # Build single matrix of all states in tip label then node number order:
  all_states <- do.call(what = cbind, args = lapply(X = lapply(X = ancestral_states[2:length(x = ancestral_states)], "[[", "matrix"), function(x) x[c(time_tree$tip.label, 1:n_nodes + n_tips), , drop = FALSE]))

  # Make vector of ordering of characters:
  ordering <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "ordering")))

  # Make vector of weights of characters:
  character_weights <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "character_weights")))

  # Make vector of minimum values:
  minimum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "minimum_values")))

  # Make vector of maximum values:
  maximum_values <- unname(unlist(x = lapply(X = cladistic_matrix[2:length(x = cladistic_matrix)], "[[", "maximum_values")))

  # Find positions in matrix with polymorphisms:
  polymorphism_positions <- grep("&", all_states)

  # Find positions in matrix with uncertainties:
  uncertainty_positions <- grep("/", all_states)

  # Find positions in matrix with inapplicables:
  inapplicable_positions <- which(x = all_states == "")

  # If polymorphisms were found:
  if (length(x = polymorphism_positions) > 0) {

    # If replacing polymorphsims with missing do so:
    if (polymorphism_state == "missing") all_states[polymorphism_positions] <- NA

    # If replacing polymorphisms with random values draw and replace:
    if (polymorphism_state == "random") all_states[polymorphism_positions] <- unlist(x = lapply(X = strsplit(all_states[polymorphism_positions], "&"), sample, size = 1))
  }

  # If uncertainties were found:
  if (length(x = uncertainty_positions) > 0) {

    # If replacing uncertainties with missing do so:
    if (uncertainty_state == "missing") all_states[uncertainty_positions] <- NA

    # If replacing uncertainties with random values draw and replace:
    if (uncertainty_state == "random") all_states[uncertainty_positions] <- unlist(x = lapply(X = strsplit(all_states[uncertainty_positions], "/"), sample, size = 1))
  }

  # If inapplicable states were found:
  if (length(x = inapplicable_positions) > 0) {

    # If replacing inapplicables with missing do so:
    if (inapplicable_state == "missing") all_states[inapplicable_positions] <- NA
  }

  # Set default converted continuous characters to FALSE:
  continuous_characters_discretized <- FALSE

  # Check for continuous characters as these will need to be modified for modelling to discrete characters:
  if (any(ordering == "cont")) {

    # Tell user this is happening:
    cat("Continuous characters found. Converting to gap-weighted discrete characters.\n")

    # Set default converted continuous characters to TRUE:
    continuous_characters_discretized <- TRUE

    # Find out which characters are continuous:
    continuous_characters_found <- which(x = ordering == "cont")

    # Rescale continous characters as zero to one values:
    list_of_continuous_values_rescaled_zero_to_one <- lapply(X = lapply(X = lapply(X = apply(all_states[, continuous_characters_found, drop = FALSE], 2, list), unlist), as.numeric), function(x) {
      x <- x - min(sort(x = x))
      x <- x / max(sort(x = x))
      return(x)
    })

    # Now discretize and store these characters (0 to 31 scale):
    all_states[, continuous_characters_found] <- do.call(what = cbind, args = lapply(X = lapply(X = lapply(X = list_of_continuous_values_rescaled_zero_to_one, function(x) as.list(x = x)), lapply, function(x) ifelse(is.na(x), NA, max(which(x = x >= (0:31) / 31)) - 1)), unlist))

    # Convert character type to ordered:
    ordering[continuous_characters_found] <- "ord"

    # Convert weights to 1/31:
    character_weights[continuous_characters_found] <- 1 / 31

    # Set minimum value to zero:
    minimum_values[continuous_characters_found] <- 0

    # Set maximum value to 31:
    maximum_values[continuous_characters_found] <- 31
  }

  # If all_weights_integers is TRUE rescale weights until they are all integers so can model appropriately with Poisson later:
  if (all_weights_integers) while (is.character(all.equal(sum(character_weights %% 1), 0))) character_weights <- (1 / (character_weights %% 1)[(character_weights %% 1) > 0])[1] * character_weights

  # Add from-to node states for each character to edge list:
  edge_list <- lapply(X = edge_list, function(x) {
    x$character_states_from_to <- matrix(all_states[x$node_number_from_to, , drop = FALSE], nrow = 2, dimnames = list(c("from", "to")))
    x
  })

  # Subfunction to define character changes:
  build_changes_matrix <- function(x) {

    # Find only comparable characters (those scored for both from and to states):
    comparable_characters <- which(x = apply(!apply(x$character_states_from_to, 2, is.na), 2, all))

    # Isolate comparable ordering:
    comparable_ordering <- ordering[comparable_characters]

    # Isolate comparable weights:
    comparable_weights <- character_weights[comparable_characters]

    # Isolate only characters that actually differ (change):
    character_differences <- which(x = x$character_states_from_to["from", comparable_characters] != x$character_states_from_to["to", comparable_characters])

    # Build character change matrix:
    character_changes <- matrix(nrow = 0, ncol = 5, dimnames = list(c(), c("character", "from", "to", "steps", "weight")))

    # If characters change then make a matrix from them:
    if (length(x = character_differences) > 0) character_changes <- rbind(character_changes, cbind(as.numeric(comparable_characters[character_differences]), as.numeric(x$character_states_from_to["from", comparable_characters[character_differences]]), as.numeric(x$character_states_from_to["to", comparable_characters[character_differences]]), ifelse(comparable_ordering[character_differences] == "unord", 1, abs(as.numeric(x$character_states_from_to["to", comparable_characters[character_differences]]) - as.numeric(x$character_states_from_to["from", comparable_characters[character_differences]]))), comparable_weights[character_differences]))

    # Store character changes as new sublist for x:
    x$character_changes <- character_changes

    # Store comparable characters as new sublist of x:
    x$comparable_characters <- comparable_characters

    # Return x:
    x
  }

  # Get character changes and comparable characters and add to edge list:
  edge_list <- lapply(X = edge_list, build_changes_matrix)

  # Check whether time bins are being compared (otherwise no need to assign character changes):
  if (!is.null(time_partitions)) {

    # Subfunction to add change times to character changes:
    add_change_times <- function(x, change_times) {

      # Isolate character changes:
      character_changes <- x$character_changes

      # If any changes involve two or more steps (requiring replacement with multiple changes):
      if (any(character_changes[, "steps"] > 1)) {

        # Get multistep character changes:
        multistep_characters <- which(x = character_changes[, "steps"] > 1)

        # For each multistep character change:
        for (i in rev(multistep_characters)) {

          # Isolate other rows:
          other_row_numbers <- setdiff(x = 1:nrow(character_changes), y = i)

          # Get unpacked changes (X:Y, e.g., 0:2 would become 0 1 2):
          unpacked_changes <- character_changes[i, "from"]:character_changes[i, "to"]

          # Update character changes with multistep changes unpacked:
          character_changes <- rbind(character_changes[other_row_numbers, ], unname(cbind(rep(character_changes[i, "character"], length.out = length(x = unpacked_changes) - 1), unpacked_changes[1:(length(x = unpacked_changes) - 1)], unpacked_changes[2:length(x = unpacked_changes)], rep(1, length.out = length(x = unpacked_changes) - 1), rep(character_changes[i, "weight"], length.out = length(x = unpacked_changes) - 1))))
        }

        # Resort character changes by character number:
        character_changes <- character_changes[order(character_changes[, "character"]), ]
      }

      # If using midpoint option set character change times as midpoint of branch:
      if (change_times == "midpoint") character_changes <- cbind(character_changes, rep(x$node_age_from_to[1] - (x$branch_duration / 2), length.out = nrow(character_changes)))

      # If using spaced then set character change times as equally spaced along branch:
      if (change_times == "spaced") character_changes <- cbind(character_changes, x$node_age_from_to[1] - (seq(from = 0, to = x$branch_duration, length.out = nrow(character_changes) + 1)[1:nrow(character_changes)] + (diff(seq(from = 0, to = x$branch_duration, length.out = nrow(character_changes) + 1))[1] / 2)))

      # If using random then set character change times as random draws from a uniform distribution:
      if (change_times == "random") character_changes <- cbind(character_changes, x$node_age_from_to[1] - stats::runif(n = nrow(character_changes), min = 0, max = x$branch_duration))

      # Add column name to change time column:
      colnames(x = character_changes)[ncol(character_changes)] <- "time"

      # Subfunction to re-sort character change times so they occur in correct order:
      sort_change_times <- function(character_changes) {

        # Sort change time for each character from oldest (first) to youngest (last) and store it:
        character_changes[, "time"] <- unname(unlist(x = lapply(X = as.list(x = unique(x = character_changes[, "character"])), function(x) sort(x = character_changes[which(x = character_changes[, "character"] == x), "time"], decreasing = TRUE))))

        # Return sorted character changes:
        character_changes
      }

      # Re-sort character change times so they occur in correct order:
      character_changes <- sort_change_times(character_changes)

      # Add bin for character change as last column:
      character_changes <- cbind(character_changes, unlist(x = lapply(X = as.list(x = character_changes[, "time"]), function(x) max(which(x = x <= time_bins)))))

      # Add column name to change time column:
      colnames(x = character_changes)[ncol(character_changes)] <- "bin"

      # Overwrite character changes with new version with changes added:
      x$character_changes <- character_changes

      # Return x:
      x
    }

    # Add character change times to edge list:
    edge_list <- lapply(X = edge_list, add_change_times, change_times = change_times)
  }

  # Subfunction to get edge sections in time bins:
  get_edge_sections_in_bins <- function(x, time_bins = time_bins) {

    # Set first appearance datum of edge:
    fad <- x$node_age_from_to[1]

    # Set last appearance datum of edge:
    lad <- x$node_age_from_to[2]

    # Get any time bin boundaries crossed (can be empty if none are):
    boundaries_crossed <- time_bins[2:(length(x = time_bins) - 1)][intersect(which(x = time_bins[2:(length(x = time_bins) - 1)] > lad), which(x = time_bins[2:(length(x = time_bins) - 1)] < fad))]

    # If boundaries are crossed:
    if (length(x = boundaries_crossed) > 0) {

      # Break up branch into binned sections as vector of FADs:
      fad <- c(fad, boundaries_crossed)

      # Break up branch into binned sections as vector of LADs:
      lad <- c(boundaries_crossed, lad)
    }

    # Build matrix of branch sections with FADs and LADs:
    branch_sections <- rbind(fad, lad)

    # Add bin number present in to column names:
    colnames(x = branch_sections) <- unlist(x = lapply(X = lapply(X = lapply(X = as.list(x = branch_sections["fad", ]), "<=", time_bins), which), max))

    # Add new list section for branch (edge) sections binned by time:
    x$binned_edge_sections <- branch_sections

    # Return output:
    x
  }

  # Get edge sections in time bins:
  edge_list <- lapply(X = edge_list, get_edge_sections_in_bins, time_bins = time_bins)

  # Add binned branch durations to edge list:
  edge_list <- lapply(X = edge_list, function(x) {
    branch_durations <- rep(0, length(x = time_bins) - 1)
    branch_durations[as.numeric(colnames(x = x$binned_edge_sections))] <- abs(apply(x$binned_edge_sections, 2, diff))
    x$binned_branch_durations <- branch_durations
    x
  })

  # Add proportional binned branch lengths to edge list:
  edge_list <- lapply(X = edge_list, function(x) {
    x$proportional_binned_edge_durations <- x$binned_branch_durations / sum(x$binned_branch_durations)
    x
  })

  # Start to build matrix of all changes with list of character changes:
  inferred_character_changes <- lapply(X = edge_list, function(x) x$character_changes)

  # Add edge number to each matrix of character changes:
  for (i in 1:length(x = inferred_character_changes)) inferred_character_changes[[i]] <- cbind(rep(i, times = nrow(inferred_character_changes[[i]])), inferred_character_changes[[i]])

  # Combine all changes into a single matrix:
  inferred_character_changes <- do.call(what = rbind, args = lapply(X = inferred_character_changes, function(x) {
    colnames(x = x)[1] <- "edge"
    x
  }))

  # Remove silly rownames from all changes:
  rownames(x = inferred_character_changes) <- NULL

  # Create NULL output variables (to be overwritten if called):
  branch_rates <- clade_rates <- time_rates <- character_rates <- NULL

  # If doing some kind of edge test (branch or clade):
  if (!is.null(branch_partitions) || !is.null(clade_partitions)) {

    # Get (weighted) number of changes on each edge:
    edge_changes <- unlist(x = lapply(X = edge_list, function(x) sum(x$character_changes[, "steps"] * x$character_changes[, "weight"])))

    # Get completeness for each edge:
    edge_completeness <- unlist(x = lapply(X = edge_list, function(x) sum(character_weights[x$comparable_characters]) / sum(character_weights)))

    # Get duration of each edge:
    edge_durations <- unlist(x = lapply(X = edge_list, function(x) x$branch_duration))

    # Set global rate:
    global_rate <- sum(edge_changes) / sum(edge_completeness * edge_durations)

    # If performing branch partition tests:
    if (!is.null(branch_partitions)) {

      # Create branch rates for output:
      branch_rates <- lapply(X = list(as.list(x = 1:nrow(time_tree$edge))), function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(edge_changes[y]), sum(edge_completeness[y] * edge_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))[[1]]

      # Add edge numbers and rates:
      branch_rates <- cbind(edge = 1:nrow(time_tree$edge), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = branch_rates[, "changes"] / branch_rates[, "completeness"])), branch_rates)

      # If using Likelihood Ratio Test:
      if (test_type == "lrt") {

        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        partitioned_data <- lapply(X = branch_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(edge_changes[y]), sum(edge_completeness[y] * edge_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))

        # Add sampled rate to paritioned data matrices:
        partitioned_data <- lapply(X = partitioned_data, function(x) {
          cbind(rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        })

        # Get LRT p-values and combine output as edge test results:
        branch_test_results <- lapply(X = partitioned_data, function(x) {
          list(rates = x[, "rate"], p_value = get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        })
      }

      # If using AIC:
      if (test_type == "aic") {

        # Build partitioned data for AIC:
        partitioned_data <- lapply(X = branch_partitions, function(x) {
          y <- cbind(partition = rep(NA, length(x = time_tree$edge.length)), rate = rep(NA, length(x = time_tree$edge.length)), changes = edge_changes, completeness = edge_completeness, duration = edge_durations)
          y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(z) rep(sum(y[z, "changes"]) / (sum(y[z, "completeness"] * y[z, "duration"])), length(x = z))))[order(unlist(x = x))]))
          y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))[order(unlist(x = x))]
          y
        })

        # Get AIC, AICc and rate results:
        branch_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
      }

      # Pack branch partitions to test into single strings for output:
      packed_branch_partitions <- pack_partitions(branch_partitions)

      # Add packed partitions to results:
      for (i in 1:length(x = branch_test_results)) branch_test_results[[i]]$partition <- packed_branch_partitions[i]

      # If not performing branch partition tests:
    } else {

      # Make empty branch partition result output:
      branch_test_results <- NULL
    }

    # If performing clade partition tests:
    if (!is.null(clade_partitions)) {

      # Create clade rates for output:
      clade_rates <- do.call(what = rbind, args = lapply(X = lapply(X = as.list(x = n_tips + c(1:time_tree$Nnode)), function(x) find_descendant_edges(x, tree = time_tree)), function(y) c(changes = sum(edge_changes[y]), completeness = sum(edge_completeness[y] * edge_durations[y]), duration = 1)))

      # Add rates and node numbers:
      clade_rates <- cbind(node = n_tips + c(1:time_tree$Nnode), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = clade_rates[, "changes"] / clade_rates[, "completeness"])), clade_rates)

      # If using Likelihood Ratio Test:
      if (test_type == "lrt") {

        # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
        partitioned_data <- lapply(X = clade_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(edge_changes[y]), sum(edge_completeness[y] * edge_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))

        # Add sampled rate to paritioned data matrices:
        partitioned_data <- lapply(X = partitioned_data, function(x) {
          cbind(rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        })

        # Get LRT p-values and combine output as edge test results:
        clade_test_results <- lapply(X = partitioned_data, function(x) {
          list(rates = x[, "rate"], p_value = get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        })
      }

      # If using AIC:
      if (test_type == "aic") {

        # Build partitioned data for AIC:
        partitioned_data <- lapply(X = clade_partitions, function(x) {
          y <- cbind(partition = rep(NA, length(x = time_tree$edge.length)), rate = rep(NA, length(x = time_tree$edge.length)), changes = edge_changes, completeness = edge_completeness, duration = edge_durations)
          y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(x) rep(sum(y[x, "changes"]) / (sum(y[x, "completeness"] * y[x, "duration"])), length(x = x))))[order(unlist(x = x))]))
          y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))[order(unlist(x = x))]
          y
        })

        # Get AIC, AICc and rate results:
        clade_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
      }

      # Pack clade partitions to test into single strings for output:
      packed_clade_partitions <- pack_partitions(clade_partitions)

      # Add packed partitions to results:
      for (i in 1:length(x = clade_test_results)) clade_test_results[[i]]$partition <- packed_clade_partitions[i]


      # If not performing clade partition tests:
    } else {

      # Make empty clade partition result output:
      clade_test_results <- NULL
    }

    # If not doing clade OR branch tests:
  } else {

    # Make empty branch partition result output:
    branch_test_results <- NULL

    # Make empty clade partition result output:
    clade_test_results <- NULL
  }

  # If performing character partition tests:
  if (!is.null(character_partitions)) {

    # Get vector of (weighted) changes for each character:
    character_changes <- unlist(x = lapply(X = as.list(x = character_numbers), function(x) {
      character_rows <- which(x = inferred_character_changes[, "character"] == x)
      sum(inferred_character_changes[character_rows, "steps"] * inferred_character_changes[character_rows, "weight"])
    }))

    # Get vector of weighted durations for each character:
    character_durations <- (character_weights / sum(character_weights)) * sum(time_tree$edge.length)

    # Get vector of completness (opportunity to observe changes) for each character:
    character_completeness <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) {
      character_presence <- rep(0, times = length(x = character_numbers))
      character_presence[x$comparable_characters] <- 1
      character_presence * x$branch_duration
    })), 2, sum) / sum(time_tree$edge.length)

    # Set global rate:
    global_rate <- sum(character_changes) / sum(character_completeness * character_durations)

    # Create character rates for output:
    character_rates <- lapply(X = list(as.list(x = 1:max(character_numbers))), function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(character_changes[y]), sum(character_completeness[y] * character_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))[[1]]

    # Add character numbers and rates:
    character_rates <- cbind(character = 1:max(character_numbers), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = character_rates[, "changes"] / character_rates[, "completeness"])), character_rates)

    # If using likelihood ratio test:
    if (test_type == "lrt") {

      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later):
      partitioned_data <- lapply(X = character_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(character_changes[y]), sum(character_completeness[y] * character_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))

      # Add sampled rate to paritioned data matrices:
      partitioned_data <- lapply(X = partitioned_data, function(x) {
        x <- cbind(as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        colnames(x = x)[1] <- "rate"
        x
      })

      # Get P-Values and combine output as edge test results:
      character_test_results <- lapply(X = partitioned_data, function(x) {
        x <- list(x[, "rate"], get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        names(x) <- c("rates", "p_value")
        x
      })
    }

    # If using AIC:
    if (test_type == "aic") {

      # Build partitioned data for AIC:
      partitioned_data <- lapply(X = character_partitions, function(x) {
        y <- cbind(partition = rep(NA, max(character_numbers)), rate = rep(NA, max(character_numbers)), changes = character_changes, completeness = character_completeness, duration = character_durations)
        y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(x) rep(sum(y[x, "changes"]) / (sum(y[x, "completeness"] * y[x, "duration"])), length(x = x))))[order(unlist(x = x))]))
        y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))[order(unlist(x = x))]
        y
      })

      # Get AIC, AICc and rate results:
      character_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
    }

    # Pack character partitions to test into single strings for output:
    packed_character_partitions <- pack_partitions(character_partitions)

    # Add packed partitions to results:
    for (i in 1:length(x = character_test_results)) character_test_results[[i]]$partition <- packed_character_partitions[i]

    # If performing branch partition tests:
  } else {

    # Make empty character partition result output:
    character_test_results <- NULL
  }

  # If performing time bin partition tests:
  if (!is.null(time_partitions)) {

    # Get weighted number of changes from each time bin:
    time_bin_changes <- unlist(x = lapply(X = as.list(x = 1:(length(x = time_bins) - 1)), function(x) {
      change_rows <- inferred_character_changes[, "bin"] == x
      sum(inferred_character_changes[change_rows, "steps"] * inferred_character_changes[change_rows, "weight"])
    }))

    # If using the Close time bin completeness approach get completeness value for each time bin:
    if (time_binning_approach == "close") time_bin_completeness <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) x$proportional_binned_edge_durations * (sum(character_weights[x$comparable_characters]) / sum(character_weights)))), 2, sum) / apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) x$proportional_binned_edge_durations)), 2, sum)

    # If using the Lloyd time bin completeness approach get completeness value for each time bin::
    if (time_binning_approach == "lloyd") time_bin_completeness <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) apply(matrix(character_weights[x$comparable_characters], ncol = 1) %*% x$binned_branch_durations, 2, sum))), 2, sum) / apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) apply(matrix(character_weights, ncol = 1) %*% x$binned_branch_durations, 2, sum))), 2, sum)

    # Get durations of edges in each time bin:
    time_bin_durations <- apply(do.call(what = rbind, args = lapply(X = edge_list, function(x) x$binned_branch_durations)), 2, sum)

    # Set global rate (NB: will differ between Close and Lloyd approaches, but Lloyd approach will match edge or character global rate):
    global_rate <- sum(time_bin_changes) / sum(time_bin_completeness * time_bin_durations)

    # Create time rates for output:
    time_rates <- lapply(X = list(as.list(x = 1:(length(x = time_bins) - 1))), function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(time_bin_changes[y]), sum(time_bin_completeness[y] * time_bin_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))[[1]]

    # Add bin numbers and rates:
    time_rates <- cbind(bin = 1:(length(x = time_bins) - 1), rate = as.numeric(gsub(pattern = NaN, replacement = 0, x = time_rates[, "changes"] / time_rates[, "completeness"])), time_rates)

    # If using Likelihood Ratio Test to compare partitions:
    if (test_type == "lrt") {

      # Build partitioned data matrices (NB: completeness and duration get combined here as they cannot be summed separately later) for LRT:
      partitioned_data <- lapply(X = time_partitions, function(x) matrix(unlist(x = lapply(X = x, function(y) c(sum(time_bin_changes[y]), sum(time_bin_completeness[y] * time_bin_durations[y]), 1))), ncol = 3, byrow = TRUE, dimnames = list(c(), c("changes", "completeness", "duration"))))

      # Add sampled rate to paritioned data matrices for LRT:
      partitioned_data <- lapply(X = partitioned_data, function(x) {
        x <- cbind(as.numeric(gsub(pattern = NaN, replacement = 0, x = c(x[, "changes"] / (x[, "completeness"] * x[, "duration"])))), x)
        colnames(x = x)[1] <- "rate"
        x
      })

      # Get P-Values and combine output as edge test results:
      time_test_results <- lapply(X = partitioned_data, function(x) {
        x <- list(x[, "rate"], get_likelihood_p(mean_rate = global_rate, sampled_rates = x[, "rate"], sampled_changes = x[, "changes"], sampled_completeness = x[, "completeness"], sampled_time = x[, "duration"]))
        names(x) <- c("rates", "p_value")
        x
      })
    }

    # If using AIC to compare partitions:
    if (test_type == "aic") {

      # Build partitioned data for AIC:
      partitioned_data <- lapply(X = time_partitions, function(x) {
        y <- cbind(partition = rep(NA, length(x = time_bin_changes)), rate = rep(NA, length(x = time_bin_changes)), changes = time_bin_changes, completeness = time_bin_completeness * time_bin_durations, duration = rep(1, length(x = time_bin_changes)))
        y[, "rate"] <- as.numeric(gsub(pattern = NaN, replacement = 0, x = unlist(x = lapply(X = x, function(x) rep(sum(y[x, "changes"]) / sum(y[x, "completeness"]), length(x = x))))))
        y[, "partition"] <- rep(1:length(x = x), unlist(x = lapply(X = x, length)))
        y
      })

      # Get AIC, AICc and rate results:
      time_test_results <- lapply(X = partitioned_data, function(x) list(rates = unname(unlist(x = lapply(X = as.list(x = unique(x = sort(x = x[, "partition"]))), function(y) x[x[, "partition"] == y, "rate"][1]))), aic = calculate_partition_aic(x), aicc = calculate_partition_aic(x, aicc = TRUE)))
    }

    # Pack time bin partitions to test into single strings for output:
    packed_time_bin_partitions <- pack_partitions(time_partitions)

    # Add packed partitions to results:
    for (i in 1:length(x = time_test_results)) time_test_results[[i]]$partition <- packed_time_bin_partitions[i]

    # If not performing time bin partition tests:
  } else {

    # Make empty time bin partition result output:
    time_test_results <- NULL
  }

  # Set global rate for output:
  global_rate <- sum(unlist(x = lapply(X = edge_list, function(x) sum(x$character_changes[, "steps"] * x$character_changes[, "weight"])))) / sum(unlist(x = lapply(X = edge_list, function(x) sum(character_weights[x$comparable_characters]) / sum(character_weights))) * unlist(x = lapply(X = edge_list, function(x) x$branch_duration)))

  # If performing Likelihood Ratio Test:
  if (test_type == "lrt") {

    # Subfunction to calculate adjusted alphas for multiple comparison corrections:
    add_cutoffs <- function(test_results, alpha, multiple_comparison_correction = multiple_comparison_correction) {

      # Get number of comparisons performed:
      n_comparisons <- length(x = test_results)

      # If using the Benjamini-Hochberg false discovery rate approach:
      if (multiple_comparison_correction == "benjaminihochberg") {

        # Set cutoff values:
        cutoff_values <- ((1:n_comparisons) / n_comparisons) * alpha

        # Get actual p-values found:
        p_values <- unlist(x = lapply(X = test_results, "[[", "p_value"))

        # Order cutoffs by p-value rank:
        cutoff_values <- cutoff_values[rank(p_values, ties.method = "random")]
      }

      # If using the Bonferroni correction set cutoff values as alpha over N:
      if (multiple_comparison_correction == "bonferroni") cutoff_values <- alpha / n_comparisons

      # Add cutoffs to output:
      for (i in 1:length(x = test_results)) test_results[[i]]$CorrectedAlpha <- cutoff_values[i]

      # Return modified test results:
      test_results
    }

    # If doing branch partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(branch_partitions)) branch_test_results <- add_cutoffs(test_results = branch_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)

    # If doing character partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(character_partitions)) character_test_results <- add_cutoffs(test_results = character_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)

    # If doing clade partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(clade_partitions)) clade_test_results <- add_cutoffs(test_results = clade_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)

    # If doing time bin partition tests then add multiple comparison alpha cutoffs:
    if (!is.null(time_partitions)) time_test_results <- add_cutoffs(test_results = time_test_results, alpha = alpha, multiple_comparison_correction = multiple_comparison_correction)
  }

  # Return compiled output:
  list(time_bins_used = time_bins, inferred_character_changes = inferred_character_changes, mean_rate = global_rate, continuous_characters_discretized = continuous_characters_discretized, branch_test_results = branch_test_results, character_test_results = character_test_results, clade_test_results = clade_test_results, time_test_results = time_test_results, branch_rates = branch_rates, character_rates = character_rates, clade_rates = clade_rates, time_rates = time_rates, time_tree = time_tree)
}

# Ages <- read.table("~/Documents/Packages/Claddis/LungfishTest/ages.txt", sep = ",")
# Matrix <- Claddis::read_nexus_matrix("~/Documents/Packages/Claddis/LungfishTest/Lloyd_etal_2012a.nex")
# time_tree <- ape::read.tree("~/Documents/Packages/Claddis/LungfishTest/Lloyd_etal_2012a.tre")
# time_bins <- c(443.8, 419.2, 358.9, 298.9, 251.9, 201.3, 145.0, 0.0)
#
# time_tree <- time_tree[sample(1:100000, 100)]
# time_tree <- lapply(X = time_tree, function(x) strap::DatePhylo(x, Ages, rlen = 2, method = "equal"))
# class(time_tree) <- "multiPhylo"
#
# time_tree <- time_tree[[1]]
# cladistic_matrix <- Matrix
# time_bins <- time_bins
# branch_partitions <- lapply(X = as.list(x = 1:nrow(time_tree$edge)), as.list)
# character_partitions <- list(list(1:91), list(Cranial = 1:81, Postcranial = 82:91))
# clade_partitions <- lapply(X = as.list(x = ape::Ntip(phy = time_tree) + (2:ape::Nnode(phy = time_tree))), as.list)
# time_partitions <- partition_time_bins(7)
# change_times = "random"
# test_type = "aic"
# alpha = 0.01
# multiple_comparison_correction = "benjaminihochberg"
# polymorphism_state = "missing"
# uncertainty_state = "missing"
# inapplicable_state = "missing"
# time_binning_approach = "lloyd"
# all_weights_integers = FALSE
# estimate_all_nodes = FALSE
# estimate_tip_values = FALSE
# inapplicables_as_missing = FALSE
# polymorphism_behaviour = "equalp"
# uncertainty_behaviour = "equalp"
# threshold = 0.01
# all_missing_allowed = FALSE
