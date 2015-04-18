DiscreteCharacterRate <- function(tree, clad.matrix, time.bins, alpha=0.01) {
  
  # Ensure time bins are in correct order:
  time.bins <- sort(time.bins, decreasing=TRUE)

  # Find the Time bin midpoints:
  Time.bin.midpoints <- (time.bins[2:length(time.bins)] + time.bins[1:(length(time.bins)-1)])/2
  
  # Get ages for each (tip and internal) node:
  node.ages <- GetNodeAges(tree)

  # Get branch ages (from and to):
  branch.ages <- cbind(node.ages[as.character(tree$edge[, 1])], node.ages[as.character(tree$edge[, 2])])

  # Remove rownmaes:
  rownames(branch.ages) <- NULL

  # Matrix to store proportion of each branch in a bin:
  prop.branch.in.bin <- matrix(0, nrow=nrow(tree$edge), ncol=length(time.bins) - 1)

  # Get branch durations:
  branch.durations <- abs(apply(branch.ages, 1, diff))

  # For each time bin:
  for(i in 2:length(time.bins)) {

    # Get first appearance dates:
    LADs <- branch.ages[, 2]
    
    # Get last appearance dates:
    FADs <- branch.ages[, 1]

    # Make last appearances after the bin the same age as the top of the bin:
    LADs[LADs < time.bins[i]] <- time.bins[i]
    
    # Make first appearances before the bin the same age as the bottom of the bin:
    FADs[FADs > time.bins[(i - 1)]] <- time.bins[(i - 1)]
    
    # Ensure difference in ages between branches not present in bin are zero:
    LADs[FADs < LADs] <- FADs[FADs < LADs]
    
    # Calculate and store the proportion of each branch ipresent in the time bin:
    prop.branch.in.bin[, (i - 1)] <- (FADs - LADs) / branch.durations
    
  }

  # Get ancestral character states:
  anc.states <- AncStateEstMatrix(clad.matrix, tree, estimate.allchars=FALSE, estimate.tips=FALSE)
  
  # Isolate tip states in tip label order:
  tip.states <- clad.matrix$matrix[tree$tip.label,]
  
  # Build matrix of all states in node number order:
  all.states <- rbind(tip.states, anc.states)
  
  # Re-label row names by node number:
  rownames(all.states) <- c(1:(Ntip(tree)+Nnode(tree)))
  
  # Create trees to store branch lengths with respect to completeness and observed and observable character changes:
  completeness.tree <- observed.tree <- observable.tree <- tree
  
  # For each branch of the phylogenetic tree:
  for(i in 1:length(tree$edge[, 1])) {
    
    # Find the characters scored for both nodes (either end of the branch):
    comp.chars <- intersect(grep(TRUE, !is.na(all.states[tree$edge[i, 1], ])), grep(TRUE, !is.na(all.states[tree$edge[i, 2], ])))
    
    # Store completeness of branch:
    completeness.tree$edge.length[i] <- sum(clad.matrix$weights[comp.chars])
    
    # List the states for both nodes for comparable characters only:
    comp.states <- all.states[tree$edge[i,], comp.chars]
    
    # Special case if there is only one character (need to re-define as matrix):
    if(length(comp.states) == 2) comp.states <- as.matrix(comp.states)
    
    # Get weightings for comparable characters only:
    weightings <- clad.matrix$weights[comp.chars]
    
    # Get orderings for comparable characters only:
    orderings <- clad.matrix$ordering[comp.chars]
    
    # Store ranges (difference between minima and maxima):
    ranges <- (clad.matrix$max.vals - clad.matrix$min.vals)[comp.chars]
    
    # Correct ranges for unordered characters:
    ranges[intersect(grep(TRUE, ranges > 1), grep(TRUE, orderings == "UNORD"))] <- 1
    
    # Find polymorphisms (if present):
    polymorphisms <- sort(unique(c(grep("&", comp.states[1,]), grep("&", comp.states[2,]))))
    
    # If polymorphisms are present:
    if(length(polymorphisms) > 0) {
      
      # For each character where at least one polymorphism is encoded:
      for(j in 1:length(polymorphisms)) {
        
        # Find state(s) at one end of the branch:
        top.states <- strsplit(comp.states[, polymorphisms[j]], "&")[[1]]
        
        # Find state(s) at the other end of the branch:
        bottom.states <- strsplit(comp.states[, polymorphisms[j]], "&")[[2]]
        
        # If states overlap:
        if(length(sort(unique(c(match(top.states, bottom.states), match(bottom.states, top.states))))) > 0) {
          
          # Set effective difference as zero by making both states zero:
          comp.states[, polymorphisms[j]] <- c("0", "0")
          
        # If states do not overlap
        } else {
          
          # Create differences matrix between each set of states:
          diff.matrix <- matrix(nrow=length(top.states), ncol=length(bottom.states))
          
          # For each state at one end of the branch:
          for(k in 1:length(top.states)) {
            
            # For each state at the other end of the branch:
            for(l in 1:length(bottom.states)) {
              
              # Record the absolute difference between individual states:
              diff.matrix[k, l] <- abs(as.numeric(top.states[k]) - as.numeric(bottom.states[l]))
              
            }
            
          }
          
          # Set first state as zero and second state as the minimum difference:
          comp.states[, polymorphisms[j]] <- c("0", as.character(min(diff.matrix)))
          
        }
        
      }
      
    }
    
    # Convert comparable states to numeric now that polymorphisms are collapsed:
    comp.states <- matrix(as.numeric(comp.states), nrow=2)
    
    # Identify differences:
    comp.diffs <- abs(comp.states[1,] - comp.states[2,])
    
    # Identify distances greater than one for unordered characters if present:
    unord.diffs <- intersect(grep(TRUE, orderings == "UNORD"), grep(TRUE, comp.diffs > 1))
    
    # If distances greater than one are found for unordered characters:
    if(length(unord.diffs) > 0) {
      
      # Recode distances greater than one for unordered characters as one:
      comp.diffs[unord.diffs] <- 1
      
    }
    
    # Update branch lengths for observed tree:
    observed.tree$edge.length[i] <- sum(comp.diffs * weightings)
    
    # Update branch lengths for observable tree:
    observable.tree$edge.length[i] <- sum(ranges * weightings)
    
  }

  # Number of characters:
  Nchar <- sum(clad.matrix$weights)

  # Get row numbers of terminal branches:
  terminal.branches <- match(1:Ntip(tree), tree$edge[, 2])
  
  # Get row numbers of internal branches:
  internal.branches <- setdiff(1:nrow(tree$edge), terminal.branches)
  
  # Vectors to store time series values:
  changes <- changes.tb <- changes.ib <- Time <- Time.tb <- Time.ib <- pctcomp <- pctcomp.tb <- pctcomp.ib <- vector(mode="numeric")
  
  # For each time bin:
  for(i in 2:length(time.bins)) {
    
    # Calculate percentage completeness for all branches:
    ifelse(sum(prop.branch.in.bin[, (i - 1)]) > 0, pctcomp[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)] * completeness.tree$edge.length) / sum(prop.branch.in.bin[, (i - 1)]) / Nchar, pctcomp[(i - 1)] <- 0)
    
    # Calculate percentage completeness for terminal branches:
    ifelse(sum(prop.branch.in.bin[, (i - 1)][terminal.branches]) > 0, pctcomp.tb[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)][terminal.branches] * completeness.tree$edge.length[terminal.branches]) / sum(prop.branch.in.bin[, (i - 1)][terminal.branches]) / Nchar, pctcomp.tb[(i - 1)] <- 0)
    
    # Calculate percentage completeness for internal branches:
    ifelse(sum(prop.branch.in.bin[, (i - 1)][internal.branches]) > 0, pctcomp.ib[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)][internal.branches] * completeness.tree$edge.length[internal.branches]) / sum(prop.branch.in.bin[, (i - 1)][internal.branches]) / Nchar, pctcomp.ib[(i - 1)] <- 0)
    
    # Calculate sum of all branches (including parts) present in bin:
    Time[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)] * tree$edge.length)

    # Calculate sum of terminal branches (including parts) present in bin:
    Time.tb[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)][terminal.branches] * tree$edge.length[terminal.branches])
    
    # Calculate sum of internal branches (including parts) present in bin:
    Time.ib[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)][internal.branches] * tree$edge.length[internal.branches])
    
    # Callculate number of all changes present in bin:
    changes[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)] * observed.tree$edge.length)
    
    # Callculate number of terminal changes present in bin:
    changes.tb[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)][terminal.branches] * observed.tree$edge.length[terminal.branches])
    
    # Callculate number of internal changes present in bin:
    changes.ib[(i - 1)] <- sum(prop.branch.in.bin[, (i - 1)][internal.branches] * observed.tree$edge.length[internal.branches])
    
  }
 
  # Get maximum likelihood numerator:
  mlenumer <- sum(changes) / sum(Time * pctcomp)
  
  # Get maximum likelihood numerator:
  mlenumer.tb <- sum(changes.tb) / sum(Time.tb * pctcomp.tb)
  
  # Get maximum likelihood numerator:
  mlenumer.ib <- sum(changes.ib) / sum(Time.ib * pctcomp.ib)
  
  # Get maximum likelihood denominator:
  mledenom <- (changes / (Time * pctcomp))[which(Time > 0)]

  # Get maximum likelihood denominator:
  mledenom.tb <- (changes.tb / (Time.tb * pctcomp.tb))[which(Time.tb > 0)]

  # Get maximum likelihood denominator:
  mledenom.ib <- (changes.ib / (Time.ib * pctcomp.ib))[which(Time.ib > 0)]
  
  
# MAYBE USE CONTINUOUS APPROXIMATION FOR POISSON INSTEAD OF ROUNDING?
  
  
  # Get log numerator:
  lognumer <- sum(log(dpois(round(changes[which(Time > 0)]), mlenumer * Time[which(Time > 0)] * pctcomp[which(Time > 0)])))

  # Get log numerator:
  lognumer.tb <- sum(log(dpois(round(changes.tb[which(Time.tb > 0)]), mlenumer.tb * Time.tb[which(Time.tb > 0)] * pctcomp.tb[which(Time.tb > 0)])))
  
  # Get log numerator:
  lognumer.ib <- sum(log(dpois(round(changes.ib[which(Time.ib > 0)]), mlenumer.ib * Time.ib[which(Time.ib > 0)] * pctcomp.ib[which(Time.ib > 0)])))
  
  # Get log denominator with zero-length branches removed:
  logdenom <- sum(log(dpois(round(changes[which(Time > 0)]), mledenom * Time[which(Time > 0)] * pctcomp[which(Time > 0)])), na.rm=TRUE)

  # Get log denominator with zero-length branches removed:
  logdenom.tb <- sum(log(dpois(round(changes.tb[which(Time.tb > 0)]), mledenom.tb * Time.tb[which(Time.tb > 0)] * pctcomp.tb[which(Time.tb > 0)])), na.rm=TRUE)
  
  # Get log denominator with zero-length branches removed:
  logdenom.ib <- sum(log(dpois(round(changes.ib[which(Time.ib > 0)]), mledenom.ib * Time.ib[which(Time.ib > 0)] * pctcomp.ib[which(Time.ib > 0)])), na.rm=TRUE)
  
  # Get test statistic:
  teststat <- -2 * (lognumer - logdenom)
  
  # Get test statistic:
  teststat.tb <- -2 * (lognumer.tb - logdenom.tb)
  
  # Get test statistic:
  teststat.ib <- -2 * (lognumer.ib - logdenom.ib)
  
  # Calculate position of test statistic in chi-square distribution to get probability (zero-length branches not calculated in df):
  chisq.p <- pchisq(teststat, length(which(Time > 0)) - 1, lower.tail=F)
  
  # Calculate position of test statistic in chi-square distribution to get probability (zero-length branches not calculated in df):
  chisq.p.tb <- pchisq(teststat.tb, length(which(Time.tb > 0)) - 1, lower.tail=F)
  
  # Calculate position of test statistic in chi-square distribution to get probability (zero-length branches not calculated in df):
  chisq.p.ib <- pchisq(teststat.ib, length(which(Time.ib > 0)) - 1, lower.tail=F)
  
# NEED TO CONSIDER WHAT SHOULD DO IF TB OR IB FAIL TEST EVEN IF POOLED DATA DOES NOT

  # Check to see if null hypothesis of equal rates across the tree can be rejected:
  if(chisq.p < alpha) {
    
    # If so then print notification and carry on:
    cat(paste("H_0 - all rates equal across time bins - is rejected at an alpha of ", alpha, " (actual p = ", chisq.p, ").\nContinuing to per-bin rate calculations.\n", sep=""))

    # Create matrix to store time bin results:
    bin.results <- matrix(0, nrow=length(time.bins) - 1, ncol=10)

    # Add column names:
    colnames(bin.results) <- c("bin", "from", "to", "in.rate", "out.rate", "ml.chisq", "ml.pval", "ml.signif", "ml.signif.hi", "ml.signif.lo")

    # Fill bin numbers (from 1:n):
    bin.results[, "bin"] <- 1:(length(time.bins) - 1)

    # Fill in from values (in Ma):
    bin.results[, "from"] <- time.bins[1:(length(time.bins) - 1)]

    # Fill in to values (in Ma):
    bin.results[, "to"] <- time.bins[2:length(time.bins)]

    # For each time bin:
    for(i in 1:(length(time.bins) - 1)) {

      # Get numbers of other time bins:
      other <- (1:(length(time.bins) - 1))[-i]

      # Create maximum likelihood denominator variable:
      mledenom <- rep(NA, nrow(bin.results))

      # Calculate maximum likelihood denominator for bin:
      mledenom[i] <- mlebin <- changes[i] / (Time[i] * pctcomp[i])

      # Calculate maximum likelihood denominator for other bins:
      mledenom[other] <- mleother <- sum(changes[other]) / sum(Time[other] * pctcomp[other])

      # Calculate log numerator (rounding non-integers):
      lognumer <- sum(log(dpois(round(changes), mlenumer * Time * pctcomp)), na.rm=TRUE)

      # Calculate log denominator:
      logdenom <- sum(log(dpois(round(changes), mledenom * Time * pctcomp)), na.rm=TRUE) 

      # Store rate for branch:
      bin.results[i, "in.rate"] <- round(mlebin, 2)

      # Store rate for branch:
      bin.results[i, "out.rate"] <- round(mleother, 2)

      # Store chi-squared value (will be zero for zero-duration branch):
      bin.results[i, "ml.chisq"] <- teststat <- -2 * (lognumer - logdenom)

      # Store probability for bin (will be 1 for zero-duration branch):
      bin.results[i, "ml.pval"] <- pchisq(teststat, 2 - 1, lower.tail=FALSE)

    }

    # Get just the branch p-values:
    bin.pvals <- bin.results[, "ml.pval"]

    # Number of zero-duration bins:
    nzero <- 0

    # Do not count zero-duration branches in sample:
    m <- nrow(bin.results) - nzero

    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
    cutoffs <- c((1:m) / m * alpha, rep(0, nzero))

    # Get indices ready for identifying significant p values:
    ifelse(length(grep(TRUE, sort(bin.pvals) <= cutoffs)) > 0, indices <- 1:max((1:m)[sort(bin.pvals) <= cutoffs]), indices <- c(1)[-1])

    # Isolate significant p values:
    signif <- order(bin.pvals)[indices]

    # Add 1 in significance column for significant p-values after FDR correction:
    bin.results[signif, "ml.signif"] <- 1

    # Indicate significantly high rate bins:
    bin.results[(bin.results[, "ml.signif"] & bin.results[, "in.rate"] > bin.results[, "out.rate"]), "ml.signif.hi"] <- 1
  
    # Indicate significantly low rate bins:
    bin.results[(bin.results[, "ml.signif"] & bin.results[, "in.rate"] < bin.results[, "out.rate"]), "ml.signif.lo"] <- 1

    # Round chi-square values in output:
    bin.results[, "ml.chisq"] <- round(bin.results[, "ml.chisq"], 2)
  
    # Round p-values in output:
    bin.results[, "ml.pval"] <- round(bin.results[, "ml.pval"], 3)
  
# WHAT TO DO WITH ZERO VALUES IN TIME SERIES? EXCLUDE?
# TIME IS KEY THING TO CHECK (IF ZERO THEN NO CHANCE TO OBSERVE ANYTHING)

  # Case if equal rates cannot be rejected:
  } else {

    # If not then print notification and stop:
    cat(paste("H_0 - all rates equal across time bins - cannot be rejected at an alpha of ", alpha, " (Actual p = ", chisq.p, ").\nCalculations of per-bin rates aborted.", sep=""))
    
    # Create NULL outputs:
    bin.results <- paste("H_0 - all rates equal across time bins - cannot be rejected at an alpha of ", alpha, " (Actual p = ", chisq.p, ").\nA single rate of ", mlenumer, "is preferred.", sep="")
    
  }

  # Get number of branches:
  n <- length(tree$edge[,1])
  
  # Get number of tips (and terminal branches):
  ntips <- n.tb <- Ntip(tree)
  
  # Get number of internal branches:
  n.ib <- n - n.tb
  
  # Get total number of character changes on tree:
  changes <- observed.tree$edge.length
  
  # Get number of character changes on terminal branches:
  changes.tb <- observed.tree$edge.length[terminal.branches]
  
  # Get number of character changes on internal branches:
  changes.ib <- observed.tree$edge.length[internal.branches]
  
  # Get list of nonterminal and nonroot nodes:
  nodes <- (Ntip(tree) + 2):(Ntip(tree) + Nnode(tree))
  
  # Percentage of characters that are observable:
  pctcomp <- completeness.tree$edge.length / Nchar
  
  # Percentage of characters that are observable on terminal branches:
  pctcomp.tb <- completeness.tree$edge.length[terminal.branches] / Nchar
  
  # Percentage of characters that are observable on internal branches:
  pctcomp.ib <- completeness.tree$edge.length[internal.branches] / Nchar
  
  # Get branch durations:
  Time <- tree$edge.length
  
  # Get terminal branch durations:
  Time.tb <- tree$edge.length[terminal.branches]
  
  # Get internal branch durations:
  Time.ib <- tree$edge.length[internal.branches]
  
  # Get number of zero length branches:
  nzero <- length(grep(TRUE, Time == 0))
  
  # Get number of zero length terminal branches:
  nzero.tb <- length(grep(TRUE, Time.tb == 0))
  
  # Get number of zero length internal branches:
  nzero.ib <- length(grep(TRUE, Time.ib == 0))
  
  # Get maximum likelihood numerator:
  mlenumer <- sum(changes) / sum(Time * pctcomp)
  
  # Get maximum likelihood numerator for terminal branches:
  mlenumer.tb <- sum(changes.tb) / sum(Time.tb * pctcomp.tb)
  
  # Get maximum likelihood numerator for internal branches:
  mlenumer.ib <- sum(changes.ib) / sum(Time.ib * pctcomp.ib)
  
  # Get maximum likelihood denominator:
  mledenom <- changes / (Time * pctcomp)
  
  # Get maximum likelihood denominator for terminal branches:
  mledenom.tb <- changes.tb / (Time.tb * pctcomp.tb)
  
  # Get maximum likelihood denominator for internal branches:
  mledenom.ib <- changes.ib / (Time.ib * pctcomp.ib)
  
  # Get log numerator:
  lognumer <- sum(log(dpois(changes, mlenumer * Time * pctcomp)))
  
  # Get log numerator for terminal branches:
  lognumer.tb <- sum(log(dpois(changes.tb, mlenumer.tb * Time.tb * pctcomp.tb)))
  
  # Get log numerator for internal branches:
  lognumer.ib <- sum(log(dpois(changes.ib, mlenumer.ib * Time.ib * pctcomp.ib)))
  
  # Get log denominator with zero-length branches removed:
  logdenom <- sum(log(dpois(changes, mledenom * Time * pctcomp)), na.rm=TRUE)
  
  # Get log denominator for terminal branches with zero-length branches removed:
  logdenom.tb <- sum(log(dpois(changes.tb, mledenom.tb * Time.tb * pctcomp.tb)), na.rm=TRUE)
  
  # Get log denominator for internal branches with zero-length branches removed:
  logdenom.ib <- sum(log(dpois(changes.ib, mledenom.ib * Time.ib * pctcomp.ib)), na.rm=TRUE)
  
  # Get test statistic:
  teststat <- -2 * (lognumer - logdenom)
  
  # Get test statistic for terminal branches:
  teststat.tb <- -2 * (lognumer.tb - logdenom.tb)
  
  # Get test statistic for internal branches:
  teststat.ib <- -2 * (lognumer.ib - logdenom.ib)
  
  # Calculate position of test statistic in chi-square distribution to get probability (zero-length branches not calculated in df):
  chisq.p <- pchisq(teststat, n - 1 - nzero, lower.tail=F)
  
  # Calculate position of test statistic in chi-square distribution for terminal branches to get probability (zero-length branches not calculated in df):
  chisq.p.tb <- pchisq(teststat.tb, n.tb - 1 - nzero.tb, lower.tail=F)
  
  # Calculate position of test statistic in chi-square distribution for internal branches to get probability (zero-length branches not calculated in df):
  chisq.p.ib <- pchisq(teststat.ib, n.ib - 1 - nzero.ib, lower.tail=F)
  
  # Check to see if null hypothesis of equal rates across the tree can be rejected:
  if(chisq.p < alpha) {
    
    # If so then print notification and carry on:
    cat(paste("H_0 - all rates equal across the tree - is rejected at an alpha of ", alpha, " (actual p = ", chisq.p, ").\nContinuing to per-branch and per-clade rate calculations.", sep=""))
    
    # Create matrix to store branch results:
    branch.results <- matrix(0, nrow=n, ncol=20)
    
    # Create matrix to store terminal branch results:
    branch.results.tb <- matrix(0, nrow=n.tb, ncol=20)
    
    # Create matrix to store internal branch results:
    branch.results.ib <- matrix(0, nrow=n.ib, ncol=20)
    
    # Create matrix to store node results:
    node.results <- matrix(0, nrow=length(nodes), ncol=29)
    
    # Add column names for branches:
    colnames(branch.results) <- colnames(branch.results.tb) <- colnames(branch.results.ib) <- c("branch", "from", "to", "in.rate", "out.rate", "ml.chisq", "ml.pval", "ml.signif.hi", "ml.signif.hi.ti", "ml.signif.lo", "ml.signif.lo.ti", "ml.signif", "ml.signif.ti", "rand.val", "rand.mean", "rand.sd", "rand.pval", "rand.signif.hi", "rand.signif.lo", "rand.signif")
    
    # Add column names for nodes:
    colnames(node.results) <- c("node", "in.rate", "in.rate.tb", "in.rate.ib", "out.rate", "out.rate.tb", "out.rate.ib", "ml.chisq", "ml.chisq.tb", "ml.chisq.ib", "ml.pval", "ml.pval.tb", "ml.pval.ib", "ml.signif.hi", "ml.signif.hi.tb", "ml.signif.hi.ib", "ml.signif.lo", "ml.signif.lo.tb", "ml.signif.lo.ib", "ml.signif", "ml.signif.tb", "ml.signif.ib", "rand.val", "rand.mean", "rand.sd", "rand.pval", "rand.signif.hi", "rand.signif.lo", "rand.signif")
    
    # Number branches 1 to N:
    branch.results[, "branch"] <- 1:n
    
    # Number terminal branches:
    branch.results.tb[, "branch"] <- terminal.branches
    
    # Number internal branches:
    branch.results.ib[, "branch"] <- internal.branches
    
    # Add from and to node numbers:
    branch.results[, c("from", "to")] <- tree$edge[, 1:2]
    
    # Add from and to node numbers for terminal branches:
    branch.results.tb[, c("from", "to")] <- tree$edge[terminal.branches, 1:2]
    
    # Add from and to node numbers for internal branches:
    branch.results.ib[, c("from", "to")] <- tree$edge[internal.branches, 1:2]
    
    # Number nodes:
    node.results[, "node"] <- nodes
    
    # For each branch:
    for (i in 1:n) {
      
      # Check if branch is terminal:
      if(length(sort(match(i, terminal.branches)))) {
        
        # If yes set as TRUE:
        branch.is.terminal <- TRUE
        
      # Branch is terminal:
      } else {
        
        # If no set as FALSE:
        branch.is.terminal <- FALSE
        
      }
      
      # Get numbers of other (not ith) branches:
      other <- (1:n)[-i]
      
      # Get numbers of other (not ith) terminal or internal branches:
      ifelse(branch.is.terminal, other.tb <- terminal.branches[-match(i, terminal.branches)], other.ib <- internal.branches[-match(i, internal.branches)])
      
      # Create maximum likelihood denominator variable:
      mledenom <- rep(NA, n)
      
      # Create maximum likelihood denominator variable for terminal or internal branches:
      ifelse(branch.is.terminal, mledenom.tb <- rep(NA, n.tb), mledenom.ib <- rep(NA, n.ib))
      
      # Calculate maximum likelihood denominator for branch:
      mledenom[i] <- mlebranch <- changes[i] / (Time[i] * pctcomp[i])
      
      # Calculate maximum likelihood denominator for terminal or internal branch:
      ifelse(branch.is.terminal, mledenom.tb[match(i, terminal.branches)] <- mlebranch.tb <- changes.tb[match(i, terminal.branches)] / (Time.tb[match(i, terminal.branches)] * pctcomp.tb[match(i, terminal.branches)]), mledenom.ib[match(i, internal.branches)] <- mlebranch.ib <- changes.ib[match(i, internal.branches)] / (Time.ib[match(i, internal.branches)] * pctcomp.ib[match(i, internal.branches)]))
      
      # Calculate maximum likelihood denominator for other branches:
      mledenom[other] <- mleother <- sum(changes[other]) / sum(Time[other] * pctcomp[other])
      
      # Calculate maximum likelihood denominator for other terminal or internal branches:
      ifelse(branch.is.terminal, mledenom.tb[-match(i, terminal.branches)] <- mleother.tb <- sum(changes[other.tb]) / sum(Time[other.tb] * pctcomp[other.tb]), mledenom.ib[-match(i, internal.branches)] <- mleother.ib <- sum(changes[other.ib]) / sum(Time[other.ib] * pctcomp[other.ib]))
      
      # Calculate log numerator:
      lognumer <- sum(log(dpois(changes, mlenumer * Time * pctcomp)), na.rm=TRUE)
      
      # Calculate log numerator for terminal or internal branches:
      ifelse(branch.is.terminal, lognumer.tb <- sum(log(dpois(changes.tb, mlenumer.tb * Time.tb * pctcomp.tb)), na.rm=TRUE), lognumer.ib <- sum(log(dpois(changes.ib, mlenumer.ib * Time.ib * pctcomp.ib)), na.rm=TRUE))
      
      # Calculate log denominator:
      logdenom <- sum(log(dpois(changes, mledenom * Time * pctcomp)), na.rm=TRUE) 
      
      # Calculate log denominator for terminal or internal branches:
      ifelse(branch.is.terminal, logdenom.tb <- sum(log(dpois(changes.tb, mledenom.tb * Time.tb * pctcomp.tb)), na.rm=TRUE), logdenom.ib <- sum(log(dpois(changes.ib, mledenom.ib * Time.ib * pctcomp.ib)), na.rm=TRUE))
      
      # Store rate for branch:
      branch.results[i, "in.rate"] <- round(mlebranch, 2)
      
      # Store rate for terminal or internal branch:
      ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "in.rate"] <- round(mlebranch.tb, 2), branch.results.ib[match(i, branch.results.ib[, "branch"]), "in.rate"] <- round(mlebranch.ib, 2))
      
      # Store rate for outside branches:
      branch.results[i, "out.rate"] <- round(mleother, 2)
      
      # Store rate for outside terminal or internal branches:
      ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "out.rate"] <- round(mleother.tb, 2), branch.results.ib[match(i, branch.results.ib[, "branch"]), "out.rate"] <- round(mleother.ib, 2))
      
      # Store chi-squared value (will be zero for zero-duration branch):
      branch.results[i, "ml.chisq"] <- teststat <- -2 * (lognumer - logdenom)
      
      # Store chi-squared value (will be zero for zero-duration branch) for internal or terminal branches:
      ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "ml.chisq"] <- teststat.tb <- -2 * (lognumer.tb - logdenom.tb), branch.results.ib[match(i, branch.results.ib[, "branch"]), "ml.chisq"] <- teststat.ib <- -2 * (lognumer.ib - logdenom.ib))
      
      # Store probability for branch (will be 1 for zero-duration branch):
      branch.results[i, "ml.pval"] <- pchisq(teststat, 2-1, lower.tail=FALSE)
      
      # Store probability for branch (will be 1 for zero-duration branch) for terminal or internal branches:
      ifelse(branch.is.terminal, branch.results.tb[match(i, branch.results.tb[, "branch"]), "ml.pval"] <- pchisq(teststat.tb, 2 - 1, lower.tail=FALSE), branch.results.ib[match(i, branch.results.ib[, "branch"]), "ml.pval"] <- pchisq(teststat.ib, 2 - 1, lower.tail=FALSE))
      
    }
    
    # Set initial row number for storing data:
    j <- 1
    
    # For each node (excluding the root):
    for (i in nodes) {
      
      # Identify branches within clade:
      clade <- GetDescendantEdges(i, tree)
      
      # Identify terminal branches within clade:
      clade.tb <- clade[sort(match(terminal.branches, clade))]
      
      # Identify internal branches within clade:
      clade.ib <- clade[sort(match(internal.branches, clade))]
      
      # Identify branches outside clade:
      nonclade <- (1:n)[-clade]
      
      # Identify terminal branches outside clade:
      nonclade.tb <- setdiff(terminal.branches, clade.tb)
      
      # Identify internal branches outside clade:
      nonclade.ib <- setdiff(internal.branches, clade.ib)
      
      # Set empty maximum likelihood denominator vector:
      mledenom <- rep(NA, n)
      
      # Set empty maximum likelihood denominator vector for terminal branches:
      mledenom.tb <- rep(NA, n.tb)
      
      # Set empty maximum likelihood denominator vector for internal branches:
      mledenom.ib <- rep(NA, n.ib)
      
      # Fill within clade values for maximum likelihood denominator:
      mledenom[clade] <- claderate <- sum(changes[clade]) / sum(Time[clade] * pctcomp[clade])
      
      # Fill within clade values for maximum likelihood denominator for terminal branches:
      mledenom.tb[match(clade.tb, terminal.branches)] <- claderate.tb <- sum(changes[clade.tb]) / sum(Time[clade.tb] * pctcomp[clade.tb])
      
      # Fill within clade values for maximum likelihood denominator for internal branches (if present):
      if(length(clade.ib) > 0) mledenom.ib[match(clade.ib, internal.branches)] <- claderate.ib <- sum(changes[clade.ib]) / sum(Time[clade.ib] * pctcomp[clade.ib])
      
      # Fill outside clade values for maximum likelihood denominator:
      mledenom[nonclade] <- noncladerate <- sum(changes[nonclade]) / sum(Time[nonclade] * pctcomp[nonclade])
      
      # Fill outside clade values for maximum likelihood denominator for terminal branches:
      mledenom.tb[match(nonclade.tb, terminal.branches)] <- noncladerate.tb <- sum(changes[nonclade.tb]) / sum(Time[nonclade.tb] * pctcomp[nonclade.tb])
      
      # Fill outside clade values for maximum likelihood denominator for internal branches (if present):
      if(length(clade.ib) > 0) mledenom.ib[match(nonclade.ib, internal.branches)] <- noncladerate.ib <- sum(changes[nonclade.ib]) / sum(Time[nonclade.ib] * pctcomp[nonclade.ib])
      
      # Set log-likelihood numerator:
      lognumer <- sum(log(dpois(changes, mlenumer * Time * pctcomp)), na.rm=T)
      
      # Set log-likelihood numerator for terminal branches:
      lognumer.tb <- sum(log(dpois(changes.tb, mlenumer.tb * Time.tb * pctcomp.tb)), na.rm=T)
      
      # Set log-likelihood numerator for internal branches (if present):
      if(length(clade.ib) > 0) lognumer.ib <- sum(log(dpois(changes.ib, mlenumer.ib * Time.ib * pctcomp.ib)), na.rm=T)
      
      # Set log-likelihood denominator:
      logdenom <- sum(log(dpois(changes, mledenom * Time * pctcomp)), na.rm=T)
      
      # Set log-likelihood denominator for terminal branches:
      logdenom.tb <- sum(log(dpois(changes.tb, mledenom.tb * Time.tb * pctcomp.tb)), na.rm=T)
      
      # Set log-likelihood denominator for internal branches (if present):
      if(length(clade.ib) > 0) logdenom.ib <- sum(log(dpois(changes.ib, mledenom.ib * Time.ib * pctcomp.ib)), na.rm=T)
      
      # Record within clade rate:
      node.results[j, "in.rate"] <- round(claderate, 2)
      
      # Record within clade rate for terminal branches:
      node.results[j, "in.rate.tb"] <- round(claderate.tb, 2)
      
      # Record within clade rate for internal branches (if present):
      ifelse(length(clade.ib) > 0, node.results[j, "in.rate.ib"] <- round(claderate.ib, 2), node.results[j, "in.rate.ib"] <- NA)
      
      # Record outside clade rate:
      node.results[j, "out.rate"] <- round(noncladerate, 2)
      
      # Record outside clade rate for terminal branches:
      node.results[j, "out.rate.tb"] <- round(noncladerate.tb, 2)
      
      # Record outside clade rate for internal branches (if present):
      ifelse(length(clade.ib) > 0, node.results[j, "out.rate.ib"] <- round(noncladerate.ib, 2), node.results[j, "out.rate.ib"] <- NA)
      
      # Record chi-square test statistic:
      node.results[j, "ml.chisq"] <- teststat <- -2 * (lognumer - logdenom)
      
      # Record chi-square test statistic for terminal branches:
      node.results[j, "ml.chisq.tb"] <- teststat.tb <- -2 * (lognumer.tb - logdenom.tb)
      
      # Record chi-square test statistic for internal branches (if present):
      ifelse(length(clade.ib) > 0, node.results[j, "ml.chisq.ib"] <- teststat.ib <- -2 * (lognumer.ib - logdenom.ib), node.results[j, "ml.chisq.ib"] <- NA)
      
      # Record p-value for chi-squared test:
      node.results[j, "ml.pval"] <- testresult <- pchisq(teststat, 2 - 1, lower.tail=F)
      
      # Record p-value for chi-squared test for terminal branches:
      node.results[j, "ml.pval.tb"] <- testresult.tb <- pchisq(teststat.tb, 2 - 1, lower.tail=F)
      
      # Record p-value for chi-squared test for terminal branches:
      ifelse(length(clade.ib) > 0, node.results[j, "ml.pval.ib"] <- testresult.ib <- pchisq(teststat.ib, 2 - 1, lower.tail=F), node.results[j, "ml.pval.ib"] <- NA)
      
      # Update row number:
      j <- j + 1
      
    }
    
    # Combine terminal and internal branch rate results and order as branch.results:
    termandintern.results <- rbind(branch.results.tb, branch.results.ib)[order(rbind(branch.results.tb, branch.results.ib)[, "branch"]), ]
    
    # Create vector to store true (1; terminal) / false (0; internal) branch type:
    branchisterminal <- rep(0, nrow(tree$edge))
    
    # Store 1 for terminals:
    branchisterminal[terminal.branches] <- 1
    
    # Combine and reduce branch type and terminal and internal branch rate results (ignores randomisations):
    termandintern.results <- cbind(branchisterminal, termandintern.results[, c("in.rate", "out.rate", "ml.chisq", "ml.pval")])
    
    # Update column names (to avoid clashes with branch.results):
    colnames(termandintern.results) <- c("is.term", "in.rate.ti", "out.rate.ti", "ml.chisq.ti", "ml.pval.ti")
    
    # Combine all branch rates versus split terminal-internal branches:
    branch.results <- cbind(branch.results, termandintern.results)
    
    # Reorder and cut down (ignores randomisation results):
    branch.results <- branch.results[, c("branch", "from", "to", "is.term", "in.rate", "in.rate.ti", "out.rate", "out.rate.ti", "ml.chisq", "ml.chisq.ti", "ml.pval", "ml.pval.ti", "ml.signif.hi", "ml.signif.hi.ti", "ml.signif.lo", "ml.signif.lo.ti", "ml.signif", "ml.signif.ti")]
    
    # Get just the branch p-values:
    branch.pvals <- branch.results[, "ml.pval"]
    
    # Get just the terminal branch p-values:
    branch.pvals.tb <- branch.results[terminal.branches, "ml.pval.ti"]
    
    # Get just the internal branch p-values:
    branch.pvals.ib <- branch.results[internal.branches, "ml.pval.ti"]
    
    # Do not count zero-duration branches in sample:
    m <- n - nzero
    
    # Do not count zero-duration terminal branches in sample:
    m.tb <- n.tb - nzero.tb
    
    # Do not count zero-duration internal branches in sample:
    m.ib <- n.ib - nzero.ib
    
    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
    cutoffs <- c((1:m) / m * alpha, rep(0, nzero))
    
    # Calculate cutoffs (significance thresholds) for terminal branches:
    cutoffs.tb <- c((1:m.tb) / m.tb * alpha, rep(0, nzero.tb))
    
    # Calculate cutoffs (significance thresholds) for internal branches:
    cutoffs.ib <- c((1:m.ib) / m.ib * alpha, rep(0, nzero.ib))
    
    # Get indices ready for identifying significant p values:
    ifelse(length(grep(TRUE, sort(branch.pvals) <= cutoffs)) > 0, indices <- 1:max((1:m)[sort(branch.pvals) <= cutoffs]), indices <- c(1)[-1])
    
    # Get indices ready for identifying significant terminal p values:
    ifelse(length(grep(TRUE, sort(branch.pvals.tb) <= cutoffs.tb)) > 0, indices.tb <- 1:max((1:m.tb)[sort(branch.pvals.tb) <= cutoffs.tb]), indices.tb <- c(1)[-1])
    
    # Get indices ready for identifying significant internal p values:
    ifelse(length(grep(TRUE, sort(branch.pvals.ib) <= cutoffs.ib)) > 0, indices.ib <- 1:max((1:m.ib)[sort(branch.pvals.ib) <= cutoffs.ib]), indices.ib <- c(1)[-1])
    
    # Isolate significant p values:
    signif <- order(branch.pvals)[indices]
    
    # Isolate significant terminal p values:
    signif.tb <- terminal.branches[order(branch.pvals.tb)[indices.tb]]
    
    # Isolate significant internal p values:
    signif.ib <- internal.branches[order(branch.pvals.ib)[indices.ib]]
    
    # Add 1 in significance column for significant p-values after FDR correction:
    branch.results[signif, "ml.signif"] <- 1
    
    # Add 1 in significance column for significant terminal p-values after FDR correction:
    branch.results[signif.tb, "ml.signif.ti"] <- 1
    
    # Add 1 in significance column for significant internal p-values after FDR correction:
    branch.results[signif.ib, "ml.signif.ti"] <- 1
    
    # Indicate significantly high rate branches:
    branch.results[(branch.results[, "ml.signif"] & branch.results[, "in.rate"] > branch.results[, "out.rate"]), "ml.signif.hi"] <- 1
    
    # Indicate significantly high rate terminal and internal branches:
    branch.results[(branch.results[, "ml.signif.ti"] & branch.results[, "in.rate.ti"] > branch.results[, "out.rate.ti"]), "ml.signif.hi.ti"] <- 1
    
    # Indicate significantly low rate branches:
    branch.results[(branch.results[, "ml.signif"] & branch.results[, "in.rate"] < branch.results[, "out.rate"]), "ml.signif.lo"] <- 1
    
    # Indicate significantly low rate terminal and internal branches:
    branch.results[(branch.results[, "ml.signif.ti"] & branch.results[, "in.rate.ti"] < branch.results[, "out.rate.ti"]), "ml.signif.lo.ti"] <- 1
    
    # Round chi-square values in output:
    branch.results[, "ml.chisq"] <- round(branch.results[, "ml.chisq"], 2)
    
    # Round chi-square values in output:
    branch.results[, "ml.chisq.ti"] <- round(branch.results[, "ml.chisq.ti"], 2)
    
    # Round p-values in output:
    branch.results[, "ml.pval"] <- round(branch.results[, "ml.pval"], 3)

    # Round p-values in output:
    branch.results[, "ml.pval.ti"] <- round(branch.results[, "ml.pval.ti"], 3)
    
    # Isolate the node p-values:
    node.pvals <- node.results[, "ml.pval"]
    
    # Isolate the node p-values for terminal branches:
    node.pvals.tb <- node.results[, "ml.pval.tb"]
    
    # Isolate the node p-values for internal branches:
    node.pvals.ib <- node.results[, "ml.pval.ib"]
    
    # Set m as number of nodes (excluding root), serves for terminal branches alone too:
    m <- length(nodes)
    
    # Set m as number of non-cherry nodes (excluding root and NAs):
    m.ib <- length(sort(node.pvals.ib))
    
    # Calculate cutoffs (significance thresholds) based on Benjamini and Hochberg 1995 JRSSB False Discovery Rate (FDR):
    cutoffs <- (1:m) / m * alpha
    
    # Calculate cutoffs for internal branches:
    cutoffs.ib <- (1:m.ib) / m.ib * alpha
    
    # Get indices ready for identifying significant p values:
    ifelse(length(grep(TRUE, sort(node.pvals) <= cutoffs)) > 0, indices <- 1:max((1:m)[sort(node.pvals) <= cutoffs]), indices <- c(1)[-1])
    
    # Get indices ready for identifying significant p values for terminal branches:
    ifelse(length(grep(TRUE, sort(node.pvals.tb) <= cutoffs)) > 0, indices.tb <- 1:max((1:m)[sort(node.pvals.tb) <= cutoffs]), indices.tb <- c(1)[-1])
    
    # Get indices ready for identifying significant p values for internal branches:
    ifelse(length(grep(TRUE, sort(node.pvals.ib) <= cutoffs.ib)) > 0, indices.ib <- 1:max((1:m.ib)[sort(node.pvals.ib) <= cutoffs.ib]), indices.ib <- c(1)[-1])
    
    # Isolate significant p-values:
    signif <- order(node.pvals)[indices]
    
    # Isolate significant p-values for terminal branches:
    signif.tb <- order(node.pvals.tb)[indices.tb]
    
    # Isolate significant p-values for internal branches:
    signif.ib <- order(node.pvals.ib)[indices.ib]
    
    # Add NAs for terminal branch only clades (cherries):
    node.results[is.na(node.results[, "ml.chisq.ib"]), c("ml.signif.hi.ib", "ml.signif.lo.ib", "ml.signif.ib")] <-  NA
    
    # Record significant clades:
    node.results[signif, "ml.signif"] <- 1
    
    # Record significant terminal branches clades:
    node.results[signif.tb, "ml.signif.tb"] <- 1
    
    # Record significant internal branches clades:
    node.results[signif.ib, "ml.signif.ib"] <- 1
    
    # Round chi-squared test statistic to 1dp:
    node.results[, "ml.chisq"] <- round(node.results[, "ml.chisq"], 1)
    
    # Round chi-squared test statistic to 1dp:
    node.results[, "ml.chisq.tb"] <- round(node.results[, "ml.chisq.tb"], 1)
    
    # Round chi-squared test statistic to 1dp:
    node.results[, "ml.chisq.ib"] <- round(node.results[, "ml.chisq.ib"], 1)
    
    # Round chi-squared p-value to 3dp:
    node.results[, "ml.pval"] <- round(node.results[, "ml.pval"], 3)
    
    # Round chi-squared p-value to 3dp:
    node.results[, "ml.pval.tb"] <- round(node.results[, "ml.pval.tb"], 3)
    
    # Round chi-squared p-value to 3dp:
    node.results[, "ml.pval.ib"] <- round(node.results[, "ml.pval.ib"], 3)
    
    # Record significantly high clade rates:
    node.results[(node.results[, "ml.signif"] & node.results[, "in.rate"] > node.results[, "out.rate"]), "ml.signif.hi"] <- 1
    
    # Record significantly high clade rates:
    node.results[(node.results[, "ml.signif.tb"] & node.results[, "in.rate.tb"] > node.results[, "out.rate.tb"]), "ml.signif.hi.tb"] <- 1
    
    # Record significantly high clade rates:
    node.results[intersect(grep(TRUE, node.results[, "ml.signif.ib"] == 1), grep(TRUE, node.results[, "in.rate.ib"] > node.results[, "out.rate.ib"])), "ml.signif.hi.ib"] <- 1
    
    # Record signficiantly low rates:
    node.results[(node.results[, "ml.signif"] & node.results[, "in.rate"] < node.results[, "out.rate"]), "ml.signif.lo"] <- 1
    
    # Record signficiantly low rates:
    node.results[(node.results[, "ml.signif.tb"] & node.results[, "in.rate.tb"] < node.results[, "out.rate.tb"]), "ml.signif.lo.tb"] <- 1
    
    # Record significantly low rates:
    node.results[intersect(grep(TRUE, node.results[, "ml.signif.ib"] == 1), grep(TRUE, node.results[, "in.rate.ib"] < node.results[, "out.rate.ib"])), "ml.signif.lo.ib"] <- 1

  # Case if equal rates cannot be rejected:
  } else {
    
    # If not then print notification and stop:
    cat(paste("H_0 - all rates equal across the tree - cannot be rejected at an alpha of ", alpha, " (Actual p = ", chisq.p, ").\nCalculations of per-branch and per-clade rates aborted.", sep=""))
    
    # Create NULL outputs:
    branch.results <- node.results <- NULL
    
  }
  
  # List output matrices:
  out <- list(node.results, branch.results, bin.results)
  
  # Add names to them:
  names(out) <- c("node.results", "branch.results", "per.bin.rates")
  
  # Return results:
  return(out)
  
}
