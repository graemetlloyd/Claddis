splitTree <- function (tree, split) {
    if (split$node > length(tree$tip.label)) {
        tr2 <- extract.clade(tree, node = split$node)
        tr2$root.edge <- tree$edge.length[which(tree$edge[, 2] ==
        split$node)] - split$bp
        tr1 <- drop.clade(tree, tr2$tip.label)
        nn <- if (!is.null(tree$node.label))
        c(tree$node.label, "NA")
        else "NA"
        tr1$tip.label[which(tr1$tip.label %in% nn)] <- "NA"
        tr1$edge.length[match(which(tr1$tip.label == "NA"), tr1$edge[,
        2])] <- split$bp
    }
    else {
        tr2 <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = tree$tip.label[split$node],
        edge.length = tree$edge.length[which(tree$edge[,
        2] == split$node)] - split$bp, Nnode = 1L)
        class(tr2) <- "phylo"
        tr1 <- tree
        tr1$edge.length[match(which(tr1$tip.label == tr2$tip.label[1]),
        tr1$edge[, 2])] <- split$bp
        tr1$tip.label[which(tr1$tip.label == tr2$tip.label[1])] <- "NA"
    }
    trees <- list(tr1, tr2)
    class(trees) <- "multiPhylo"
    trees
}

reroot <- function (tree, node.number, position) {
    if (class(tree) != "phylo")
    stop("tree object must be of class 'phylo.'")
    tt <- splitTree(tree, list(node = node.number, bp = position))
    p <- tt[[1]]
    d <- tt[[2]]
    tip <- if (length(which(p$tip.label == "NA")) > 0)
    "NA"
    else p$tip.label[which(p$tip.label %in% tree$node.label)]
    p <- root(p, outgroup = tip, resolve.root = T)
    bb <- which(p$tip.label == tip)
    p$tip.label[bb] <- "NA"
    ee <- p$edge.length[which(p$edge[, 2] == bb)]
    p$edge.length[which(p$edge[, 2] == bb)] <- 0
    cc <- p$edge[which(p$edge[, 2] == bb), 1]
    dd <- setdiff(p$edge[which(p$edge[, 1] == cc), 2], bb)
    p$edge.length[which(p$edge[, 2] == dd)] <- p$edge.length[which(p$edge[,
    2] == dd)] + ee
    tt <- paste.tree(p, d)
    return(tt)
}

rerootingMethod <- function (tree, x, model = c("ER", "SYM"), ...) {
    if (hasArg(tips))
    tips <- list(...)$tips
    else tips <- NULL
    if (!is.matrix(model))
    model <- model[1]
    n <- length(tree$tip.label)
    if (!is.matrix(x)) {
        yy <- to.matrix(x, sort(unique(x)))
        if (is.null(tips))
        tips <- FALSE
    }
    else {
        if (is.null(tips))
        tips <- TRUE
        yy <- x
    }
    yy <- yy[tree$tip.label, ]
    yy <- yy/rowSums(yy)
    XX <- matrix(NA, tree$Nnode + n, ncol(yy))
    rownames(XX) <- 1:(tree$Nnode + n)
    colnames(XX) <- colnames(yy)
    YY <- apeAce(tree, yy, model = model)
    XX[n + 1, ] <- YY$lik.anc[1, ]
    Q <- matrix(c(0, YY$rates)[YY$index.matrix + 1], ncol(XX),
    ncol(XX), dimnames = list(colnames(XX), colnames(XX)))
    diag(Q) <- -colSums(Q, na.rm = TRUE)
    for (i in 1:(tree$Nnode + n)) {
        if (i != (n + 1)) {
            if (i > n || tips) {
                tt <- reroot(tree, node.number = i, position = tree$edge.length[which(tree$edge[,
                2] == i)])
                XX[i, ] <- apeAce(tt, yy, model = model, fixedQ = Q)$lik.anc[1,
                ]
            }
        }
    }
    rownames(XX)[1:n] <- tree$tip.label
    XX <- if (tips)
    XX
    else XX[1:tree$Nnode + n, ]
    return(list(loglik = YY$loglik, Q = Q, marginal.anc = XX))
}
