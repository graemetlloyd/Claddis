MatrixExp <- function (mat, t = 1, method = NULL, ...) {
    if (!is.matrix(mat) || (nrow(mat) != ncol(mat)))
    stop("\"mat\" must be a square matrix")
    qmodel <- if (is.qmatrix(mat) && !is.null(method) && method ==
    "analytic")
    msm.form.qmodel(mat)
    else list(iso = 0, perm = 0, qperm = 0)
    if (!is.null(method) && method == "analytic") {
        if (!is.qmatrix(mat))
        warning("Analytic method not available since matrix is not a Markov model intensity matrix. Using \"pade\".")
        else if (qmodel$iso == 0)
        warning("Analytic method not available for this Markov model structure. Using \"pade\".")
    }
    if (length(t) > 1)
    res <- array(dim = c(dim(mat), length(t)))
    for (i in seq(along = t)) {
        if (is.null(method) || !(method %in% c("pade", "series",
        "analytic"))) {
            if (is.null(method))
            method <- eval(formals(expm::expm)$method)
            resi <- expm::expm(t[i] * mat, method = method, ...)
        }
        else {
            ccall <- .C("MatrixExpR", as.double(mat), as.integer(nrow(mat)),
            res = double(length(mat)), as.double(t[i]), as.integer(match(method,
            c("pade", "series"))), as.integer(qmodel$iso),
            as.integer(qmodel$perm), as.integer(qmodel$qperm),
            as.integer(0), NAOK = TRUE)
            resi <- matrix(ccall$res, nrow = nrow(mat))
        }
        if (length(t) == 1)
        res <- resi
        else res[, , i] <- resi
    }
    res
}

expm <- function (Y) {
    Z <- MatrixExp(Y)
    dimnames(Z) <- dimnames(Y)
    return(Z)
}

rstate <- function (y) {
    if (length(y) == 1)
    return(names(y)[1])
    else return(names(which(rmultinom(1, 1, y/sum(y))[, 1] ==
    1)))
}

getPars<-function(bt,xx,model,Q,tree,tol,m,liks=TRUE,pi="equal"){
    obj<-fitMk(bt,xx,model,fixedQ=Q,output.liks=liks,pi=pi)
    N<-length(bt$tip.label)
    II<-obj$index.matrix+1
    lvls<-obj$states
    if(liks){
        L<-obj$lik.anc
        rownames(L)<-N+1:nrow(L)
        if(!is.binary.tree(tree)){
            ancNames<-matchNodes(tree,bt)
            L<-L[as.character(ancNames[,2]),]
            rownames(L)<-ancNames[,1]
        }
        L<-rbind(xx,L)
        rownames(L)[1:N]<-1:N
    } else L<-NULL
    Q<-matrix(c(0,obj$rates)[II],m,m,dimnames=list(lvls,lvls))
    if(any(rowSums(Q,na.rm=TRUE)<tol)){
        message(paste("\nWarning: some rows of Q not numerically distinct from 0; setting to",tol,"\n"))
        ii<-which(rowSums(Q,na.rm=TRUE)<tol)
        for(i in 1:length(ii)) Q[ii[i],setdiff(1:ncol(Q),ii[i])]<-tol/(ncol(Q)-1)
    }
    diag(Q)<--rowSums(Q,na.rm=TRUE)
    return(list(Q=Q,L=L,loglik=logLik(obj)))
}

# convert vector of x to binary matrix
# written by Liam J. Revell 2012
to.matrix<-function(x,seq){
    X<-matrix(0,length(x),length(seq),dimnames=list(names(x),seq))
    for(i in 1:length(seq)) X[x==seq[i],i]<-1
    return(X)
}

# function does the stochastic mapping, conditioned on our model & given the conditional likelihoods
# written by Liam J. Revell 2013
smap<-function(tree,x,N,m,root,L,Q,pi,logL){
    # create the map tree object
    mtree<-tree; mtree$maps<-list()
    mtree$mapped.edge<-matrix(0,nrow(tree$edge),m,dimnames=list(paste(tree$edge[,1],",",tree$edge[,2],sep=""),colnames(L)))
    # now we want to simulate the node states & histories by pre-order traversal
    NN<-matrix(NA,nrow(tree$edge),2) # our node values
    NN[which(tree$edge[,1]==root),1]<-rstate(L[as.character(root),]*pi/sum(L[as.character(root),]*pi)) # assign root
    for(j in 1:nrow(tree$edge)){
        # conditioned on the start value, assign end value of node (if internal)
        p<-expm(Q*tree$edge.length[j])[NN[j,1],]*L[as.character(tree$edge[j,2]),]
        NN[j,2]<-rstate(p/sum(p))
        NN[which(tree$edge[,1]==tree$edge[j,2]),1]<-NN[j,2]
        # now simulate on the branches
        accept<-FALSE
        while(!accept){
            map<-sch(NN[j,1],tree$edge.length[j],Q)
            if(names(map)[length(map)]==NN[j,2]) accept=TRUE
        }
        mtree$maps[[j]]<-map
        for(k in 1:length(mtree$maps[[j]]))
        mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]<-mtree$mapped.edge[j,names(mtree$maps[[j]])[k]]+mtree$maps[[j]][k]
    }
    mtree$Q<-Q
    mtree$logL<-logL
    if(!inherits(mtree,"simmap")) class(mtree)<-c("simmap",setdiff(class(mtree),"simmap"))
    return(mtree)
}

# function generates a history along a branch
# written by Liam J. Revell 2013
sch<-function(start,t,Q){
    tol<-t*1e-12
    dt<-setNames(0,start)
    while(sum(dt)<(t-tol)){
        s<-names(dt)[length(dt)]
        dt[length(dt)]<-if(-Q[s,s]>0) rexp(n=1,rate=-Q[s,s]) else t-sum(dt)
        if(sum(dt)<(t-tol)){
            dt<-c(dt,0)
            if(sum(Q[s,][-match(s,colnames(Q))])>0)
            names(dt)[length(dt)]<-rstate(Q[s,][-match(s,colnames(Q))]/sum(Q[s,][-match(s,colnames(Q))]))
            else names(dt)[length(dt)]<-s
        } else dt[length(dt)]<-dt[length(dt)]-sum(dt)+t
    }
    return(dt)
}

# function uses numerical optimization to solve for the stationary distribution
# written by Liam J. Revell 2013
statdist<-function(Q){
    foo<-function(theta,Q){
        Pi<-c(theta[1:(nrow(Q)-1)],1-sum(theta[1:(nrow(Q)-1)]))
        sum((Pi%*%Q)^2)
    }
    k<-nrow(Q)
    if(nrow(Q)>2){
        fit<-optim(rep(1/k,k-1),foo,Q=Q,control=list(reltol=1e-16))
        return(setNames(c(fit$par[1:(k-1)],1-sum(fit$par[1:(k-1)])),rownames(Q)))
    } else {
        fit<-optimize(foo,interval=c(0,1),Q=Q)
        return(setNames(c(fit$minimum,1-fit$minimum),rownames(Q)))
    }
}

## S3 print method for objects of class "simmap" & multiSimmap
## based on print.phylo in ape
print.simmap<-function(x,printlen=6,...){
    N<-Ntip(x)
    M<-x$Nnode
    cat(paste("\nPhylogenetic tree with",N,"tips and",M,"internal nodes.\n\n"))
    cat("Tip labels:\n")
    if(N>printlen) cat(paste("\t",paste(x$tip.label[1:printlen],collapse=", "),", ...\n",sep=""))
    else print(x$tip.label)
    ss<-sort(unique(c(getStates(x,"tips"),getStates(x,"nodes"))))
    cat(paste("\nThe tree includes a mapped, ",length(ss),"-state discrete character with states:\n",
    sep=""))
    if(length(ss)>printlen) cat(paste("\t",paste(ss[1:printlen],collapse=", "),", ...\n",sep=""))
    else cat(paste("\t",paste(ss,collapse=", "),"\n",sep=""))
    rlab<-if(is.rooted(x)) "Rooted" else "Unrooted"
    cat("\n",rlab,"; includes branch lengths.\n",sep="")
}
print.multiSimmap<-function(x,details=FALSE,...){
    N<-length(x)
    cat(N,"phylogenetic trees with mapped discrete characters\n")
    if(details){
        n<-sapply(x,Ntip)
        s<-sapply(x,function(x) length(unique(c(getStates(x,"tips"),getStates(x,"nodes")))))
        for(i in 1:N) cat("tree",i,":",n[i],"tips,",s[i],"mapped states\n")
    }
}

## S3 summary method for objects of class "simmap" & "multiSimmap"
summary.simmap<-function(object,...) describe.simmap(object,...)
summary.multiSimmap<-function(object,...) describe.simmap(object,...)

## for backward compatibility with any function using apeAce internally
apeAce<-function(tree,x,model,fixedQ=NULL,...){
    if(hasArg(output.liks)){ 
        output.liks<-list(...)$output.liks
        return(fitMk(tree,x,model,fixedQ,...))
    } else { 
        output.liks<-TRUE
        return(fitMk(tree,x,model,fixedQ,output.liks=TRUE,...))
    }
}

make.simmap <- function (tree, x, model = "SYM", nsim = 1, ...) {
    if (class(tree) == "multiPhylo") {
        ff <- function(yy, x, model, nsim, ...) {
            zz <- make.simmap(yy, x, model, nsim, ...)
            if (nsim > 1)
            class(zz) <- NULL
            return(zz)
        }
        if (nsim > 1)
        mtrees <- unlist(sapply(tree, ff, x, model, nsim,
        ..., simplify = FALSE), recursive = FALSE)
        else mtrees <- sapply(tree, ff, x, model, nsim, ...,
        simplify = FALSE)
        class(mtrees) <- "multiPhylo"
    }
    else {
        if (hasArg(pi))
        pi <- list(...)$pi
        else pi <- "equal"
        if (hasArg(message))
        pm <- list(...)$message
        else pm <- TRUE
        if (hasArg(tol))
        tol <- list(...)$tol
        else tol <- 0
        if (hasArg(Q))
        Q <- list(...)$Q
        else Q <- "empirical"
        if (hasArg(burnin))
        burnin <- list(...)$burnin
        else burnin <- 1000
        if (hasArg(samplefreq))
        samplefreq <- list(...)$samplefreq
        else samplefreq <- 100
        if (hasArg(vQ))
        vQ <- list(...)$vQ
        else vQ <- 0.1
        prior <- list(alpha = 1, beta = 1, use.empirical = FALSE)
        if (hasArg(prior)) {
            pr <- list(...)$prior
            prior[names(pr)] <- pr
        }
        if (class(tree) != "phylo")
        stop("'tree' should be an object of class 'phylo'")
        if (!is.matrix(x))
        xx <- to.matrix(x, sort(unique(x)))
        else xx <- x
        xx <- xx[tree$tip.label, ]
        xx <- xx/rowSums(xx)
        tree <- bt <- reorder.phylo(tree, "cladewise")
        if (!is.binary.tree(bt))
        bt <- multi2di(bt)
        N <- length(tree$tip)
        m <- ncol(xx)
        root <- N + 1
        if (is.character(Q) && Q == "empirical") {
            XX <- getPars(bt, xx, model, Q = NULL, tree, tol,
            m)
            L <- XX$L
            Q <- XX$Q
            logL <- XX$loglik
            if (pi[1] == "equal")
            pi <- setNames(rep(1/m, m), colnames(L))
            else if (pi[1] == "estimated")
            pi <- statdist(Q)
            else pi <- pi/sum(pi)
            if (pm)
            printmessage(Q, pi, method = "empirical")
            mtrees <- replicate(nsim, smap(tree, x, N, m, root,
            L, Q, pi, logL), simplify = FALSE)
        }
        else if (is.character(Q) && Q == "mcmc") {
            if (prior$use.empirical) {
                qq <- apeAce(bt, xx, model)$rates
                prior$alpha <- qq * prior$beta
            }
            XX <- mcmcQ(bt, xx, model, tree, tol, m, burnin,
            samplefreq, nsim, vQ, prior)
            L <- lapply(XX, function(x) x$L)
            Q <- lapply(XX, function(x) x$Q)
            logL <- lapply(XX, function(x) x$loglik)
            if (pi[1] == "equal") {
                pi <- setNames(rep(1/m, m), colnames(L))
                pi <- lapply(1:nsim, function(x, y) y, y = pi)
            }
            else if (pi[1] == "estimated") {
                pi <- lapply(Q, statdist)
            }
            else {
                pi <- pi/sum(pi)
                pi <- lapply(1:nsim, function(x, y) y, y = pi)
            }
            if (pm)
            printmessage(Reduce("+", Q)/length(Q), pi, method = "mcmc")
            mtrees <- mapply(smap, L = L, Q = Q, pi = pi, logL = logL,
            MoreArgs = list(tree = tree, x = x, N = N, m = m,
            root = root), SIMPLIFY = FALSE)
        }
        else if (is.matrix(Q)) {
            XX <- getPars(bt, xx, model, Q = Q, tree, tol, m)
            L <- XX$L
            logL <- XX$loglik
            if (pi[1] == "equal")
            pi <- setNames(rep(1/m, m), colnames(L))
            else if (pi[1] == "estimated")
            pi <- statdist(Q)
            else pi <- pi/sum(pi)
            if (pm)
            printmessage(Q, pi, method = "fixed")
            mtrees <- replicate(nsim, smap(tree, x, N, m, root,
            L, Q, pi, logL), simplify = FALSE)
        }
        if (length(mtrees) == 1)
        mtrees <- mtrees[[1]]
        else class(mtrees) <- "multiPhylo"
    }
    (if (hasArg(message))
    list(...)$message
    else TRUE)
    if ((if (hasArg(message))
    list(...)$message
    else TRUE) && class(tree) == "phylo")
    message("Done.")
    return(mtrees)
}

drop.clade <- function (tree, tip) {
    nn <- if (!is.null(tree$node.label))
    c(tree$node.label, "NA")
    else "NA"
    tree <- drop.tip(tree, tip, trim.internal = FALSE)
    while (sum(tree$tip.label %in% nn) > 1) tree <- drop.tip(tree,
    tree$tip.label[tree$tip.label %in% nn], trim.internal = FALSE)
    tree
}

paste.tree <- function (tr1, tr2) {
    if (length(tr2$tip) > 1) {
        temp <- tr2$root.edge
        tr2$root.edge <- NULL
        tr1$edge.length[match(which(tr1$tip.label == "NA"), tr1$edge[,
        2])] <- tr1$edge.length[match(which(tr1$tip.label ==
        "NA"), tr1$edge[, 2])] + temp
    }
    tr.bound <- bind.tree(tr1, tr2, where = which(tr1$tip.label ==
    "NA"))
    return(tr.bound)
}

fitMk<-function(tree,x,model="SYM",fixedQ=NULL,...){
    if(hasArg(output.liks)) output.liks<-list(...)$output.liks
    else output.liks<-FALSE
    N<-Ntip(tree)
    M<-tree$Nnode
    if(is.matrix(x)){
        x<-x[tree$tip.label,]
        m<-ncol(x)
        states<-colnames(x)
    } else {
        x<-to.matrix(x,sort(unique(x)))
        m<-ncol(x)
        states<-colnames(x)
    }
    if(hasArg(pi)) pi<-list(...)$pi
    else pi<-"equal"
    if(pi[1]=="equal") pi<-setNames(rep(1/m,m),states)
    else if(pi[1]=="estimated"){
        pi<-if(!is.null(fixedQ)) statdist(fixedQ) else statdist(summary(fitMk(tree,x,model),quiet=TRUE)$Q)
        cat("Using pi estimated from the stationary distribution of Q assuming a flat prior.\npi =\n")
        print(round(pi,6))
        cat("\n")
    }
    else pi<-pi/sum(pi)
    if(is.null(fixedQ)){
        if(is.character(model)){
            rate<-matrix(NA,m,m)
            if(model=="ER") k<-rate[]<-1
            else if(model=="ARD"){
                k<-m*(m-1)
                rate[col(rate)!=row(rate)]<-1:k
            } else if(model=="SYM"){
                k<-m*(m-1)/2
                ii<-col(rate)<row(rate)
                rate[ii]<-1:k
                rate<-t(rate)
                rate[ii]<-1:k
            }
        } else {
            if(ncol(model)!=nrow(model))
            stop("model is not a square matrix")
            if(ncol(model)!=ncol(x))
            stop("model does not have the right number of columns")
            rate<-model
            k<-max(rate)
        }
        Q<-matrix(0,m,m)
    } else {
        rate<-matrix(NA,m,m)
        k<-m*(m-1)
        rate[col(rate)!=row(rate)]<-1:k
        Q<-fixedQ
    }
    index.matrix<-rate
    tmp<-cbind(1:m,1:m)
    rate[tmp]<-0
    rate[rate==0]<-k+1
    liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
    pw<-reorder(tree,"pruningwise")
    lik<-function(pp,output.liks=FALSE,pi){
        if(any(is.nan(pp))||any(is.infinite(pp))) return(1e50)
        comp<-vector(length=N+M,mode="numeric")
        Q[]<-c(pp,0)[rate]
        diag(Q)<--rowSums(Q)
        parents<-unique(pw$edge[,1])
        root<-min(parents)
        for(i in 1:length(parents)){
            anc<-parents[i]
            ii<-which(pw$edge[,1]==parents[i])
            desc<-pw$edge[ii,2]
            el<-pw$edge.length[ii]
            v<-vector(length=length(desc),mode="list")
            for(j in 1:length(v))
            v[[j]]<-matexpo(Q*el[j])%*%liks[desc[j],]
            vv<-if(anc==root) Reduce('*',v)[,1]*pi else Reduce('*',v)[,1]
            comp[anc]<-sum(vv)
            liks[anc,]<-vv/comp[anc]
        }
        if(output.liks) return(liks[1:M+N,])
        logL<--sum(log(comp[1:M+N]))
        return(if(is.na(logL)) Inf else logL)
    }
    if(is.null(fixedQ)){
        fit<-nlminb(rep(0.1,k),function(p) lik(p,pi=pi),lower=rep(0,k),upper=rep(1e50,k))
        obj<-list(logLik=-fit$objective,
        rates=fit$par,
        index.matrix=index.matrix,
        states=states,
        pi=pi)
        if(output.liks) obj$lik.anc<-lik(obj$rates,TRUE,pi=pi)
    } else {
        fit<-lik(Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],pi=pi)
        obj<-list(logLik=-fit,
        rates=Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],
        index.matrix=index.matrix,
        states=states,
        pi=pi)
        if(output.liks) obj$lik.anc<-lik(obj$rates,TRUE,pi=pi)
    }
    class(obj)<-"fitMk"
    return(obj)
}

## print method for objects of class "fitMk"
print.fitMk<-function(x,digits=6,...){
    cat("Object of class \"fitMk\".\n\n")
    cat("Fitted (or set) value of Q:\n")
    Q<-matrix(NA,length(x$states),length(x$states))
    Q[]<-x$rates[x$index.matrix]
    diag(Q)<-0
    diag(Q)<--rowSums(Q)
    colnames(Q)<-rownames(Q)<-x$states
    print(round(Q,digits))
    cat("\nFitted (or set) value of pi:\n")
    print(x$pi)
    cat(paste("\nLog-likelihood:",round(x$logLik,digits),"\n\n"))
}

## summary method for objects of class "fitMk"
summary.fitMk<-function(object,...){
    if(hasArg(digits)) digits<-list(...)$digits
    else digits<-6
    if(hasArg(quiet)) quiet<-list(...)$quiet
    else quiet<-FALSE
    if(!quiet) cat("Fitted (or set) value of Q:\n")
    Q<-matrix(NA,length(object$states),length(object$states))
    Q[]<-object$rates[object$index.matrix]
    diag(Q)<-0
    diag(Q)<--rowSums(Q)
    colnames(Q)<-rownames(Q)<-object$states
    if(!quiet) print(round(Q,digits))
    if(!quiet) cat(paste("\nLog-likelihood:",round(object$logLik,digits),"\n\n"))
    invisible(list(Q=Q,logLik=object$logLik))
}

## logLik method for objects of class "fitMk"
logLik.fitMk<-function(object,...) object$logLik

## AIC method
AIC.fitMk<-function(object,...,k=2){
    np<-length(object$rates)
    -2*logLik(object)+np*k
}

apeAce<-function(tree,x,model,fixedQ=NULL,...){
    if(hasArg(output.liks)){
        output.liks<-list(...)$output.liks
        return(fitMk(tree,x,model,fixedQ,...))
    } else {
        output.liks<-TRUE
        return(fitMk(tree,x,model,fixedQ,output.liks=TRUE,...))
    }
}

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
