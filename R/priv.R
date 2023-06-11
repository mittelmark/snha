# private methods
# just inside inside the snha package
# x = adjacency matrix, factor strength
# create data using a precision matrix
Snha_pmatrix <- function (x,factor=1,mu=100,n=100) {
    Adj2cov = function (x,factor=1) {
        x[x!=0]=x[x!=0]+factor
        diag(x)=apply(x,1,sum)+factor
        return(x)
    }
    Cov2data = function (x,mu=rep(0,nrow(x)),n=100) {
        data=MASS::mvrnorm(n,mu=mu,Sigma=x)
        colnames(data)=colnames(x)
        return(data)
    }
    
    Pcor2cor = function (x, tol = sqrt(.Machine$double.eps)) {
        x = -x
        rn=rownames(x)
        diag(x) = -diag(x)
        x = MASS::ginv(x, tol = tol)
        rownames(x)=colnames(x)=rn
        return(stats::cov2cor(x))
    }
    y=Snha_d2u(x)
    COV=Adj2cov(y,factor=factor)
    COR=Pcor2cor(COV)
    data=Cov2data(COR,mu=rep(mu,nrow(x)),n=n)
    return(t(data))
}

# private functions

# convert a directed to an undirected graph
Snha_d2u <- function (g) {
    g[lower.tri(g)]=g[lower.tri(g)]+t(g)[lower.tri(g)]
    g[upper.tri(g)]=g[upper.tri(g)]+t(g)[upper.tri(g)]    
    g[g>0]=1
    g[g<0]=-1
    return(g)
}

# remove duplicated chains
ReduceChains = function (g) {
    ichains=c()
    chs=sort(names(g$chains))
    for (cho in chs) {
        if (length(g$chains[[cho]]) == 1) {
            next
        }
        co=g$chains[[cho]]
        for (chi in chs) {
            if (length(g$chains[[chi]])==1) {
                next
                
            }
            if (length(g$chains[[cho]]) == 1) {
                next
            }

            if (chi == cho) { next }
            ci=g$chains[[chi]]
            if (length(co)> length(ci)) { next }
            if (length(co) == length(ci)) {
                if (all(co==ci) | all(co==rev(ci))) {
                    ichains=c(ichains,chi)
                    g$chains[[chi]] = ''
                }
            } else {
                cop=paste(co,collapse='')
                cipo=paste(ci,collapse='')
                cipr=paste(rev(ci),collapse='')
                if (grepl(cop,cipo,fixed=TRUE)) {
                    ichains=c(ichains,cho)
                    g$chains[[cho]] = ''

                } else if (grepl(cop,cipr,fixed=TRUE)) {
                    ichains=c(ichains,cho)
                    g$chains[[cho]] = ''
                }
            }
            
        }
    }
    for (ch in ichains) {
        g$chains[[ch]]=NULL
    }
    return(g)
}


## return the component ids for an adjacency matrix
Components = function (A) {
    A=as.matrix(A)
    A=A+t(A)
    A[A>0]=1
    comp=c()
    P=Snha_shortest_paths(A)
    nodes=rownames(A)
    x=1
    while (length(nodes) > 0) {
        n=nodes[1]
        idx=which(P[n,] < Inf)
        ncomp=rep(x,length(idx))
        names(ncomp)=rownames(P)[idx]
        comp=c(comp,ncomp)
        nodes=setdiff(nodes,rownames(P)[idx])
        x=x+1
    }
    return(comp[rownames(A)])
}

# add edges to connect components of a graph
# used for layout mechanism to not 
# have infinite path lengths
ConnectComponents = function (A) {
    A=as.matrix(A)
    A=A+t(A)
    A[A>0]=1
    P=Snha_shortest_paths(A)
    if (!any(P==Inf)) {
        return(A)
    }
    comp=Components(A)
    nodes=c()
    tab=table(comp)
    for (n in names(tab)) {
        c=names(which(comp==n))
        if (tab[[n]] > 2) {
            Am=A[c,c]
            # todo min
            deg=apply(Am,1,sum)
            idx=which(deg>0)
            minval=min(deg[idx])
            idx=which(deg == minval)[1]
            node=c[idx]
        } else {
            node = c[1]
        }
        nodes=c(nodes,node)
    }
    A[nodes,nodes]=1
    diag(A)=0
    return(A)
}

# simple shortest path  bfs search
Snha_shortest_paths = function (A,mode="directed") {
    if (class(A)[1] == "snha") {
        A=A$theta
    }
    if (mode == "undirected") {
        A=A+t(A)
        A[A!=0]=1
    }
    S=A
    S[]=Inf
    diag(S)=0
    x=1
    S[A > 0 & A < Inf]=1
    while (TRUE) { 
        flag = FALSE 
        for (m in 1:nrow(S)) {
            ns=which(S[m,] == x)
            for (n in ns) {
                for (o in which(A[n,]==1)) {
                    if (o != m) {
                        flag = TRUE
                        if (S[m,o] > x + 1) {
                            S[m,o]=x+1
                            if (mode == "undirected") {
                                S[o,m]=x+1
                            }
                        }
                    }
                }
            }
        }
        if (!flag) {
            break
        }
        x=x+1
    }
    return(S)
}

# creating a simple correlation plot for snha object

# moved to snha.R as public function
