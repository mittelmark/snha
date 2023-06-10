#' \name{mgraph}
#' \alias{mgraph}
#' \title{create a mgraph object, an adjacency matrix of a specific type}
#' \description{
#'     This is a function to create a few different graph types such as barabasi, random, cluster etc to
#'     test the snha algorithms.
#' }
#' \usage{mgraph(
#'   x=NULL,
#'   type="random",
#'   mode="directed",
#'   nodes=10,
#'   edges=12,
#'   m=1,
#'   k=3,
#'   p=NULL,
#'   power=1)
#' }
#' \arguments{
#' \item{x}{either a adjacency matrix or an adjacency list, if given type is not used}
#' \item{type}{graph type, one of 'angie', 'band', 'barabasi', 'circle', 'cluster',
#'       'hubs', 'random', 'regular' or 'werner', default: 'random'}
#' \item{mode}{either 'undirected' or 'directed', default: 'directed'}
#' \item{nodes}{number of nodes for a given type, default: 10}
#' \item{edges}{number of edges for a given type, not used for type "barabasi", default: 12}
#' \item{m}{number of edges added in each iteration of type is 'barabasi', default: 1}
#' \item{k}{the degree for a regular graph, or the number of clusters for a cluster graph or the number of hubs for a hubs graph, default: 2}
#' \item{p}{the probabilty for edges if type is 'gnp',default=NULL}
#' \item{power}{the power for preferential attachment if type is 'barabasi', 
#'       1 is linear preferential attachment, below one is sublinear, smaller hubs, 0 is no hubs,
#'       above 1 is super linear with super hubs, default: 1}
#' }
#' \keyword{network}
#' \keyword{graph}
#' \value{mgraph object, which contains an adjacency matrix}
#' \examples{
#' M <- matrix(0,nrow=7,ncol=7)
#' rownames(M)=colnames(M)=LETTERS[1:7]
#' M[c('A','B'),'C']=1
#' M['C','D']=1
#' M['D',c('E','F')]=1
#' M['E','F']=1
#' G=mgraph(M)
#' set.seed(125)
#' R = mgraph(type="random",nodes=8,edges=9)
#' summary(R)
#' A = mgraph(type="angie",nodes=8,edges=9)
#' A
#' C = mgraph(type="circle",nodes=8)
#' B = mgraph(type="band",nodes=8)
#' W=mgraph(type='werner')
#' W
#' B=mgraph(type='barabasi')
#' B
#' }
#' \seealso{ \link[snha:snha]{snha}, \link[snha:plot.snha]{plot.snha}}
#' 
mgraph <- function (x=NULL,type="random",mode="directed",nodes=10,edges=12,m=1,k=3,p=NULL,power=1) {
    types=c("angie","band","barabasi","circle","cluster","gnp","hubs","random","regular","werner")
    i=pmatch(type,types)
    if (is.na(i)) {
        stop(paste("Unkown type: Known types are: '",paste(types,collapse="', '"),"'!",sep=""))
    }
    type=types[i]
    if (!is.null(x[1])) {
        if (is.matrix(x) & nrow(x)==ncol(x) & all(rownames(x)== colnames(x))) {
            class(x)="mgraph"
        }
        if (is.null(rownames(x)[1])) {
            nms=mgraph_autonames(nrow(x))
            rownames(x)=colnames(x)=nms
        }
    } else if (is.null(x[1])) {
        if (is.null(type)) {
            stop("Error: Either a matrix or a type must be given for mgraph!")
        }
        x=matrix(0,nrow=nodes,ncol=nodes)
        nms=mgraph_autonames(nodes)
        rownames(x)=colnames(x)=nms
        if (type == "random") {
            v=1:length(x[upper.tri(x)])
            idx=sample(v,edges)
            x[upper.tri(x)][idx]=1
        } else if (type %in% c("band","circle")) {
            for (i in 1:((nrow(x)-1))) {
                x[i,i+1]=1
            }
            if (type == "circle") {
                x[nrow(x),1]=1
            }
        } else if (type == "barabasi") {
            x[2, 1]=1
            for (n in 3:ncol(x)) {
                if (m==n) {
                    sel=n-1
                } else {
                    sel=m
                }
                deg=mgraph_degree(x)^power
                # preferential attachment to nodes with higher degree
                idx=sample(1:(n-1),sel,prob=deg[1:(n-1)]/sum(deg[1:(n-1)]))
                x[idx,n]=1
            }
        } else if (type == "gnp") {
            if (is.null(p[1])) {
                stop("Error: For graphs of type 'gnp' you must give the edge probabilty 'p'!")
            } 
            if (p < 0 | p > 1) {
                stop("Error: p must be within 0 and 1!")
            }
            x[]=stats::rbinom(length(x),1,p=p)
            diag(x)=0
        } else if (type == "angie") {
            nds=rownames(x)
            done=nds[1]
            nds=nds[-1]
            while (length(nds)>0) {
                tar=sample(done,1)
                x[tar,nds[1]]=1
                done=c(done,nds[1])
                nds=nds[-1]
            }
            while (sum(x) < edges) {
                idx=sample(which(x[upper.tri(x)]==0),1)
                x[upper.tri(x)][idx]=1
            }
        } else if (type == "werner") {
            x=x[1:6,1:6]
            x[c(1,2),3]=1
            x[3,4]=1
            x[4,5:6]=1
            x[5,6]=1
        } else if (type %in% c("cluster","hubs")) {
            size = nodes %/% k
            rem  = nodes %%  k
            esize= edges %/% k
            ree  = edges %%  k
            pos=0
            for (i in 1:k) {
                s=size
                e=esize
                if (rem>0) {
                    s=size+1; rem=rem-1
                }
                if (ree>0) {
                    e=e+1; ree=ree-1
                }
                if (type == "cluster") {
                    g=mgraph(type="angie",nodes=s,edges=e)
                } else {
                    g=matrix(c(0,rep(1,s-1),rep(0,s*s-s)),ncol=s,byrow=TRUE)
                }
                x[(pos+1):(pos+nrow(g)),(pos+1):(pos+nrow(g))]=g 
                pos=pos+s

            }
        } else if (type == "regular") {
            if (k %in% c(1,3)) {
                if (nodes %% 2 != 0) {
                    stop("regular graphs with k = 1 or k = 3 must have an even number of nodes")
                }
            }
            if (k == 1) {
                for (i in seq(1,nrow(x)-1,by=2)) {
                    x[i,i+1]=1
                }
            } else if (k <= 3) {
                x=mgraph(type="circle",nodes=nodes)
                if (k == 3) {
                    for (i in 1:(nrow(x)/2)) {
                        x[i,i+(nrow(x)/2)]=1
                     }
                 }
            } else {
                stop("Error: Only values of k <= 3 are implemented for regular graphs")
            }
        }
    } else {
        stop("either a matrix or a type must be given")
    }
    if (mode=="undirected") {
        x=mgraph_d2u(x)
    }
    class(x)="mgraph"
    return(x)
    
}
plot.mgraph=plot.snha

#' \name{mgraph_autonames}
#' \alias{mgraph_autonames}
#' \title{create names for nodes and other data structures}
#' \description{
#'     This function aids in creating standard node labels for graphs.
#' }
#' \usage{mgraph_autonames(
#'   n,
#'   prefix=NULL)
#' }
#' \arguments{
#' \item{n}{how many labels}
#' \item{prefix}{the node prefix, per defaultif not given LETTERS with numbers will be used, default: NULL}
#  }
#' \examples{
#' mgraph_autonames(12)
#' mgraph_autonames(12,LETTERS[1:4])
#' mgraph_autonames(12,"R")
#' }
#' 

mgraph_autonames <- function (n,prefix=NULL) {
    if (is.null(prefix[1])) {
        if (n<=26) {
            nms=LETTERS[1:n]
        } else if (n <= 26*9) {
            nms=mgraph_autonames(n,prefix=LETTERS)
        } else {
            nms=mgraph_autonames(n,prefix="N")
        }
        return(nms)
    } else {
        nms=prefix
        ln <- ceiling(n / length(nms))
        nms_tmp <- rep(nms, ln)[1:n]
        ptf <- rep(1:ln, each=length(nms))[1:n]
        frmt <- formatC(ptf, flag="0", width=log10(ln) + 1)
        return(paste0(nms_tmp, frmt))
    }
}

#' \name{mgraph_degree}
#' \alias{mgraph_degree}
#' \title{return the number of edges for each node}
#' \description{
#'     This function returns the number of edges adjacent to every node,
#'     for directed graphs it will as well return the in and out going edges.
#' }
#' \usage{mgraph_degree(
#'   g,
#'   mode='undirected')
#' }
#' \arguments{
#' \item{g}{a mgraph object}
#' \item{mode}{how should the graph been analyzed, either 'undirected', 'out' or 'in', default: 'undirected'}
#  }
#' \examples{
#' A=mgraph(type="regular",nodes=8,k=3)
#' mgraph_degree(A)
#' mgraph_degree(A,mode="out")
#' mgraph_degree(A,mode="in")
#' }
#' 

mgraph_degree <- function (g,mode="undirected") {
    if (mode == "undirected") {
        g=mgraph_d2u(g)
    } else if (mode == "in") {
        g=t(g)
    } 
    g[abs(g)!=0]=1
    d=apply(g,1,sum)
    names(d)=rownames(g)
    return(d)
}

#' \name{mgraph_d2u}
#' \alias{mgraph_d2u}
#' \title{create an undirected graph out of a directed graph}
#' \description{
#'   This function gets an directed graph and convertes all edges from directed ones to undirected ones.
#'   The number of edges should stay the same, the edge sign (+ or -) stays the same.
#' }
#' \usage{mgraph_d2u(g)}
#' \arguments{
#' \item{g}{a mgraph object or an adjacency matrix}
#  }
#' \examples{
#' A=mgraph(type="angie",nodes=4,edges=4)
#' A
#' U=mgraph_d2u(A)
#' U
#' }
#' 

mgraph_d2u <- function (g) {
    g[lower.tri(g)]=g[lower.tri(g)]+t(g)[lower.tri(g)]
    g[upper.tri(g)]=g[upper.tri(g)]+t(g)[upper.tri(g)]    
    g[g>0]=1
    g[g<0]=-1
    return(g)
}

#' \name{mgraph_nodeColors}
#' \alias{mgraph_nodeColors}
#' \title{create node colors for directed graphs}
#' \description{
#'   This function simplifies automatic color coding of nodes for directed graphs.
#'   Nodes will be colored based on their degree properties, based
#'   on their incoming and outcoming edges.
#' }
#' \usage{mgraph_nodeColors(g,
#'     col=c("skyblue","grey80","salmon")
#' )}
#' \arguments{
#' \item{g}{a mgraph object or an adjacency matrix or an adjacency list}
#' \item{col}{default colors for nodes with only incoming, in- and outgoing and only outgoing edges, default: c("skyblue","grey80","salmon")}
#  }
#' \examples{
#' par(mfrow=c(1,2),mai=rep(0.1,4))
#' A=mgraph(type="random",nodes=6,edges=8)
#' cols=mgraph_nodeColors(A)
#' mgraph_degree(A,mode="in")
#' mgraph_degree(A,mode="out")
#' cols
#' plot(A, layout="star")
#' plot(A, layout="star",vertex.color=cols) 
#' }
#'

mgraph_nodeColors <- function (g,col=c("skyblue","grey80","salmon")) {
    colors = rep(col[2],nrow(g))
    out    = mgraph_degree(g,mode="out")
    inc     = mgraph_degree(g,mode="in")    
    colors[out >  0 & inc == 0] = col[3]
    colors[out == 0 & inc >  0] = col[1]    
    return(colors)
}

#' \name{mgraph_u2d}
#' \alias{mgraph_u2d}
#' \title{reate a directed graph out of an undirected graph}
#' \description{
#'   This function creates a directed graph from an undirected one by the given input nodes. Input nodes can be chosen by names or a number for random selection of input nodes will be given.
#' Input nodes will have at least shortest path distance to other input nodes of pathlength two.
#' Selected input nodes will draw in each iteration outgoing edges to other nodes in the nth iteration neighborhood. The input
#' nodes will alternatively select the next edges on the path to not visited nodes. All edges will be only visited onces.
#' }
#' \usage{mgraph_u2d(
#'     g,input=2,negative=0.0,shuffle=FALSE
#' )}
#' \arguments{
#' \item{g}{mgraph object created with mgraph}
#' \item{input}{number or names of input nodes in the graph, if number of input nodes is smaller than number of components, for each component one input node is automatically created}
#' \item{negative}{proportion of inhibitive associations in the network value between 0 and 1 are acceptable, Default 0.0}
#' \item{shuffle}{should just the edge directions beeing shuffled, if TRUE the graph will be very random without a real structure or chains of associations, default: FALSE}
#  }
#' \examples{
#' par(mfrow=c(2,2),mai=rep(0.1,4))
#' G=mgraph(type="angie",nodes=7,edges=9)
#' mgraph_degree(G,mode='in')
#' U=mgraph_d2u(G)
#' H=mgraph_u2d(U,input=c("G","E"))
#' identical(G,H)
#' I=mgraph_u2d(U,shuffle=TRUE)
#' identical(G,I)
#' lay=snha_layout(G)
#' plot(G,layout=lay,vertex.color=mgraph_nodeColors(G))
#' plot(U,layout=lay)
#' plot(H,layout=lay,vertex.color=mgraph_nodeColors(H))
#' plot(I,layout=lay,vertex.color=mgraph_nodeColors(I))
#' }
#' 

mgraph_u2d <- function (g,input=2,negative=0.0,shuffle=FALSE) {
    if (shuffle) {
        h=g
        idx=sample(rownames(h))
        h=h[idx,idx]
        h[lower.tri(h)]=0
        h=h[rownames(g),rownames(g)]
        class(h)="mgraph"
        return(h)          
    } else {
        A=g
        # creates from undirected a directed one
        # A must be an symmetric adjacency matrix
        # MN The new check is 5x faster and takes less memory
        #if (!identical(A[lower.tri(A)], t(A)[lower.tri(A)])) {
        if (!all(A == t(A))) {
            stop("adjacency matrix must be symmetric")
        }
        if (negative < 0 || negative > 1) {
            stop("negative proportions must be within 0 and 1")
        }
        # undirected matrix
        U=A
        # future directed matrix
        D=U
        D[]=0
        neighbours=c()
        if (class(input)=='numeric') {
            nodes=c()
            
            comps=Mgraph_components(A)
            if (max(comps)>1) {
                for (i in 1:max(comps)) {
                    if (length(which(comps==i)) > 1) {
                        node=sample(names(which(comps==i)),1)
                        neighbours=c(neighbours,rownames(U)[which(U[node,]==1)],node)
                        nodes=c(nodes,node)
                    }
                }
            }
            n=input-length(nodes)
            while (n>0) {
                rnames=setdiff(rownames(U),c(neighbours))
                node=sample(rnames,1)
                n=n-1
                neighbours=c(neighbours,rownames(U)[which(U[node,]==1)],node)
                nodes=c(nodes,node)
            }
        } else {
            nodes=input
            n=0
        }
        visits=list()
        while (length(nodes)>0) {
            node=nodes[1]
            edges=which(U[node,]==1)
            newnodes=colnames(U)[edges]
            if (length(nodes)==1) {
                nodes=newnodes
            } else {
                nodes=c(nodes[2:length(nodes)],newnodes)
            }
            D[node,edges]=1
            U[node,]=0
            U[,node]=0
        }
        A[]=D
        if (negative>0) {
            idx=which(A==1)
            n=floor(length(idx)*negative)
            if(n>0) {
                min=sample(idx,n)
                A[min]=-1
            }
        }
        class(A)="mgraph"
        return(A)
    }
}


Mgraph_components <- function (g) {
    A=g
    A=as.matrix(A)
    A=A+t(A)
    A[A>0]=1
    comp=c()
    P=Mgraph_shortest.paths(A)
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

Mgraph_shortest.paths <- function (g,mode="directed") {
    A=g
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
