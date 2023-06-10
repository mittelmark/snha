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

