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
#'       'hubs', 'nicolas','random', 'regular' or 'werner', default: 'random'}
#' \item{mode}{either 'undirected' or 'directed', default: 'directed'}
#' \item{nodes}{number of nodes for a given type, default: 10}
#' \item{edges}{number of edges for a given type, not used for type "barabasi", default: 12}
#' \item{m}{number of edges added in each iteration of type is 'barabasi', can be values such as 1,2 or fractionla numbers such as 1.3 which means 70 percent single chain addition and 30 percent  double chain additions, default: 1}
#' \item{k}{the degree for a regular graph, or the number of clusters for a cluster graph or the number of hubs for a hubs graph, default: 2}
#' \item{p}{the probabilty for edges if type is 'gnp',default=NULL}
#' \item{power}{the power for preferential attachment if type is 'barabasi', 
#'       1 is linear preferential attachment, below one is sublinear, smaller hubs, 0 is no hubs,
#'       above 1 is super linear with super hubs, default: 1}
#' }
#' \keyword{network}
#' \keyword{graph}
#' \value{mgraph object, which contains an adjacency matrix and optional an attribute for a layout}
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
#' D=mgraph(type='barabasi',m=1.3)
#' D
#' N=mgraph(type='nicolas')
#' ## Nicolas graph brings its own layout
#' plot(N,layout=attr(N,'layout'))
#' }
#' \seealso{ \link[snha:snha]{snha}, \link[snha:plot.snha]{plot.snha}}
#' 
mgraph <- function (x=NULL,type="random",mode="directed",nodes=10,edges=12,m=1,k=3,p=NULL,power=1) {
    types=c("angie","band","barabasi","circle","cluster","gnp","hubs","nicolas","random","regular","werner")
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
                if (floor(m)==n) {
                    sel=n-1
                } else {
                    if (floor(m)==m) {
                        sel=m
                    } else {
                        low=floor(m)
                        high=ceiling(m)
                        probs=c((high-m),(m-low))
                        sel=sample(c(low,high),1,prob=probs)
                    }
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
        } else  if (type == "nicolas") {
            mt=matrix(0,nrow=5,ncol=5)
            rownames(mt)=colnames(mt)=LETTERS[1:5]
            mt['A',c('B','E')]=1
            mt['B',c('C','E')]=1
            mt['C','D']=1
            mt['D',c('A','B')]=1
            mt['E','D']=1
            lay=matrix(c(0,0, 0,1, 0.5,1.75, 1,1, 1,0),
                       ncol=2,byrow=TRUE)
            rownames(lay)=rownames(mt)
            colnames(lay)=c("x","y")
            attr(mt,'layout')=lay
            x=mt
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

#' \name{mgraph_harmonic_centrality}
#' \alias{mgraph_harmonic_centrality}
#' \title{return the harmonic centraility for each node of a graph}
#' \description{
#'     This function returns the harmonic centrality for every node, of a graph
#'     directed graphs are converted to unidirecred graphs before and then the following
#'     formula is applied: \code{h = sum(1/d) / (n-1)} where d are the shortest paths
#'     of a node to all other nodes of the graph and n ist the number of nodes in a graph.
#' }
#' \usage{mgraph_harmonic_centrality(
#'   g,
#'   norm=TRUE)
#' }
#' \arguments{
#' \item{g}{a mgraph object}
#' \item{norm}{should the harmonic centrality normalized over the number of nodes, default: TRUE}
#  }
#' \references{
#' \itemize{
#' \item Marchiori, M., & Latora, V. (2000). Harmony in the small-world. Physica A: Statistical Mechanics and its Applications, 285(3-4), 539-546.
#' }
#' }
#' \examples{
#'  A=mgraph(type="regular",nodes=8,k=3)
#'  plot(A,layout="sam")
#'  mgraph_degree(A)
#'  mgraph_harmonic_centrality(A)
#'  B=mgraph(type="angie",nodes=20,edges=30)
#'  plot(B,layout="sam")
#'  mgraph_degree(B)
#'  mgraph_harmonic_centrality(B)
#'  cols=mgraph_nodeColors(B,type="harmonic")
#'  plot(B,layout="sam",vertex.color=cols)
#'  lg=mgraph_nodeColors(B,type="harmonic",legend=TRUE)
#'  legend("bottom",fill=lg$col, legend=lg$legend,ncol=5)
#' }
#' 

mgraph_harmonic_centrality <- function (g,norm=TRUE) {
    sp=Snha_shortest_paths(g,mode="undirected")
    diag(sp)=Inf
    if (norm) {
        res=apply(sp,1,function (x) { (sum(1/x))/(length(x)-1) })
    } else {
        res=apply(sp,1,function (x) { (sum(1/x)) })
    }        
    
    names(res)=rownames(g)
    return(res)
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
#'     col=c("skyblue","grey80","salmon"),
#'     type="inout",
#'     legend=FALSE
#' )}
#' \arguments{
#' \item{g}{a mgraph object or an adjacency matrix or an adjacency list}
#' \item{col}{default colors for nodes with only incoming, in- and outgoing and only outgoing edges, default: c("skyblue","grey80","salmon")}
#' \item{type}{type of color selection, either "inout" or "harmonic" (centrality), default: "inout"}
#' \item{legend}{should be colors and legend labels returned instead of colors, default: FALSE}
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
#' cols=mgraph_nodeColors(A,type="harmonic")
#' plot(A, layout="star",vertex.color=cols)
#' lg=mgraph_nodeColors(A,type="harmonic",legend=TRUE)
#' legend("bottom",fill=lg$col, legend=lg$legend)
#' }
#'

mgraph_nodeColors <- function (g,col=c("skyblue","grey80","salmon"),type="inout",legend=FALSE) {
    if (type == "inout") {
        if (legend) {
            return(list(legend=c("out","inout","in"),
                 col=col))
        }
        colors = rep(col[2],nrow(g))
        out    = mgraph_degree(g,mode="out")
        inc     = mgraph_degree(g,mode="in")    
        colors[out >  0 & inc == 0] = col[3]
        colors[out == 0 & inc >  0] = col[1]    
        return(colors)
    } else if (type == "harmonic") {
        cols=c("#eeeeee","#eebbbb","#ee9999","#ee7777","#ee5555","#ee3333")
        h = mgraph_harmonic_centrality(g)
        c=cut(h,breaks=c(0,0.01,0.2,0.4,0.6,0.8,1),include.lowest=TRUE)
        if (legend) {
            return(list(legend=levels(c),
                   col=cols))
        }
        lvls=as.character(c)
        c=as.numeric(c)
        cs = cols[c]
        names(cs)=lvls
        return(cs)
    } else {
        stop("Error: type must be either 'inout' or 'harmonic'!")
    }
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
            
            comps=Components(A)
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


#' \name{mgraph_accuracy}
#' \alias{mgraph_accuracy}
#' \title{quality measures for a predicted graph}
#' \description{
#'    The function `mgraph_accuracy` measures the prediction
#'    quality of a predicted graph in comparison to a known true one.
#'    The graph comparisons are done on the basis of undirected graphs only.
#'    Directed graphs are converted to undirected internally.
#' 
#'    The most used accuracy measure to balance the different measures is the F-score where:
#' 
#' \itemize{
#'    \item the F0.5 measure (beta=0.5) gives more weight on precision, less weight on sensitivity/recall.
#'    \item the F1 measure (beta=1.0) balances the weight on precision and sensitivity/recall.
#'    \item the F2 measure (beta=2.0) gives less weight on precision, more weight on sensitivity/recall
#'   }
#' }
#' \usage{mgraph_accuracy(g.true,g.pred)}
#' \arguments{
#' \item{g.true}{snha, mgraph or or adjacency matrix of a the true graph}
#' \item{g.pred}{snha, mgraph or or adjacency matrix of a the predicted graph}
#' }
#' \value{returns list object with various accuracy measures, such as Sens(itivity), Spec(ificity), Prec(ision), Acc(uracy), Non-information-rate (NIR), balanced classification rate (BCR), F1 measure and Mathhews correlation coefficient (MCC)}
#' \examples{
#' ang=mgraph(type="angie",nodes=12,edges=18)
#' data=snha_graph2data(ang)
#' pred=snha(t(data))$theta
#' round(unlist(mgraph_accuracy(ang,pred)),2)
#' #  using vectors
#' tvec=c(1,1,1,0,0,0)
#' pvec=c(1,1,0,1,0,0)
#' table(tvec,pvec)
#' mgraph_accuracy(tvec,pvec)
#' }
#' 

mgraph_accuracy <- function (g.true,g.pred) {
    # Convert to undirected
    fbeta = function (prec,sens,beta=1) {
        return(((1+beta^2)*prec*sens)/(beta^2*prec+sens))
    }
    if (!is.vector(g.true)) {
        g.true <- mgraph_d2u(g.true)
        g.pred <- mgraph_d2u(g.pred)
        g.true <- g.true[upper.tri(g.true)]
        g.pred <- g.pred[upper.tri(g.pred)]
    }
    # Convert explictly to 'double' to avoid integer overflow
    TP=as.double(length(which(g.true != 0 & g.pred != 0)))
    FP=as.double(length(which(g.true == 0 & g.pred != 0)))
    TN=as.double(length(which(g.true == 0 & g.pred == 0)))
    FN=as.double(length(which(g.true != 0 & g.pred == 0)))
    # Included NA check
    Acc  = if ( (TP+TN+FP+FN) != 0 ) { (TP+TN)/(TP+TN+FP+FN) } else { NaN }
    Sens = if ( (TP+FN) != 0 ) { TP / (TP+FN) } else { NaN }
    Spec = if ( (TN+FP) != 0 ) { TN / (TN+FP) } else { NaN }
    BCR  =  (Sens+Spec)/2
    Prec = if ( (TP+FP) != 0 ) { TP / (TP+FP) } else { NaN }
    NIR  = max(prop.table(table(g.true)))
    pt   = stats::prop.test(sum(diag(table(g.true,g.pred))),length(g.true),p=NIR,alternative="greater")
    Acc.conf.int=I=pt$conf.int
    Acc.p.value=pt$p.value
    F1= fbeta(Prec,Sens,beta=1.0)
    F2= fbeta(Prec,Sens,beta=2.0)
    F0.5= fbeta(Prec,Sens,beta=0.5)
    MCC=if ( sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) != 0) { (TP*TN-FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN)) } else { NaN }
    norm_MCC=(MCC+1)/2
    return(list(TP=TP,FP=FP,TN=TN,FN=FN,Sens=Sens,Spec=Spec, Prec=Prec,
                BCR=BCR,F1=F1,F2=F2,F0.5=F0.5,
                MCC=MCC, norm_MCC=norm_MCC, 
                Acc=Acc, NIR=NIR, Acc.p.value=pt$p.value))
}

#' \name{mgraph_nd}
#' \alias{mgraph_nd}
#' \title{network deconvolved data matrix using the algorithm of Feizi et. al. 2013}
#' \description{
#'    This function is a R implementation for "A General Method to Distinguish Direct Dependencies over Networks" by Feizi et. al (2013). For the Matlab code look here:
#'   \url{http://compbio.mit.edu/nd/code/ND.m}
#' }
#' \usage{mgraph_nd(x,beta=0.99,alpha=1,control=FALSE)}
#' \arguments{
#' \item{x}{symmetric relevance matrix, where high similarity between nodes is expressed with high values such as in an correlation matrix}
#' \item{beta}{scaling paramater, the largest absolute eigenvalue is mapped to beta, values should be between 0 and 1, default: 0.99}
#' \item{alpha}{fraction of edges of the observed dependency matrix to be kept, should be between 0 (none) and 1 (all), default: 1}
#' \item{control}{if FALSE only direct weights are displayed, if TRUE also non-observed interactions are displayed, default: FALSE}
#' }
#' \value{returns matrix with deconvoluted relevance matrix}
#' \examples{
#' W=mgraph(type="werner")
#' W
#' data=snha_graph2data(W)
#' C=cor(t(data))
#' round(C,2)
#' round(mgraph_nd(C,beta=0.9,alpha=0.3),2)
#' # manually prepare the deconvoluted matrix
#' D=mgraph_nd(C,beta=0.9,alpha=0.3)
#' diag(D)=1
#' par(mfrow=c(2,2))
#' plot(W)
#' plot(snha(t(data)))
#' plot(snha(t(data),nd=TRUE))
#' plot(snha(D))
#' }
#' \references{
#' \itemize{
#' \item Feizi, S., Marbach, D., Medard, M., & Kellis, M. (2013). Network deconvolution as a general method to distinguish direct dependencies in networks.
#'     Nature biotechnology, 31(8), 726-733. DOI 10.1038/nbt.2635
#' }
#' }
#' \author{
#' \itemize{
#' \item @2013 KELLIS-LAB, Soheil Feizi, Matlab code
#' \item @2021 Detlef Groth, University of Potsdam, R code
#' }
#' }
#' 

mgraph_nd = function (x,beta=0.99,alpha=1,control=FALSE) {
   n = ncol(x)
   diag(x)=0
   y=stats::quantile(x,1-alpha,type=5)
   mat_th=x
   mat_th[mat_th<y]=0
   mat_th=(mat_th+t(mat_th))/2
   E=eigen(mat_th)
   lam_n=abs(min(E$values))
   lam_p=abs(max(E$values))
   m1=lam_p*(1-beta)/beta;
   m2=lam_n*(1+beta)/beta;
   m=max(m1,m2);
   E$values=E$values/(m+E$values)
   EV=matrix(0,nrow=n,ncol=n)
   diag(EV)=E$values
   mat_new1=E$vectors %*% EV %*% t(E$vectors)
   if (!control) {
       ind_edges = (mat_th>0)*1.0;
       ind_nonedges = (mat_th==0)*1.0;
       m1 = max(x * ind_nonedges); # yes * as checked by Octave Op: .* is like * in R
       m2 = min(mat_new1);
       mat_new2 = (mat_new1+max(m1-m2)) * ind_edges + (x * ind_nonedges);
   } else {
       m2=min(mat_new1)
       mat_new2 = (mat_new1+max(-m2,0));
   }
   m1 = min(mat_new2);
   m2 = max(mat_new2);
   mat_nd = (mat_new2-m1) / (m2-m1);
   rownames(mat_nd)=colnames(mat_nd)=colnames(x)
   return(mat_nd)
}

#' \name{mgraph_lmc}
#' \alias{mgraph_lmc}
#' \title{using linear models to check the create snha graph}
#' \description{
#'    This function checks a given snha graph against its own data using linear
#'    models. Edges are removed if they do not add  more than 2 percent of the 
#'    adjusted  R-square value. The function might be only called directly if the user would like to
#'    change the default values for del and add. 
#'    The function snha uses defaults of 0.02 for del(eting) and 0.05 for add(ing) values.
#' }
#' \usage{mgraph_lmc(x,del=0.02,add=NULL,method='pearson')}
#' \arguments{
#'     \item{x}{snha graph object}
#'     \item{del}{r-square threshold to delete edges, edges not adding more than this value to the linear model to the target node will be deleted, default: 0.02}
#'     \item{add}{r-square threshold to add edges, if NULL no edges will be added, recommended value is 0.05, default: NULL}
#'     \item{method}{The method which was used to generate a graph, if spearman or mutual information were given data are rank transformed, default: 'pearson'}
#' }
#' \value{returns adjacency matrix where edges migh be removed if they are not adding explanation to the model}
#' \examples{
#' set.seed(123)
#' B=mgraph(type="barabasi",m=2,nodes=6)
#' data=t(snha_graph2data(B))
#' aa=snha(data)
#' aa$theta
#' unlist(mgraph_accuracy(B,aa$theta))
#' ab=mgraph_lmc(aa)
#' unlist(mgraph_accuracy(B,ab))
#' ab
#' aa$theta-ab
#' ac=mgraph_lmc(aa,add=0.02)
#' unlist(mgraph_accuracy(B,ac))
#' }
#' 

# linear model check for a snha graph object
# adding or deleting edges if linear model suggests it
#  TODO: rank based approach in case method of snha is spearman

mgraph_lmc <- function (x,del=0.02,add=NULL,method="pearson") {
    g=x
    if (method == "pearson") {
        data=g$data
    } else {
        data=apply(data,2,rank)
    }
    fixGraph <- function (A) {
        for (i in 1:(ncol(A)-1)) {
            for (j in i:ncol(A)) {
                if (A[i,j]<del & A[j,i] < del) {
                    A[i,j]=A[j,i]=0
                }
            }
        }
        A[A>0]=1
        return(A)
    }
    T=g$theta
    S=g$sigma
    U=T
    cnames=colnames(T)
    # addition
    if (!is.null(add)) {
        for (i in 1:ncol(T)) {
            idx = which(T[i,] > 0 | abs(S[i,]) > 0.1)
            idx = setdiff(idx,i)
            if (length(idx)>1) {
                L=mgraph_lms(data[,c(i,idx)])
                # it is slower due to measure more linear models
                # than in deletion mode
                for(j in 1:length(idx)) {
                    if (L[1,j+1] > add & T[i,idx[j]] == 0) {
                        U[i,idx[j]] = L[1,j+1]
                    }
                }
            }
        }
        U[U>0]=1
        U[upper.tri(U)]=mapply(max,U[upper.tri(U)],t(U)[upper.tri(U)])
        U[lower.tri(U)]=t(U)[lower.tri(U)]
        T=U
    }
    for (i in 1:ncol(T)) {
        #idx=which(S[i,]^2 > add | T[i,] > 0)
        idx=which(T[i,] > 0)
        if (length(idx)>1) {
            L=mgraph_lms(data[,c(i,idx)])
            for(j in 1:length(idx)) {
                U[i,idx[j]]= L[1,j+1]
            }
        } else if (length(idx) == 1) {
            U[i,idx]=g$sigma[i,idx]
        }
    }
    X=U
    X[T>0 & X==0] = 0.001
    U=fixGraph(X)
    return(U)
}


#' \name{mgraph_lms}
#' \alias{mgraph_lms}
#' \title{linear model backward variable selection}
#' \description{
#'    This function is a simple version of linear model building using backward
#'    variable selection. At each step the variable with the lowest decrease value
#'    in the adjusted R-square value is removed from the model and the r-square value is stored in the return matrix.
#' }
#' \usage{mgraph_lms(x)}
#' \arguments{
#' \item{x}{a data matrix or data frame}
#' }
#' \value{returns matrix with the adjusted r-square values attributed by the variables}
#' \examples{
#' data(swiss)
#' round(mgraph_lms(swiss)*100,2)
#' W = mgraph(type='werner')
#' data=t(snha_graph2data(W))
#' round(mgraph_lms(data)*100,2)
#' }
#' 

mgraph_lms <- function (x) {
    L=matrix(0,nrow=ncol(x),ncol=ncol(x))
    rownames(L)=colnames(L)=colnames(x)
    cnames=colnames(x)
    for (i in 1:ncol(x)) {
        cn=setdiff(cnames,cnames[i])
        while (length(cn)>1) {
            frm=stats::formula(paste(cnames[i],"~",paste(cn,collapse="+")))
            adjrs.all=summary(lm(frm,data=as.data.frame(x)))$adj.r.squared
            min=1
            min.cn=cn[1]
            for (c in cn) {
                frm=stats::formula(paste(cnames[i],"~",paste(setdiff(cn,c),collapse="+")))
                adjrs.cn=summary(lm(frm,data=as.data.frame(x)))$adj.r.squared
                if (adjrs.all-adjrs.cn<min) {
                    min=adjrs.all-adjrs.cn
                    min.cn=c
                }   
            }    
            cn=setdiff(cn,min.cn)
            L[cnames[i],min.cn]=min
        }
        frm=stats::formula(paste(cnames[i],"~",cn[1]))
        L[cnames[i],cn[1]]=summary(lm(frm,data=as.data.frame(x)))$adj.r.squared
    }
    L[L<0]=0
    return(L)
}


#' \name{mgraph_trf}
#' \alias{mgraph_trf}
#' \title{checking possible chains with 3 nodes for triad structures}
#' \description{
#'    Experimental:
#' 
#'    This function checks a given snha graph for possible triad structures
#'    which were not recognized by the original algorithm.
#'    The method checks all nodes within path length 2 for a possible edge
#'    by comparing the ratio of the correlation coefficients.
#'    The function snha uses defaults of 0.02 for del(eting) and 0.05 for add(ing) values.
#' }
#' \usage{mgraph_trf(x,frac=1.2)}
#' \arguments{
#' \item{x}{snha graph object}
#' \item{frac}{amount which the absolute correlation exceeds in compariuson of the product of the other two correlations in the possible triad, value should be larger than 1, default: 1.2}
#' }
#' \value{returns snha graph object with added triads to chains and to the adjacency matrix}
#' \examples{
#' set.seed(124)
#' W=mgraph(type="werner")
#' data=t(snha_graph2data(W))
#' aa=snha(data)
#' aa$theta
#' unlist(mgraph_accuracy(W,aa$theta))
#' # remove false positive by linear model
#' ab=mgraph_lmc(aa)
#' aa$theta=ab
#' aa=mgraph_trf(aa)
#' unlist(mgraph_accuracy(W,aa$theta))
#' }
#' 

mgraph_trf <- function (x,frac=1.2) {
    C=abs(x$sigma)
    sp=Snha_shortest_paths(x$theta)
    cn=colnames(sp)
    for (i in 1:nrow(sp)) {
        p2=which(sp[i,]==2) 
        for (p in p2) {
            sns=intersect(which(sp[p,]==1),which(sp[i,]==1))
            if (length(sns)>0) {
                for (s in sns) {
                    f=C[i,p]/(C[i,s]*C[p,s])
                    if (C[i,p]>0.1 & x$p.values[i,p]<0.1 & f > frac) {
                        x$theta[i,p]=1
                        x$theta[p,i]=1
                        key=paste('t-chain',cn[i],cn[p],cn[s],sep="-")
                        x$chains[[key]]=paste(cn[i],cn[p],cn[s],sep="-")
                    }
                }
            }
        }
    }
    return(x)
}
