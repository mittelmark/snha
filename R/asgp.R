# These are private functions which do the real work of creating the 
# Association chains and then later the graph
Asgp = new.env()


Asgp$getChains = function (cor.mt, square=TRUE,top=10, 
                           threshold=0.01, maxl=3) {
  
  getMiddleChain = function (forchain,mt2) {
    res=c()
    for (fi in 2:(length(forchain)-1)) {
      fl=forchain[fi]
      found=FALSE
      for (si in (fi+1):length(forchain)) {
        sl=forchain[si]
        mt2=mt2[order(mt2[,sl],decreasing = TRUE),]
        schain=rev(rownames(mt2))
        mt2=mt2[order(mt2[,fl],decreasing = TRUE),]                    
        fchain=rownames(mt2)
        if (identical(fchain,schain)) {
          # collect results
          res=fchain        
          found = TRUE
          break
        }
      }
      if (found) { break }
    }
    return(res)
  }
  # remove non-correlated values
  #cor.mt[abs(cor.mt)<0.05]=0
  cor.mt=abs(cor.mt)
  if (square) {
    cor.mt=cor.mt^2
  } else {
    threshold=0.1
  }
  if (top>nrow(cor.mt)) {
    top=nrow(cor.mt)
  }
  results=list()
  chained=list()
  for (i in 1:ncol(cor.mt)) {
    # for each node create a r-ordered list of nodes
    node=colnames(cor.mt)[i]
    mt=cor.mt[order(cor.mt[,i],decreasing = TRUE),][1:top,]
    # ignore nodes without any correlation to any other node
    # might be ignored
    # initial full chain
    chain=rownames(mt)
    cor.mt2=cor.mt[chain,chain]
    j = length(chain)
    lchain=chain[j]
    mt=cor.mt2[order(cor.mt2[,lchain],decreasing = TRUE),]
    # reverse chain
    revchain=rev(rownames(mt))
    forchain=chain
    # as long reversed revchain and chain are different
    # shorten  until the same chain ist found in both
    # or it is to short
    # direct chains
    l=length(forchain)
    while (l>maxl) {
      l=l-1
      lchain=chain[l]
      forchain=forchain[1:l]
      mt2=cor.mt2[forchain,forchain]
      mt2=mt2[order(mt2[,lchain],decreasing = TRUE),]
      if (abs(mt2[nrow(mt2),lchain])<threshold) {
        next
      }
      # revchain is shorter automatically
      # as we have reduced mt2
      revchain=rev(rownames(mt2))
      if (identical(revchain,forchain)) {
        # collect results
        key=paste("a-chain-",node,sep="")
        results[[key]]=forchain        
        break
      }
      # otherwise
      # new: try all chains with node inside, not only 
      # in the beginning as above
      res=getMiddleChain(forchain,mt2)
      if (length(res) > 0) {
        key=paste("m-chain-",node,sep="")
        results[[key]]=res
        break ;# todo remove break??
        }
    }
  }
  return(results);
}

Asgp$chains2edgelists = function (chainlist) {
    self=Asgp
    edgelists=c()
    for (name in names(chainlist)) {
        for (i in seq(1,length(chainlist[[name]])-1,by=1)) {
            edgelists=c(edgelists,paste(chainlist[[name]][i],chainlist[[name]][i+1],sep="--"))
        }
    }
    return(edgelists)
}
Asgp$data2chainGraph = function (data,method='pearson',
                                square=TRUE,
                                threshold=0.01,maxl=3,top=10,
                                p.adjust='none',
                                alpha=0.01,
                                cor.p.value=NULL) {
    self=Asgp
    if (nrow(data)==ncol(data) &&  
        length(which(colnames(data)==rownames(data))) == nrow(data)) {
        if (max(data)>1) {
            stop("symmetric matrices must be correlation matrices")
        }
        cormt=data
        if (length(cor.p.value)==0) {
            # pseudo p-values
            # will add 1's later
            cor.p.value=matrix(0,nrow=nrow(data),ncol=nrow(data))
        }
        
    } else {
        cormt=cor(data,method=method,use='pairwise.complete.obs')
        if (method == "kendall") {
            square=TRUE
            cormt=cormt
        }
        cor.p.value=self$corTest(data,method=method,p.adjust=p.adjust)$p.value 
    }
    # if all pairs are NA
    cor.p.value[is.na(cormt)]=1
    cormt[is.na(cormt)]=0
    chains=self$getChains(cormt,square=square,threshold=threshold,
                          maxl=maxl,top=top)
    edgelist=self$chains2edgelists(chains)
    A=matrix(0,nrow=ncol(data),ncol=ncol(data))
    colnames(A)=rownames(A)=colnames(data)
    for (edge in edgelist) {
        nodes=strsplit(edge,"--")[[1]]
        A[nodes[1],nodes[2]]=A[nodes[2],nodes[1]]=1
    }
    if (max(cor.p.value)==0) {
        # cor matrix was given
        # set all non edges to not significant
        # as we can't calculate p-values
        cor.p.value[A==0]=1
    }
    if (alpha < 1) {
        # remove non signif edges
        A=self$removeNonsignifGraphEdges(A,cor.p.value,alpha=alpha)
    }
    asgm=list(theta=A,p.values=cor.p.value,data=data,sigma=cormt,
              method=method,threshold=threshold,alpha=alpha,chains=chains,data=data)
    class(asgm)="snha"
    return(asgm)
}

Asgp$corTest = function (data,method='pearson',p.adjust='none') {
    self=Asgp
        
    cor.p.values <- function(r, n) {
        df <- n - 2
        STATISTIC <- c(sqrt(df) * r / sqrt(1 - r^2))
        p <- pt(STATISTIC, df)
        return(2 * pmin(p, 1 - p))
    }
    
    
    pmatrix = function (M,method="pearson") {
        P=matrix(0,ncol=ncol(M),nrow=ncol(M))
        
        rownames(P)=colnames(P)=colnames(M)
        if (method=="spearman") {
            M=apply(M,2, rank,na.last="keep")
        }
        r=cor(M,use="pairwise.complete.obs")
        ncd=ncol(M)
        #N=unlist(sapply(1:(ncd-1), function(i) sapply((i+1):ncd, function(j) nrow(na.omit(M[,c(i,j)])))) )
        N=c()
        for (i in 1:(ncd-1)) {
            for (j in (i+1):ncd) {
                N=c(N,nrow(na.omit(M[,c(i,j)])))
            }
        }
        P[upper.tri(P)]=apply(cbind(r=r[upper.tri(r)],N=N),1, function (x) cor.p.values(x[1],x[2]))
        P[lower.tri(P)]=t(P)[lower.tri(P)]
        return(list(r=r,P=P))
    }

    if (method=="pearson" | method == "spearman") {
        res=pmatrix(data,method=method)
        cor.mt=res$r
        p.value=res$P
        diag(p.value)=0
        diag(cor.mt)=1
    } else {
        cor.mt=cor(data,method=method,use='pairwise.complete.obs')
        p.value=matrix(0,nrow=ncol(data),ncol=ncol(data))
        ncd=ncol(data)
        for (i in 1:(ncd-1)) {
            for (j in (i+1):ncd) {
                p.value[i,j]=p.value[j,i]=cor.test(data[,i],data[,j],method=method,exact=FALSE)$p.value
            }
        }
    }
    rownames(p.value)=colnames(p.value)=colnames(data)
    if (p.adjust != "none") {
        p.adj=p.adjust(p.value[upper.tri(p.value)],method=p.adjust)
        p.value[upper.tri(p.value)]=p.adj
        p.value[lower.tri(p.value)]=t(p.value)[lower.tri(p.value)]
        
    }
    
    return(list(r=cor.mt,p.value=p.value,p.adjust=p.adjust))
}

Asgp$removeNonsignifGraphEdges = function (A,p.value,alpha=0.05,...) {
    self=Asgp
    if (sum(A)==0) {
        return(A)
    }
    dels=c()
    for (i in 1:(ncol(A)-1)) {
        for (j in i:(ncol(A))) {
            if (p.value[i,j] > alpha) {
                A[i,j]=A[j,i]=0
            }
        } 
    }
    return(A)
}
# TODO: - really required?
Asgp$new = function (data) {
    self=Asgp
}

