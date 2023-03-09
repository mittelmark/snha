
LL.lmvnorm = function (x, method="pearson") {
    x=scale(x)
    mean=rep(0,ncol(x))
    sigma=cor(x,method=method,use="pairwise.complete.obs")
    p =   ncol(x)
    d =   chol(sigma)
    b = backsolve(d, t(x) - mean, transpose = TRUE)
    csum = colSums(b^2)
    logretval = -sum(log(diag(d))) - 0.5 * p * log(2 * pi) - 0.5 * csum
    names(logretval) = rownames(x)
    return(logretval)
}
LL.getLL = function (data,cor.method='pearson',cols=c()) {
    #ttll=sum(dmvnorm(scale(data),
    #                  sigma=cor(data,method=cor.method),
    #                  mean=rep(0,ncol(data)),log=TRUE))
    ttll=sum(LL.lmvnorm(data))
    return(ttll)
}

LL.getChainLL = function (data,chain) {
    if (is.character(chain)) {
        # given as colnames
        # extract indices and keep order
        cnames=1:ncol(data)
        names(cnames)=colnames(data)
        chain=cnames[chain]
    }
    rest=setdiff(1:ncol(data),chain) 
    # calculate p-value for chain in matrix
    ll.tot=LL.getLL(data)
    ll.chain=LL.getLL(data[,chain])
    ll.rest=LL.getLL(data[,rest])
    ll.sum=ll.rest+ll.chain

    df.ch=(length(chain)^2-length(chain))/2
    df.rs=(length(rest)^2-length(rest))/2
    diff=abs(ll.tot-ll.sum)*2
    # lower.tail ?
    p.value=stats::pchisq(diff, df=df.ch+df.rs,lower.tail=FALSE)
    # calculate if block is 
    # block loop
    ll.block=0
    
    for (i in 1:(length(chain)-1)) {
        ll.block=ll.block+LL.getLL(data[,c(chain[i],chain[i+1])])-
                                       LL.getLL(data[,chain[i+1],drop=FALSE])
    }
    ll.block=ll.block+LL.getLL(data[,chain[length(chain)],drop=FALSE])
    df.block=(length(chain)^2-length(chain))/2-(length(chain)-1)
    ch.block=abs(ll.chain-ll.block)*2
    block.p.value=stats::pchisq(ch.block,df=df.block,lower.tail=FALSE)
    # no typo ll.chain <=> ll.block as due to
    # Christian ll.block has to be larger than ll.chain
    # he was calculatiing the same values using Matlab (Mail from February/March 2020
    return(list(ll.total=ll.tot,ll.chain=ll.block,
                ll.rest=ll.rest,ll.block=ll.chain,df=df.ch+df.rs,
                chisq=diff,p.value=p.value,
                block.df=df.block,block.ch=ch.block,block.p.value=block.p.value))
}

LL.makeDF = function (data,chains) {
    resDF = NULL
    cormt=cor(data)
    diag(cormt)=0
    total.rs=sum(cormt)/2
    for (nm in names(chains)) {
        r2sum=0
        for (i in 1:(length(chains[[nm]])-1)) {
            r2sum=r2sum+cormt[chains[[nm]][i],chains[[nm]][i+1]]^2
        }
        r2per=(r2sum/total.rs)*100
        newDF=data.frame(chain=nm,
                         members=paste(chains[[nm]],collapse="-"),
                         r2sum=round(r2sum,3),r2per=round(r2per,2))
        llres=LL.getChainLL(data,chains[[nm]])
        for (ln in names(llres)) {
            newDF=cbind(newDF,llres[[ln]])
            colnames(newDF)[ncol(newDF)]=ln
        }
        if (is.null(resDF)) {
            resDF=newDF 
        } else {
            resDF=rbind(resDF,newDF)
        }
    }
    return(resDF)
}
