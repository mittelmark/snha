library(snha)
library(lattice)
benchmark <-  function (iter=3) {
    huge.predict <- function (x,method="mb",criterion="stars") {
        if (!requireNamespace("huge", quietly = TRUE)) {
            stop("Package \"huge\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
        library(huge)
        out.mb = huge(x,method=method,verbose=FALSE)
        M=as.matrix(huge.select(out.mb,criterion=criterion,verbose=FALSE)$refit)
        rownames(M)=colnames(M)=colnames(x)
        G=list()
        G$theta=M
        return(G)
    }   
    qgraph.predict <- function (x,type="ebic") {
        if (!requireNamespace("qgraph", quietly = TRUE)) {
            stop("Package \"qgraph\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
        require("qgraph")
        C=cor(x)
        if (type=="ebic") {
            res=EBICglasso(C,n=nrow(x),0.5,threshold=TRUE)
        } else if (type=="bic") {
            res=EBICglasso(C,n=nrow(x),0,threshold=TRUE)
        } else {
            stop("Error: type must be either 'ebic' or 'bic'")
        }
        res[res!=0]=1
        res[res==1 & C < 0]=-1
        
        return(list(theta=res))
    }
    glmnet.predict <- function (x,type="lasso",rs=0.04,alpha=1) {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
            stop("Package \"glment\" needed for this function to work. Please install it.",
                 call. = FALSE)
        }
        if (is.data.frame(x)) {
            x=as.matrix(x)
        }
        A=matrix(0,nrow=ncol(x),ncol=ncol(x))
        rownames(A)=colnames(A)=colnames(x)
        Coef=A
        rsqs=c()
        for (i in 1:ncol(x)) {
            xi=x[,-i]
            yi=x[,i]
            CV=glmnet::cv.glmnet(as.matrix(xi),yi,alpha=alpha,family="gaussian")
            lasso.model <- glmnet::glmnet(as.matrix(xi),yi, family = "gaussian",
                                          alpha=alpha,  lambda = CV$lambda.1se )
            rsqs=c(rsqs,lasso.model$dev.ratio)
            # normalize coefficients Agresti method(?)
            sds=apply(as.matrix(xi),2,sd)
            cs = as.matrix(coef(lasso.model, s = "lambda.1se"))
            std.coef = cs[-1, 1] * sds
            Coef[i,names(std.coef)]=std.coef
            c=coef(CV)[,1]
            c=c[2:length(c)]
            nm=names(which(c!=0))
            A[i,nm]=c[nm]
        }
        names(rsqs)=colnames(x)
        rsq=(abs(Coef)/apply(abs(Coef),1,sum))*rsqs
        C=A
        A[rsq<rs]=0
        A[A<0]=-1
        A[A>0]=1
        if (alpha == 1) {
            type="lasso"
        } else if (alpha == 0) {
            type="ridge"
        } else {
            type="elnet"
        }
        G=list()
        G$theta=A
        G$type=type
        G$coef=C
        G$std.coef=Coef
        G$r.squared=rsq
        return(G)
    }
    
    performance <-  function (graph,data,FUN, ...) {
        t1=Sys.time()
        pred=FUN(data,...) 
        t2=Sys.time()
        tdiff=difftime(t2,t1,units="secs")
        acc=mgraph_accuracy(graph,pred$theta)
        return(list(Sens=acc$Sens,Spec=acc$Spec,Prec=acc$Prec,MCC=acc$MCC,Time=tdiff))
    }
    bench <-  function (graph,g.data=NULL) {
        if (!is.matrix(g.data)) {
            g.data = t(snha_graph2data(graph,iter=200,prop=0.02))
        }
        res=performance(graph,g.data,snha)
        data=data.frame(method="1def",t(unlist(res)))
        res=performance(graph,g.data,snha,lmc=TRUE)
        data=rbind(data,
                   data.frame(method="2lmc",t(unlist(res))))
        res=performance(graph,g.data,snha,lmc=TRUE,lma=TRUE)
        data=rbind(data,
                   data.frame(method="3lma",t(unlist(res))))
        res=performance(graph,g.data,snha,prob=TRUE)
        data=rbind(data,
                   data.frame(method="4boot",t(unlist(res))))
        res=performance(graph,g.data,snha,lmc=TRUE,prob=TRUE)
        data=rbind(data,
                   data.frame(method="5lmb",t(unlist(res))))
        # glmnet
        res=performance(graph,g.data,glmnet.predict)
        data=rbind(data,
                   data.frame(method="g.ls",t(unlist(res))))
        # huge criterion ric does  not work!!
        res=performance(graph,g.data,huge.predict,method="glasso",criterion="stars")
        data=rbind(data,
                   data.frame(method="h.ls.s",t(unlist(res))))
        #res=performance(graph,g.data,huge.predict,method="mb",criterion="stars")
        #data=rbind(data,
        #           data.frame(method="h.mb.s",t(unlist(res))))
        res=performance(graph,g.data,huge.predict,method="ct",criterion="stars")
        data=rbind(data,
                   data.frame(method="h.ct.s",t(unlist(res))))
        res=performance(graph,g.data,huge.predict,method="glasso",criterion="ebic")
        data=rbind(data,
                   data.frame(method="h.ls.e",t(unlist(res))))
        res=performance(graph,g.data,huge.predict,method="glasso",criterion="ric")
        data=rbind(data,
                   data.frame(method="h.ls.r",t(unlist(res))))
        res=performance(graph,g.data,qgraph.predict,type="ebic")
        data=rbind(data,
                   data.frame(method="q.ls.e",t(unlist(res))))
        res=performance(graph,g.data,qgraph.predict,type="bic")
        data=rbind(data,
                   data.frame(method="q.ls.b",t(unlist(res))))
        return(data)
    }
    W=mgraph(type="werner")
    B0=mgraph(type="band",nodes=10)    
    B10=mgraph(type="barabasi",m=1)
    B13=mgraph(type="barabasi",m=1.3)
    B20=mgraph(type="barabasi",m=2)

    for (i in 1:iter) {
        cat("starting iter",i,"\n")
        acc=bench(W)
        res=data.frame(gtype=rep("W",nrow(acc)),acc)
        if (i == 1) {
            df=res
        } else {
            df=rbind(df,res)
        }
        acc=bench(B0)
        df=rbind(df,data.frame(gtype=rep("B00",nrow(acc)),acc))
        acc=bench(B10)
        df=rbind(df,data.frame(gtype=rep("B10",nrow(acc)),acc))
        acc=bench(B13)
        df=rbind(df,data.frame(gtype=rep("B13",nrow(acc)),acc))
        acc=bench(B20)
        df=rbind(df,data.frame(gtype=rep("B20",nrow(acc)),acc))
        SF=huge.generator(graph="scale-free",d=20)
        acc=bench(SF$theta,g.data=SF$data)
        df=rbind(df,data.frame(gtype=rep("SF",nrow(acc)),acc))
        cat("ending iter",i,"\n")
    }
    return(df)
}
print(Sys.time())
results=benchmark()
print(Sys.time())
write.table(results,file="results.tab",sep="\t",quote=FALSE)
agg=aggregate(results[,c('Sens','Spec','Prec','MCC','Time')],by=list(results[,'method'],results[,'gtype']),mean)
print(agg)
pdf("results.pdf",width=12,height=14)
pset=list(layout.heights=list(top.padding=-3, bottom.padding=-2),grid.layout=c(1,4))
p1=bwplot(Sens ~ method | gtype,data=results,ylab="Sens",par.settings=pset,ylim=c(0,1.1))
p2=bwplot(Spec ~ method | gtype,data=results,ylab="Spec",par.settings=pset,ylim=c(0,1.1))
p3=bwplot(Prec ~ method | gtype,data=results,ylab="Prec",par.settings=pset,ylim=c(0,1.1))
p4=bwplot(MCC ~ method | gtype,data=results,ylab="MCC",par.settings=pset,ylim=c(0,1.1))
print(p1[c(2,3,4,1)], position=c(0, .75, 1, 1), more=TRUE)
print(p2[c(2,3,4,1)], position=c(0, .5, 1, 0.75), more=TRUE)
print(p3[c(2,3,4,1)], position=c(0, .25, 1, 0.5), more=TRUE)
print(p4[c(2,3,4,1)], position=c(0, .00, 1, 0.25))
bwplot(Time ~ method | gtype,data=results)
dev.off()
