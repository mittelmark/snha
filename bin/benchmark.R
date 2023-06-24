library(snha)
library(lattice)
benchmark <-  function (iter=3) {
    performance <-  function (graph,data,FUN, ...) {
        t1=Sys.time()
        pred=FUN(data,...) 
        t2=Sys.time()
        tdiff=difftime(t2,t1,units="secs")
        acc=mgraph_accuracy(graph,pred$theta)
        return(list(Sens=acc$Sens,Spec=acc$Spec,Prec=acc$Prec,MCC=acc$MCC,Time=tdiff))
    }
    bench <-  function (graph) {
        g.data = t(snha_graph2data(graph,iter=200,prop=0.02))
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
        return(data)
    }
    W=mgraph(type="werner")
    B0=mgraph(type="band",nodes=10)    
    B1=mgraph(type="barabasi",m=1)
    B2=mgraph(type="barabasi",m=2)
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
        df=rbind(df,data.frame(gtype=rep("B0",nrow(acc)),acc))
        acc=bench(B1)
        df=rbind(df,data.frame(gtype=rep("B1",nrow(acc)),acc))
        acc=bench(B2)
        df=rbind(df,data.frame(gtype=rep("B2",nrow(acc)),acc))
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
