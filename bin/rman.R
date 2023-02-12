#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

for (file in args) {
    if (!file.exists(file)) {
        stop(paste("Error: File",file,"does not exists"))
    }
    fin  = file(file, "r")
    flag = FALSE
    while(length((line = readLines(fin,n=1)))>0) {
        if (grepl("^#' *\\\\name\\{.+\\}",line)) {
                  name = gsub("\\$","_",gsub("^#' *\\\\name\\{(.+)\\}.*","\\1",line))
                              print(paste("name = '",name,"'",sep=''))
            flag = TRUE
            mandir = sub("R/.+","man/",file)
            fout = file(paste(mandir,name,".Rd",sep=""),"w")
            
        } else if (flag & (grepl('^[a-zA-Z" ].*',line) | grepl("^$",line))) {
            flag = FALSE
            close(fout) 
            next
        } 
        if (flag) {
            line = sub("^#' ?","",line)
            cat(line,"\n",file=fout)
        }
    }
    close(fin)
}
