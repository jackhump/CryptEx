#Rscript to read in the support frame and remove any all-NA columns.
args = commandArgs(trailingOnly=TRUE)
support_frame <- args[1]
support <- read.table(support_frame,header=T)
support <- support[,apply(X=support,MARGIN=2,FUN=function(x) !(sum(is.na(x))==length(x)))]
write.table(support, support_frame,quote=F,row.names=F,sep="\t")

