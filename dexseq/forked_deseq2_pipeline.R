
library(DESeq2)
library(optparse)

option_list <- list(
    make_option(c('--support.frame'), help=''),
    make_option(c('--code'), help=''),
    make_option(c('--gff'), help=''),
    make_option(c('--iFolder'), help=''),
    make_option(c('--annotation.file'), help=''),
    make_option(c('--keep.dups'), help='', default=FALSE),
    make_option(c('--keep.sex'), help='', default=FALSE) 
)

### Just for debugging 
#support.frame <- "/cluster/project8/vyp/Tabrizi_Huntington_RNASeq/support/htt_support2.txt"
#code <- "htt"
#keep.dups <- FALSE 
#annotation.file <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human/biomart/biomart_annotations_human.tab"
#iFolder <- "/scratch2/vyp-scratch2/Tabrizi_Huntington_RNASeq/processed/Nov2014/"

## vincent debugging
#support.frame <- "data/TDP43_m323k.tab"
#code <- "m323k"
#annotation.file <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/mouse/biomart/biomart_annotations_mouse.tab"
#iFolder <- "/scratch2/vyp-scratch2/IoN_RNASeq/Fratta_RNASeq/brain/m323k"
#keep.dups <- FALSE
#keep.sex <- FALSE

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.frame <- opt$support.frame
code <- opt$code
gff <- opt$gff
iFolder <- opt$iFolder
annotation.file <- opt$annotation.file
keep.dups <- opt$keep.dups
keep.sex <- opt$keep.sex


########################## read arguments

extra.plots <- TRUE 
remove.hb   <- FALSE 

message("Should I keep the sex chromosomes? Option set is ", keep.sex)

###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)
list.covars <- grep(names(support), pattern = '^covar.*', value  = TRUE)


annotation <- read.table(annotation.file, header = TRUE, sep = '\t', na.string = c('', 'NA'), quote = "" )
names(annotation) <- ifelse (names(annotation) == "external_gene_name", "external_gene_id", names(annotation)) # trying to agree on the column names


### deseq output folders and files
deseq2.folder <- paste(iFolder, '/deseq2', sep = '')
for (folder in c(deseq2.folder)) {
  if (! file.exists(folder)) dir.create(folder)
}


########## load the count data

if (keep.dups) print( deseq.counts <- paste(deseq2.folder, '/deseq_counts_', code, '_keep_dups.RData', sep = '') )
if (!keep.dups) print( deseq.counts <- paste(deseq2.folder, '/deseq_counts_', code, '.RData', sep = '') )
load(deseq.counts)


### Remove the sex chromosome genes 
if (!keep.sex) {
  genes.on.XY <- as.character(subset(annotation, chromosome_name %in% c('X' ,'Y', 'chrX', 'chrY'), 'EnsemblID', drop = TRUE))
  message('Prior to removing chr XY probes: ', nrow(genes.counts))
  genes.counts <- genes.counts[ ! dimnames(genes.counts)[[1]] %in% genes.on.XY, ]                 
  message('After removing chr XY probes: ', nrow(genes.counts))
}


### Remove the hemoglobin genes (for whole blood only) 
if(remove.hb) { 
   message("Removing HBB, HBA1 and HBA2") 
   hb.genes <- annotation[which(annotation$external_gene_id %in% c("HBB", "HBA1", "HBA2")), 'EnsemblID'] 
   genes.counts <- genes.counts[ ! dimnames(genes.counts)[[1]] %in% hb.genes, ]                 
   message('After removing hb genes: ', nrow(genes.counts))
} 


### Loop over all proposed conditions
for (condition in list.conditions) {
  num.cond <- FALSE 
  message("Processing for ", condition) 
  samples.to.use   <- !is.na(support[,condition]) 
  
  genes.counts.loc <- genes.counts[, samples.to.use ]

  support.loc <-  support[  samples.to.use, ]

## determines if the condition be treated as a factor or a numeric variable 
  if (substr(condition, 10, 12) == "Num") { 
     num.cond <- TRUE 
     support.loc$condition <- as.numeric(support.loc[, condition])
     loc.code <-  condition
  } else { 
     support.loc$condition <- factor(support.loc[, condition])
     loc.code <-  paste(unique(support.loc$condition), collapse = '_')
  } 

################### create the appropriate folders and specify output file
  loc.deseq2.folder <- paste(iFolder, '/deseq2/', loc.code, sep = '')
  deseq2.figs <- paste(loc.deseq2.folder, '/figs', sep = '')

  for (folder in c(loc.deseq2.folder, deseq2.figs)) {
    if (! file.exists(folder)) dir.create(folder)
  }
  
  if (keep.dups) output.file <- paste(loc.deseq2.folder, '/deseq_', code, '_differential_expression_keep_dups.tab', sep = '')
  if (!keep.dups) output.file <- paste(loc.deseq2.folder, '/deseq_', code, '_differential_expression.tab', sep = '')

###################
  use.covars <- FALSE
  if (length(list.covars) > 0) {
     use.covars <- TRUE
  }
  
  if (use.covars) {
    message("Using ", length(list.covars), " covariates") 
    formula1 <- paste(" ~ ", paste(list.covars, collapse = "+"), " + ", condition, sep = "") 
    formula1 <- as.formula(formula1) 
    formula0 <- paste(" ~ ", paste(list.covars, collapse = "+"), sep = "") 
    formula0 <- as.formula(formula0) 
   
    design.deseq <- support.loc[, which(names(support.loc) %in% c(condition, list.covars))]
  } else {
    formula1 <- as.formula(paste0("~ ", condition))
    formula0 <- as.formula("~ 1")
    design.deseq <- support.loc[, c(condition), drop = FALSE]
  }

  
  
  CDS <- DESeqDataSetFromMatrix(countData = genes.counts.loc, colData = design.deseq, design = formula1)

#################### Do the actual model fitting 
  CDS <- DESeq(CDS, test = "LRT", reduced = formula0, 
               minReplicatesForReplace = 5 ) 
  deseq.res <- results(CDS)  

  ### output the design information
  row.names(design.deseq) <- support.loc$sample
  write.table(x = design.deseq, file = paste(loc.deseq2.folder, '/design.tab', sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')
  
############# Make the results table into a sensible format 
  deseq.res.df <- data.frame(deseq.res) 
  print(head(deseq.res.df)) 
  deseq.res.df$EnsemblID <- row.names( deseq.res.df)
  deseq.res.df <- merge(deseq.res.df, annotation, by = 'EnsemblID', all.x = TRUE)
  deseq.res.df <- deseq.res.df[order(deseq.res.df$pvalue),]
  print(head(deseq.res.df)) 

#################### now write the output to a file 
  write.table(deseq.res.df, file = output.file, quote = FALSE, row.names = FALSE, sep = "\t" ) 
  save(list = c("CDS"), file = paste0(loc.deseq2.folder, "/deseq2_object.RData")) 
  
######### Now add a PCA for the subset of individuals being considered
  if(extra.plots) { 
     if (keep.dups) { 
        output.pca <- paste(deseq2.figs, '/', loc.code, '_pca_keepdups.pdf', sep = '')
     } else { 
        output.pca <- paste(deseq2.figs, '/', loc.code, '_pca.pdf', sep = '')
     } 

     rld <- rlog(CDS)
     pdf(output.pca)
     plotPCA(rld, intgroup = condition) 
     dev.off() 

## Visualise the counts versus condition for the genes with best p-values 
     if (keep.dups) { 
        output.sig.genes <- paste(deseq2.figs, '/', loc.code, "_siggenes_keepdups.pdf", sep = '') 
     } else { 
        output.sig.genes <- paste(deseq2.figs, '/', loc.code, "_siggenes.pdf", sep = '') 
     }

     if (!num.cond) { 
     pdf(output.sig.genes) 
     sig.genes <- which(deseq.res.df$padj < 0.1) 
     for (i in sig.genes) { 
        plotCounts(CDS, gene = deseq.res.df$Ensembl[i], intgroup = condition, 
                   transform = TRUE, main = deseq.res.df$external_gene_id[i])
        mtext(paste("pval: ", deseq.res.df$padj[i], sep = ""), 3)  
     } 
     dev.off() 
     } 

     if (keep.dups) { 
        disp.plot <- paste(deseq2.figs, '/', loc.code, "_disp_keepdups.pdf", sep = '') 
     } else { 
        disp.plot <- paste(deseq2.figs, '/', loc.code, "_disp.pdf", sep = '') 
     } 

     pdf(disp.plot) 
     plotDispEsts(CDS)  
     dev.off() 
     
  } # End of extra plots  

} # End looping over conditions 

warnings()


sessionInfo()
