#!/usr/bin/env Rscript

library(DEXSeq)
library(DESeq2) 
library(optparse)

########################## read arguments

option_list <- list(
    make_option(c('--support.frame'), help='', default='data/RNASeq_AD_Tc1J20.tab'),
    make_option(c('--code'), help=''),
    make_option(c('--gff'), help='', default="/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/GTF/Gallus_gallus.Galgal4.78.gff"),
    make_option(c('--iFolder'), help='', default="/scratch2/vyp-scratch2/Daudet_RNASeq/processed"),
    make_option(c('--annotation.file'), help=''),
    make_option(c('--keep.dups'), help='', default=FALSE)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

support.frame <- opt$support.frame
code <- opt$code
gff <- opt$gff
iFolder <- opt$iFolder
annotation.file <- opt$annotation.file
keep.dups <- opt$keep.dups

#gff <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/Tc1_mouse/GTF/Tc1.gff"
#annotation.file <- '/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/Tc1_mouse/tc1_annotations.tab'
#iFolder <- '/scratch2/vyp-scratch2/IoN_RNASeq/Frances/processed'
#code <- 'Zanda_AD_Tc1J20'


#gff <- annotation.file <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/chicken/biomart/biomart_annotations_chicken.tab"
#support.frame <- "support/Daudet_RNASeq.tab"
#code <- "Daudet"

annotations <- read.table(annotation.file, header = TRUE, sep = '\t', quote = "")
names(annotations) <- ifelse (names(annotations) == "external_gene_name", "external_gene_id", names(annotations)) # trying to agree on the column names

###check input files and data frame
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)


#support <- subset(support, !is.na(condition))  #I think we want all the count data here
my.ids <- support$sample


if (!keep.dups) print(count.files <- paste(iFolder, '/dexseq/', my.ids, '_dexseq_counts.txt', sep = ''))
if (keep.dups) print(count.files <- paste(iFolder, '/dexseq/', my.ids, '_dexseq_counts_keep_dups.txt', sep = ''))


if (sum(!file.exists(count.files)) > 0) {
  print(subset(count.files, ! file.exists(count.files)))
  stop('Some input files are missing')
}


### dexseq output folders
fig.folder <- paste(iFolder, '/fig', sep = '')
dexseq.folder <- paste(iFolder, '/dexseq', sep = '')
deseq.folder <- paste(iFolder, '/deseq2', sep = '')
deseq.figs <- paste(iFolder, '/deseq2/figs', sep = '')

for (folder in c(dexseq.folder, fig.folder, deseq.folder)) {
  if (! file.exists(folder)) dir.create(folder)
}

deseq.geneTable <- paste(deseq.folder, '/gene_info.csv', sep = '')
dexseq.exonTable <- paste(dexseq.folder, '/exon_info.csv', sep = '')

if (!keep.dups) {
  rpkm.file <- paste(deseq.folder, '/rpkm_values.csv', sep = '')
  sizeFactors.file <- paste(deseq.folder, '/size_factors_deseq.csv', sep = '')
  
  deseq.counts <- paste(deseq.folder, '/deseq_counts_', code, '.RData', sep = '')  ##this contains the key data
  dexseq.counts <- paste(dexseq.folder, '/dexseq_counts_', code, '.RData', sep = '')  ##this contains the key data
  
  dexseq.counts.norm <- paste(dexseq.folder, '/dexseq_counts_norm_', code, '.tab', sep = '')  ##this contains the key data
}


if (keep.dups) {
  rpkm.file <- paste(deseq.folder, '/rpkm_values_keep_dups.csv', sep = '')
  sizeFactors.file <- paste(deseq.folder, '/size_factors_deseq_keep_dups.csv', sep = '')
  
  deseq.counts <- paste(deseq.folder, '/deseq_counts_', code, '_keep_dups.RData', sep = '')  ##this contains the key data
  dexseq.counts <- paste(dexseq.folder, '/dexseq_counts_', code, '_keep_dups.RData', sep = '')  ##this contains the key data
  
  dexseq.counts.norm <- paste(dexseq.folder, '/dexseq_counts_norm_', code, 'keep_dups.tab', sep = '')  ##this contains the key data
}




#### Now compute things

message('Reading the count data')
support$condition <- 1 ## needs a dummy condition for the function below to work
start_t<- proc.time() 

DexSeqExons <- DEXSeqDataSetFromHTSeq(count.files, sampleData = support, flattenedfile = gff)

# converting counts to integer mode
# the design formula contains a numeric variable with integer values,
# specifying a model with increasing fold change for higher values.
# did you mean for this to be a factor? if so, first convert
# this variable to a factor using the factor() function
# Error: all(object@modelFrameBM$sample %in% colData(object)$sample) is not TRUE
# In addition: Warning message:
# In DESeqDataSet(rse, design, ignoreRank = TRUE) :
# some variables in design formula are characters, converting to factors

message("Time used in dexseqdatasetfromHTseq") 
message(proc.time()-start_t)

##I am puzzled by this "twice the column number" thing
my.counts <- counts( DexSeqExons)[, 1:length(count.files) ]  
colnames(my.counts) <- colData(DexSeqExons)$sample.1[ 1:length(count.files) ]
# create gene-based counts from the DEXSeq exon counts
gene.names <- gsub(rownames(my.counts), pattern = ':.*', replacement = '')
genes.counts <- aggregate(x = my.counts, by = list(gene.names), FUN = sum)
my.genes <- genes.counts$Group.1
# drop the first column
genes.counts <- as.matrix( genes.counts[, -1])
row.names(genes.counts) <- my.genes
save(list = 'genes.counts', file = deseq.counts)



###### Now compute the RPKM data if the feature length information is available
feature.length.file <- gsub( pattern = '.gff', replacement = '_length_features.tab', gff)
if (!file.exists(feature.length.file)) {
  message('There is no feature length file ', feature.length.file)
}
if (file.exists(feature.length.file)) {
  message('There is a feature length file so I will now compute the RPKM values using DESeq normalization')
  feature.lengths <- read.table( feature.length.file, header = TRUE)
  rpkms <- genes.counts

  support$condition.dummy <- "nocondition"
  support$condition.dummy [1:floor(nrow(support)/2)] <- "somecondition"  ##silly hack to make the function below work
  
  dds <- DESeqDataSetFromMatrix(countData = genes.counts,
                                colData = support[, 'condition.dummy', drop = FALSE],
                                design = ~ condition.dummy)
  dds <- estimateSizeFactors(dds)
  sizeFactors <- sizeFactors(dds)

  write.csv(x = data.frame( sample =  dimnames(rpkms)[[2]], sizeFactors = sizeFactors),
            file = sizeFactors.file, row.names = FALSE, quote = TRUE)
  
  if (sum(is.na(sizeFactors)) > 0) stop('Some of the size factors are equal to NA')
  average.depth <- sum(as.numeric(rpkms))/(10^6*ncol(rpkms))

########## Now compute the length of each feature
  compute.lengths <- data.frame( id = dimnames(rpkms)[[1]], length = NA)
  compute.lengths$length <- feature.lengths$Length[ match( compute.lengths$id, table = feature.lengths$GeneID) ]

##### the composite features now
  problematic <- subset(compute.lengths, is.na(length))
  problematic$length <- sapply (as.character(problematic$id),
                                FUN = function(x) {sum(subset(feature.lengths, GeneID %in% strsplit(x, split = '\\+')[[1]], Length, drop = TRUE))})
  
  compute.lengths$length [ as.numeric(row.names(problematic)) ] <- problematic$length

###### Now normalize for everything
  for (i in 1:ncol(rpkms)) {
    rpkms[, i ] <- rpkms[, i ] /((compute.lengths$length/1000)*average.depth * sizeFactors [ i ])
  }
  
  my.median <- apply(rpkms, MAR = 1, FUN = median)
  rpkms <- as.data.frame(rpkms [ order( my.median, decreasing = TRUE), ]) ##reorder

  samples.names <- names(rpkms)

####### Now fix the gene names
  rpkms$external_gene_id <- as.character(annotations$external_gene_id [ match(row.names(rpkms), annotations$EnsemblID) ])

  problematic <- subset(rpkms, is.na(external_gene_id))
  problematic$external_gene_id <- sapply (as.character(row.names(problematic)),
                                          FUN = function(x) {paste(
                                            subset(annotations, EnsemblID %in% strsplit(x, split = '\\+')[[1]], external_gene_id, drop = TRUE),
                                            collapse = '+')})
  problematic$external_gene_id <- as.character(problematic$external_gene_id)
  rpkms[ row.names(problematic), 'external_gene_id' ] <- problematic$external_gene_id


########## finalize the computation
  rpkms$ensemblID <- row.names( rpkms ) 
  rpkms <- rpkms[, c('external_gene_id', 'ensemblID', samples.names) ]

  write.csv( x = rpkms, file = rpkm.file, row.names = FALSE, quote = TRUE)
  message("Done printing RPKM values")
}



################# process the gff table
gff.table <- read.table(file = gff, sep = '\t')[, c('V1', 'V4', 'V5', 'V7', 'V9')]
names(gff.table) <- c('chromosome', 'start', 'end', 'strand', 'long.name')
gff.table$EnsemblID <-  gsub(pattern = '.*gene_id ', replacement = '',  gff.table$long.name)

gff.table.genes <- subset(gff.table, grepl(pattern = '^gene_id', gff.table$long.name))
gff.table.exons <- subset(gff.table,!grepl(pattern = '^gene_id', gff.table$long.name))

gff.table.genes <- merge(gff.table.genes, annotations, by = 'EnsemblID', all.x = TRUE)
write.csv( x = gff.table.genes, file = deseq.geneTable, row.names = FALSE)
write.csv( x = gff.table.exons, file = dexseq.exonTable, row.names = FALSE)



if (file.exists(feature.length.file)) {  ##if a RPKM file has been created
  

######################## basic PCA analysis
  my.conditions <- names(support)[grepl(names(support), pattern = '^condition.*')]  ##which condition should we use?
  if (length(my.conditions) > 0) {
    my.legend <- TRUE
    if (length(my.conditions) == 1) {my.cond <- my.conditions;} else {
      message('Multiple conditions, I need to select one for the plot')
      my.numb <- sapply(X = support[, my.conditions], FUN = function(x) {sum(!is.na(x))})
      my.cond <- names(sort(my.numb, decreasing = TRUE)[1])
      message('My condition for PCA analysis: ', my.cond)
    }
  } else {
    my.legend <- FALSE ##no legend in this case
  }
  

  rpkms.num <- as.matrix(rpkms[, samples.names])
  my.sd <- apply(rpkms.num, MAR = 1, FUN = sd)
  mat.for.pca <- t(rpkms.num[my.sd > median(my.sd), ])
  pca.data <- prcomp(mat.for.pca, scale = TRUE)
  
  output.pca <- paste(fig.folder, '/', code, '_pca.pdf', sep = '')
  pdf(output.pca)
  

  for (i in 1:2) {
    message('PCA ', i)
    col <- 'black'
    if (my.legend) {
      col = as.numeric(factor(support[, my.cond]))
      col <- ifelse (is.na(col), 'grey', col)
    }
    
    plot(x = pca.data$x[,2*i-1],
         y = pca.data$x[,2*i],
         xlab = ifelse (i == 1, 'PC1', 'PC3'),
       ylab = ifelse (i == 1, 'PC2', 'PC4'),
         col = col,
         pch = '+')
    
    text(x = pca.data$x[,2*i -1],
         y = pca.data$x[,2*i],
         labels = as.character(support$sample),
         pos = 3)
    
    if (my.legend) {
      my.levels <- levels(factor(support[, my.cond]))
      legend(col = 1:length(my.levels),
           legend = my.levels,
             pch = '+',
             x = 'bottomright')
    }
    
  }
  print(output.pca)
  dev.off()
}


##########
#d <- dist(mat.for.pca)   # find distance matrix
#hc <- hclust(d)                # apply hirarchical clustering
#hc$labels <- paste(gsub(pattern = '_dexseq_counts.txt', replacement = '', basename (hc$labels)), support[, my.cond], sep= '-')
#dimnames(genes.counts)[[2]] <- basename(dimnames(genes.counts)[[2]])

#output.hclust <- paste(fig.folder, '/', code, '_hclust.pdf', sep = '')

#pdf(output.hclust)
#heatmap(genes.counts)
#plot(hc)  
#dev.off()
