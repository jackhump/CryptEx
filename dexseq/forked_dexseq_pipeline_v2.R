#!/usr/bin/env Rscript

#FOR DEBUGGING
gff <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//TDP-43_patient_human/Prudencio/GFF/TDP-43_patient_human_Prudencio.strict.500.total.cryptics.gff"
keep.sex <- TRUE 
keep.dups <- FALSE 
support.frame <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/support/Prudencio_ALS_brain_support.tab"
code <- "Prudencio" 
annotation.file <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human_hg38/biomart/biomart_annotations_human.tab" 
iFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//TDP-43_patient_human/Prudencio/strict_500"
cryptic <- TRUE


options(echo=T) 

library(DEXSeq)
library(BiocParallel)
library(optparse)


option_list <- list(
    make_option(c('--support.frame'), help=''),
    make_option(c('--code'), help=''),
    make_option(c('--gff'), help=''),
    make_option(c('--iFolder'), help=''),
    make_option(c('--annotation.file'), help=''),
    make_option(c('--keep.dups'), help='', default=FALSE),
    make_option(c('--keep.sex'), help='', default=FALSE),
    make_option(c('--cryptic'), help='', default=FALSE) 
)

########################## read arguments
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


#annotation.file <- '/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/Tc1_mouse/tc1_annotations.tab'
#iFolder <- '/scratch2/vyp-scratch2/IoN_RNASeq/Frances/processed'
#support.frame <- 'data/RNASeq_AD_Tc1J20.tab'
#code <- 'Zanda_AD_Tc1J20'
#gff <- '/cluster/project8/vyp/vincent/Software/pipeline/RNASeq/bundle/Tc1_mouse/GTF/Tc1.gff'
#keep.dups <- FALSE
#keep.sex <- FALSE

support.frame <- opt$support.frame
code <- opt$code
iFolder <- opt$iFolder
annotation.file <- opt$annotation.file
gff <- opt$gff
keep.dups <- opt$keep.dups
keep.sex <- opt$keep.sex
cryptic <-  opt$cryptic

dexseq.compute <- TRUE 


###check input files and data frame
message('gff file is ', gff)
message('Now reading ', support.frame)
support <- read.table(support.frame, header = TRUE, stringsAsFactors = FALSE)
#remove any columns that are just NA values - occurs in special cases
support <- support[,apply(X=support,MARGIN=2,FUN=function(x) !(sum(is.na(x))==length(x)))]
my.ids <- support$sample
list.conditions <- grep(names(support), pattern = '^condition.*', value  = TRUE)

annotation <- read.table(annotation.file, header = TRUE, sep = '\t', na.string = c('NA', ''), quote = "")
names(annotation) <- ifelse (names(annotation) == "external_gene_name", "external_gene_id", names(annotation)) # trying to agree on the column names

BPPARAM = MulticoreParam(workers=4)
# if more than 8 samples are present then the job is run on 12 cores instead so make sure to use them all!
if(nrow(support) > 8)
  BPPARAM = MulticoreParam(workers=12)
}

if( cryptic ){
    (files <- paste(iFolder, '/counts/',my.ids, '_dexseq_counts.txt', sep = '') )
} else {
    ( files <- paste(iFolder, '/', my.ids, '/dexseq/', my.ids, '_dexseq_counts.txt', sep = '') )
}


if (sum(!file.exists(files)) > 0) {
  message(files [ !file.exists(files) ])
  stop('Some input files are missing')
}


### dexseq output folders
dexseq.folder <- paste(iFolder, '/dexseq', sep = '')
dexseq.counts <- paste(dexseq.folder, '/dexseq_counts_', code, '.RData', sep = '')  ##this contains the key data
if (!file.exists(dexseq.folder)) dir.create(dexseq.folder)



my.ids <- support$sample
if (!keep.dups) countFiles <- paste(iFolder, '/', my.ids, '/dexseq/', my.ids, '_dexseq_counts.txt', sep = '')
if (keep.dups) countFiles <- paste(iFolder, '/', my.ids, '/dexseq/', my.ids, '_dexseq_counts_keep_dups.txt', sep = '')
if ( cryptic ) countFiles <- paste(iFolder, '/counts/', my.ids, '_dexseq_counts.txt', sep = '')

countFiles

for (condition in list.conditions){
  message('Condition ', condition)
  support.loc <- support

  ##handle the type variable
  support.loc$condition <- factor(support[, condition])
  loc.countFiles <- countFiles[ !is.na(support.loc$condition) ]
  support.loc <-  support.loc[ !is.na(support.loc$condition), ]


  ##handle the type variable
  type.loc <- gsub(x = condition, pattern = 'condition', replacement = 'type')
  if ( (! type.loc %in% names(support.loc)) & ('type' %in% names(support.loc))) {type.loc <- 'type'}  ##if only "type" is present, use it
  if (type.loc %in% names(support)) {
    support.loc$type <- factor(support.loc[, type.loc])
    support.loc <- subset(support.loc, !is.na(type))
  }

  loc.code <-  paste(unique(support.loc$condition), collapse = '_')
  message('Support data frame', loc.code)
  print(support.loc)


  ################### create the appropriate folders
  loc.dexseq.folder <- paste(iFolder, '/dexseq/', loc.code, sep = '')
  dexseq.figs <- paste(loc.dexseq.folder, '/figs', sep = '')
  dexseq.data <- paste(loc.dexseq.folder, '/dexseq_', code, '_', loc.code, '.RData', sep = '')  ##this will contain the final output of dexseq

  for (folder in c(loc.dexseq.folder, dexseq.figs)) {
    if (! file.exists(folder)) dir.create(folder)
  }
  message("we got this far!")

  if (dexseq.compute) {  
    #load(dexseq.counts)  ##object mycounts is key
    #DexSeqExons <- subset(DexSeqExons, c(rep(TRUE, 300), rep(FALSE, nrow(counts(DexSeqExons)) - 300))) 
    use.covariate <- FALSE
    if ('type' %in% names(support.loc)) {
      if (length(unique(as.character(support.loc$type))) > 1){
        use.covariate <- TRUE
      }
    }

    if (use.covariate) {
      formuladispersion <- ~ sample + (condition + type) * exon
      formula0 <-  ~ sample + type * exon + condition
      formula1 <-  ~ sample + type * exon + condition * exon
      my.design <- support.loc[, c('type', 'condition')]
      my.design$type <- factor(my.design$type) ## probably not needed
      my.design$condition <- factor(my.design$condition)  ## probably not needed
      my.design.loc <- my.design  ##just to print basically
    } 
    else {
      formuladispersion <-  ~ sample + condition * exon
      formula0 <-  ~ sample + exon + condition
      formula1 <-  ~ sample + exon + condition * exon
      my.design <- factor(support.loc[, c('condition')])
      my.design.loc <- support.loc[, c('condition'), drop = FALSE]  ##just to print basically
    }

    row.names(my.design.loc) <- factor(support.loc$sample)
    message("here we go again!")
    DexSeqExons.loc <- DEXSeqDataSetFromHTSeq(loc.countFiles,
                                              sampleData = my.design.loc,
                                              design = formula1,
                                              flattenedfile = gff)
    
    write.table(x = my.design.loc, file = paste(loc.dexseq.folder, '/design.tab', sep = ''), row.names = TRUE, quote = FALSE, sep = '\t')

    
    
    #message('Updating the dexseq object')
    #DexSeqExons.loc <- DEXSeqDataSet(countData= featureCounts(mycounts),
    #                                 sampleData = my.design.loc,
    #                                 design= formula1,
    #                                 groupID=geneIDs(mycounts),
    #                                 featureID=exonIDs(mycounts),
    #                                 transcripts = DexSeqExons@featureData@data$transcripts,
    #                                 featureRanges = GRanges(DexSeqExons@featureData@data$chr, IRanges (start = DexSeqExons@featureData@data$start, end = DexSeqExons@featureData@data$end)) )


    #DexSeqExons.loc <- DexSeqExons.loc[1:20,]  ##VP
    
    message('Starting the computations')
    DexSeqExons.loc <- estimateSizeFactors(DexSeqExons.loc)

    
    message('Here is the part that takes a lot of time')
    DexSeqExons.loc <- DEXSeq::estimateDispersions(DexSeqExons.loc, BPPARAM=BPPARAM)
    #fData(DexSeqExons.loc)$dispersion <- fData(DexSeqExons.loc)$dispBeforeSharing
    message('Done with estimateDispersions')    
    DexSeqExons.loc <- DEXSeq::testForDEU(DexSeqExons.loc, BPPARAM=BPPARAM)
    message('Done with testDEU')    
    DexSeqExons.loc <- DEXSeq::estimateExonFoldChanges(DexSeqExons.loc, BPPARAM=BPPARAM)
    message('Done with estimateFoldChange')
    
  ######################### output basic table
    res <- DEXSeq::DEXSeqResults (DexSeqExons.loc)
    logname <- grep(names(res), pattern = 'log2fold', value = TRUE)
    res.clean <- as(res[, c('groupID', 'featureID', 'exonBaseMean', logname, 'dispersion', 'stat', 'pvalue')], 'data.frame')
    names(res.clean)<- c("EnsemblID", "exonID", "meanBase", "log2FoldChange", "dispersion", "stat", "pvalue")
    
    res.clean$FDR <- p.adjust(res.clean$pvalue, method = 'fdr')    
    res.clean$chromosome <- as.character(seqnames( res$genomicData))
    res.clean$exon.start <- start(res$genomicData)
    res.clean$exon.end <- end(res$genomicData)


    res.clean$external_gene_id <- annotation$external_gene_id[ match(res.clean$EnsemblID, table = annotation$EnsemblID) ]
    res.clean <- res.clean[, c('external_gene_id', "EnsemblID", "exonID", "meanBase", "log2FoldChange", "dispersion", "stat", "pvalue", "FDR", "chromosome", "exon.start", "exon.end")]  ### reorder the names nicely

    if ('strand' %in% names(annotation)) res.clean$strand <- annotation$strand[ match(res.clean$EnsemblID, table = annotation$EnsemblID) ] ## add strand if available
    res.clean <- res.clean[ order(res.clean$pvalue),]  ##reorder the rows
    
    write.csv(x = res.clean,
              file=paste(loc.dexseq.folder, "/", code, "_", loc.code, "_SignificantExons.csv", sep = ''),
              row.names = FALSE)
    if (cryptic) {
      res.clean.cryptics <- subset(res.clean,res.clean$FDR < 0.01 & grepl("i",res.clean$exonID))
      res.clean.cryptics.up <- subset(res.clean.cryptics, res.clean.cryptics$log2FoldChange > 0)
      res.clean.cryptics.down <- subset(res.clean.cryptics, res.clean.cryptics$log2FoldChange < 0)
      res.clean.cryptics.NA <- subset(res.clean.cryptics, is.na(res.clean.cryptics$log2FoldChange))
      write.table(x = res.clean.cryptics,
           file=paste(loc.dexseq.folder, "/", code, "_", loc.code, "_CrypticExons.tab", sep = ''),
            row.names = FALSE)
      codes <- c("Significant Cryptic events (FDR < 0.01):", "Up-going:", "Down-going:", "NA (possible error):")
      counts <- c(dim(res.clean.cryptics)[1], dim(res.clean.cryptics.up)[1], dim(res.clean.cryptics.down)[1], dim(res.clean.cryptics.NA)[1] ) 
      report <- data.frame(codes,counts)
      write.table(x = report,
          file=paste0(loc.dexseq.folder, "/", code, "_", loc.code, "Cryptic_Report.tab"),
          row.names = F, quote = F)
      write.table(x = res.clean.cryptics.up[,c(10:12,1,3)],
          file=paste0(loc.dexseq.folder, "/", code, "_", loc.code, "Cryptics_UP.bed"),
          quote=F, row.names=F, col.names=F, sep="\t")
      write.table(x = res.clean.cryptics.down[,c(10:12,1,3)],
          file=paste0(loc.dexseq.folder, "/", code, "_", loc.code, "Cryptics_DOWN.bed"),
          quote=F, row.names=F, col.names=F, sep="\t")    
      message('Saving results in ', dexseq.data)
      save(list = c('res.clean', 'DexSeqExons.loc'), file = dexseq.data)
    }
  } 
  else {
    message("looky here!")
    load(dexseq.data)
  }

  message("ok still going after avoiding the DEXSeq lifting")

  ########################## Now plot a subset
  file.remove(list.files(dexseq.figs, pattern = 'DEXSeq*', full.names = TRUE)) ##remove the old plots
  #Altered so only plots significant CRYPTIC events
  n.sig <- sum(res.clean$FDR < 0.01, na.rm = TRUE)
  if (n.sig <= 50) {
    res.cleanSigs <- subset(res.clean, FDR<0.01 & grepl("i",res.clean$exonID))
  } else res.cleanSigs <- subset(res.clean, grepl("i",res.clean$exonID))[1:50,]


  genes.to.plot <- unique(res.cleanSigs$EnsemblID)
  pretty.gene.names <- as.character(annotation$external_gene_id[ match(genes.to.plot, table = annotation$EnsemblID) ])

  for (i in 1:length(genes.to.plot)) {
    gene <- as.character(genes.to.plot[i])

    if (!is.na(pretty.gene.names[ i ])) {
      gene.pretty <- as.character(pretty.gene.names[ i ])
      
      message(i, ' ', gene, ' ', gene.pretty)
      
      output.pdf <- paste(dexseq.figs, '/DEXSeq-', gene.pretty, '.pdf', sep = '')
      pdf(output.pdf, width = 16, height = 9.8)
      plotDEXSeq(res,
                 geneID = gene,  ##I suspect it has to be gene, otherwise it crashes
                 cex.axis = 1.2,
                 cex=1.3,
                 lwd=2,
                 legend=TRUE,
                 displayTranscripts = TRUE,
                 names = TRUE,
                 main = gene.pretty)
      dev.off()
      print(output.pdf)
    }
  }

  message("apparently printed some genes")


  #################################### plot some graphs

  pdf(file=paste(dexseq.figs, '/DEXSeq-MeanVsDispPoints.pdf', sep = ''))
  plotDispEsts (DexSeqExons.loc)
  dev.off()

    
  pdf(file=paste(dexseq.figs, '/DEXSeq-MeanVsDispCircles.png', sep = ''))
  plotMA(data.frame(baseMean = res.clean[,6],log2FoldChange = res.clean[,7], padj = res.clean[,5] < 0.1),
         ylim=c(-4,4), cex=0.8)
  dev.off()

  message('Done with ', condition)
  rm(list = c('DexSeqExons.loc', 'res', 'res.clean'))
  gc()
}

warnings()

sessionInfo()