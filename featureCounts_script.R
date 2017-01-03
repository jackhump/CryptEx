
library(Rsubread)
library(stringr)
library(data.table)
library(optparse)
options(echo=T)
########################## read arguments

option_list <- list(
	make_option(c('--bamFile',"-b"), help = '', default = '/SAN/vyplab/HuRNASeq/misc/mouse_chr8_onegene_test.bam'),
    make_option(c('--paired',"-p"), help = '', default = 'yes'),
    make_option(c('--countStrand'), help='', default='no'),
    make_option(c('--gff', '-g'), help='', default = "/SAN/vyplab/HuRNASeq/misc/mouse_chr8_test.gtf"),
    make_option(c('--output', '-o'), help='', default="/SAN/vyplab/HuRNASeq/misc/mouse_chr8_test_counts.tab")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

bamFile <- opt$bamFile
paired <- opt$paired
countStrand <- opt$countStrand
GFF <- opt$gff
outFile <- opt$output

# get strandedness correct
strandedness <- 0
if(countStrand == "yes"){
	strandedness <- 1
}
if(countStrand == "reverse"){
	strandedness <- 2
}

isPaired <- FALSE
if(paired == "yes"){
	isPaired <- TRUE
}

counts <- featureCounts(bamFile,
	reportReads=FALSE,
	annot.ext = GFF,
	isGTFAnnotationFile=TRUE,
	GTF.featureType = "exonic_part",
	strandSpecific = strandedness,
	countMultiMappingReads=FALSE,
	useMetaFeatures=FALSE,
	allowMultiOverlap=TRUE,
	minOverlap=1,
	fracOverlap=0.1,
	fraction=FALSE,
	ignoreDup=TRUE,
	isPairedEnd=isPaired,
	tmpDir=dirname(outFile)
)
gff <- fread(GFF)

exonID <- str_split_fixed(str_split_fixed(gff$V9, ";", 3)[,2], "\"", 3)[,2]
row.names(counts$counts) <- paste(row.names(counts$counts), exonID, sep = ":")

write.table(counts$counts, file = outFile, sep = "\t", col.names = FALSE, row.names = TRUE, quote=F)



