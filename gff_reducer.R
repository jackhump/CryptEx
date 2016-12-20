# R script to reduce cryptic GFF file
library(data.table)
library(optparse)

option_list <- list(
	make_option(c('--resFolder'), help = '', default = '/SAN/vyplab/HuRNASeq/misc/mouse_chr8_onegene_test.bam'),
    make_option(c('--output'), help = '', default = 'yes'),
    make_option(c('--exon_GFF'), help='', default='no'),
    make_option(c('--outFile'), help='', default='no'),
)  

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

resFolder <- opt$resFolder
output <- opt$output
exon_GFF <- opt$exon_GFF
outFile <- opt$outFile

cryptic_gff <- fread( paste0(output, ".cryptics.gff" ) )
genes <- str_split_fixed(cryptic_gff$V9, "gene_id ", 2)[,2]
genes <- str_split_fixed(genes, '\\"', 3 )[,2]

gff <- fread( exon_GFF )
gff$EnsemblID <- str_split_fixed(gff$V9, "gene_id ", 2)[,2]
gff$EnsemblID <- str_split_fixed(gff$EnsemblID, '\\"', 3 )[,2]

# only keep the genes that have cryptic exons
reduced <- subset(gff, EnsemblID %in% genes)

gff$EnsemblID <- NULL

total <- rbind(gff, cryptic_gff)

write.table(total, paste0( resFolder,"/",outFile,".unsorted.reduced.total.gff" ), quote =F, sep = "\t" )