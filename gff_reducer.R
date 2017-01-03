# R script to reduce cryptic GFF file
library(data.table)
library(optparse)
library(stringr)

option_list <- list(
    make_option( c('--resFolder') ),
    make_option( c('--output') ),
    make_option( c('--exon_GFF') ),
    make_option( c('--outFile') )
)  

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

resFolder <- opt$resFolder
output <- opt$output
exon_GFF <- opt$exon_GFF
outFile <- opt$outFile

cryptic_gff <- fread( paste0(output, ".cryptics.gff" ) )
cryptic_gff$EnsemblID <- str_split_fixed(cryptic_gff$V9, "gene_id ", 2)[,2]
cryptic_gff$EnsemblID <- str_split_fixed(cryptic_gff$EnsemblID, '\\"', 3 )[,2]
#stupid bug in featureCounts screws up when EnsemblID is longer than 259 characters

print(length(cryptic_gff$EnsemblID) )
print(length( nchar(cryptic_gff$EnsemblID) < 260 ) )
print(dim(cryptic_gff) )

cryptic_gff <- cryptic_gff[ nchar(cryptic_gff$EnsemblID) < 260 ,] 

gff <- fread( exon_GFF )
gff$EnsemblID <- str_split_fixed(gff$V9, "gene_id ", 2)[,2]
gff$EnsemblID <- str_split_fixed(gff$EnsemblID, '\\"', 3 )[,2]

# only keep the genes that have cryptic exons
reduced <- subset(gff, gff$EnsemblID %in% cryptic_gff$EnsemblID)

total <- rbind(reduced, cryptic_gff)
total$EnsemblID <- NULL


###save.image( paste0( resFolder, "/", "test.Rdata" ) )

write.table(total, paste0( resFolder,"/",outFile,".unsorted.reduced.total.gff" ), quote =F, sep = "\t", col.names=F, row.names=F )

