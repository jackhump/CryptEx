# dinucleotide and trinucleotide enrichment
library(data.table)
library(stringr)
library(optparse)
library(ggplot2)

set.seed(1234)

files.exist <- function(files){
	for(file in files){
		if( !file.exists(file) ){
			stop( paste0(file, " does not exist!"))
		}
}

option_list <- list(
    make_option(c('--motifFolder'), help=''),
    make_option(c('--code'), help=''),
)

########################## read arguments
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

if(length(opt) > 1){
	code <- opt$code
	condition.names <- opt$condition.names
	outFolder <- opt$motifFolder
	species <- opt$species
}

cryptic.exons.fasta <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/Motif_Finding/dataset_1_cryptic.exons.fasta"
adjacent.introns.fasta <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/Motif_Finding/dataset_1_control.exons.fasta"

outFolder <- "/SAN/vyplab/HuRNASeq/ENCODE/HNRNPK/HepG2_ENCSR853ZJS/cryptex/Motif_Finding"
cryptic.exons.fasta <- "/SAN/vyplab/HuRNASeq/ENCODE/HNRNPK/HepG2_ENCSR853ZJS/cryptex/Motif_Finding/HepG2_ENCSR853ZJS_cryptic.exons.fasta"
control.exons.fasta <-"/SAN/vyplab/HuRNASeq/ENCODE/HNRNPK/HepG2_ENCSR853ZJS/cryptex/Motif_Finding/HepG2_ENCSR853ZJS_control.exons.fasta"

cryptic.exons.fasta <- paste0( outFolder, "/", code, "_cryptic.exons.fasta" )
control.exons.fasta <- paste0( outFolder, "/", code, "_control.exons.fasta")


files.exist( c( cryptic.exons.fasta, control.exons.fasta) )

# compare cryptic exons with adjacent intron sequences

cryptic.exons <- as.data.frame(fread(cryptic.exons.fasta,header=F,stringsAsFactors=F))
cryptic.exons <- data.frame(gene.id =  cryptic.exons[seq(1,nrow(cryptic.exons) - 1,2),],
							sequence = cryptic.exons[seq(2,nrow(cryptic.exons),2),])

control.exons <- as.data.frame(fread(control.exons.fasta,header=F,stringsAsFactors=F))
control.exons <- data.frame(gene.id =  control.exons[seq(1,nrow(control.exons) - 1,2),],
							sequence = control.exons[seq(2,nrow(control.exons),2),])

# random sample the control exons to compare lists of same length - save time!
control.exons <- control.exons[ sample( nrow( control.exons ),nrow(cryptic.exons ) ), ]

exons.list <- list(cryptic.exons, control.exons)
exons.names <- c("cryptic.exons","control.exons")



plots <- list()

for( n in 1:5){

	for(exon.type in 1:length(exons.list)){

		subseq.list <- list()
		subseq <- c()
		exons <- exons.list[[exon.type]]

		exons$sequence <- toupper(exons$sequence)

		for( exon in 1:length(exons$sequence)){
			for(i in 1:( str_length(exons$sequence[exon]) - n ) ){
				subseq[i] <- str_sub(exons$sequence[exon], i, (i+n) )	
			}
			subseq.list[[exon]] <- subseq
		}

		subseq.total <- do.call(what = c , args = subseq.list)
		# solution without dplyr - the updated version has changed things somewhat
		d <- as.data.frame( table(subseq.total) )
		d$prop <- d$Freq / sum(d$Freq)
		names(d) <- paste0(exons.names[exon.type],".",names(d))
		assign(paste0(exons.names[exon.type],".results"), d)
	}

	#apply(d, MAR = 1, FUN = function(X) prop.test(x = X[2], n = N, p = (1/16) ) )
	results <- cbind(control.exons.results, cryptic.exons.results)
	control.N <- sum(results$control.exons.Freq)
	cryptic.N <- sum(results$cryptic.exons.Freq)

	for(i in 1:nrow(results)){
		results$prop.test.pvalue[i] <- prop.test(x = c(results$cryptic.exons.Freq[i], results$control.exons.Freq[i]), 
							n = c(cryptic.N,control.N), p = NULL)$p.value
		results$log2FoldChange <- log2(results$cryptic.exons.prop / results$control.exons.prop )
	}
	results$prop.test.padj <- p.adjust(results$prop.test.pvalue, method = "bonferroni")
	results$control.exons.subseq.total <- gsub("T","U",results$control.exons.subseq.total)

	titles <- c("Dinucleotide", "Trinucleotide", "Quadnucleotide", "5mer")

	p <- ggplot(results, aes(x = log2FoldChange, y = -log10(prop.test.padj), label = control.exons.subseq.total) ) + 
		geom_text(size = 10, ) + 
		xlim(c(-1,1.5)) +
		ylim(c(0,250)) +
		ggtitle(paste0(gsub("_"," ",code)," (",species,")\n", titles[n], " enrichment between flanked cryptic exons\nand adjacent intronic sequence")) + 
		xlab("log2(cryptic vs null proportions)") + 
		ylab("-log10(adjusted p)") + 
		#ylim(c(0,(max.y+10))) + 
		#xlim(c(-1,1)) + 
		theme_bw()

	write.table(results, paste0(outFolder,"/", tolower( titles[n] ), 
		"dinucleotide_enrichment.tab"),quote=F, row.names=F,sep="\t")

plots[[n]] <- p

}

# plot <- ggplot(results, aes(x = log2FoldChange, y = 0, label = control.exons.x)) + 
# 		geom_text(size=7, position = position_jitter(height = 0.3)) +
# 		xlab("log2(cryptic vs null proportions)") +
# 		xlim(c(-1.5,1.5)) +
# 		ylim(c(-1,1)) +
# 		ylab("") + 
# 		ggtitle(paste0(gsub("_"," ",code)," (",species,")\nDinucleotide enrichment between flanked cryptic exons\nand adjacent intronic sequence")) + 
# 		theme_bw()

dinucleotide.pdf <- paste0(outFolder,"/dinucleotide_enrichment.pdf")
ggsave(dinucleotide.pdf)

save.image(paste0(RM.outFolder,code,"_Feature_Enrichment.RData") )
