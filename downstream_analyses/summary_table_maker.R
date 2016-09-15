## putting together the summary cryptic.tables
library(dplyr)
setwd("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets")

# for human
if(species = "human"){
	initial <- "human/both_datasets_unison_splicing_analysis.tab"
	expression <- c("human/cryptic_gene_expression/mRNA_cryptic_gene_expression.tab",
					"human/cryptic_gene_expression/total_cryptic_gene_expression.tab")
	expression.datasets <- c("mRNA","total")
	conservation <- "human/conservation/both_ENCODE_cryptic_exon_conservation.tab"
	prediction <- "human/toxic_exons/both_ENCODE_protein_prediction.tab" 
	splice.scoring <- "human/SJ_scoring/human_union_SJ_scores.tab"
	individual.datasets <- c(c("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_control_TDP_splicing_analysis.tab",
	"/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/splice_junction_analysis/dataset_2_control_TDP_splicing_analysis.tab"))	
}

# all the initial details about position, PSI and classification
cryptic.table <- read.table(initial,header=T, sep = " ", stringsAsFactors = F)
cryptic.table <- select(cryptic.table, gene.name = fix.gene.names, 
					   EnsemblID, 
					   exonID, 
					   chromosome:exon.end, 
					   strand = fix.strand, 
					   fiveprime_delta_PSI = upstream_delta_psi, 
					   threeprime_delta_PSI = downstream_delta_psi, 
					   classification = PSI.class, 
					   seen_by_Ling = in.Ling )
# add which datasets found in
cryptic.table$seen_by_Ling[is.na(cryptic.table$seen_by_Ling)] <- ""
cryptic.table$seen_by_Ling[nchar(as.character(cryptic.table$seen_by_Ling))!=0 & !is.na(cryptic.table$seen_by_Ling)] <- "Ling"

cryptic.table$overlap.one <- ""
cryptic.table$overlap.two <- ""
d <- read.table(individual.datasets[1],header=T)
overlap.num <- match(paste(cryptic.table$EnsemblID, cryptic.table$exonID), paste(d$EnsemblID, d$exonID))
cryptic.table$overlap.one[!is.na(overlap.num)] <- expression.datasets[1]
d <- read.table(individual.datasets[2],header=T)
overlap.num <- match(paste(cryptic.table$EnsemblID, cryptic.table$exonID), paste(d$EnsemblID, d$exonID))
cryptic.table$overlap.two[!is.na(overlap.num)] <- expression.datasets[2]
cryptic.table <- rename(cryptic.table, expression.data[1] = overlap.one, expression.data[2] = overlap.two )

cryptic.table <- mutate(cryptic.table, observed_in = paste(seen_by_Ling,overlap.one, overlap.two, sep=";"))
cryptic.table$observed_in <- gsub(";$|^;|;;", "", cryptic.table$observed_in)
cryptic.table <- select(cryptic.table, -(seen_by_Ling:overlap.two))

# log2 FC for each dataset
expression.data <- read.table(expression[1],header=T,sep="\t") 
cryptic.table$expression.one <- expression.data$log2FoldChange[match(cryptic.table$EnsemblID, expression.data$EnsemblID)]
expression.data <- read.table(expression[2],header=T,sep="\t") 
cryptic.table$expression.two <- expression.data$log2FoldChange[match(cryptic.table$EnsemblID, expression.data$EnsemblID)]
expression.names <- (ncol(cryptic.table) - 1)
names(cryptic.table)[c(expression.names, expression.names + 1)] <- paste0(expression.datasets,"_gene_expression")
# per exon PhyloP
conservation.table <- read.table(conservation,header=T)
cryptic.table$PhyloP_conservation <- conservation.table$phyloP.score[match(paste(cryptic.table$gene.name, cryptic.table$exonID), paste(conservation.table$gene.name, conservation.table$exonID))]
# inclusion prediction
prediction.table <- read.table(prediction, header=T,sep="\t")
cryptic.table$inclusion_prediction <- paste(gsub("nontoxic","benign",prediction.table$verdict.stop.codon), prediction.table$verdict.frameshift, sep = " + ")[ match(paste(cryptic.table$EnsemblID, cryptic.table$exonID), paste(prediction.table$EnsemblID, prediction.table$exonID))]
# if human add 3' and 5' splice sites
if(species == "human"){
	splice.table <- read.table(splice.scoring, header=T, sep = "\t")
	cryptic.table$fiveprime_splice_site <- splice.table$five_prime_splice_site[match(paste(cryptic.table$EnsemblID, cryptic.table$exonID), paste(splice.table$EnsemblID, splice.table$exonID))]	
}


# for mouse

# all the initial details about position, PSI and classification


# add which datasets found in

# log2 FC for each dataset

# per exon PhyloP

# inclusion prediction

# if human add 3' and 5' splice sites