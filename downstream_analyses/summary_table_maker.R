## putting together the summary cryptic.tables
library(dplyr)
setwd("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets")

outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/summary_tables/"
if(!file.exists(outFolder)){dir.create(outFolder)}
species <- c("human","mouse")

for(i in 1:length(species)){
# for human
	if(species[i] == "human"){
		initial <- "human/both_datasets_unison_splicing_analysis.tab"
		expression <- c("human/cryptic_gene_expression/mRNA_cryptic_gene_expression.tab",
						"human/cryptic_gene_expression/total_cryptic_gene_expression.tab")
		expression.datasets <- c("mRNA","total")
		conservation <- "human/conservation/both_ENCODE_cryptic_exon_conservation.tab"
		prediction <- "human/toxic_exons/both_ENCODE_protein_prediction.tab" 
		splice.scoring <- "human/SJ_scoring/human_union_SJ_scores.tab"
		individual.datasets <- c("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_control_TDP_splicing_analysis.tab",
		"/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/splice_junction_analysis/dataset_2_control_TDP_splicing_analysis.tab"
		)	
	}

	if(species[i] == "mouse"){
		initial <- "mouse/both_datasets_unison_splicing_analysis.tab"
		expression <- c("mouse/cryptic_gene_expression/chiang_cryptic_gene_expression.tab",
						"mouse/cryptic_gene_expression/cleveland_cryptic_gene_expression.tab")
		expression.datasets <- c("ES_cell","adult_brain")
		conservation <- "mouse/conservation/both_mouse_cryptic_exon_conservation.tab"
		prediction <- "mouse/toxic_exons/both_mouse_protein_prediction.tab"
		individual.datasets <- c("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Chiang_processed/splice_junction_analysis/Chiang_processed_CTL_TDP_splicing_analysis.tab",
			"/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/splice_junction_analysis/Cleveland_TDP43_CTL_TDP_splicing_analysis.tab")
	}

	if(!all(file.exists(c(initial,expression,conservation,prediction,individual.datasets,splice.scoring)))){
		message("files are missing")
		break}

	# all the initial details about position, PSI and classification
	cryptic.table <- read.table(initial,header=T, sep = " ", stringsAsFactors = F)
	cryptic.table <- select(cryptic.table, gene_name = fix.gene.names, 
						   EnsemblID, 
						   exon_ID = exonID, 
						   chromosome,
						   exon_start = exon.start,
						   exon_end = exon.end, 
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
	overlap.num <- match(paste(cryptic.table$EnsemblID, cryptic.table$exon_ID), paste(d$EnsemblID, d$exonID))
	cryptic.table$overlap.one[!is.na(overlap.num)] <- expression.datasets[1]
	d <- read.table(individual.datasets[2],header=T)
	overlap.num <- match(paste(cryptic.table$EnsemblID, cryptic.table$exon_ID), paste(d$EnsemblID, d$exonID))
	cryptic.table$overlap.two[!is.na(overlap.num)] <- expression.datasets[2]

	# merge to create one observed in column
	cryptic.table <- mutate(cryptic.table, observed_in = paste(seen_by_Ling,overlap.one, overlap.two, sep=";"))
	cryptic.table$observed_in <- gsub(";$|^;|^;;", "", cryptic.table$observed_in)
	cryptic.table$observed_in <- gsub(";;", ";", cryptic.table$observed_in)
	cryptic.table <- select(cryptic.table, -(seen_by_Ling:overlap.two))

	# log2 FC for each dataset
	expression.data <- read.table(expression[1],header=T,sep="\t") 
	cryptic.table$expression.one <- expression.data$log2FoldChange[match(cryptic.table$EnsemblID, expression.data$EnsemblID)]
	expression.data <- read.table(expression[2],header=T,sep="\t") 
	cryptic.table$expression.two <- expression.data$log2FoldChange[match(cryptic.table$EnsemblID, expression.data$EnsemblID)]
	expression.names <- (ncol(cryptic.table) - 1)
	names(cryptic.table)[c(expression.names, expression.names + 1)] <- paste0(expression.datasets,"_log2FoldChange")

	# per exon PhyloP
	conservation.table <- read.table(conservation,header=T)
	cryptic.table$per_exon_phyloP_conservation <- conservation.table$phyloP.score[match(paste(cryptic.table$gene_name, cryptic.table$exon_ID), paste(conservation.table$gene.name, conservation.table$exonID))]

	# inclusion protein prediction
	prediction.table <- read.table(prediction, header=T,sep="\t")
	cryptic.table$inclusion_prediction <- paste(gsub("nontoxic","benign",prediction.table$verdict.stop.codon), prediction.table$verdict.frameshift, sep = "/")[ match(paste(cryptic.table$EnsemblID, cryptic.table$exon_ID), paste(prediction.table$EnsemblID, prediction.table$exonID))]

	# if human add 3' and 5' splice sites
	if(species[i] == "human"){
		splice.table <- read.table(splice.scoring, header=T, sep = "\t")
		cryptic.table$maxEnt_fiveprime_score <- splice.table$five_prime_splice_score[match(paste(cryptic.table$gene_name, cryptic.table$exon_ID), paste(splice.table$gene.name, splice.table$exonID))]	
		cryptic.table$maxEnt_threeprime_score <- splice.table$three_prime_splice_score[match(paste(cryptic.table$gene_name, cryptic.table$exon_ID), paste(splice.table$gene.name, splice.table$exonID))]	
	}
	out.table <- paste0(outFolder,species[i],"_summary_table.tsv")

	write.table(cryptic.table, out.table, sep="\t", quote=T, na = "NA", row.names = F, col.names = T)

	# create a cryptic exon BED file to load as an IGV track for pictures
	exon.bed <- select(cryptic.table, chromosome, exon_start, exon_end, gene_name, exon_ID, strand)
	exon.bed.out <- paste0(outFolder,species[i],"_cryptic_exons.bed")
	write.table(exon.bed, exon.bed.out, sep="\t", quote=F, row.names=F,col.names = F)

	# create a BED file of the canonical intron coordinates for IGV.

	IGV.table <- read.table(initial,header=T, sep = " ", stringsAsFactors = F)
	IGV.table <- select(IGV.table, 
						   chr = chromosome,
						   start = canonical.start, 
						   end = canonical.end,
						   gene.id = fix.gene.names,
						   score = exonID,
						   strand = fix.strand ) %>% 
				mutate(start = start - 200, end = (end + 200) ) 
	IGV.out.table <- paste0(outFolder,species[i],"_IGV_table.bed")
	write.table(IGV.table,IGV.out.table, col.names=F,row.names=F,sep="\t",quote=F)
						

}

#