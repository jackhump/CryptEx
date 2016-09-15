# Toxic Exons script
# Jack Humphrey
library(data.table)
library(stringr)
library(GenomicRanges)
library(Biostrings)
library(ggplot2)
library(dplyr)

options(echo=T)
# # TESTING
# code <- "Cleveland_TDP43"
# species <- "mouse"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/"
# # Cryptic Exon results - the result of sending the original dexseq output through the Splicing Analysis pipeline:
# method <- "cryptic_exons"
# dexseq.res <- "paper_TDP_mouse/Cleveland_TDP43/splice_junction_analysis/Cleveland_TDP43_CTL_TDP_splicing_analysis.tab"

# # Annotated exons - the output from the DEXSeq pipeline. 
# method <- "annotated_exons"
# dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/strict_500/dexseq/CTL_TDP/Cleveland_TDP43_CTL_TDP_SignificantExons.csv"
# control.SJs <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//paper_TDP_mouse/Cleveland_TDP43/splice_junction_analysis/Cleveland_TDP43_SJs_control.tab"
# case.SJs <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//paper_TDP_mouse/Cleveland_TDP43/splice_junction_analysis/Cleveland_TDP43_SJs_case.tab"

# # ENCODE dataset 1
# code <- "dataset_1"
# species <- "human"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1"
# method <- "cryptic_exons"
# dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_control_TDP_splicing_analysis.tab"

# # ENCODE dataset 2
# code <- "dataset_2"
# species <- "human"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2"
# method <- "cryptic_exons"
# dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/splice_junction_analysis/dataset_2_control_TDP_splicing_analysis.tab"


# Both Human TDP43 ENCODE for paper
code <- "both_ENCODE"
species <- "human"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/human/"
dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/human/both_datasets_unison_splicing_analysis.tab"
method <- "cryptic_exons"

# Both Mouse TDP43 for paper
code <- "both_mouse"
species <- "mouse"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/mouse/"
dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/mouse/both_datasets_unison_splicing_analysis.tab"
method <- "cryptic_exons"

# This script should use either the splice junction supported cryptic exons or the annotated exons as the input. 
# It then should determine the canonical intron that the exon sits within and then attributes upstream and downstream exon coordinates.
# It should then generate DNA sequences for the upstream, central and downstream exon
# and create transcripts with the central exon spliced in and out.
# the two transcripts are then translated in the correct reading frame (determined by the codon phase)
# and assessed for stop codons.

# canonical_junction_detector <- function(SJ.summary,results.df,mode="discovery"){
# 	GRanges_object <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T)
# 	if(mode == "discovery"){
# 		junctions.list <- apply(results.df, MAR=1,FUN=function(x) canonical_junction_query(x[10],x[11],x[12], GRanges_object))
# 		}
# 	if(mode == "replication"){
# 		junctions.list <- apply(results.df, MAR=1,FUN=function(x) canonical_junction_replication(x[10],x[15],x[16], GRanges_object))	
# 	}
# 	#output is a list of GRange objects - unuseable.
# 	junctions.list <- unlist(GRangesList(junctions.list))
# 	#convert into a dataframe, extracting the relevent information from the GRanges object.
# 	#names(GRanges) is a vector of rownames, confusingly.
# 	canonical.df <- data.frame(row.names=names(junctions.list),
# 			canonical.chr=seqnames(junctions.list),
# 			canonical.start=start(junctions.list),
# 			canonical.end=end(junctions.list),
# 			canonical.unique.count = score(junctions.list),
# 			canonical.strand = strand(junctions.list),
# 			intron.motif = mcols(junctions.list)[3])
# 	return(canonical.df)
# }

# this fixes any wayward gene names
# fix.gene.names <- function(results.df, annotation){
# 	# This is a function to find the most likely gene name for a given genomic range.
# 	# If multiple genes overlap the query range then output the names appended together with "+".
# 	# Get the strand information as well. If multiple genes are returned and their strands agree then output that strand.
# 	# If multiple genes are returned and the strands do not agree then output NA.
# 	gene.names.query <- function(chromosome, canonical.start, canonical.end, anno.GRange){
# 		gene.name <- anno.GRange[seqnames(anno.GRange) == chromosome & start(anno.GRange) <= as.numeric(canonical.start) & end(anno.GRange) >= as.numeric(canonical.end)]
# 		#if multiple gene names come up then collapse all the gene names into one string separated by "+".
# 		if(length(gene.name) > 1){
# 			mcols(gene.name)[,2] <- paste(mcols(gene.name)[,2],collapse="+")	
# 		# and strand as well!
# 		# if the different genes have differing strands then assign strand as NA.
# 			if(length(unique(as.list(strand(gene.name)))) != 1){
# 				strand(gene.name) <- NA
# 			}
# 			gene.name <- gene.name[1]
# 		}
# 		return(gene.name)
# 	}
# 	annotation <- as.data.frame(fread(annotation,header=T,stringsAsFactors=F))
# 	#sort out annotation and turn into a GRanges object. remove gm genes
# 	names(annotation)[4:5] <- c("start","end")
# 	anno.GRange <-  makeGRangesFromDataFrame(annotation,keep.extra.columns=T)
# 	fixed.gene.names <- apply(results.df, MAR=1,FUN=function(x) gene.names.query(x[10],x[11],x[12],anno.GRange)) 
# 	fixed.gene.names <- unlist(GRangesList(fixed.gene.names))
# 	fixed.gene.names <- data.frame(row.names=names(fixed.gene.names),gene <- mcols(fixed.gene.names)[,2], fixed.strand=strand(fixed.gene.names))
# 	names(fixed.gene.names)[1] <- "fixed.gene.id"
# 	#fixed.gene.names$fixed.strand <- gsub("+","1",fixed.gene.names$fixed.strand,fixed=T)
# 	#fixed.gene.names$fixed.strand <- gsub("-","-1",fixed.gene.names$fixed.strand,fixed=T)
# 	return(fixed.gene.names)
# }

translate_toxic <- function(seq.spliced.out, seq.spliced.in, fix.strand,upstream.codon.phase,downstream.codon.phase,upstream.fasta,downstream.fasta){
	seq.out <- DNAString(seq.spliced.out)
	codon.phase <- as.numeric(upstream.codon.phase)
	if(fix.strand == "-"){
		seq.out <- reverseComplement(seq.out)
		codon.phase <- as.numeric(downstream.codon.phase) 
	}
# Get the three different reading frames
	out.orf <- lapply(1:3, function(pos) subseq(seq.out, start=pos))
	out.orf <- suppressWarnings(lapply(out.orf, translate))
# One of the three reading frames is correct. One should have either no stop codons or one at the start so ignore the first letter when counting.
# in the absence of codon phase information then guess which is the correct reading frame by choosing the frame with the lowest or no stop codons
	readframe <- NA
	if( is.na(codon.phase) ){
		readframe <- lapply(out.orf, function(x) countPattern("*",x))
	}
	correct.frame <- ifelse(test = is.na(codon.phase),
		yes = c(1,2,3)[readframe < 1],
		no = codon.phase + 1 )
	
	out.correct <- out.orf[correct.frame]
	if( is.na(correct.frame) ){
		out.correct <- "no reading frame found"
	}
	spliced.out.final <- as.character(out.correct[[1]])
	
	if( is.na(correct.frame) ){
		spliced.in.final <- "no reading frame found"
		verdict.stop.codons <- "special case"
		verdict.frameshit <- "special case"
		return(c(correct.frame, spliced.out.final,spliced.in.final, verdict.stop.codons, verdict.frameshift))
	}

	
	seq.in <- DNAString(seq.spliced.in)
	if(fix.strand == "-") seq.in <- reverseComplement(seq.in)
	seq.in.orf <- suppressWarnings(translate(subseq(seq.in,start = correct.frame)))
	spliced.in.final <- as.character(seq.in.orf)

	# if a frameshift has occured then the tail end of the sequences will differ. 
	# first work out the number of amino acids coming from the downstream exon.
	tail.seq.length <- floor(str_length(downstream.fasta) / 3 )
	if(fix.strand == "-"){
		tail.seq.length <- floor(str_length(upstream.fasta) / 3)
	}
	# create substrings of tail end of the two sequences for comparisons. 
	spliced.in.tail <- subseq(spliced.in.final, start = length(spliced.in.final) - tail.seq.length)
	spliced.out.tail <- subseq(spliced.out.final, start = length(spliced.out.final) - tail.seq.length)

	
	# take out of function, do with the completed dataframe
	verdict.stop.codons <- ifelse(test = grepl("*", spliced.in.final,fixed=T),
		yes = "PTC",
		no = "nontoxic")
	verdict.frameshift <- ifelse(spliced.in.tail == spliced.out.tail, 
		yes = "frame conserved",
		no = "frame shifted")

	
	
return(c(correct.frame, spliced.out.final,spliced.in.final, verdict.stop.codons, verdict.frameshift))
#return(verdict)
}

# can I vectorise? Take instead the vectors of the numbers I need
translate_toxic_permutations_vectorized <- function(permute.seq.spliced.in,fix.strand,correct.frame, spliced.out){
	seq.in <- DNAStringSet(permute.seq.spliced.in)
	seq.in[fix.strand == "-"] <- reverseComplement(seq.in[fix.strand == "-"])
	seq.in.orf <- suppressWarnings(translate(subseq(seq.in,start = as.numeric(correct.frame))))
	spliced.in.final <- as.character(seq.in.orf)

	spliced.out.substring <- subseq(as.character(spliced.out), start = str_length(spliced.out) - 3)
	spliced.in.substring <- subseq(spliced.in.final, start = str_length(spliced.in.final) - 3)

	verdict.frameshift <- ifelse(test = spliced.out.substring == spliced.in.substring,
			yes = "frame conserved", no = "frame shifted")
	# take out of function, do with the completed dataframe
	verdict.stop.codon <- ifelse(test = grepl("*", spliced.in.final,fixed=T),
		yes = "PTC",
		no = "nontoxic")

	verdict.both <- data.frame(verdict.stop.codon = verdict.stop.codon,verdict.frameshift = verdict.frameshift)

	return(verdict.both)	
}




if(species=="mouse"){
	annotation <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/mouse/biomart/biomart_annotations_mouse.tab"
	ling.list <- c("2610507B11Rik",    "A230046K03Rik",	"Adipor2",	"Adnp2",	"Ahnak",	"Atraid",	"Cluh",	"Edem2",	"Ercc6",	"Fam21",	"Fam73a",	"Flnb",	"Ggct",	"Gsta4",	"Gtf2e2",	"Hace1",	"Hgsnat",	"Ift81",	"Lnp",	"Mettl6",	"Mib1",	"Mier1",	"Necap1",	"Nme6",	"Pir",	"Pno1",	"Ppp6c",	"Ptcd2",	"Pycr2",	"Sars",	"Smg5",	"Smg5",	"Snapc3",	"Spata13",	"Spata7",	"Spcs2",	"Sptbn4",	"Sulf1",	"Synj2bp",	"Tecpr1",	"Tnfaip1",	"Tnks2",	"Tnnt1",	"Trim8",	"Uggt2",	"Usp15",	"Usp15",	"Wbscr22",	"Zfp13",	"Zfp809")
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/mm10.fa"
	exon.gff <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gff"
    exon.gtf <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gtf"
	annotation <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/mouse/biomart/biomart_annotations_mouse.tab"
}

if(species=="human"){
	annotation <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human_hg38/biomart/biomart_annotations_human.tab"
	ling.list <- c("EPB41L4A",    "CEP72",	"INSR",	"FAM114A2",	"PFKP",	"ST5",	"RNFT2",	"RNFT2",	"ALX1",	"AGRN",	"AGRN",	"ATG4B",	"AGRN",	"AGRN",	"ST5",	"SETD5",	"KDELC2",	"MUC16",	"PKN1",	"IRF9",	"UPF2",	"GPSM2",	"XPO4",	"RASA4",	"RASA4B",	"PARP6",	"KRT7",	"TRAPPC12",	"RANBP1",	"HERC6",	"BLZF1",	"ZFP91",	"HDGFRP2",	"MAP3K8",	"SSFA2",	"CENPK",	"ITPR3",	"KYNU",	"IRF9",	"COL4A6",	"KYNU")
    genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa "
    annotation <- "/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human_hg38/biomart/biomart_annotations_human.tab"
    exon.gff <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gff"
    exon.gtf <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gtf"
}

# Annotated exons only. Load the DEXSeq results. Create the merged splice junction lists and extract the canonical intron.


toxic.outFolder <- paste(outFolder,"toxic_exons", sep = "/")
fasta.outFolder <- paste(toxic.outFolder, "DNAfasta",sep="/")
for(folder in c(toxic.outFolder, fasta.outFolder)){
if (! file.exists(folder)) dir.create(folder,recursive=T)
}

# load in a DEXSeq results file. Needs to have a list of exons with log2FCs and P-values


# merge together the splice junction files from the case and control bams separately. 


# query these splice junction lists to find the canonical intron. These coordinates are the upstream.end and the downstream.start


# OPTIONAL: if cryptic exon rather than annotated, query the splice junction lists to find the cryptic exon start and end.


# ALL THE BELOW CAN BE DONE WITH THE CRYPTIC.RES.CLASSIFIED TABLE


# query the GFF file to find the upstream.start and the downstream.end.
# canonical.start will match V5, canonical.end will match V4
#d <- read.table("Downloads/Cleveland_TDP43_CTL_TDP_crypt_res.tab",header=T)
# if(method == "annotated_exons"){
# 	d <- as.data.frame(fread(dexseq.res,header=T))
# 	# only interested in significantly altered annotated exons
# 	d <- subset(d, d$FDR < 0.05 & !grepl("i",d$exonID))
# 	# If splice junction lists have been compiled before in the splice junction lists then load them in:
# 	control_total_SJ_counts <- as.data.frame(fread(control.SJs))

# 	# for developing:
# 	#d_test <- head(d,250)

# 	canonical_results_control <- canonical_junction_detector(SJ.summary = control_total_SJ_counts, 
# 															 results.df = d, 
# 															 mode = "discovery")
# 	# remove any exons with length = 0 because why the bloody hell are they there?!
# 	d <- d[d$exon.end - d$exon.start != 0,] 
# 	# also remove any oddball "exons" formed from overlapping splice junctions (SJs don't make any sense)
# 	d <- d[d$downstream.case.cryptic.3prime - d$upstream.case.cryptic.5prime > 0,]

# 	d$canonical.chr <- canonical_results_control$canonical.chr[match(rownames(d),rownames(canonical_results_control))]
# 	d$canonical.start <- canonical_results_control$canonical.start[match(rownames(d),rownames(canonical_results_control))]
# 	d$canonical.end <- canonical_results_control$canonical.end[match(rownames(d),rownames(canonical_results_control))]
# 	d$canonical.strand <- canonical_results_control$canonical.strand[match(rownames(d),rownames(canonical_results_control))]
# 	d$canonical.intron.motif <- canonical_results_control$intron.motif[match(rownames(d),rownames(canonical_results_control))]
# 	# sort gene name and strand
# 	fixed.gene.names <- fix.gene.names(d,annotation)
# 	d$fix.gene.names <- as.character(fixed.gene.names$fixed.gene.id[match(row.names(d),row.names(fixed.gene.names))])
# 	d$fix.gene.names <- ifelse(test = is.na(d$fix.gene.names),yes = d$external_gene_id,no = d$fix.gene.names)
# 	d$fix.strand <- as.character(fixed.gene.names$fixed.strand[match(row.names(d),row.names(fixed.gene.names))])
# 	d$fix.strand <- ifelse(test=is.na(d$fix.strand),yes=d$strand,no=as.character(d$fix.strand))
# }

if(method == "cryptic_exons"){
	d <- read.table(dexseq.res,header=T)
	d <- subset(d,class == "SJ.SUPPORTED.UP")

	gtf <- as.data.frame(fread(exon.gtf))
	gff <- as.data.frame(fread(exon.gff))
	cds <- subset(gtf, V3 == "CDS")
	#rm(gtf)

	# Match in the upstream start and downstream end coordinates. First try the CDS gtf file as this will give us the codon phase.
	d$upstream.start <- cds$V4[ match( paste(d$canonical.chr, (d$canonical.start - 1) ), paste(cds$V1, cds$V5) )]
	d$downstream.end <- cds$V5[match( paste(d$canonical.chr, (d$canonical.end + 1) ), paste(cds$V1, cds$V4) )]
	
	# check that all the upstream and downstream exons are definitely coding.
	gtf.no.exons <- filter(gtf, V3 != "exon")
	d$upstream.exon.class <- gtf.no.exons$V3[ match( paste(d$canonical.chr, (d$canonical.start - 1) ), paste(gtf.no.exons$V1, gtf.no.exons$V5) )] 
	d$downstream.exon.class <- gtf.no.exons$V3[match( paste(d$canonical.chr, (d$canonical.end + 1) ), paste(gtf.no.exons$V1, gtf.no.exons$V4) )]

	# any coordinates not filled by the CDS (because they fall inside a UTR) are then matched using the GFF.
	d[is.na(d$upstream.start),]$upstream.start <- gff$V4[ match( paste(d[is.na(d$upstream.start),]$canonical.chr, (d[is.na(d$upstream.start),]$canonical.start - 1) ), paste(gff$V1, gff$V5) )]
	d[is.na(d$downstream.end),]$downstream.end <- gff$V5[match( paste(d[is.na(d$downstream.end),]$canonical.chr, (d[is.na(d$downstream.end),]$canonical.end + 1) ), paste(gff$V1, gff$V4) )]

	# any recalcitrant coordinates refusing to be matched should be ditched entirely.
	#d <- d[!is.na(d$upstream.start),]
	#d <- d[!is.na(d$downstream.end),]

	# match on codon phase for the upstream and downstream exons from the CDS gtf
	d$upstream.codon.phase <- cds$V8[match(
		paste(d$chr, d$upstream.start, d$canonical.start - 1),
		paste(cds$V1, cds$V4, cds$V5))]

	d$downstream.codon.phase <- cds$V8[match(
		paste(d$chr, d$canonical.end + 1, d$downstream.end),
		paste(cds$V1, cds$V4, cds$V5))]

	# create a reference table
	reference <- d
	# clean any missing or upstream or downstream values
	d <- filter(d, upstream.exon.class == "CDS" & downstream.exon.class == "CDS")
}

# this is for the cryptic exons that have splice junctions on either side. 
# regulated (annotated or cryptic) exon

if(method == "cryptic_exons"){
# If there is no strong splice junction from one end of the cryptic exon to the next exon then assume some kind of run on.
# In that case the boundary of the central exon at that side is continuous with the adjacent exon.

		
# need absolute coordinates for cryptic exon start and end. If an extension then assume that the extension is continuous with the adjacent exon.
d$cryptic.exon.start <- ifelse(test = (d$PSI.class == "CASSETTE.LIKE" | (d$PSI.class == "FIVEPRIME.BIAS" & d$fix.strand == "+") | (d$PSI.class == "THREEPRIME.BIAS" & d$fix.strand == "-") ), 
							yes = d$upstream.case.cryptic.5prime, no = d$canonical.start ) 
d$cryptic.exon.end <- ifelse(test = (d$PSI.class == "CASSETTE.LIKE" | (d$PSI.class == "FIVEPRIME.BIAS" & d$fix.strand == "-") | (d$PSI.class == "THREEPRIME.BIAS" & d$fix.strand == "+") ), 
							yes = d$downstream.case.cryptic.3prime, no = d$canonical.end ) 

central.bed <- data.frame(chr = d$canonical.chr,
				start = d$cryptic.exon.start,
				end = d$cryptic.exon.end - 1,
				name = d$fix.gene.names,
				score = ".",
				strand = d$fix.strand )
	}

# if(method == "annotated_exons"){
# 	central.bed <- data.frame(chr = d$canonical.chr,
# 				start = d$exon.start,
# 				end = d$exon.end,
# 				name = d$fix.gene.names,
# 				score = ".",
# 				strand = d$fix.strand )
# 	}

central.bed.out <- paste0(fasta.outFolder, "/",code,"_central.bed")
	# downstream exon

upstream.bed <- data.frame(chr = d$canonical.chr,
				start = d$upstream.start -1,
				end = d$canonical.start -1,
				name = d$fix.gene.names,
				score = ".",
				strand = d$fix.strand )
upstream.bed.out <- paste0(fasta.outFolder, "/",code,"_upstream.bed")

downstream.bed <- data.frame(chr = d$canonical.chr,
			start = d$canonical.end,
			end = d$downstream.end,
			name = d$fix.gene.names,
			score = ".",
			strand = d$fix.strand )
downstream.bed.out <- paste0(fasta.outFolder, "/",code,"_downstream.bed")

write.table(upstream.bed, upstream.bed.out, quote=F, col.names=F,row.names=F,sep="\t")
write.table(downstream.bed, downstream.bed.out, quote=F,col.names=F,row.names=F,sep="\t")
write.table(central.bed, central.bed.out, quote=F, col.names=F,row.names=F,sep="\t")


#bed.list <- c(upstream.bed, central.bed,downstream.bed)
#name.list <- c("upstream","central","downstream")
#bed.df <- data.frame(bed.list,name.list)
fasta.df <- data.frame(bed.out = c(upstream.bed.out, central.bed.out, downstream.bed.out),file.name = c("upstream","central","downstream"))


for(i in 1:nrow(fasta.df)) {
# write each bed file to the toxic exon folder
	fasta.name <- paste0(fasta.df[i,2],".fasta")
	fasta.out <- paste0( dirname(as.character(fasta.df[i,1])),"/", fasta.name)
# create bedtools command
	fasta.cmd <- paste0("bedtools getfasta -name -fi ",genome.fa," -bed ",fasta.df[i,1] ," -fo ",fasta.out)

# execute bedtools getfasta and read in the output	
	system(fasta.cmd)
	fasta <- as.data.frame(fread(fasta.out))
# remove all the ">geneid" rows
	fasta <- data.frame(sequence =  fasta[seq(1,nrow(fasta),2),])
	names(fasta) <- fasta.name
	assign(fasta.name,fasta)
}

# use BedTools to get the DNA sequence for each bed file, respecting strand. Add these all to the dataframe
d.fasta <- cbind(d, upstream.fasta,central.fasta,downstream.fasta)

# splice together either upstream + downstream (spliced out) or upstream + cryptic + downstream (spliced in)
d.fasta$seq.spliced.out <- paste0(d.fasta$upstream.fasta,d.fasta$downstream.fasta)
d.fasta$seq.spliced.in <- paste0(d.fasta$upstream.fasta,d.fasta$central.fasta, d.fasta$downstream.fasta)

d.fasta <- d.fasta[d.fasta$seq.spliced.in != "",]
d.fasta <- d.fasta[d.fasta$seq.spliced.out != "",]
# Looking at strand, reassign the upstream and downstream labels so that the sequences can be spliced together.


# Splice together the upstream and downstream sequence and translate. Which reading frame gives a continous sequence free of stop codons?




# Splice together the upstream, central and downstream exon and translate in the same reading frame.

if(method == "cryptic_exons"){
# for use with Cryptic exon results
# translate_toxic <- function(seq.spliced.out, seq.spliced.in, fix.strand,upstream.codon.phase,downstream.codon.phase,upstream.fasta,downstream.fasta)
toxic.seq <- apply(d.fasta, MAR=1, FUN = function(x) translate_toxic(x[64],x[65],x[46],x[57],x[58], x[61], x[63]) )
toxic.df <- data.frame(correct.frame = toxic.seq[1,], spliced.out = toxic.seq[2,], spliced.in = toxic.seq[3,], verdict.stop.codon = toxic.seq[4,], verdict.frameshift = toxic.seq[5,]) 
}


d.toxic <- cbind(d.fasta,toxic.df)

d.toxic.out <- paste0(toxic.outFolder,"/",code,"_protein_prediction.tab")
write.table(d.toxic, d.toxic.out, row.names=F, quote=F,sep='\t')
# is there a frame shift or not?
# determine length of downstream exon and determine whether the spliced in peptide contains the downstream exon's peptide sequence
# if so then there's no frameshift.

# permute the central sequence to see how many exons come up as non-toxic by chance. 

#d.toxic <- d.toxic[!is.na(d.toxic$correct.frame),]

#nontoxic <- sum(d.toxic$verdict.stop.codons == "nontoxic" & d.toxic$verdict.frameshift == "frame conserved" ) / nrow(reference)
#<- sample( seq(1:n), size = n, replace =F )
n <- length(d.toxic$central.fasta)
set.seed(1234)

permutation.list <- list()
for(go in 1:1000){
	#permute order of central exons
central.permuted <- d.toxic$central.fasta[ sample( seq(1:n), size = n, replace = F ) ]
d.toxic$permute.seq.spliced.in <- paste0(d.toxic$upstream.fasta, central.permuted, d.toxic$downstream.fasta)
#translate_toxic_permutations(permute.seq.spliced.in, fix.strand, correct.frame)
# use a modified translation function to assess how many of the permuted central exons cause premature stop codons/frameshifting
#translate_toxic_permutations_vectorized <- function(permute.seq.spliced.in,fix.strand,correct.frame, spliced.out)
percent.toxic <- translate_toxic_permutations_vectorized(d.toxic$permute.seq.spliced.in, d.toxic$fix.strand, d.toxic$correct.frame, d.toxic$spliced.out)
permutation.list[[go]] <- percent.toxic
}


# create one big dataframe of the permuted results and add a compound column
permutation <- do.call(what = rbind, args = permutation.list)
permutation <- mutate(permutation,  compound.verdict = paste(verdict.stop.codon, verdict.frameshift, sep =" / "), exon.type = "permuted exons")

# Plot the proportions of each verdict including those where the cryptic exon does not fall in the CDS.
#reference$verdict.stop.codon <- d.toxic$verdict.stop.codon[match(paste(reference$EnsemblID,reference$exonID), paste(d.toxic$EnsemblID, d.toxic$exonID))]
# reference$verdict.frameshift <- d.toxic$verdict.frameshift[match(paste(reference$EnsemblID,reference$exonID), paste(d.toxic$EnsemblID, d.toxic$exonID))]

# reference <- mutate(reference, compound.verdict = paste(verdict.stop.codon, verdict.frameshift, sep =" / ") ) %>% 
# 			 mutate(compound.verdict = gsub("NA / NA", "not in CDS", compound.verdict))


comparison <-  mutate(d.toxic, compound.verdict = paste(verdict.stop.codon, verdict.frameshift, sep =" / ") ) %>%
			   select(verdict.stop.codon, verdict.frameshift, compound.verdict ) %>%
			   mutate(exon.type = "cryptic exons")

prediction <- rbind(comparison, permutation)			   

prediction$compound.verdict <- factor(prediction$compound.verdict, levels = rev(c("PTC / frame shifted", "PTC / frame conserved", "nontoxic / frame shifted", "nontoxic / frame conserved" )))

exon.number <- as.numeric(table(prediction$exon.type))

exon.names <- c("cryptic exons", "permuted exons")

exon.labels <- paste0(exon.names, " (",exon.number,")")



prediction.pdf <- paste0(toxic.outFolder,"/",code,"_prediction_plot.pdf")

pdf(prediction.pdf)
prediction.plot <- ggplot(prediction, aes(x = exon.type, fill = compound.verdict)) + 
		geom_bar(position = "fill") + 
		scale_y_continuous(name = "", labels = scales::percent, limits = c(0,1)) + 
		ggtitle(paste(code, " ", species)) +
		scale_x_discrete(name="Exon type",
					limits=exon.names,
                    labels=exon.labels) + 
		scale_fill_manual(name = "in silico prediction", 
			values = c("chartreuse2","darkgoldenrod1","darkorange2","firebrick2"),
			breaks = c("PTC / frame shifted", "PTC / frame conserved", "nontoxic / frame shifted", "nontoxic / frame conserved" )) +
        xlab("") +
        theme_bw() + 
		theme(panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
          axis.line = element_line(colour = "black"),
          axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1) 
          )
# plot the distribution of the percentage of permuted exons that are clean too
#permutation <- data.frame(clean = percent.clean, class = "permutation")
#permutation.plot <- ggplot(permutation, aes(y = clean, x = class)) + geom_boxplot() +
#		scale_y_continuous(labels = scales::percent, limits = c(0,1))
print(prediction.plot)
# print(permutation.plot)
# grid.arrange(prediction.plot, permutation.plot, ncol = 2)
dev.off()

toxic.image <- paste0(toxic.outFolder,"/",code,"_toxic_prediction.Rdata")
save.image(toxic.image)


quit()




# perm.plot <- paste0(toxic.outFolder,"/",code,"_permutation_plot.pdf")

# v.df <- data.frame(proportion = verdicts)

# # Create simple barplot with confidence intervals of the binomial test
# # use PlyR to summarise the resulting dataframes - cryptic and permuted.



# #hist(verdicts, col="red")
# #abline(v=nontoxic)
# pdf(perm.plot)
# true_value <- paste0("true value: ",round(nontoxic,3) )
# bin_width <- 2 / nrow(d.toxic)
# p <- ggplot(v.df, aes(x=proportion)) +
# 	geom_histogram(fill = "red",colour = "black", alpha = 0.8,binwidth=bin_width) + 
# 	geom_vline(xintercept = nontoxic, colour = "blue" ) +
# 	xlab("proportion of non-toxic exons") +
# 	geom_label(label = true_value, aes(x= nontoxic, y = 500)) + 
# 	ggtitle(paste0( gsub("_"," ",code), " \ncentral exon permuted 1000 times") )
# print(p)
# dev.off()

# toxic.out <- paste0(toxic.outFolder,"/",code,"_toxic.tab")
# write.table(d.toxic,toxic.out,quote=F,sep="\t",row.names=F, col.names=T)

# exit
# # gtf exons for Dock4
# # 1272787 chr12 ensembl CDS 40834621 40834730  .  +  2
# # 1272789 chr12 ensembl CDS 40836630 40836702  .  +  0

# # This is a 0-based vs 1-based coordinate problem
# # upstream exon (for fasta): chr12  40834620  40834730  
# # downstream exon : 		   chr12  40836629  40836702 

