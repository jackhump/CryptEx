# max-ent splice site scoring

library(data.table)
library(stringr)
library(gridExtra)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(tidyr)

# what do I need? A list of cryptic introns!

# code <- "dataset_1"
# species <- "human"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1"
# splicing_analysis.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_control_TDP_splicing_analysis.tab"
# cryptic.exons.bed <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/RepeatMasker/cryptic_exons.exons.bed"
# cryptic.nulls.bed <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/RepeatMasker/cryptic_nulls.exons.bed"
# case.SJs <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_SJs_case.tab"



# code <- "Cleveland_TDP43"
# species <- "mouse"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/"
# splicing_analysis.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/splice_junction_analysis/Adult_mouse_brain_TDP43_knockdown_CTL_TDP_splicing_analysis.tab"
# cryptic.exons.bed <- paste0(outFolder,"RepeatMasker/cryptic_exons.exons.bed")
# cryptic.nulls.bed <- paste0(outFolder,"RepeatMasker/cryptic_nulls.exons.bed")
# case.SJs <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//paper_TDP_mouse/Cleveland_TDP43/splice_junction_analysis/Cleveland_TDP43_SJs_case.tab"

code <- "human_union"
species <- "human"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/human/"
splicing_analysis.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/human/both_datasets_unison_splicing_analysis.tab"
cryptic.exons.bed <- paste0(outFolder, "Feature_Enrichment/cryptic_exons.exons.bed")
cryptic.nulls.bed <- paste0(outFolder, "Feature_Enrichment/cryptic_nulls.exons.bed")
case.SJs <- paste0(outFolder,"case_SJs_merge.tab")

SJ.scoring.outFolder <- paste0(outFolder,"/SJ_scoring/")
if (! file.exists(SJ.scoring.outFolder)) dir.create(SJ.scoring.outFolder)

# For working out the coordinates of the adjacent annotated exon:
adjacent_exon_detector <- function(SJ.summary,d){
	exon.spacer <- 500
	adjacent_intron_query <- function(canonical.chr, canonical.start, canonical.end, SJ.GRange){
		canonical.start <- as.numeric(canonical.start)
		canonical.end <- as.numeric(canonical.end)
		junction <- list()
		junction <- SJ.GRange[
						#( width(SJ.GRange) > cryptic.flanked.length &
						seqnames(SJ.GRange) == canonical.chr  &
						# for capturing left adjacent introns 
						 ( ( end(SJ.GRange) < canonical.start ) & ( end(SJ.GRange) > (canonical.start - exon.spacer) ) |
						  ( start(SJ.GRange) > canonical.end  & start(SJ.GRange) < (canonical.end + exon.spacer) ) ) ]
		# insist that the intron must be longer than the length of the cryptic exon
		junction <- head(junction[order(score(junction),decreasing=T)],1)
	return(junction)
	}
	SJ.GRange <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T)
	junctions.list <- apply(d, MAR=1,FUN=function(x) adjacent_intron_query(x[10],x[15],x[16], SJ.GRange))	
	#output is a list of GRange objects - unuseable.
	junctions.list <- unlist(GRangesList(junctions.list))
	#convert into a dataframe, extracting the relevent information from the GRanges object.
	#names(GRanges) is a vector of rownames, confusingly.
	adjacent.df <- data.frame(row.names=names(junctions.list),
			chr=seqnames(junctions.list),
			adjacent.intron.start=start(junctions.list),
			adjacent.intron.end=end(junctions.list),
			adjacent.unique.count = score(junctions.list),
			strand = strand(junctions.list),
			intron.motif = mcols(junctions.list)[3])
	return(adjacent.df)
}


if(species == "mouse"){
	#iCLIP.peaks <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/iCLIP/mouse/peaks_id81678_rnd100_flank15_fdr0.05_group_5533_F210I-WT-rep-2-3_sum_G_mm10--ensembl75_from_5447-5448_bedGraph-cDNA-hits-in-genome.bed.gz_lowFDR_stranded.bed" 
	iCLIP.peaks <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/iCLIP/mouse/peaks_id81678_rnd100_flank15_fdr0.05_group_5533_F210I-WT-rep-2-3_sum_G_mm10--ensembl75_from_5447-5448_bedGraph-cDNA-hits-in-genome.bed.gz_lowFDR_stranded.bed" 
	RepeatMasker.bed <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/RepeatMasker/mm10_repeat_masker.bed"
	#phyloP.bw <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/phyloP/mm10.60way.phyloP60way.bw"
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/mm10.fa"
}

if(species == "human"){
	RepeatMasker.bed <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/RepeatMasker/hg38_repeat_masker.bed"
	iCLIP.peaks <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_human/iCLIP/human_ES_cells/peaks_id71003_rnd100_flank15_fdr0.05_regionAsOne_20100222_LUjt3_4_hg38_ensembl59_G_iCLIP_TDP-43_Embrionic-Stem-Cells_bedGraph-cDNA-hits-in-genome.bed_lowFDR_stranded.bed"
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa "
}

d <- read.table(splicing_analysis.res,header=T,stringsAsFactors = F)
d$fix.strand <- gsub("-1","-",d$fix.strand,fixed=T)
d$fix.strand <- gsub("1","+",d$fix.strand,fixed=T)
d$canonical.start <- as.numeric(d$canonical.start)
d$canonical.end <- as.numeric(d$canonical.end)

# subset just the cryptic exons

d <- filter(d, class == 'SJ.SUPPORTED.UP')

d$cryptic.start <- ifelse(d$upstream.case.cryptic.5prime > 0,
					yes = d$upstream.case.cryptic.5prime,
					no = d$exon.start)
d$cryptic.end <- ifelse(d$downstream.case.cryptic.3prime > 0,
					yes = d$downstream.case.cryptic.3prime,
					no = d$exon.end)
good.exons <- d$cryptic.end > d$cryptic.start
d <- d[good.exons,]

# add on extra flanking sequence to capture the canonical splice sites too!
cryptic.introns.bed <- data.frame(chr = d$chr, 
	intron.start = d$canonical.start - 4, 
	intron.end = d$canonical.end + 3, 
	gene.id = paste0(d$fix.gene.names,"_",d$exonID),
	exonID = ".", 
	strand = d$fix.strand, 
	cryptic.start = d$cryptic.start,
	cryptic.end = d$cryptic.end,
	PSI.class = d$PSI.class)

# use BedTools to retrieve the sequence of these introns.
cryptic.introns.bed.out <- paste0(SJ.scoring.outFolder,"/",code,".cryptic.introns.bed")
write.table(cryptic.introns.bed, cryptic.introns.bed.out, row.names=F,sep="\t",col.names=F,quote=F )
fasta.out <-  paste0(SJ.scoring.outFolder,"/",code,".cryptic.introns.fasta")
MOTIF.cmd <- paste0("bedtools getfasta -s -name -fi ",genome.fa," -bed ",cryptic.introns.bed.out," -fo ",fasta.out)
system(MOTIF.cmd)

introns.fasta <- as.data.frame(fread(fasta.out,header=F,stringsAsFactors=F))
#
introns.fasta <- data.frame(gene.id =  introns.fasta[seq(1,nrow(introns.fasta) - 1,2),],
							sequence = introns.fasta[seq(2,nrow(introns.fasta),2),])
introns.fasta$gene.id <- gsub(">","",introns.fasta$gene.id)

cryptic.introns.bed$sequence <- as.character(introns.fasta$sequence[match(cryptic.introns.bed$gene.id, introns.fasta$gene.id)])

cryptic.intron.length <- cryptic.introns.bed$intron.end - cryptic.introns.bed$intron.start
#cryptic.introns.bed <- head(cryptic.introns.bed)

# for each intron create a table of all the possible ninemers.
ninemers.list <- list()
for(i in 1:nrow(cryptic.introns.bed)){	
	intron.length <- cryptic.intron.length[i]
	ninemers <- data.frame(gene.id = cryptic.introns.bed$gene.id[i], start = seq(1,intron.length -8,1))
	ninemers$subsequence <- str_sub(cryptic.introns.bed$sequence[i], ninemers$start,ninemers$start + 8)
	ninemers$splice_site <- "five.prime"
	ninemers$start <- ninemers$start - 4
	ninemers.list[[i]] <- ninemers
} 
ninemers.unique <- do.call(what = rbind, args = ninemers.list)
ninemers.unique <- unique(ninemers.unique$subsequence)

ninemers.unique.out <- paste0(SJ.scoring.outFolder,"/",code,"_unique_ninemers.txt")
write.table(ninemers.unique, ninemers.unique.out,quote=F,row.names=F,col.names=F)

setwd("/SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/")
cmd <- paste0("perl /SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/score5.pl ",ninemers.unique.out)
res <- fread(cmd)
names(res) <- c("subsequence", "score")

for(i in 1:length(ninemers.list)){
	ninemers.list[[i]]$score <- res$score[match(ninemers.list[[i]]$subsequence, res$subsequence)]
}

# scoring the 3' sites - 23mers
twentythreemers.list <- list()
for(i in 1:nrow(cryptic.introns.bed)){	
	intron.length <- cryptic.intron.length[i]
	twentythreemers <- data.frame(gene.id = cryptic.introns.bed$gene.id[i], start = seq(1,intron.length -22,1))
	twentythreemers$subsequence <- str_sub(cryptic.introns.bed$sequence[i], twentythreemers$start,twentythreemers$start + 22)
	twentythreemers$splice_site <- "three.prime"
	twentythreemers.list[[i]] <- twentythreemers
} 
twentythreemers.unique <- do.call(what = rbind, args = twentythreemers.list)
twentythreemers.unique <- unique(twentythreemers.unique$subsequence)

twentythreemers.unique.out <- paste0(SJ.scoring.outFolder,"/",code,"_unique_twentythreemers.txt")
write.table(twentythreemers.unique, twentythreemers.unique.out,quote=F,row.names=F,col.names=F)

setwd("/SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/")
cmd <- paste0("perl /SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/score3.pl ",twentythreemers.unique.out)
res <- fread(cmd)

names(res) <- c("subsequence", "score")
for(i in 1:length(twentythreemers.list)){
	twentythreemers.list[[i]]$score <- res$score[match(twentythreemers.list[[i]]$subsequence, res$subsequence)]
}

all.scores.fiveprime <- do.call(what=rbind, args = ninemers.list )
all.scores.threeprime <- do.call(what=rbind, args = twentythreemers.list )


# read in flanked beds - why can't I just do this from the initial unflanked bed?
# cryptic.exons and cryptic.introns are just what you expect. cryptic.nulls are randomly selected intronic chunks for comparison. Seems pointless to me as you will just get a uniform distribution of splice scores.

#cryptic.introns <- as.data.frame(fread(cryptic,header=F))
cryptic.exons <- as.data.frame(fread(cryptic.exons.bed,header=F))
cryptic.nulls <- as.data.frame(fread(cryptic.nulls.bed,header=F))
cryptic.introns <- cryptic.introns.bed
#all.scores <- as.data.frame(fread("~/Documents/Cryptic_Exons/SJ_scoring/dataset_1_all_scores.tab",header=T))



names(cryptic.exons) <- c("chr","start","end","gene.name","exonID","strand","iteration")
names(cryptic.nulls) <- c("chr","start","end","gene.name","exonID","strand","iteration")
#cryptic.introns <- cryptic.introns[,1:6]
#names(cryptic.introns) <- c("chr","start","end","gene.name","exonID","strand")

cryptic.introns$exonID <- str_split_fixed(cryptic.introns$gene.id,pattern="_",n=2)[,2]
cryptic.introns$gene.name <- str_split_fixed(cryptic.introns$gene.id,pattern="_",n=2)[,1]

# the beds I'm using are the 100bp flanked beds. Remove the flanks!
cryptic.exons$start <- cryptic.exons$start + 100
cryptic.exons$end <- cryptic.exons$end - 101
cryptic.nulls$start <- cryptic.nulls$start + 100
cryptic.nulls$end <- cryptic.nulls$end - 101

write.table(cryptic.exons, paste0(outFolder,"/SJ_scoring/cryptic.exons.bed"),sep="\t",col.names=F,row.names=F,quote=F)

# calculate the length of each intron and convert the cryptic exon coordinates into the relative position along the intron.
# match on the intron start for each cryptic exon.
cryptic.exons$intron.start <- cryptic.introns$intron.start[match(paste(cryptic.exons$gene.name, cryptic.exons$exonID), paste(cryptic.introns$gene.name,cryptic.introns$exonID))] 
cryptic.exons$intron.end <- cryptic.introns$intron.end[match(paste(cryptic.exons$gene.name, cryptic.exons$exonID), paste(cryptic.introns$gene.name,cryptic.introns$exonID))] 

# only match 5'ss for cryptic exons that have them - ie 3' extension or cassette like.

cryptic.exons$PSI.class <- cryptic.introns$PSI.class[match(paste(cryptic.exons$gene.name, cryptic.exons$exonID), paste(cryptic.introns$gene.name,cryptic.introns$exonID))] 

cryptic.exons$five_prime_splice_position <- ifelse(cryptic.exons$strand == "+", 
						yes = (cryptic.exons$end - cryptic.exons$intron.start) - 6,
						no = (cryptic.exons$intron.end - (cryptic.exons$start + 6)) )

cryptic.exons$five_prime_splice_score <- ifelse(cryptic.exons$PSI.class != "FIVEPRIME.BIAS", 
	yes = all.scores.fiveprime$score[match(
			paste0(cryptic.exons$gene.name,"_",cryptic.exons$exonID,"_",cryptic.exons$five_prime_splice_position),
			paste0(all.scores.fiveprime$gene.id,"_",all.scores.fiveprime$start))],
	no = NA)

cryptic.exons$five_prime_splice_seq <- all.scores.fiveprime$subsequence[match(
			paste0(cryptic.exons$gene.name,"_",cryptic.exons$exonID,"_",cryptic.exons$five_prime_splice_position),
			paste0(all.scores.fiveprime$gene.id,"_",all.scores.fiveprime$start))]

# THREE PRIME SPLICE SITE

cryptic.exons$three_prime_splice_position <- as.numeric(ifelse(cryptic.exons$strand == "+", 
						yes = (cryptic.exons$start - cryptic.exons$intron.start) - 19,
						no = (cryptic.exons$intron.end - (cryptic.exons$end + 19))))

cryptic.exons$three_prime_splice_score <- ifelse(cryptic.exons$PSI.class != "THREEPRIME.BIAS", 
	yes = all.scores.threeprime$score[match(
			paste0(cryptic.exons$gene.name,"_",cryptic.exons$exonID,"_",cryptic.exons$three_prime_splice_position),
			paste0(all.scores.threeprime$gene.id,"_",all.scores.threeprime$start))],
	no = NA)

cryptic.exons$three_prime_splice_seq <- all.scores.threeprime$subsequence[match(
			paste0(cryptic.exons$gene.name,"_",cryptic.exons$exonID,"_",cryptic.exons$three_prime_splice_position),
			paste0(all.scores.threeprime$gene.id,"_",all.scores.threeprime$start))]

cryptic.exons$exon.type <- "cryptic"
# select just the relevent columns for graphing
cryptic.exons.slim <- data.frame(exon.type = cryptic.exons$exon.type,
		gene.id = paste(cryptic.exons$gene.name, cryptic.exons$exonID,sep="_"), 
		five_prime_splice_score = cryptic.exons$five_prime_splice_score, 
		three_prime_splice_score = cryptic.exons$three_prime_splice_score)

# CRYPTIC NULLS

cryptic.nulls <- filter(cryptic.nulls, iteration == 45)
cryptic.nulls$intron.start <- cryptic.introns$intron.start[match(paste(cryptic.nulls$gene.name, cryptic.nulls$exonID), paste(cryptic.introns$gene.name,cryptic.introns$exonID))] 
cryptic.nulls$intron.end <- cryptic.introns$intron.end[match(paste(cryptic.nulls$gene.name, cryptic.nulls$exonID), paste(cryptic.introns$gene.name,cryptic.introns$exonID))] 

cryptic.nulls$PSI.class <- cryptic.introns$PSI.class[match(paste(cryptic.nulls$gene.name, cryptic.nulls$exonID), paste(cryptic.introns$gene.name,cryptic.introns$exonID))] 



cryptic.nulls$five_prime_splice_position <- ifelse(cryptic.nulls$strand == "+", 
						yes = (cryptic.nulls$end - cryptic.nulls$intron.start) - 6,
						no = (cryptic.nulls$intron.end - (cryptic.nulls$start + 6)) )

cryptic.nulls$five_prime_splice_score <- ifelse(cryptic.nulls$PSI.class != "FIVEPRIME.BIAS", 
		yes = all.scores.fiveprime$score[match(
			paste0(cryptic.nulls$gene.name,"_",cryptic.nulls$exonID,"_",cryptic.nulls$five_prime_splice_position),
			paste0(all.scores.fiveprime$gene.id,"_",all.scores.fiveprime$start))],
		no = NA)

cryptic.nulls$five_prime_splice_seq <- all.scores.fiveprime$subsequence[match(
			paste0(cryptic.nulls$gene.name,"_",cryptic.nulls$exonID,"_",cryptic.nulls$five_prime_splice_position),
			paste0(all.scores.fiveprime$gene.id,"_",all.scores.fiveprime$start))]

# THREE PRIME SPLICE SITE

cryptic.nulls$three_prime_splice_position <- as.numeric(ifelse(cryptic.nulls$strand == "+", 
						yes = (cryptic.nulls$start - cryptic.nulls$intron.start) - 19,
						no = (cryptic.nulls$intron.end - (cryptic.nulls$end + 19))))

cryptic.nulls$three_prime_splice_score <- ifelse(cryptic.nulls$PSI.class != "THREEPRIME.BIAS", 
		yes = all.scores.threeprime$score[match(
			paste0(cryptic.nulls$gene.name,"_",cryptic.nulls$exonID,"_",cryptic.nulls$three_prime_splice_position),
			paste0(all.scores.threeprime$gene.id,"_",all.scores.threeprime$start))],
		no = NA)

cryptic.nulls$three_prime_splice_seq <- all.scores.threeprime$subsequence[match(
			paste0(cryptic.nulls$gene.name,"_",cryptic.nulls$exonID,"_",cryptic.nulls$three_prime_splice_position),
			paste0(all.scores.threeprime$gene.id,"_",all.scores.threeprime$start))]

cryptic.nulls$exon.type <- "null"
# select just the relevent columns for graphing
cryptic.nulls.slim <- data.frame(exon.type = cryptic.nulls$exon.type,
		gene.id = paste(cryptic.nulls$gene.name, cryptic.nulls$exonID,sep="_"), 
		five_prime_splice_score = cryptic.nulls$five_prime_splice_score, 
		three_prime_splice_score = cryptic.nulls$three_prime_splice_score)

# CANONICAL SPLICE SITES

canonical <- cryptic.exons[,1:10]
# FIVE PRIME - this is just -3 on the position
canonical$five_prime_splice_score <- all.scores.fiveprime$score[match(
	paste0(canonical$gene.name,"_",canonical$exonID,"_","-3"),
	paste0(all.scores.fiveprime$gene.id,"_", all.scores.fiveprime$start))]

# THREE PRIME - this is the final position for each gene.
# use dplyr to output a list of final positions for each cryptic exon.
canonical.threeprime.pos <- group_by(all.scores.threeprime, gene.id) %>% summarise(position = max(start))

canonical$three_prime_splice_position <- canonical.threeprime.pos$position[match(
	paste0(canonical$gene.name,"_",canonical$exonID),
	paste0(canonical.threeprime.pos$gene.id))]

canonical$three_prime_splice_score <- all.scores.threeprime$score[match(
	paste0(canonical$gene.name,"_",canonical$exonID,"_",canonical$three_prime_splice_position),
	paste0(all.scores.threeprime$gene.id,"_", all.scores.threeprime$start))]
canonical.slim <- data.frame(exon.type = "canonical",
		gene.id = paste(canonical$gene.name, canonical$exonID,sep="_"), 
		five_prime_splice_score = canonical$five_prime_splice_score, 
		three_prime_splice_score = canonical$three_prime_splice_score)

## ANNOTATED ADJACENT EXONS

SJ.summary <- as.data.frame(fread(case.SJs))

adjacent.introns.junctions <- adjacent_exon_detector(SJ.summary,d)

adjacent.introns.junctions$gene.name <- d$fix.gene.names[match(rownames(adjacent.introns.junctions),rownames(d))]
adjacent.introns.junctions$exonID <- d$exonID[match(rownames(adjacent.introns.junctions),rownames(d))]
adjacent.introns.junctions$gene.strand <- d$fix.strand[match(rownames(adjacent.introns.junctions),rownames(d))]


adjacent.introns.junctions$canonical.start <- d$canonical.start[match(rownames(adjacent.introns.junctions),rownames(d))]
adjacent.introns.junctions$canonical.end <- d$canonical.end[match(rownames(adjacent.introns.junctions),rownames(d))]


# work out the adjacent exon coordinates
adjacent.exons <- data.frame(chr = adjacent.introns.junctions$chr,
	adjacent.exon.start = ifelse(adjacent.introns.junctions$adjacent.intron.start < adjacent.introns.junctions$canonical.end,
								yes = adjacent.introns.junctions$adjacent.intron.end,
								no = adjacent.introns.junctions$canonical.end ),
	adjacent.exon.end = ifelse(adjacent.introns.junctions$adjacent.intron.start < adjacent.introns.junctions$canonical.end,
								yes = adjacent.introns.junctions$canonical.start - 1,
								no = adjacent.introns.junctions$adjacent.intron.start - 1),
	strand = adjacent.introns.junctions$gene.strand,
	gene.id = paste0(adjacent.introns.junctions$gene.name,"_",adjacent.introns.junctions$exonID)
	)

adjacent.exons.threeprime.bed <- data.frame(chr = adjacent.exons$chr,
	start = ifelse(adjacent.introns.junctions$strand == "+",
				yes = adjacent.exons$adjacent.exon.start - 20,
				no = adjacent.exons$adjacent.exon.end  - 3),
	end = ifelse(adjacent.introns.junctions$strand == "+",
				yes = adjacent.exons$adjacent.exon.start + 3,
				no = adjacent.exons$adjacent.exon.end + 20),
	gene.id = adjacent.exons$gene.id,
	blank = ".",
	strand = adjacent.exons$strand
	)

adjacent.exons.fiveprime.bed <- data.frame(chr = adjacent.exons$chr,
	start = ifelse(adjacent.introns.junctions$strand == "+",
				yes = adjacent.exons$adjacent.exon.end - 3,
				no = adjacent.exons$adjacent.exon.start  - 6),
	end = ifelse(adjacent.introns.junctions$strand == "+",
				yes = adjacent.exons$adjacent.exon.end + 6,
				no = adjacent.exons$adjacent.exon.start + 3),
	gene.id = adjacent.exons$gene.id,
	blank = ".",
	strand = adjacent.exons$strand
	)

adjacent.bed.names <- c("adjacent.exons.threeprime.bed","adjacent.exons.fiveprime.bed")
adjacent.bed.list <- list(adjacent.exons.threeprime.bed, adjacent.exons.fiveprime.bed)
for(i in 1:2){
	# create bed file and fasta file
	adjacent.bed.out <- paste0(SJ.scoring.outFolder,adjacent.bed.names[i])
	adjacent.fasta.out <- paste0(SJ.scoring.outFolder,adjacent.bed.names[i],".fasta")
	# save bed to disk
	write.table(adjacent.bed.list[i],adjacent.bed.out,quote=F, row.names =F ,col.names=F,sep="\t")
	# create getfasta command and execute with system()
	fasta.cmd <- paste0("bedtools getfasta -s -name -fi ",genome.fa," -bed ",adjacent.bed.out," -fo ",adjacent.fasta.out)
	system(fasta.cmd)
}


# THREE PRIME
threeprime.fasta <- paste0(paste0(SJ.scoring.outFolder,adjacent.bed.names[1],".fasta"))
setwd("/SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/")
cmd <- paste0("perl /SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/score3.pl ",threeprime.fasta)
adjacent.exons.threeprime.res <- fread(cmd)
adjacent.exons$three_prime_splice_score <- adjacent.exons.threeprime.res$V2[match(rownames(adjacent.exons),rownames(adjacent.exons.threeprime.res))]
# FIVE PRIME
fiveprime.fasta <- paste0(paste0(SJ.scoring.outFolder,adjacent.bed.names[2],".fasta"))
setwd("/SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/")
cmd <- paste0("perl /SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/score5.pl ",fiveprime.fasta)
adjacent.exons.threeprime.res <- fread(cmd)
adjacent.exons$five_prime_splice_score <- adjacent.exons.threeprime.res$V2[match(rownames(adjacent.exons),rownames(adjacent.exons.threeprime.res))]

adjacent.exons$exon.type <- "adjacent"

adjacent.exons.slim <- data.frame(exon.type = adjacent.exons$exon.type,
	gene.id = adjacent.exons$gene.id, 
		five_prime_splice_score = adjacent.exons$five_prime_splice_score, 
		three_prime_splice_score = adjacent.exons$three_prime_splice_score)

# RANDOM SEQUENCE WITH INVARIANT AG/GT 
	# are the cryptic splice sites more complicated than just having the AG and GT?
# FIVE PRIME SPLICE SITES
	# create a vector of 1000 random ninemers with GT at positions 4:5
	# DNA_BASES is a built in vector of A G C T
DNA_BASES <- c("A","C","G","T")
ninemers.random <- c()
for(i in 1:1000){
	s <- sample(DNA_BASES, 9, replace = T)
	s[4] <- "G"
	s[5] <- "T"
	s <- paste(s, collapse = "")
	ninemers.random[i] <- s
}
# create table and write it out
ninemers.random.table <- data.frame(seq = ninemers.random)
ninemers.random.out <- paste0(SJ.scoring.outFolder,"random.ninemers.fasta")
write.table(ninemers.random.table, ninemers.random.out, col.names=F,row.names=F,quote=F,sep="\t")
# run maxEnt 5' splice site analysis
setwd("/SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/")
cmd <- paste0("perl /SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/score5.pl ",ninemers.random.out)
ninemers.random.res <- fread(cmd)
ninemers.random.table$five_prime_splice_score <- ninemers.random.res$V2[match(rownames(ninemers.random.table),rownames(ninemers.random.res))]

# THREE PRIME SPLICE SITES
twentythreemers.random <- c()
for(i in 1:1000){
        s <- sample(DNA_BASES, 23, replace = T)
        s[19] <- "A"
        s[20] <- "G"
        s <- paste(s, collapse = "")
        twentythreemers.random[i] <- s
}
# create table and write it out
twentythreemers.random.table <- data.frame(seq = twentythreemers.random)
twentythreemers.random.out <- paste0(paste0(SJ.scoring.outFolder,"random.twentythreemers.fasta"))
write.table(twentythreemers.random.table, twentythreemers.random.out, col.names=F,row.names=F,quote=F,sep="\t")
# run maxEnt 5' splice site analysis
setwd("/SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/")
cmd <- paste0("perl /SAN/biomed/biomed5/biomed5/GEUV_exon_annotaiton/exon_annotation/score3.pl ",twentythreemers.random.out)
twentythreemers.random.res <- fread(cmd)
twentythreemers.random.table$three_prime_splice_score <- twentythreemers.random.res$V2[match(rownames(twentythreemers.random.table),rownames(twentythreemers.random.res))]

#
random.sequence.slim <- data.frame(exon.type = "random sequence",
	gene.id = "random",
                five_prime_splice_score = ninemers.random.table$five_prime_splice_score,
                three_prime_splice_score = twentythreemers.random.table$three_prime_splice_score)

#### OUTPUT FOR SUMMARY TABLE
cryptic.exons.out <- cryptic.exons.slim
cryptic.exons.out$gene.name <- str_split_fixed(cryptic.exons.out$gene.id, "_",2)[,1] 
cryptic.exons.out$exonID <- str_split_fixed(cryptic.exons.out$gene.id, "_",2)[,2] 

cryptic.exons.out.table <- paste0(SJ.scoring.outFolder,"/",code,"_SJ_scores.tab")
write.table(cryptic.exons.out, cryptic.exons.out.table, sep="\t", quote=F, row.names=F)

### GRAPHING

all.exons <- rbind(cryptic.exons.slim,canonical.slim, random.sequence.slim)
# melt table for plotting
all.exons.gather <-  gather(all.exons, splice.site, score, five_prime_splice_score:three_prime_splice_score)
all.exons.gather$splice.site <- gsub("_"," ", gsub("_splice_score","",all.exons.gather$splice.site)) 
all.exons.gather$PSI.class <- cryptic.introns$PSI.class[match(all.exons.gather$gene.id, cryptic.introns$gene.id)]

# plot splice scores for cassette exons
cassette.exons.gather <- filter(all.exons.gather, (PSI.class == "CASSETTE.LIKE" | exon.type == "random sequence"))
cassette.exons.plot <- ggplot(cassette.exons.gather, aes(x = exon.type, y = score, 
	#group = gene.id, 
	fill = exon.type)) + 
			#geom_boxplot() +
			#geom_point() +
			geom_boxplot(notch=T) +
			#geom_line(colour="grey") +
			xlab("") +
			ylab("") +
			ylim(-35,15) +
			theme_bw() +
			#ylim(-20,20) +
			#scale_fill_manual(values = c("firebrick3","skyblue","skyblue")) +
			theme(panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
          axis.line = element_line(colour = "black"),
          legend.position="none") +
			ggtitle("Cassette-like cryptic exons") +
			facet_wrap(~splice.site, ncol = 2)


# test the differences between cryptic and canonical splice site scores in the 5' and 3' splice sites of cassette exons.
categories <- c("five prime", "three prime")
cassette.results <- list()
cassette.differences <- list()
for( i in  1:2 ){
	x <- filter(cassette.exons.gather, splice.site == categories[i] & exon.type == "cryptic")$score
	y <- filter(cassette.exons.gather, splice.site == categories[i] & exon.type == "canonical")$score
	cassette.results[[i]] <- t.test(y,x, paired =T)
	cassette.differences[[i]] <- y -x
}
cassette.results[[3]] <- t.test(cassette.differences[[1]], cassette.differences[[2]])





all.extensions <- filter(all.exons.gather, PSI.class != "CASSETTE.LIKE" & !is.na(score))
all.extensions <- filter(all.extensions, (splice.site == "five prime" & PSI.class == "THREEPRIME.BIAS") | (splice.site == "three prime" & PSI.class == "FIVEPRIME.BIAS" ) )
all.extensions$extension.type <- ifelse(all.extensions$PSI.class == "FIVEPRIME.BIAS",
	yes = "5\' extension",
	no = "3\' extension")

extensions.plot <- ggplot(all.extensions, aes(x = exon.type, y = score, fill = exon.type, label = gene.id)) + 
			#geom_boxplot() +
			#geom_line(colour="skyblue") +
			geom_boxplot(notch=T) +
			#geom_jitter(width=0.05) +
			#geom_text() +
			xlab("") +
			ylab("") +
			ylim(-35,15) +
			theme_bw() +
			theme(panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
          axis.line = element_line(colour = "black"),
          legend.position="none") +
			ggtitle("Extension-like cryptic exons") +
			facet_wrap(~extension.type, ncol = 2)



# Paired t-test of the two sets of extension scores
categories <- c("3\' extension","5\' extension")
extensions.results <- list()
extensions.differences <- list()
for(i in 1:2){
	x <- filter(all.extensions, extension.type == categories[i] & exon.type == "cryptic")$score
	y <- filter(all.extensions, extension.type == categories[i] & exon.type == "canonical")$score
	extensions.results[[i]] <- t.test(y, x, paired = T)
	extensions.differences[[i]] <- y - x
}			
# two sample t test of the differences between 5' and 3' extension scores
extensions.results[[3]] <- t.test(extensions.differences[[1]], extensions.differences[[2]])

sink(paste0(SJ.scoring.outFolder, code, "_SJ_scores_all_results.tab"))
cat("### Counts of the different classes of cryptic exon")
table(d$PSI.class)
cat("### Five prime splice sites in cassette-like cryptic exons")
cassette.results[[1]]
cat("### Three prime splice sites in cassette-like cryptic exons")
cassette.results[[2]]
cat("### Comparison of differences between the three prime and five prime splice site scores in cassette exons")
cassette.results[[3]]
cat("### 3\' extension exons - 5\' splice sites")
extensions.results[[1]]
cat("### 5\' extension exons - 3\' splice sites")
extensions.results[[2]]
cat("### Comparison of differences between the 3\' and 5\' extension splice site scores")
extensions.results[[3]]
sink()

#pdf(score.graph)
cassette.graph <- paste0(SJ.scoring.outFolder,code,"_cassette_scores.pdf")
extensions.graph <- paste0(SJ.scoring.outFolder,code,"_extensions_scores.pdf")

ggsave(file = cassette.graph, cassette.exons.plot)

ggsave(file = extensions.graph, extensions.plot)

SJ_scores.RData <- paste0(SJ.scoring.outFolder,code,"_SJ_scores.RData")
save.image(SJ_scores.RData)

# each gene has a five prime and three prime delta PSI based on the ratios of splice junctions. Does this correlate with maxEnt scores for those splice sites? Or with the differences between canonical and cryptic splice sites?

scores <- spread(all.exons.gather, exon.type,score) %>% filter(complete.cases(scores)) 
scores$PSI <- ifelse(scores$splice.site == "five prime", 
	yes = d$downstream_delta_psi[match(scores$gene.id, paste0(d$fix.gene.names, "_", d$exonID))],
	no = d$upstream_delta_psi[match(scores$gene.id, paste0(d$fix.gene.names, "_", d$exonID)) ])



