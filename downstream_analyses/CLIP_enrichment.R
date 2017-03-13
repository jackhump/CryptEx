#e/iCLIP enrichment

library(optparse)
suppressMessages(library(data.table))
library(ggplot2)
library(stringr)
suppressMessages(library(GenomicRanges))
library(plyr)
library(RColorBrewer)
suppressMessages(library(dplyr))
library(gplots)

options(echo=T) 
# create bed files of the cryptic exons and annotated exons.
PSI.threshold <- 0.05
min.canonical.control.SJs <- 5



files.exist <- function(files.list){
	for(my.file in files.list){
		if(!file.exists(my.file)){
			stop(paste(my.file,"doesn't exist!"))
		}
	}
}
files.are.empty <- function(files.list){
	if(length(files.list) == 0){
		stop(paste(files.list,"is empty!"))
	}
	
}


option_list <- list(
    make_option(c('--support.frame'), help=''),
    make_option(c('--code'), help=''),
    make_option(c('--condition.names'), help=''),
    make_option(c('--outFolder'), help=''),
    make_option(c('--species'), help=''),
    make_option(c('--CLIP'), help=''),
    make_option(c('--CLIPname')) 
)

########################## read arguments
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

if(length(opt) > 1){
	support.frame <- opt$support.frame
	code <- opt$code
	condition.names <- opt$condition.names
	outFolder <- opt$outFolder
	species <- opt$species
	CLIPlist <- opt$CLIP
}

#if CLIP list contains multiple arguments then split
if( grepl(",", CLIPlist)){
	CLIPlist <- strsplit(x = CLIPlist, split = ",")[[1]]
}

# code is CELL_EXPERIMENT_TARGET
# all previous steps used code = CELL_EXPERIMENT
# stupid hack
# only need the TARGET for the graph titles
title <- code
#code <- paste( str_split_fixed(code , '_', 3 )[,1:2], collapse = "_" )

print(title)
print(code)

splicing_analysis.res <- paste0(outFolder,"/splice_junction_analysis/",code,"_",condition.names,"_splicing_analysis.tab")
files.exist(splicing_analysis.res)

case.SJs <- paste0(outFolder,"/splice_junction_analysis/",code,"_SJs_case.tab")
files.exist(case.SJs)
# START HERE

CLIP.outFolder <- paste0(outFolder,"/CLIP_enrichment/")
if (! file.exists(CLIP.outFolder)) dir.create(CLIP.outFolder,recursive = T)

if(species == "mouse"){
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/mm10.fa"
	}

if(species == "human"){
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa "
	}
# if(species == "human"){
# 	RepeatMasker.bed <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/RepeatMasker/stranded/hg38_repeat_masker.bed"
# 	iCLIP.peaks <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_human/iCLIP/human_ES_cells/peaks_id71003_rnd100_flank15_fdr0.05_regionAsOne_20100222_LUjt3_4_hg38_ensembl59_G_iCLIP_TDP-43_Embrionic-Stem-Cells_bedGraph-cDNA-hits-in-genome.bed_lowFDR_stranded.bed"
# 	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa "
# 	iCLIP.folder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/iCLIP/human/"
	
# 	#CLIP.file.list <- c(TDP.iCLIP.A,TDP.iCLIP.B,TDP.iCLIP.C)
# 	#CLIP.name.list <- c("TDP43.iCLIP.A","TDP43.iCLIP.B","TDP43.iCLIP.C")
# }
############
#KEY are 
#CLIP.file.list  - a vector of strings of the paths to each CLIP file
#CLIP.name.list  - a vector of strings for the name of each CLIP file



adjacent_intron_detector <- function(SJ.summary,d){
	flank.col <- which(names(d) == "cryptic.flanked.length")
	exon.spacer <- 500
	adjacent_intron_query <- function(canonical.chr, canonical.start, canonical.end, cryptic.flanked.length,SJ.GRange){
		canonical.start <- as.numeric(canonical.start)
		canonical.end <- as.numeric(canonical.end)
		cryptic.flanked.length <- as.numeric(cryptic.flanked.length)
		junction <- list()
		junction <- SJ.GRange[
						#( width(SJ.GRange) > cryptic.flanked.length &
						seqnames(SJ.GRange) == canonical.chr  &
						# for capturing left adjacent introns 
						 ( ( end(SJ.GRange) < canonical.start ) & ( end(SJ.GRange) > (canonical.start - exon.spacer) ) |
						  ( start(SJ.GRange) > canonical.end  & start(SJ.GRange) < (canonical.end + exon.spacer) ) ) ]
		# insist that the intron must be longer than the length of the cryptic exon
		junction <- junction[width(junction) > cryptic.flanked.length]
		junction <- head(junction[order(score(junction),decreasing=T)],1)
	return(junction)
	}
	SJ.GRange <-  makeGRangesFromDataFrame(SJ.summary,keep.extra.columns=T)
	junctions.list <- apply(d, MAR=1,FUN=function(x) adjacent_intron_query(x[10],x[15],x[16], x[flank.col], SJ.GRange))	
	#output is a list of GRange objects - unuseable.
	junctions.list <- unlist(GRangesList(junctions.list))
	#convert into a dataframe, extracting the relevent information from the GRanges object.
	#names(GRanges) is a vector of rownames, confusingly.
	adjacent.df <- data.frame(row.names=names(junctions.list),
			chr=seqnames(junctions.list),
			adjacent.start=start(junctions.list),
			adjacent.end=end(junctions.list),
			adjacent.unique.count = score(junctions.list),
			strand = strand(junctions.list),
			intron.motif = mcols(junctions.list)[3])
	return(adjacent.df)
}

# random.start <- apply(adjacent.introns,MAR=1, FUN=function(x) 
# 	sample(seq(from = as.numeric(x[2]), 
# 				to = as.numeric(x[3]) - as.numeric(x[4]),
# 				by = 1) ))

adjacent.random.sampler <- function(chr, intron.start, intron.end, gene.name, exonID, strand, intron.length, cryptic.flanked.length){
	intron.start <- as.numeric(intron.start)
	intron.length <- as.numeric(intron.length)
	cryptic.flanked.length <- as.numeric(cryptic.flanked.length)

	v <- seq(intron.start, intron.start + (intron.length - cryptic.flanked.length), 1)
	random.start <- sample(v,1)
	end <- random.start + cryptic.flanked.length
	piece.length <- end - random.start
	return(c(chr,random.start, end, gene.name, exonID, strand,piece.length))
}

# need a function that generates random coordinates within an intron but deliberately avoids starting within the cryptic exon.
cryptic.intron.random.sampler <- function(chr, intron.start, intron.end, gene.name, exonID, strand, cryptic.start, cryptic.flanked.length, intron.length){
	intron.start <- as.numeric(intron.start)
	intron.length <- as.numeric(intron.length)
	cryptic.start <- as.numeric(cryptic.start)
	cryptic.flanked.length <- as.numeric(cryptic.flanked.length)
	
	v <- seq(intron.start, intron.start + (intron.length - cryptic.flanked.length), 1)
	random.start <- sample(v,1)
	#message(random.start)
	# keep generating random start coordinates until one is found that doesn't overlap the flanked cryptic exon coordinates.
	while( random.start == cryptic.start | # if there's an exact overlap
	( (random.start > cryptic.start ) & random.start < (cryptic.start + cryptic.flanked.length) ) | #  or if the start of the random sample falls within the flanked cryptic exon
		(	(random.start < cryptic.start) & (random.start + cryptic.flanked.length) > cryptic.start ) ){ # or if the end of the random sample falls within the flanked cryptic exon
		#message(random.start)
	 	random.start <- sample(v,1) 
	 }
	sample.end <- random.start + cryptic.flanked.length
	return(c(chr,random.start, sample.end, gene.name, exonID, strand))
}


d <- read.table(splicing_analysis.res,header=T,stringsAsFactors = F)

d$fix.strand <- gsub("-1","-",d$fix.strand,fixed=T)
d$fix.strand <- gsub("1","+",d$fix.strand,fixed=T)
d$canonical.start <- as.numeric(d$canonical.start)
d$canonical.end <- as.numeric(d$canonical.end)

# subset just the supported cryptic exons
d <- d[d$class == 'SJ.SUPPORTED.UP',]

flank <- 100
d$cryptic.start <- ifelse(d$upstream.case.cryptic.5prime > 0,
					yes = d$upstream.case.cryptic.5prime,
					no = d$exon.start)
d$cryptic.end <- ifelse(d$downstream.case.cryptic.3prime > 0,
					yes = d$downstream.case.cryptic.3prime,
					no = d$exon.end)
d$cryptic.start <- as.numeric(d$cryptic.start)
d$cryptic.end <- as.numeric(d$cryptic.end)
naughty.exons <- d$cryptic.end > d$cryptic.start
d <- d[naughty.exons,]
d$cryptic.length <- d$cryptic.end - d$cryptic.start
d$intron.length <- d$canonical.end - d$canonical.start
d$cryptic.flanked.length <- d$cryptic.length + (2 * flank)
cryptic.exons.flank.bed <- data.frame(chr = d$chr,
					 start = d$cryptic.start - flank , 
					 end = (d$cryptic.end + flank), 
					 gene.name = d$fix.gene.names, 
					 exonID = d$exonID, 
					 strand = d$fix.strand,
					 cryptic.flanked.length = d$cryptic.flanked.length, 
					 intron.length = d$intron.length )
# remove any entries where the introns are shorter than the flanked cryptic exon.
cryptic.exons.flank.bed <- cryptic.exons.flank.bed[cryptic.exons.flank.bed$cryptic.flanked.length < cryptic.exons.flank.bed$intron.length,]

# Create the beds of the cryptic and adjacent introns. exclude any introns that are not larger than 2 * length of flanked cryptic exon. 
# reduce interval by 50bp from both ends to avoid sampling the immediate 3' and 5' end.
cryptic.introns.bed <- data.frame(chr = d$chr, 
	intron.start = d$canonical.start + 50, 
	intron.end = d$canonical.end - 50, 
	gene.name = d$fix.gene.names,
	exonID = d$exonID, 
	strand = d$fix.strand, 
	cryptic.start = d$cryptic.start - flank,
	cryptic.flanked.length = d$cryptic.flanked.length, 
	intron.length = d$intron.length)
# remove any introns that are not at least longer than thrice the length of the flanked cryptic exon
cryptic.introns.bed <- cryptic.introns.bed[ (cryptic.introns.bed$cryptic.flanked.length * 3) < cryptic.introns.bed$intron.length, ]


# query the list of splice junctions to find the coordinates of introns adjacent to the cryptic containing intron.
#message(case.SJs)
SJ.summary <- as.data.frame(fread(case.SJs))
# find the adjacent introns from the merged list of splice junctions generated in the previous analysis.
adjacent.introns.bed <- adjacent_intron_detector(SJ.summary,d)

# avoid sampling the 3' or 5' splice sites!
adjacent.introns.bed$adjacent.start <- adjacent.introns.bed$adjacent.start + 50
adjacent.introns.bed$adjacent.end <- adjacent.introns.bed$adjacent.end - 50

adjacent.introns.bed$gene.name <- d$fix.gene.names[match(rownames(adjacent.introns.bed),rownames(d))]
adjacent.introns.bed$exonID <- d$exonID[match(rownames(adjacent.introns.bed),rownames(d))]
adjacent.introns.bed$intron.length <- adjacent.introns.bed$adjacent.end - adjacent.introns.bed$adjacent.start
adjacent.introns.bed$cryptic.flanked.length <- d$cryptic.flanked.length[match(rownames(adjacent.introns.bed),rownames(d))]
adjacent.introns.bed <- adjacent.introns.bed[,c(1,2,3,7,8,5,9,10)]
# make sure adjacent introns are at least as long as the flanked cryptic exon
adjacent.introns.bed <- adjacent.introns.bed[ adjacent.introns.bed$intron.length > adjacent.introns.bed$cryptic.flanked.length,]

# randomly sample a position within each adjacent intron that is the same length as the flanked cryptic exon

# produce a list of randomly selected pieces of adjacent introns.
iterations <- 100
time <- proc.time()
random.adjacent.introns.exons.list <- list()
for(i in 1:iterations){
	r <- apply(adjacent.introns.bed,MAR=1, FUN = function(x) adjacent.random.sampler(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8]))
	r <- as.data.frame(t(r))[,1:6]
	r$iteration <- i
	random.adjacent.introns.exons.list[[i]] <-  r
}


adjacent.introns.random.samples.bed <- do.call(what = rbind, args = random.adjacent.introns.exons.list)
names(adjacent.introns.random.samples.bed) <- c("chr","start","end","ExonID","gene.name","strand","iteration")

proc.time() - time

# For the cryptic introns, randomly pick a position within that intron. If it overlaps the cryptic exon coordinates then pick again.
iterations <- 100
random.cryptic.introns.piece.list <- list()
for(i in 1:iterations){
	r <- apply(cryptic.introns.bed, MAR=1, FUN = function(x) cryptic.intron.random.sampler(x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9]))
	r <- as.data.frame(t(r))[,1:6]
	r$iteration <- i
	random.cryptic.introns.piece.list[[i]] <- as.data.frame(r)
}
cryptic.introns.random.samples.bed <- do.call(what = rbind, args = random.cryptic.introns.piece.list)
names(cryptic.introns.random.samples.bed)[2:3] <- c("start","end")

# write out the bed files of the cryptic introns and adjacent introns for the Conservation analysis later
write.table(cryptic.introns.bed, paste0(CLIP.outFolder,"cryptic.introns.bed"), row.names=F,sep="\t",col.names=T,quote=F )
write.table(adjacent.introns.bed, paste0(CLIP.outFolder,"adjacent.introns.bed"), row.names=F,sep="\t",col.names=T,quote=F )

# Ready the bed files - the different RepeatMasker files and the iCLIP files.

# Intersect each set of introns with each bed file

bed.names <- c("cryptic_exons","cryptic_nulls","adjacent_nulls")
bed.files <- list(cryptic.exons.flank.bed, cryptic.introns.random.samples.bed, adjacent.introns.random.samples.bed)

# for testing
#bed.names <- c("cryptic_exons","cryptic_nulls")
#bed.files <- list(cryptic.exons.flank.bed, cryptic.introns.random.samples.bed)

# RepeatMasker overlap - how many of the cryptic segments overlap a repeat region compared to the null exons?
RM.files <- list()
#rep.type <- c("ALL","ALU+","ALU-","LINE+","LINE-","SINE+","SINE-","LC+","LC-","SR+","SR-")
# create names for the CLIP data
CLIP.name.list <- basename( CLIPlist)
CLIP.name.list <- gsub( ".bed", "", CLIP.name.list, fixed=T)

feature.type <- c(CLIP.name.list)
all.feature.bed.list <- c(CLIPlist)
bed.num <- 0

# Loop over for each bed file
for(i in 1:length(bed.names)){
	bed <- bed.files[[i]]
	bed <- bed[,1:7]
	if(names(bed)[7] != "iteration"){
		bed <- bed[,1:6]
		bed$iteration <- 1
	}
	bed.out <- paste0(CLIP.outFolder,bed.names[i],".exons.bed")
	write.table(bed, bed.out, row.names=F,sep="\t",col.names=F,quote=F )
	
	for( feature in 1:length(all.feature.bed.list) ){
		rep.bed <- all.feature.bed.list[feature]

		# clause for iCLIP intersection 
		if( grepl("peaks", feature.type[feature]) ) {
			rep.bed <- paste(rep.bed," -s")
		}
		RM.cmd <- paste0("bedtools intersect -a ",bed.out, " -b ", rep.bed, " -c")
		RM.bed <- fread(RM.cmd)
		names(RM.bed) <- c("chr","start","end","gene.name","exonID","strand","iteration","number.overlap.repeats")
		RM.bed$exon.type <- bed.names[i]
		RM.bed$feature.type <- feature.type[feature]
		bed.num <- bed.num + 1
		RM.files[[bed.num]] <- RM.bed
	}
}
RM.merge <- as.data.frame(do.call(what = rbind, args = RM.files))

RM.merge$overlap <- ifelse(RM.merge$number.overlap.repeats > 0, 1, 0)

# make heatmaps
CLIP.heatmap <- subset(RM.merge, exon.type == "cryptic_exons" & feature.type %in% CLIP.name.list)
CLIP.heatmap <- as.matrix(xtabs(formula = number.overlap.repeats ~ feature.type + gene.name, data = CLIP.heatmap))

CLIP.heatmap.pdf <- paste0(CLIP.outFolder,"/",title,"_CLIP_heatmap.pdf")

hmcol <- brewer.pal(6,"Blues")

pdf(CLIP.heatmap.pdf)
heatmap(t(CLIP.heatmap[,colSums(CLIP.heatmap) != 0]), margins = c(5,10),cex.lab = 0.5,cex.axis = 0.5,cexRow=0.5,cexCol=0.5, col=(hmcol))
dev.off()


summarised <- ddply(RM.merge,c("feature.type","exon.type"), summarise, 
	N = length(overlap), 
	num.overlapping = sum(overlap > 0), 
	proportion = num.overlapping / N,
	ci_minus = binom.test(x = num.overlapping, n = N)$conf.int[1],
	ci_plus = binom.test(x = num.overlapping, n = N)$conf.int[2] )
summarised$exon.type <- factor(x=summarised$exon.type, 
	levels = c("cryptic_exons","cryptic_nulls", "adjacent_nulls") )

summarised$feature.type <- gsub("."," ", summarised$feature.type,fixed=T)
#RM.type <- gsub(".","\n",RM.type,fixed=T)

for(exon in c("cryptic_exons","cryptic_nulls","adjacent_nulls")){
	s <- subset(summarised, exon.type == exon)
	names(s) <- paste(exon,names(s),sep=".")
	file.name <- paste0("summarised_",exon)
	assign(file.name, s)
}
all.exons <- cbind(summarised_cryptic_exons,summarised_cryptic_nulls,summarised_adjacent_nulls)
# for each row

all.exons$cryptic_vs_null.p.value <- 1E-10
all.exons$cryptic_vs_adjacent.p.value <- 1E-10

for(i in 1:nrow(all.exons)){
	all.exons$cryptic_vs_null.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,11]), n = c(all.exons[i,3],all.exons[i,10]), p = NULL, alternative = "greater")$p.value
	all.exons$cryptic_vs_adjacent.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,18]), n = c(all.exons[i,3],all.exons[i,17]), p = NULL,alternative = "greater")$p.value
}
all.exons$cryptic_vs_null.padj <- p.adjust(all.exons$cryptic_vs_null.p.value, method = "bonferroni")
all.exons$cryptic_vs_adjacent.padj <- p.adjust(all.exons$cryptic_vs_adjacent.p.value, method = "bonferroni")

p.values <-  all.exons[, c(1,(ncol(all.exons) - 1), ncol(all.exons)) ]
names(p.values) <- c("feature.type","cryptic_nulls","adjacent_nulls")
p.melt <- melt(p.values, id.vars="feature.type", variable.name = "exon.type", value.name= "p.value")
summarised$p.value <- p.melt$p.value[match(paste(summarised$feature.type, summarised$exon.type), paste(p.melt$feature.type, p.melt$exon.type))]

summarised$exon.type <- gsub("_"," ", summarised$exon.type)
summarised$exon.type <- factor(x=summarised$exon.type, 
	levels = c("cryptic exons","cryptic nulls", "adjacent nulls") )

summarised$sig.stars <- ifelse(summarised$p.value >= 0.05, NA, ifelse(summarised$p.value >= 0.001, "*", ifelse(summarised$p.value >= 1e-16, "**", "***")))

summarised.out <- paste0(CLIP.outFolder,"/",title,"_CLIP_enrichment_results.tab")
write.table(summarised, summarised.out, sep="\t", quote=F, row.names = F)


 
# plot graphs for iCLIP and RepeatMasker
 
iCLIP.pdf <- paste0(CLIP.outFolder,"/",title,"_CLIP_enrichment_plot.pdf")

iCLIP.summarised <- summarised

pdf(iCLIP.pdf)
p <- ggplot(iCLIP.summarised, aes(x = feature.type, y = proportion, fill = exon.type, group = exon.type, label = sig.stars) ) + 
	geom_bar(stat = "identity",position = "dodge") +
	geom_errorbar(aes(ymin=ci_minus, ymax=ci_plus), 
		width=.2,                    # Width of the error bars
        position=position_dodge(.9)) +
	scale_x_discrete("Feature type") +
	#geom_text(aes(group = exon.type)) + 
	scale_fill_discrete(name="exon type", labels=c("cryptic exon", "cryptic intron null","adjacent intron null")) +
	ylab("Proportion of overlapping exons") + 
	theme_bw() + 
	theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1) ) +
	scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
	ggtitle(paste0(gsub("_"," ",title)," (",species,")\niCLIP peaks\nN = ",nrow(cryptic.exons.flank.bed)) )
print(p)
dev.off()

