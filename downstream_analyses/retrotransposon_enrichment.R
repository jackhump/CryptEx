# retrotransposon enrichment
# strip out all default options and anything to do with iCLIP

library(optparse)
library(data.table)
library(ggplot2)
library(stringr)
library(GenomicRanges)
library(plyr)
library(RColorBrewer)
library(dplyr)
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
    make_option(c('--species'), help='') 
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
}

# stupid hack
title <- code
#code <- paste( str_split_fixed(code , '_', 3 )[,1:2], collapse = "_" )

print(title)
print(code)

splicing_analysis.res <- paste0(outFolder,"/splice_junction_analysis/",code,"_",condition.names,"_splicing_analysis.tab")
case.SJs <- paste0(outFolder,"/splice_junction_analysis/",code,"_SJs_case.tab")

# START HERE

RM.outFolder <- paste0(outFolder,"/retrotransposon_enrichment/")
if (! file.exists(RM.outFolder)) dir.create(RM.outFolder,recursive = T)

if(species == "mouse"){
	RepeatMasker.bed <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/RepeatMasker/stranded/mm10_repeat_masker.bed"
	#phyloP.bw <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/phyloP/mm10.60way.phyloP60way.bw"
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/mm10.fa"
	}

if(species == "human"){
	RepeatMasker.bed <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/RepeatMasker/stranded/hg38_repeat_masker.bed"
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa "
	}


# What about the entire lengths of the introns? Are they enriched for particular repeat elements and this is what gives such odd human results?
# To do this properly I need control introns. This could be for example the adjacent intron to the canonical intron. 

# I'm going to retool some of the functions I wrote before that discover splice junctions.

#this function looks for SJs that span from the upstream end of the canonical intron to the 5' end of the cryptic exon. Only the most abundant SJ will be reported.
#does this have to be canonical.start or canonical.start +1?
#applying the above function to each cryptic exon and cleaning up the result


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



# First create the cryptic exon bed with 100bp flanking either side. 


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
message(case.SJs)
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
write.table(cryptic.introns.bed, paste0(RM.outFolder,"cryptic.introns.bed"), row.names=F,sep="\t",col.names=T,quote=F )
write.table(adjacent.introns.bed, paste0(RM.outFolder,"adjacent.introns.bed"), row.names=F,sep="\t",col.names=T,quote=F )

# Ready the bed files - the different RepeatMasker files and the iCLIP files.

# Intersect each set of introns with each bed file

bed.names <- c("cryptic_exons","cryptic_nulls","adjacent_nulls")
bed.files <- list(cryptic.exons.flank.bed, cryptic.introns.random.samples.bed, adjacent.introns.random.samples.bed)

# for testing
#bed.names <- c("cryptic_exons","cryptic_nulls")
#bed.files <- list(cryptic.exons.flank.bed, cryptic.introns.random.samples.bed)

# RepeatMasker overlap - how many of the cryptic segments overlap a repeat region compared to the null exons?
RM.files <- list()
rep.type <- c("ALL","ALU+","ALU-","LINE+","LINE-","SINE+","SINE-","LC+","LC-","SR+","SR-","LTR+","LTR-","DNA+","DNA-" )
feature.type <- rep.type
rep.files <- gsub(pattern = "\\+|\\-|ALL",replacement = "",rep.type,fixed=F)
rep.bed.files <- paste0(str_split_fixed(RepeatMasker.bed, "[.]",2)[,1],"_", rep.files,".bed")


rep.bed.files <- sub("_.bed",".bed", rep.bed.files)

file.exists(rep.bed.files)
all.feature.bed.list <- rep.bed.files
bed.num <- 0

# Loop over for each bed file
for(i in 1:length(bed.names)){
	bed <- bed.files[[i]]
	bed <- bed[,1:7]
	if(names(bed)[7] != "iteration"){
		bed <- bed[,1:6]
		bed$iteration <- 1
	}
	bed.out <- paste0(RM.outFolder,bed.names[i],".exons.bed")
	write.table(bed, bed.out, row.names=F,sep="\t",col.names=F,quote=F )
	for( feature in 1:length(all.feature.bed.list) ){
		rep.bed <- all.feature.bed.list[feature]
		if( grepl("\\+$", feature.type[feature])) {
			rep.bed <- paste(rep.bed, " -s") # sense repeats
		}
		if( grepl("\\-$", feature.type[feature]) ){
			rep.bed <- paste(rep.bed," -S") # antisense repeats
		}
		RM.cmd <- paste0("bedtools intersect -a ",bed.out, " -b ", rep.bed, " -c")
		RM.bed <- fread(RM.cmd)
		names(RM.bed) <- c("chr","start","end","gene.name","exonID","strand","iteration","number.overlap.repeats")
		RM.bed$exon.type <- bed.names[i]
		RM.bed$feature.type <- feature.type[feature]
		bed.num <- bed.num + 1
		RM.files[[bed.num]] <- RM.bed
		# write out the All repeats file
		if( "ALL" %in% feature){
			all_repeats.out <- paste0(RM.outFolder,"/All_repeats_intersect.tab")
			allrep.cmd <- paste0("bedtools intersect -a ",bed.out, " -b ", rep.bed, " -wa -wb")
			all_repeats.intersect <- fread(allrep.cmd) 
			write.table(all_repeats.intersect, all_repeats.out, sep = '\t', quote = F, col.names = T, row.names = F)
		}
	}
}
RM.merge <- do.call(what = rbind, args = RM.files)

RM.merge$overlap <- ifelse(RM.merge$number.overlap.repeats > 0, 1, 0)

# make heatmaps
RM.heatmap <- subset(RM.merge, exon.type == "cryptic_exons" & feature.type %in% rep.type)
RM.heatmap <- as.matrix(xtabs(formula = number.overlap.repeats ~ feature.type + gene.name, data = RM.heatmap))

RM.heatmap.pdf <- paste0(RM.outFolder,"/",title,"_repeatmasker_heatmap.pdf")

hmcol <- brewer.pal(6,"Blues")

pdf(RM.heatmap.pdf)
heatmap(t(RM.heatmap[,colSums(RM.heatmap) != 0]), margins = c(5,10),cex.lab = 0.5,cex.axis = 0.5,cexRow=0.5,cexCol=0.5, col=(hmcol))
dev.off()

# plot the overlaps

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
summarised$feature.type <- gsub("\\+"," sense",summarised$feature.type,fixed=F)
summarised$feature.type <- gsub("\\-"," antisense",summarised$feature.type,fixed=F)
summarised$feature.type <- gsub("SR","Simple repeat",summarised$feature.type,fixed=F)
summarised$feature.type <- gsub("LC","Low complexity",summarised$feature.type,fixed=F)
summarised$feature.type <- gsub("ALU","Alu",summarised$feature.type,fixed=F)
summarised$feature.type <- gsub("ALL","All repeats",summarised$feature.type,fixed=F)
summarised$feature.type <- gsub("LTR","Long terminal repeats",summarised$feature.type,fixed=F)
summarised$feature.type <- gsub("DNA","DNA repeats",summarised$feature.type,fixed=F)

for(exon in c("cryptic_exons","cryptic_nulls","adjacent_nulls")){
	s <- subset(summarised, exon.type == exon)
	names(s) <- paste(exon,names(s),sep=".")
	file.name <- paste0("summarised_",exon)
	assign(file.name, s)
}
all.exons <- cbind(summarised_cryptic_exons,summarised_cryptic_nulls,summarised_adjacent_nulls)
# for each row

# statistical tests

all.exons$cryptic_vs_null.p.value <- 1E-10
all.exons$cryptic_vs_adjacent.p.value <- 1E-10

for(i in 1:nrow(all.exons)){
	all.exons$cryptic_vs_null.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,11]), n = c(all.exons[i,3],all.exons[i,10]), p = NULL, alternative = "greater")$p.value
	all.exons$cryptic_vs_adjacent.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,18]), n = c(all.exons[i,3],all.exons[i,17]), p = NULL,alternative = "greater")$p.value
}
all.exons$cryptic_vs_null.padj <- p.adjust(all.exons$cryptic_vs_null.p.value, method = "bonferroni")
all.exons$cryptic_vs_adjacent.padj <- p.adjust(all.exons$cryptic_vs_adjacent.p.value, method = "bonferroni")
# keep only the adjusted p values?
p.values <-  all.exons[, c(1,(ncol(all.exons) - 1), ncol(all.exons)) ]
names(p.values) <- c("feature.type","cryptic_nulls","adjacent_nulls")
p.melt <- melt(p.values, id.vars="feature.type", variable.name = "exon.type", value.name= "p.value")
summarised$p.value <- p.melt$p.value[match(paste(summarised$feature.type, summarised$exon.type), paste(p.melt$feature.type, p.melt$exon.type))]

summarised$exon.type <- gsub("_"," ", summarised$exon.type)
summarised$exon.type <- factor(x=summarised$exon.type, 
	levels = c("cryptic exons","cryptic nulls", "adjacent nulls") )

summarised$sig.stars <- ifelse(summarised$p.value >= 0.05, NA, ifelse(summarised$p.value >= 0.001, "*", ifelse(summarised$p.value >= 1e-16, "**", "***")))

summarised.out <- paste0(RM.outFolder,"/",title,"_retrotransposon_enrichment_results.tab")
write.table(summarised, summarised.out, sep="\t", quote=F, row.names = F)


 
# plot graphs 	
 
RM.pdf <- paste0(RM.outFolder,"/",title,"_retrotransposon_enrichment_plot.pdf")

summarised$feature.type <- factor(x = summarised$feature.type,
	levels = c("All repeats","LINE sense","LINE antisense","SINE sense","SINE antisense","Alu sense","Alu antisense","Low complexity sense","Low complexity antisense","Simple repeat sense","Simple repeat antisense", "Long terminal repeats sense", "Long terminal repeats antisense", "DNA repeats sense", "DNA repeats antisense" ))

pdf(RM.pdf)
p <- ggplot(summarised, aes(x = feature.type, y = proportion, fill = exon.type, group = exon.type) ) + 
	geom_bar(stat = "identity",position = "dodge") +
	geom_errorbar(aes(ymin=ci_minus, ymax=ci_plus), 
		width=.2,                    # Width of the error bars
                  position=position_dodge(.9)) +
	ylab("Proportion of overlapping exons") + 
	theme_bw() +
	theme(axis.text.x  = element_text(angle=45, vjust = 1, hjust = 1) ) +
	scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
	ggtitle(paste0(gsub("_"," ",title)," (",species,")\ncryptic exons overlapping repetitive elements\nN = ",nrow(cryptic.exons.flank.bed)))
print(p)
dev.off()

save.image( paste0(RM.outFolder, "/", title, "_retrotransposon_enrichment.Rdata"))






# for use with MEME


# output fasta sequence of the cryptic exons and adjacent exons (for use with MEME's discriminative mode)
MOTIF.outFolder <- paste(outFolder,"Motif_Finding", sep = "/")
if(species == "mouse"){	
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/mm10.fa"
}
if(species == "human"){
	genome.fa <- "/cluster/scratch3/vyp-scratch2/reference_datasets/human_reference_sequence/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fa "
}
if (! file.exists(MOTIF.outFolder)) dir.create(MOTIF.outFolder)

control.exons.bed <- subset(adjacent.introns.random.samples.bed, iteration <= 100)
bed.files <- list(cryptic.exons.flank.bed, control.exons.bed, adjacent.introns.bed)
bed.names <- c('cryptic.exons','control.exons', "adjacent.introns")
for(i in 1:length(bed.files)){
	MOTIF.bed <- bed.files[[i]]
	MOTIF.bed[,4] <- paste(MOTIF.bed[,4], MOTIF.bed[,5],sep = "_")
	if(!is.null(MOTIF.bed$iteration)){
		MOTIF.bed[,4] <- paste(MOTIF.bed[,4],MOTIF.bed$iteration,sep="_")
	}
	MOTIF.bed[,5] <- "."
	MOTIF.bed.name <- paste0(code,"_", bed.names[i])
	MOTIF.fasta.name <- paste0(MOTIF.outFolder,"/", code,"_",bed.names[i],".fasta")
	MOTIF.bed.out <- paste0(MOTIF.outFolder,"/",MOTIF.bed.name)
	write.table(MOTIF.bed, MOTIF.bed.out, quote=F,col.names=F,row.names=F,sep="\t")
	MOTIF.cmd <- paste0("bedtools getfasta -s -name -fi ",genome.fa," -bed ",MOTIF.bed.out," -fo ",MOTIF.fasta.name)
	system(MOTIF.cmd)
	assign(paste0(bed.names[i],".fasta"),MOTIF.fasta.name)
}
homer.outFolder <- paste(MOTIF.outFolder,"homer",sep="/")
if(! file.exists(homer.outFolder) ){ dir.create(homer.outFolder)}
homer.command <- paste0(homer.outFolder,"/homer_command.sh")

clusterFolder <- paste0(outFolder,"/cluster/")

homer_clusterArgs <- c("#$ -S /bin/bash",
"#$ -l h_vmem=4G",
"#$ -l tmem=4G",
"#$ -l h_rt=36:00:00",
"#$ -pe smp 1",
"#$ -R y",
paste0("#$ -o ", clusterFolder, "/out"),
paste0("#$ -e ", clusterFolder, "/error"),
paste0("#$ -N ", title, "_homer"),
paste0("#$ -wd ", outFolder),
"#$ -l h_rt=24:00:00" )

HOMER.cmd <- c(homer_clusterArgs, 
	paste("findMotifs.pl ",
		cryptic.exons.fasta,
		" fasta ", homer.outFolder, 
		" -fastaBg ", control.exons.fasta,
		" -rna  -nofacts -p 4 -S 20 -len 5,6,7,8,9 -noconvert -nogo" ) )

cat(HOMER.cmd, file = homer.command, sep = "\n", append = F)
# qsub it
system( paste0("qsub ", homer.command) )

save.image(paste0(RM.outFolder,code,"_Feature_Enrichment.RData") )

# dinucleotide and trinucleotide enrichment

#cryptic.exons.fasta <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/Motif_Finding/dataset_1_cryptic.exons.fasta"
#control.exons.fasta <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/Motif_Finding/dataset_1_control.exons.fasta"

# compare cryptic exons with adjacent intron sequences

cryptic.exons <- as.data.frame(fread(cryptic.exons.fasta,header=F,stringsAsFactors=F))
cryptic.exons <- data.frame(gene.id =  cryptic.exons[seq(1,nrow(cryptic.exons) - 1,2),],
							sequence = cryptic.exons[seq(2,nrow(cryptic.exons),2),])

control.exons <- as.data.frame(fread(adjacent.introns.fasta,header=F,stringsAsFactors=F))
control.exons <- data.frame(gene.id =  control.exons[seq(1,nrow(control.exons) - 1,2),],
							sequence = control.exons[seq(2,nrow(control.exons),2),])

exons.list <- list(cryptic.exons, control.exons)
exons.names <- c("cryptic.exons","control.exons")

for(exon.type in 1:length(exons.list)){

	subseq.list <- list()
	subseq <- c()
	exons <- exons.list[[exon.type]]

	exons$sequence <- toupper(exons$sequence)

	for( exon in 1:length(exons$sequence)){
		for(i in 1:( str_length(exons$sequence[exon]) -1 ) ){
			subseq[i] <- str_sub(exons$sequence[exon], i, (i+1) )	
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

p <- ggplot(results, aes(x = log2FoldChange, y = -log10(prop.test.padj), label = control.exons.subseq.total) ) + 
	geom_text(size = 10, ) + 
	xlim(c(-1,1.5)) +
	ylim(c(0,250)) +
	ggtitle(paste0(gsub("_"," ",code)," (",species,")\nDinucleotide enrichment between flanked cryptic exons\nand adjacent intronic sequence")) + 
	xlab("log2(cryptic vs null proportions)") + 
	ylab("-log10(adjusted p)") + 
	#ylim(c(0,(max.y+10))) + 
	#xlim(c(-1,1)) + 
	theme_bw()

write.table(results, paste0(MOTIF.outFolder,"/dinucleotide_enrichment.tab"),quote=F, row.names=F,sep="\t")
# plot <- ggplot(results, aes(x = log2FoldChange, y = 0, label = control.exons.x)) + 
# 		geom_text(size=7, position = position_jitter(height = 0.3)) +
# 		xlab("log2(cryptic vs null proportions)") +
# 		xlim(c(-1.5,1.5)) +
# 		ylim(c(-1,1)) +
# 		ylab("") + 
# 		ggtitle(paste0(gsub("_"," ",code)," (",species,")\nDinucleotide enrichment between flanked cryptic exons\nand adjacent intronic sequence")) + 
# 		theme_bw()

dinucleotide.pdf <- paste0(MOTIF.outFolder,"/dinucleotide_enrichment.pdf")
ggsave(dinucleotide.pdf)

save.image(paste0(RM.outFolder,code,"_Feature_Enrichment.RData") )







quit()







# # analyse all the eCLIP data from ENCODE

# if(species == "human") {

# bed.names <- c("cryptic_exons","cryptic_nulls","adjacent_nulls")
# bed.files <- list(cryptic.exons.flank.bed, cryptic.introns.random.samples.bed, adjacent.introns.random.samples.bed)

# for(cell in c("K562","HepG2")) {
# 	ENCODE.eCLIP.files <- paste0("/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/eCLIP/ALL_ENCODE/",cell,"/processed")
# 	ENCODE.eCLIP.list <- list.files(ENCODE.eCLIP.files,full.names=T)
# 	ENCODE.eCLIP.names <- str_split_fixed(list.files(ENCODE.eCLIP.files,full.names=F),"_",2)[,1]


# 	# for testing
# 	#bed.names <- c("cryptic_exons","cryptic_nulls")
# 	#bed.files <- list(cryptic.exons.flank.bed, cryptic.introns.random.samples.bed)

# 	# RepeatMasker overlap - how many of the cryptic segments overlap a repeat region compared to the null exons?
# 	eCLIP.files <- list()
# 	eCLIP.type <- c(ENCODE.eCLIP.names)
# 	eCLIP.outFolder <- RM.outFolder
# 	bed.num <- 0

# 	# Loop over for each bed file
# 	for(i in 1:length(bed.names)){
# 		bed <- bed.files[[i]]
# 		bed <- bed[,1:7]
# 		if(names(bed)[7] != "iteration"){
# 			bed <- bed[,1:6]
# 			bed$iteration <- 1
# 		}
# 		bed.out <- paste0(eCLIP.outFolder,bed.names[i],".exons.bed")
# 		write.table(bed, bed.out, row.names=F,sep="\t",col.names=F,quote=F )
# 		for( type in 1:length(eCLIP.type) ){
# 			rep.bed <- ENCODE.eCLIP.list[type]
# 			eCLIP.cmd <- paste0("bedtools intersect -a ",bed.out, " -b ", rep.bed, " -c -s")
# 			eCLIP.bed <- fread(eCLIP.cmd)
# 			names(eCLIP.bed) <- c("chr","start","end","gene.name","exonID","strand","iteration","number.overlap.repeats")
# 			eCLIP.bed$exon.type <- bed.names[i]
# 			eCLIP.bed$feature.type <- eCLIP.type[type]
# 			bed.num <- bed.num + 1
# 			eCLIP.files[[bed.num]] <- eCLIP.bed
# 		}
# 	}


# 	eCLIP.merge <- do.call(what = rbind, args = eCLIP.files)

# 	#eCLIP.sum <- lapply(seq_along(eCLIP.files), function(i) cbind(eCLIP.files[[i]][1],sum(eCLIP.files[[i]]$number.overlap.repeats > 0) ))
# 	eCLIP.merge$overlap <- ifelse(eCLIP.merge$number.overlap.repeats > 0, 1, 0)

# 	eCLIP.heatmap <- subset(eCLIP.merge, exon.type == "cryptic_exons")
# 	eCLIP.heatmap <- as.matrix(xtabs(formula = number.overlap.repeats ~ feature.type + gene.name, data = eCLIP.heatmap))

# 	eCLIP.heatmap.pdf <- paste0(eCLIP.outFolder,"/",code,"_ENCODE_",cell,"_heatmap_total.pdf")

# 	hmcol <- brewer.pal(6,"Greens")

# 	pdf(eCLIP.heatmap.pdf)
# 	heatmap(t(eCLIP.heatmap[rowSums(eCLIP.heatmap) != 0 ,colSums(eCLIP.heatmap) != 0]), margins = c(5,0),cex.lab = 0.5,cex.axis = 0.5,cexRow=0.5,cexCol=0.5, col=(hmcol))
# 	dev.off()



# 	eCLIP.summarised <- ddply(eCLIP.merge,c("feature.type","exon.type"), summarise, 
# 		N = length(overlap), 
# 		num.overlapping = sum(overlap > 0), 
# 		proportion = num.overlapping / N,
# 		ci_minus = binom.test(x = num.overlapping, n = N)$conf.int[1],
# 		ci_plus = binom.test(x = num.overlapping, n = N)$conf.int[2] )
# 	eCLIP.summarised$exon.type <- factor(x=eCLIP.summarised$exon.type, 
# 		levels = c("cryptic_exons","cryptic_nulls", "adjacent_nulls") )
# 	eCLIP.summarised$feature.type <- gsub("."," ", eCLIP.summarised$feature.type,fixed=T)
# 	eCLIP.type <- gsub(".","\n",eCLIP.type,fixed=T)

# 	eCLIP.summarised$feature.type <- gsub("-", " ", eCLIP.summarised$feature.type)
# 	eCLIP.summarised$feature.type <- gsub("TARDBP","TDP-43", eCLIP.summarised$feature.type)

# 	# do proportion tests
# 	for(exon in c("cryptic_exons","cryptic_nulls","adjacent_nulls")){
# 		s <- subset(eCLIP.summarised, exon.type == exon)
# 		names(s) <- paste(exon,names(s),sep=".")
# 		file.name <- paste0("summarised_",exon)
# 		assign(file.name, s)
# 	}
# 	all.exons <- cbind(summarised_cryptic_exons,summarised_cryptic_nulls,summarised_adjacent_nulls)
# 	# for each row

# 	all.exons$cryptic_vs_null.p.value <- 1E-10
# 	all.exons$cryptic_vs_adjacent.p.value <- 1E-10

# 	for(i in 1:nrow(all.exons)){
# 		all.exons$cryptic_vs_null.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,11]), n = c(all.exons[i,3],all.exons[i,10]), p = NULL, alternative = "greater")$p.value
# 		all.exons$cryptic_vs_adjacent.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,18]), n = c(all.exons[i,3],all.exons[i,10]), p = NULL,alternative = "greater")$p.value
# 	}
# 	all.exons$cryptic_vs_null.padj <- p.adjust(all.exons$cryptic_vs_null.p.value, method = "bonferroni", n = length(all.exons$cryptic_vs_null.p.value) * 2) # correct for 
# 	all.exons$cryptic_vs_adjacent.padj <- p.adjust(all.exons$cryptic_vs_adjacent.p.value, method = "bonferroni", n = length(all.exons$cryptic_vs_adjacent.p.value) * 2)

# 	p.values <-  all.exons[, c(1,(ncol(all.exons) - 1), ncol(all.exons)) ]
# 	names(p.values) <- c("feature.type","cryptic_nulls","adjacent_nulls")
# 	p.melt <- melt(p.values, id.vars="feature.type", variable.name = "exon.type", value.name= "p.value")
# 	eCLIP.summarised$p.value <- p.melt$p.value[match(paste(eCLIP.summarised$feature.type, eCLIP.summarised$exon.type), paste(p.melt$feature.type, p.melt$exon.type))]

# 	table.name <- paste0(cell,"_all.exons")
# 	assign(table.name, all.exons)
# 	# plot graphs

# 	eCLIP.pdf <- paste0(eCLIP.outFolder,"/",code,"_ENCODE_",cell,"_eCLIP_enrichment_plot.pdf")

# 	# all RBPs
# 	p <- ggplot(eCLIP.summarised, aes(x = feature.type, y = proportion, fill = exon.type, group = exon.type) ) + 
# 		geom_bar(stat = "identity",position = "dodge") +
# 		#geom_point(position=position_dodge(width=0.5)) + 
# 		geom_errorbar(aes(ymin=ci_minus, ymax=ci_plus), 
# 			width=.2,                    # Width of the error bars
# 	                  position=position_dodge(.9)) +
# 		scale_x_discrete("Feature type") +
# 		scale_fill_discrete(name="exon type", labels=c("cryptic exon", "cryptic intron null","adjacent intron null")) +
# 		ylab("Proportion of overlapping exons") + 
# 		ylim(0,1) +
# 		theme(axis.text.x  = element_text(angle=45, vjust = 0.5) ) +
# 		ggtitle(paste0(gsub("_"," ",code)," (",species,")\n",cell," eCLIP peaks\nN = ",nrow(cryptic.exons.flank.bed),"\nPSI.threshold = ",PSI.threshold,"\ncanonical SJs min = ",min.canonical.control.SJs) )

# 	ggsave(filename = eCLIP.pdf, plot = p, width = 40, height = 10, units = "in")

# 	ranked <- select(all.exons, cryptic_exons.feature.type, cryptic_vs_null.p.value, cryptic_vs_adjacent.p.value) %>% # select out wanted columns
# 		 separate(cryptic_exons.feature.type, into = c("RBP","replicate"), sep = " ") %>% # split dataset name into gene and replicate number
# 		 gather(variable,value, -(RBP:replicate)) %>%  # convert into long form
# 		 unite(temp, variable, replicate) %>% # concatanate p value type with replicate number
# 		 spread(temp,value) %>% # convert to wide form
# 		 mutate(max.p = pmax(cryptic_vs_adjacent.p.value_1, cryptic_vs_adjacent.p.value_2, cryptic_vs_null.p.value_1, cryptic_vs_null.p.value_2) ) %>% # add new column which contains the maximum p value of the four tests
# 		 arrange(max.p) %>% # arrange by ascending p value
# 		 #mutate(max.p.adjust = max.p) # version with no multiple testing correction - arguing w/ Vincent
# 		 mutate(max.p.adjust = p.adjust(max.p, method = "bonferroni", n = (length(max.p) ) ) ) # correct for multiple testing

# 	ranked.hits <- filter(ranked, max.p.adjust < 0.05)
# 	write.table(ranked.hits, paste0(eCLIP.outFolder,"/",code, "_ENCODE_", cell, "_eCLIP_enrichment_hits.tab"))
# 	eCLIP.summarised.hits <- separate(eCLIP.summarised, feature.type, into = c("RBP","replicate"), sep = " ", remove = F) %>%
# 			filter(RBP %in% ranked.hits$RBP)

# 	p.hits <- ggplot(eCLIP.summarised.hits, aes(x = feature.type, y = proportion, fill = exon.type, group = exon.type) ) + 
# 		geom_bar(stat = "identity",position = "dodge") +
# 		#geom_point(position=position_dodge(width=0.5)) + 
# 		geom_errorbar(aes(ymin=ci_minus, ymax=ci_plus), 
# 			width=.2,                    # Width of the error bars
# 	                  position=position_dodge(.9)) +
# 		scale_x_discrete("Feature type") +
# 		scale_fill_discrete(name="exon type", labels=c("cryptic exon", "cryptic intron null","adjacent intron null")) +
# 		ylab("Proportion of overlapping exons") + 
# 		scale_y_continuous(labels = scales::percent, limits = c(0,0.5)) +
# 		theme(axis.text.x  = element_text(angle=45, vjust = 0.5) ) +
# 		ggtitle(paste0(gsub("_"," ",code)," (",species,")\n",cell," eCLIP peaks\nN = ",nrow(cryptic.exons.flank.bed),"\nPSI.threshold = ",PSI.threshold,"\ncanonical SJs min = ",min.canonical.control.SJs) ) +
# 		theme(legend.position="none") 

# 	eCLIP.hits.pdf <- paste0(eCLIP.outFolder,"/",code,"_ENCODE_",cell,"_eCLIP_enrichment_hits.pdf")

# 	ggsave(filename = eCLIP.hits.pdf, plot = p.hits, width = 7, height = 7, units = "in")

# 	summarised.hits.name <- paste0(cell,".eCLIP.summarised.hits")
# 	assign(summarised.hits.name, eCLIP.summarised.hits)

# 	eCLIP.merge$feature.type <- paste(eCLIP.merge$feature.type,cell,sep="_")

# 	merge.name <- paste0(cell,".eCLIP.merge")
# 	assign(merge.name,eCLIP.merge)

# }



# # TARDBP to TDP-43
# K562.eCLIP.merge$feature.type <- gsub("TARDBP","TDP-43",K562.eCLIP.merge$feature.type)

# K562.eCLIP.merge.hits <- filter(K562.eCLIP.merge, grepl(paste(K562.eCLIP.summarised.hits$RBP, collapse="|"),feature.type) & exon.type == "cryptic_exons")
# HepG2.eCLIP.merge.hits <- filter(HepG2.eCLIP.merge, grepl(paste(HepG2.eCLIP.summarised.hits$RBP, collapse="|"),feature.type) & exon.type == "cryptic_exons")

# both.cells <- as.matrix( xtabs(number.overlap.repeats ~ gene.name + feature.type, data = rbind(K562.eCLIP.merge.hits,HepG2.eCLIP.merge.hits) ) )

# eCLIP.heatmap.pdf <- paste0(eCLIP.outFolder,"/",code,"_ENCODE_both_cell_types_heatmap.pdf")

# hmcol <- brewer.pal(6,"Greens")



# pdf(eCLIP.heatmap.pdf,width = 7,        # 5 x 300 pixels
# height = 7,           # 300 pixels per inch
# pointsize = 3)
#     # turn off column clustering
# my_palette <- c("white","skyblue3")
# col_breaks <- c(seq(0,1,length = 3) ) # for red


# heatmap.2((both.cells[rowSums(both.cells) != 0 ,colSums(both.cells) != 0]),
# # cellnote = (both.cells.replicates[rowSums(both.cells.replicates) != 0 ,colSums(both.cells.replicates) != 0]),  # same data set for cell labels
#   main = "Both cell types ENCODE eCLIP", # heat map title
#   notecol="black",      # change font color of cell labels to black
#   density.info="none",  # turns off density plot inside color legend
#   trace="none",         # turns off trace lines inside the heat map
#   margins =c(12,9),     # widens margins around plot
#   col=my_palette,       # use on color palette defined earlier
#   breaks=col_breaks,    # enable color transition at specified limits
#   dendrogram="column",
#   key=F,
#   colsep=seq(1,ncol(both.cells)),
#   rowsep=seq(1,nrow(both.cells)),
#   sepcolor = "lightgray",
#   sepwidth=c(0.005,0.005)     # only draw a row dendrogram
#   #Colv="NA"
#   )            # turn off column clustering


# dev.off()



# eCLIP.RData = paste0(eCLIP.outFolder,"/",code,"_Feature_Enrichment.RData")
# save.image(eCLIP.RData)


# # Replication - try to replicate eCLIP hits using iCLIP data

# replication.outFolder <- paste0(RM.outFolder,"replication/")
# if(!file.exists(replication.outFolder)){
# 	dir.create(replication.outFolder)
# }
# bed.num <- 0
# replication.overlap.files <- list()
# # Loop over for each bed file
# for(i in 1:length(bed.names)){
# 	bed <- bed.files[[i]]
# 	bed <- bed[,1:7]
# 	if(names(bed)[7] != "iteration"){
# 		bed <- bed[,1:6]
# 		bed$iteration <- 1
# 	}
# 	bed.out <- paste0(replication.outFolder,bed.names[i],".exons.bed")
# 	write.table(bed, bed.out, row.names=F,sep="\t",col.names=F,quote=F )
# 	for( type in 1:length(replication.file.list) ){
# 		rep.bed <- replication.files[type]
# 		replication.cmd <- paste0("bedtools intersect -a ",bed.out, " -b ", rep.bed, " -c -s")
# 		replication.bed <- fread(replication.cmd)
# 		names(replication.bed) <- c("chr","start","end","gene.name","exonID","strand","iteration","number.overlap.repeats")
# 		replication.bed$exon.type <- bed.names[i]
# 		replication.bed$feature.type <- replication.file.list[type]
# 		bed.num <- bed.num + 1
# 		replication.overlap.files[[bed.num]] <- replication.bed
# 	}
# }
# replication.merge <- do.call(what = rbind, args = replication.overlap.files)

# replication.merge$overlap <- ifelse(replication.merge$number.overlap.repeats > 0, 1, 0)

# replication.summarised <- ddply(replication.merge,c("feature.type","exon.type"), summarise, 
# 	N = length(overlap), 
# 	num.overlapping = sum(overlap > 0), 
# 	proportion = num.overlapping / N,
# 	ci_minus = binom.test(x = num.overlapping, n = N)$conf.int[1],
# 	ci_plus = binom.test(x = num.overlapping, n = N)$conf.int[2] )
# replication.summarised$exon.type <- factor(x=replication.summarised$exon.type, levels = c("cryptic_exons","cryptic_nulls", "adjacent_nulls") )
# replication.summarised$feature.type <- gsub("."," ", replication.summarised$feature.type,fixed=T)
# #replication.type <- gsub(".","\n",replication.type,fixed=T)

# replication.summarised$feature.type <- gsub("-", " ", replication.summarised$feature.type)
# replication.summarised$feature.type <- gsub("TARDBP","TDP-43", replication.summarised$feature.type)

# # do proportion tests
# for(exon in c("cryptic_exons","cryptic_nulls","adjacent_nulls")){
# 	s <- subset(replication.summarised, exon.type == exon)
# 	names(s) <- paste(exon,names(s),sep=".")
# 	file.name <- paste0("summarised_",exon)
# 	assign(file.name, s)
# }
# all.exons <- cbind(summarised_cryptic_exons,summarised_cryptic_nulls,summarised_adjacent_nulls)
# # for each row

# all.exons$cryptic_vs_null.p.value <- 1E-10
# all.exons$cryptic_vs_adjacent.p.value <- 1E-10

# for(i in 1:nrow(all.exons)){
# 	all.exons$cryptic_vs_null.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,11]), n = c(all.exons[i,3],all.exons[i,10]), p = NULL, alternative = "greater")$p.value
# 	all.exons$cryptic_vs_adjacent.p.value[i] <- prop.test(x = c(all.exons[i,4], all.exons[i,18]), n = c(all.exons[i,3],all.exons[i,10]), p = NULL,alternative = "greater")$p.value
# }
# all.exons$cryptic_vs_null.padj <- p.adjust(all.exons$cryptic_vs_null.p.value, method = "bonferroni", n = length(all.exons$cryptic_vs_null.p.value) * 2) # correct for 
# all.exons$cryptic_vs_adjacent.padj <- p.adjust(all.exons$cryptic_vs_adjacent.p.value, method = "bonferroni", n = length(all.exons$cryptic_vs_adjacent.p.value) * 2)

# p.values <-  all.exons[, c(1,(ncol(all.exons) - 1), ncol(all.exons)) ]
# names(p.values) <- c("feature.type","cryptic_nulls","adjacent_nulls")
# p.melt <- melt(p.values, id.vars="feature.type", variable.name = "exon.type", value.name= "p.value")
# replication.summarised$p.value <- p.melt$p.value[match(paste(replication.summarised$feature.type, replication.summarised$exon.type), paste(p.melt$feature.type, p.melt$exon.type))]

# table.name <- paste0(cell,"_all.exons")
# assign(table.name, all.exons)
# # plot graphs

# replication.pdf <- paste0(replication.outFolder,"/",code,"_replication_enrichment_plot.pdf")

# # all RBPs
# p <- ggplot(replication.summarised, aes(x = feature.type, y = proportion, fill = exon.type, group = exon.type) ) + 
# 	geom_bar(stat = "identity",position = "dodge") +
# 	#geom_point(position=position_dodge(width=0.5)) + 
# 	geom_errorbar(aes(ymin=ci_minus, ymax=ci_plus), 
# 		width=.2,                    # Width of the error bars
#                   position=position_dodge(.9)) +
# 	scale_x_discrete("Feature type") +
# 	scale_fill_discrete(name="exon type", labels=c("cryptic exon", "cryptic intron null","adjacent intron null")) +
# 	ylab("Proportion of overlapping exons") + 
# 	theme(axis.text.x  = element_text(angle=45, vjust = 0.5) ) +
# 	scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
# 	ggtitle(paste0(gsub("_"," ",code)," (",species,")\n",cell," iCLIP peaks\nN = ",nrow(cryptic.exons.flank.bed),"\nPSI.threshold = ",PSI.threshold,"\ncanonical SJs min = ",min.canonical.control.SJs) )

# ggsave(filename = replication.pdf, plot = p, width = 7, height = 7, units = "in")

# both.cells.replicates <- as.matrix( xtabs(number.overlap.repeats ~ gene.name + feature.type, data = rbind(K562.eCLIP.merge.hits,HepG2.eCLIP.merge.hits, replication.merge[replication.merge$exon.type == "cryptic_exons"]) ) )

# replication.heatmap.pdf <- paste0(replication.outFolder,"/",code,"_all_plus_replicates_heatmap.pdf")


# pdf(replication.heatmap.pdf,width = 7,        # 5 x 300 pixels
# height = 7,           # 300 pixels per inch
# pointsize = 3)
# my_palette <- c("white","skyblue3")
# col_breaks <- c(seq(0,1,length = 3) ) # for red

 
# #both.cells.replicates[both.cells.replicates == "0"] <- NA

# #heatmap((both.cells[rowSums(both.cells) != 0 ,colSums(both.cells) != 0]), margins = c(5,0),cex.lab = 0.5,cex.axis = 0.5,cexRow=0.5,cexCol=0.5, col=(hmcol))
# heatmap.2((both.cells.replicates[rowSums(both.cells.replicates) != 0 ,colSums(both.cells.replicates) != 0]),
# # cellnote = (both.cells.replicates[rowSums(both.cells.replicates) != 0 ,colSums(both.cells.replicates) != 0]),  # same data set for cell labels
#   main = "Both cell types ENCODE eCLIP + iCLIP replicates", # heat map title
#   notecol="black",      # change font color of cell labels to black
#   density.info="none",  # turns off density plot inside color legend
#   trace="none",         # turns off trace lines inside the heat map
#   margins =c(12,9),     # widens margins around plot
#   col=my_palette,       # use on color palette defined earlier
#   breaks=col_breaks,    # enable color transition at specified limits
#   dendrogram="both",
#   key=F,
#   colsep=seq(1,ncol(both.cells.replicates)),
#   rowsep=seq(1,nrow(both.cells.replicates)),
#   sepcolor = "lightgray",
#   sepwidth=c(0.005,0.005)     # only draw a row dendrogram
#   #Colv="NA"
#   )            # turn off column clustering

# dev.off()







# quit()
# quit

