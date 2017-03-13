# Conservation 
library(data.table)
library(ggplot2)
library(dplyr)
library(optparse)
library(stringr)

options(echo =T)
# create bed files of the cryptic exons and annotated exons.
# use these bed files to find the mean PhyloP conservation score for all the nucledotides within each exon.

# ##################
# ## CLEVELAND TDP43
# ##################
# species <- "mouse"
# code <- "Cleveland_TDP43"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43"
# splicing_analysis.res <- paste0(outFolder,"/splice_junction_analysis/Cleveland_TDP43_CTL_TDP_splicing_analysis.tab")
# dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/strict_500/dexseq/CTL_TDP/Cleveland_TDP43_CTL_TDP_SignificantExons.csv"

# ###################
# ## ENCODE DATASET1 
# ###################
# code <- "dataset_1"
# species <- "human"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1"
# splicing_analysis.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_control_TDP_splicing_analysis.tab"
# case.SJs <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_SJs_case.tab"
# dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/strict_500/dexseq/control_TDP/dataset_1_control_TDP_SignificantExons.csv" 


# ###################
# ## ENCODE DATASET 2 
# ###################
# code <- "dataset_2"
# species <- "human"
# outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2"
# splicing_analysis.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/splice_junction_analysis/dataset_2_control_TDP_splicing_analysis.tab"
# case.SJs <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting//TDP_ENCODE_human/dataset_2/splice_junction_analysis/dataset_2_SJs_case.tab"
# dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/strict_500/dexseq/control_TDP/dataset_2_control_TDP_SignificantExons.csv" 

# Both Human TDP43 ENCODE for paper
code <- "both_ENCODE"
species <- "human"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/human/"
dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/human/both_datasets_unison_splicing_analysis.tab"


# Both Mouse TDP43 for paper
code <- "both_mouse"
species <- "mouse"
outFolder <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/mouse/"
dexseq.res <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/mouse/both_datasets_unison_splicing_analysis.tab"

option_list <- list(
    make_option(c('--code'), help=''),
    make_option(c('--outFolder'), help=''),
    make_option(c('--species'), help=''),
    make_option(c('--condition.names'))
)

########################## read arguments
option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

if(length(opt) > 1){
	code <- opt$code
	outFolder <- opt$outFolder
	species <- opt$species
	condition.names <- opt$condition.names
}

title <- code
code <- paste( str_split_fixed(code , '_', 3 )[,1:2], collapse = "_" )

dexseq.res <- paste0(outFolder,"/splice_junction_analysis/",code,"_",condition.names,"_splicing_analysis.tab")



conservation.outFolder <- paste0(outFolder,"/conservation/")
if (! file.exists(conservation.outFolder)) dir.create(conservation.outFolder)

if(species == "mouse"){
	phyloP.bw <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/phyloP/mm10.60way.phyloP60way.bw"
}
if(species == "human"){
	phyloP.bw <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/reference/phyloP/hg38.phyloP100way.bw"
}
# find cryptic and adjacent introns
allfiles <- list.files(outFolder,recursive = T, full.names = T)
# for ENCODE this was done twice so just take the top one as they are identical
cryptic.introns <- allfiles[ grepl( "cryptic.introns.bed", allfiles)][1]
adjacent.introns <- allfiles[ grepl( "adjacent.introns.bed", allfiles)][1]
#cryptic.introns <- paste0(outFolder,"/Feature_Enrichment/cryptic.introns.bed")
#adjacent.introns <- paste0(outFolder,"/Feature_Enrichment/adjacent.introns.bed")

cryptic.introns.bed <- as.data.frame(fread(cryptic.introns))
adjacent.introns.bed <- as.data.frame(fread(adjacent.introns))

# work out adjacent canonical exons and grab their conservation scores!
# read in the adjacent and cryptic intron bed files from the Feature Enrichment analysis. 

#They have been chopped at each end by 50bp so add back in!
cryptic.introns.bed$intron.start <- cryptic.introns.bed$intron.start - 50
cryptic.introns.bed$intron.end <- cryptic.introns.bed$intron.end + 50
adjacent.introns.bed$adjacent.start <- adjacent.introns.bed$adjacent.start - 50
adjacent.introns.bed$adjacent.end <- adjacent.introns.bed$adjacent.end + 50

# make sure the adjacent introns match the cryptic introns!
adjacent.exons <- cryptic.introns.bed
adjacent.exons$adjacent.intron.start <- adjacent.introns.bed$adjacent.start[match(paste(adjacent.exons$gene.name, adjacent.exons$exonID),
											paste(adjacent.introns.bed$gene.name, adjacent.introns.bed$exonID)) ]

adjacent.exons$adjacent.intron.end <- adjacent.introns.bed$adjacent.end[match(paste(adjacent.exons$gene.name, adjacent.exons$exonID),
											paste(adjacent.introns.bed$gene.name, adjacent.introns.bed$exonID))]

adjacent.exons$exon.start <- with(adjacent.exons, ifelse(intron.start > adjacent.intron.start, yes = adjacent.intron.end + 1, no = intron.end +1))
adjacent.exons$exon.end <-   with(adjacent.exons, ifelse(intron.start > adjacent.intron.start, yes = intron.start - 1, no = adjacent.intron.start - 2))

adjacent.exons.bed <- filter(adjacent.exons, complete.cases(adjacent.exons)) %>% select(chr, exon.start, exon.end, gene.name:strand)

#stupid hack
bedfiles.folder <- dirname(cryptic.introns)

write.table(adjacent.exons.bed, paste0(bedfiles.folder,"/adjacent_exons.exons.bed"),sep ="\t", row.names= F, col.names =F, quote =F)


# # pick x non-significant exons at random
# random_sample <- 5000


# a <- as.data.frame(fread(dexseq.res))
# a <- a[!grepl("i",a$exonID),]
# a$strand <- gsub("-1","-",a$strand,fixed=T)
# a$strand <- gsub("1","+",a$strand,fixed=T)

# annotated.nonsig <- subset(a, FDR > 0.5 )
# random_non_sig <- sample(x = rownames(annotated.nonsig), size = random_sample, replace = F)
# annotated.nonsig <- annotated.nonsig[random_non_sig,]
# annotated.nonsig.bed <- data.frame(chr = annotated.nonsig$chr, start = annotated.nonsig$exon.start, end = annotated.nonsig$exon.end,gene.name = annotated.nonsig$external_gene_id, exonID = annotated.nonsig$exonID, strand = annotated.nonsig$strand)
# annotated.nonsig.bed.out <- paste0(conservation.outFolder,"annotated_nonsig_random_exons.bed")
# write.table(annotated.nonsig.bed, annotated.nonsig.bed.out,sep="\t","quote"=F,col.names=F,row.names=F)




bed.names <- c("cryptic_exons","adjacent_exons","cryptic_nulls")
bed.files.list <- paste0(bedfiles.folder,"/",bed.names,".exons.bed")

file.exists(bed.files.list)
#bed.names <- c("annotated exons",bed.names)
#bed.files <- list(cryptic.exons.bed, control.exons.bed,annotated.sig.bed,annotated.nonsig.bed)
#bed.files.list <- c(annotated.nonsig.bed.out, bed.files.list)

# PhyloP conservation!
exon.conservation <- list()
for(i in 1:length(bed.names)){
	bed.out <- bed.files.list[i]
	#bed <- read.table(bed.out,header=T)
	cmd <- paste0("cat ",bed.out," | while read chr start end name exonID strand; do bigWigSummary ", phyloP.bw, " $chr $start $end 1; done  2>&1 | awk \' $0 ~ \"data\"{print \"NA\"}$0 !~ \"data\" {print $0}\' | sed \'/^$/d\'")
	conservation <- fread(cmd)
	names(conservation)[1] <- "phyloP.score"
	#bed$phyloP.score <- conservation$phyloP.score[match]
	conservation$exon.type <- bed.names[i]
	# there is still missing data in PhyloP so not all of my exons are covered. Grrr!
	#bed$conservation.score <- conservation
	conservation.out <- paste0(bed.names[i],".conservation")
	#assign(conservation.out,bed)
	#assign(conservation.out,conservation)
	#row.names(bed)
	assign(conservation.out,conservation)
	exon.conservation[[i]] <- conservation
}

# read in the cryptic exon list
cryptic.exons.flank.bed <- as.data.frame(fread(paste0(bedfiles.folder,"/",bed.names[1],".exons.bed")))
names(cryptic.exons.flank.bed) <- c("chr","start","end","gene.name","exonID","strand","iteration")

# write out table of conservation scores for each exon
cryptic_exons.conservation$gene.name <- cryptic.exons.flank.bed$gene.name[match(rownames(cryptic_exons.conservation), rownames(cryptic.exons.flank.bed))]
cryptic_exons.conservation$exonID <- cryptic.exons.flank.bed$exonID[match(rownames(cryptic_exons.conservation), rownames(cryptic.exons.flank.bed))]
cryptic.conservation.table <- paste0(conservation.outFolder,code,"_cryptic_exon_conservation.tab")
write.table(cryptic_exons.conservation, cryptic.conservation.table, sep="\t",row.names=F,quote=F)

# prepare tables for graphing
exon.conservation.merge <- as.data.frame(do.call(args = exon.conservation, what = rbind))
conservation.graph <- paste0(conservation.outFolder,code,"_exon_conservation_plot.pdf")

if(grepl("Cleveland_TDP43", code)){
    code <- "Adult_mouse_brain_TDP43_knockdown"
  }
  if(grepl("Chiang",code)){
    code <- "Mouse_ES_Cell_TDP43_deletion"
  }
  if(grepl("dataset",code)){
    code <- paste0(gsub("dataset","Human K562 Cell TDP43 knockdown", code))
  }
  graph_title <- paste0( gsub("_"," ", title) )
  

# get order of exons the right way round!
exon.conservation.merge$exon.type <- factor(exon.conservation.merge$exon.type, levels = bed.names)
exon.n <- as.numeric(table(exon.conservation.merge$exon.type))
bed.labels <- paste0(gsub("_"," ",bed.names), "\n(",exon.n,")")

pdf(conservation.graph)
p <- ggplot(exon.conservation.merge, aes(x = gsub(".","\n",exon.type,fixed =T), y = phyloP.score, fill = exon.type, colour = exon.type)) + 
	geom_boxplot(colour = "black",notch = T ) +
	xlab("") + 
	ylab("average phyloP score") + 
	ylim(c(-1,6)) +
	scale_fill_manual(values = c("skyblue","firebrick3","skyblue","skyblue")) +
	scale_x_discrete(name="Exon type",
					limits=bed.names,
                    labels=bed.labels) +
	ggtitle(paste0(graph_title,"\naverage phyloP conservation score per exon")) + 
	theme_bw() + 
	theme(panel.background = element_rect(fill = "white"),
          plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
          axis.line = element_line(colour = "black"),
          legend.position="none"
          )
print(p)
dev.off()

# two-sample t test of cryptic vs adjacent
test <- t.test(cryptic_exons.conservation$phyloP.score, cryptic_nulls.conservation$phyloP.score)
test.results <- paste0(conservation.outFolder,code,"_cryptic_vs_null_t_test.tab")
capture.output(test, file = test.results )


