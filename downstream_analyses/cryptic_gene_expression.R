# check for expression trends in the cryptic exons
# load each dataset and create a vector of cryptic exon gene names
# then load the DESeq gene expression tables and filter just the differential expression of the cryptic exon genes
# create boxplots of differential expression of cryptic exon containing genes in the 4 different datasets.
library(data.table)
library(ggplot2)
library(dplyr)

#mouse.ling <- c("2610507B11Rik",    "A230046K03Rik",     "Adipor2",      "Adnp2",        "Ahnak",        "Atraid",       "Cluh", "Edem2",        "Ercc6",        "Fam21",        "Fam73a",       "Flnb", "Ggct", "Gsta4",        "Gtf2e2",       "Hace1",        "Hgsnat",       "Ift81",        "Lnp",  "Mettl6",       "Mib1", "Mier1",        "Necap1",       "Nme6", "Pir",  "Pno1", "Ppp6c",        "Ptcd2",        "Pycr2",        "Sars", "Smg5", "Smg5", "Snapc3",       "Spata13",      "Spata7",       "Spcs2",        "Sptbn4",       "Sulf1",        "Synj2bp",      "Tecpr1",       "Tnfaip1",      "Tnks2",        "Tnnt1",        "Trim8",        "Uggt2",        "Usp15",        "Usp15",        "Wbscr22",      "Zfp13",        "Zfp809")
#human.ling <- c("EPB41L4A",    "CEP72",  "INSR", "FAM114A2",     "PFKP", "ST5",  "RNFT2",        "RNFT2",        "ALX1", "AGRN", "AGRN", "ATG4B",        "AGRN", "AGRN", "ST5",  "SETD5",        "KDELC2",       "MUC16",        "PKN1", "IRF9", "UPF2", "GPSM2",        "XPO4", "RASA4",        "RASA4B",       "PARP6",        "KRT7", "TRAPPC12",     "RANBP1",       "HERC6",        "BLZF1",        "ZFP91",        "HDGFRP2",      "MAP3K8",       "SSFA2",        "CENPK",        "ITPR3",        "KYNU", "IRF9", "COL4A6",       "KYNU") 

options(echo=T)
#setwd("/Users/Jack/project/")
#outFolder <- "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1"

#splicing_analysis.res <- "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_control_TDP_splicing_analysis.tab"


# read in the splicing analysis and create a vector of cryptic exon ensemblIDs
dataset.names <- c("Mouse adult brain","Mouse embryonic stem cell","K562 mRNA","K562 total RNA")
setwd("/cluster/project8/vyp/")
splicing.analysis.res  <- c("Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/splice_junction_analysis/Cleveland_TDP43_CTL_TDP_splicing_analysis.tab",
                "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Chiang_processed/splice_junction_analysis/Chiang_processed_CTL_TDP_splicing_analysis.tab",
                "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/splice_junction_analysis/dataset_1_control_TDP_splicing_analysis.tab",
				"Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/splice_junction_analysis/dataset_2_control_TDP_splicing_analysis.tab")

deseq.res <- c("Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/expression/deseq2/CTL_TDP/deseq_Cleveland_TDP43_differential_expression.tab",
				"Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Chiang_processed/expression/deseq2/CTL_TDP/deseq_Chiang_processed_differential_expression.tab",
				"Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/expression/deseq2/control_TDP/deseq_dataset_1_differential_expression.tab",
				"Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/expression/deseq2/control_TDP/deseq_dataset_2_differential_expression.tab")
 
outFolder <- "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/"

# create a cryptic_exon_gene_expression folder for both species
for(i in 1:length(unique(species.code))){
	expression.folder <- paste0(outFolder,unique(species.code)[i],"/cryptic_gene_expression")
	if(!file.exists(expression.folder)){
	dir.create(expression.folder)}
}


dataset.code <- c("cleveland","chiang","mRNA","total")
species.code <- c("mouse","mouse","human","human")
# I think this is only still required to generate named lists of cryptic exons 
#deseq.res.list <- list()
for(i in 1:length(dataset.names)){
	# load in cryptic exon analysis and subset just the upward cryptics
	crypt.res <- as.data.frame(fread(splicing.analysis.res[i]))
	crypt.res$dataset <- dataset.names[i]
	crypt.res <- subset(crypt.res, class == "SJ.SUPPORTED.UP")
	assign(dataset.code[i],crypt.res)
}
#	# load in differential expression analysis
#	deseq <- as.data.frame(fread(deseq.res[i]))
#	cryptic.ids <- d$EnsemblID
#	# subset just the differentially expressed cryptic exon genes
#	deseq <- deseq[deseq$EnsemblID %in% cryptic.ids & deseq$padj < 0.1,]
#	deseq$dataset <- paste0(dataset.names[i],"\n(",nrow(deseq),")")
#	# create a reduced table of the significant cryptic exon genes
#	deseq.res.list[[i]] <- data.frame(ID = deseq$EnsemblID, FC = deseq$log2FoldChange, FDR = deseq$padj, dataset = deseq$dataset)
#}
#
#deseq <- do.call(what = rbind, args = deseq.res.list)
#
#p <- ggplot(deseq, aes(y = FC, x = dataset, fill = dataset)) + 
#	geom_hline(yintercept = 0, linetype = 1, colour = "gray") + 
#	geom_boxplot(notch=F) + 
#	ylab("log2FoldChange") + 
#	ggtitle(paste0("Cryptic exon gene expression (FDR < 0.1)")) +
#	theme(legend.position="none") 
#
#pdf(paste0(outFolder,"supp_cryptic_exon_gene_expression.pdf"))
#print(p)
#dev.off()


# slightly different. Take the the union lists of cryptic exons for the unions of the two mouse and two human.
# does gene expression trend down for all the genes in the union in each dataset? Despite not all the cryptic exons being shared across both?

# create union list of cryptic exons for each species
mouse.union <- rbind(cleveland,chiang)
human.union <- rbind(mRNA,total)
# for each dataset calculate gene expression signature for each cryptic exon gene
deseq.union.res.list <- list()
for( i in 1:4){
	if(i <= 2){
		union <- mouse.union
		ling <- mouse.ling }
        if(i > 2){
		union <- human.union
		ling <- human.ling }
	# read in differential expression data and subset just the significant changes
	deseq <- as.data.frame(fread(deseq.res[i]))
    deseq$dataset <- dataset.names[i]
    deseq.sig <- deseq[deseq$padj < 0.1,]
	# read in the cryptic exon analysis
	specific <- as.data.frame(fread(splicing.analysis.res[i]))
	specific.cryptic.ids <- specific[specific$class == "SJ.SUPPORTED.UP",]$EnsemblID
	specific.cryptic.ids <- unique(specific.cryptic.ids)
	# retain just the differential expression for the cryptic exon genes
	specific.deseq <- filter(deseq.sig, grepl( paste(specific.cryptic.ids, collapse = "|" ), EnsemblID) )
	specific.deseq$group <- "cryptic exon genes"
	min.cryptic.baseMean <- min(specific.deseq$baseMean)
	# work out the identity of the cryptic exon genes found in the opposite dataset
	union.cryptic.ids <- union$EnsemblID
	union.cryptic.ids <- unique(union.cryptic.ids)
	other.cryptic.ids <- setdiff(union.cryptic.ids, specific.cryptic.ids)
	# retain the differential expresion results for only those genes
	other.deseq <- filter(deseq.sig, grepl( paste(other.cryptic.ids, collapse = "|" ), EnsemblID) )
	other.deseq$group <- "other cryptic exon genes"
	
	# write out the gene expression for the specific and other cryptic exon genes for each dataset.
	genes.out.table <- rbind(other.deseq, specific.deseq)
	genes.out <- paste0(outFolder, species.code[i],"/cryptic_gene_expression/",dataset.code[i],"_cryptic_gene_expression.tab")
	write.table(genes.out.table, genes.out, quote=F, sep="\t", row.names=F)
	# use all remaining genes as a control. BUT keep only those that have a meanbase above the lowest meanbase of all the cryptic exon genes.
 
	all.genes <- deseq.sig
	all.genes <- filter(all.genes, baseMean >= min.cryptic.baseMean)
	all.genes$group <- "all genes"
 
	together <- rbind(all.genes,specific.deseq, other.deseq)
	deseq.union.res.list[[i]] <- data.frame( ID = together$EnsemblID, FC = together$log2FoldChange, FDR = together$padj, dataset = together$dataset, group = together$group )

}

deseq.union <- do.call(what = rbind, args = deseq.union.res.list)
deseq.union$direction <- ifelse(deseq.union$FDR >= 0.1, NA, ifelse(deseq.union$FC > 0, "UP","DOWN") )
deseq.table <- table(deseq.union$dataset, deseq.union$direction)

deseq.union <- deseq.union[!is.na(deseq.union$dataset),]




# even better idea! each dataset's DESeq output should be mined with four lists of genes:
# * all genes (FDR < 0.1) <- just keep the median value for this
# * cryptic exon genes discovered in that dataset
# * cryptic exon genes discovered in both datasets of that species
# * cryptic exon genes discovered in both datasets PLUS the cryptic exon genes claimed by Ling


rpkms <- c("Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Cleveland_TDP43/expression/deseq2/rpkm_values.csv",
				"Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/paper_TDP_mouse/Chiang_processed/expression/deseq2/rpkm_values.csv",
				"Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_1/expression/deseq2/rpkm_values.csv",
				"Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP_ENCODE_human/dataset_2/expression/deseq2/rpkm_values.csv")
 
 rpkm.union.res.list <- list()
 for(i in 1:length(rpkms)){
 	if(i <= 2){
		union <- mouse.union
		ling <- mouse.ling }
        if(i > 2){
		union <- human.union
		ling <- human.ling }

 	r <- as.data.frame(fread(rpkms[i],na.strings = ""))
 	r <- mutate(r, control.mean.rpkm = rowMeans(select(r, starts_with("c",ignore.case=T) ) ) )
 	all.genes <- r
 	all.genes$group <- "all genes"
 	median <- data.frame(ID = "median",control.mean.rpkm = median(all.genes$control.mean.rpkm), dataset = dataset.names[i], group = "control.median", is.zero = F)

 	specific <- as.data.frame(fread(splicing.analysis.res[i],na.strings = ""))
	specific.cryptic.ids <- specific[specific$class == "SJ.SUPPORTED.UP",]$EnsemblID
	specific.cryptic.ids <- unique(specific.cryptic.ids)
	specific.rpkm <- filter(r, grepl( paste(specific.cryptic.ids, collapse = "|" ), ensemblID) )
	specific.rpkm$group <- "cryptic exon genes"

	union.cryptic.ids <- union$EnsemblID
	union.cryptic.ids <- unique(union.cryptic.ids)
	other.cryptic.ids <- setdiff(union.cryptic.ids, specific.cryptic.ids)
	other.rpkm <- filter(r, grepl( paste(other.cryptic.ids, collapse = "|" ), ensemblID) )
	other.rpkm$group <- "other cryptic exon genes"

	together <- rbind(all.genes,specific.rpkm, other.rpkm)
	together$dataset <- dataset.names[i]
	together <- data.frame( ID = together$ensemblID, control.mean.rpkm = together$control.mean.rpkm, dataset = together$dataset, group = together$group )
	together$is.zero <- ifelse(together$control.mean.rpkm < 1, T,F)

	rpkm.union.res.list[[i]] <- rbind(median, together)
 }
rpkm.union <- do.call(what = rbind, args = rpkm.union.res.list)
rpkm.table <- table(rpkm.union$dataset, rpkm.union$group)
rpkm.union.graph <- subset(rpkm.union, group != "all genes")
#
#p3 <- ggplot(rpkm.union.graph, aes(x = group, y = control.mean.rpkm, fill = dataset) ) +
#		geom_boxplot(colour = "black") + 
#		#geom_jitter(width=0.4, alpha = 0.8) +
#		#geom_dotplot(binaxis = "y", stackdir = "center") +
#		#geom_bar() +
#		theme(legend.position="none",
#		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) +
#		facet_wrap(~dataset,ncol=2) +
#		ylab("log10(mean RPKM per gene across control samples)") + 
#		ylim(c(0,100)) +
#		xlab("")
#
#
#p2 <- ggplot(deseq.union, aes(y = FC, x = group, fill = dataset, label = length(group) )) + 
#    geom_hline(yintercept = 0, linetype = 1, colour = "gray") + 
#	geom_boxplot(outlier.shape = NA) + 
#	#geom_text(position = position_dodge(0.9)) +
#	ylim(c(-3,3) ) +
#        ylab("log2FoldChange") + 
#        ggtitle(paste0("Differential expression control vs depletion(FDR < 0.1)")) +
#        theme(legend.position="none",
#		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) +
#	facet_wrap(~dataset, ncol = 2)
#
#
#
#
#expression.pdf <- (paste0(outFolder,"supp_cryptic_exon_gene_expression.pdf"))
#pdf(expression.pdf)
#print(p3)
#print(p2)
#grid.arrange(tableGrob(rpkm.table),tableGrob(deseq.table), nrow = 2)
#dev.off()
#

d.sum <- group_by(deseq.union, dataset, group, direction ) %>% summarise(N = length(group) )

#d.sum <- filter(d.sum, !grepl("other",group))

rpkm.filter <- filter(rpkm.union, !grepl("median",group ) )
rpkm.group <- group_by(rpkm.filter, dataset, group )
rpkm.summarise <- summarise(rpkm.group, N = length(group))

# stupid way just to get total numbers of genes
d.sum$total <- rpkm.summarise$N[match(paste(d.sum$dataset, d.sum$group), paste(rpkm.summarise$dataset, rpkm.summarise$group))]
d.sum <- mutate(d.sum, prop = (N/total) )
d.sum$CI.plus <- 1
d.sum$CI.minus <- 1
for(i in 1:nrow(d.sum) ){
	d.sum$CI.minus[i] <- binom.test(x = d.sum$N[i], n = d.sum$total[i])$conf.int[1]
	d.sum$CI.plus[i] <- binom.test(x = d.sum$N[i], n = d.sum$total[i])$conf.int[2]
}

d.sum.all.genes <- filter(d.sum, group == "all genes")
d.sum.cryptics <- filter(d.sum, group == "cryptic exon genes")
d.sum.other <- filter(d.sum, group == "other cryptic exon genes")
names(d.sum.cryptics) <- paste0("cryptic.", names(d.sum.cryptics))
names(d.sum.other) <- paste0("other.",names(d.sum.other))

d.sum.all.genes$cryptic.N <- d.sum.cryptics$cryptic.N[match(paste(d.sum.all.genes$dataset,d.sum.all.genes$direction), paste(d.sum.cryptics$cryptic.dataset, d.sum.cryptics$cryptic.direction))]
d.sum.all.genes$cryptic.total <- d.sum.cryptics$cryptic.total[match(paste(d.sum.all.genes$dataset,d.sum.all.genes$direction), paste(d.sum.cryptics$cryptic.dataset, d.sum.cryptics$cryptic.direction))]
d.sum.all.genes$cryptic.prop <- d.sum.cryptics$cryptic.prop[match(paste(d.sum.all.genes$dataset,d.sum.all.genes$direction), paste(d.sum.cryptics$cryptic.dataset, d.sum.cryptics$cryptic.direction))]

# other dataset's cryptic exon genes
d.sum.all.genes$other.N <- d.sum.other$other.N[match(paste(d.sum.all.genes$dataset,d.sum.all.genes$direction), paste(d.sum.other$other.dataset, d.sum.other$other.direction))]
d.sum.all.genes$other.total <- d.sum.other$other.total[match(paste(d.sum.all.genes$dataset,d.sum.all.genes$direction), paste(d.sum.other$other.dataset, d.sum.other$other.direction))]
d.sum.all.genes$other.prop <- d.sum.other$other.prop[match(paste(d.sum.all.genes$dataset,d.sum.all.genes$direction), paste(d.sum.other$other.dataset, d.sum.other$other.direction))]

d.cast <- d.sum.all.genes
d.cast$cryptic.phyper <- 1
for(i in 1:nrow(d.cast)){
	d.cast$cryptic.phyper[i] <- phyper(q = d.cast$cryptic.N[i] - 1, m = d.cast$N[i], n = ( d.cast$total[i] - d.cast$N[i] ), k = ( d.cast$cryptic.total[i] - d.cast$cryptic.N[i]), lower.tail = F )
}}
d.cast$cryptic.p.adjust <- p.adjust(p = d.cast$cryptic.phyper, method = "bonferroni")
d.cast$cryptic.sig.stars <- ifelse(test = d.cast$cryptic.p.adjust >= 0.05, yes = NA, no = ifelse(test = d.cast$cryptic.p.adjust < 0.05 & d.cast$cryptic.p.adjust > 0.001, yes = "*", no = ifelse(d.cast$cryptic.p.adjust < 0.001 & d.cast$cryptic.p.adjust > 1e-16, yes = "**", no ="***")))


d.cast$other.phyper <- 1
for(i in 1:nrow(d.cast)){
	d.cast$other.phyper[i] <- phyper(q = d.cast$other.N[i] - 1, m = d.cast$N[i], n = ( d.cast$total[i] - d.cast$N[i] ), k = ( d.cast$other.total[i] - d.cast$other.N[i]), lower.tail = F )
}}

all.p.values <- c(d.cast$cryptic.phyper,d.cast$other.phyper)
all.p.adjust <- p.adjust(p = all.p.values, method = "bonferroni")

d.cast$cryptic.p.adjust <- all.p.adjust[1:length(d.cast$cryptic.phyper)]
d.cast$cryptic.sig.stars <- ifelse(test = d.cast$cryptic.p.adjust >= 0.05, yes = NA, no = ifelse(test = d.cast$cryptic.p.adjust < 0.05 & d.cast$cryptic.p.adjust > 0.001, yes = "*", no = ifelse(d.cast$cryptic.p.adjust < 0.001 & d.cast$cryptic.p.adjust > 1e-16, yes = "**", no ="***")))

d.cast$other.p.adjust <- all.p.adjust[(length(d.cast$cryptic.phyper) + 1):length(all.p.adjust)]
d.cast$other.sig.stars <- ifelse(test = d.cast$other.p.adjust >= 0.05, yes = NA, no = ifelse(test = d.cast$other.p.adjust < 0.05 & d.cast$other.p.adjust > 0.001, yes = "*", no = ifelse(d.cast$other.p.adjust < 0.001 & d.cast$other.p.adjust > 1e-16, yes = "**", no ="***")))



write.table(d.cast, "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Figure_2_differential_expression.tab", sep="\t",quote=F, row.names=F)

p.de.other <- ggplot(d.sum, aes(x = paste(group,direction), fill = dataset, label = N ) ) +
		geom_bar(aes(y = prop), stat = "identity",colour = "black") + 
		geom_errorbar(aes(ymin=CI.minus, ymax=CI.plus), width = 0.2 ) +
		#geom_text(aes(y = 1) ) +
		facet_wrap(~dataset,ncol=2) + 
		theme(legend.position="none",
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) + 
		scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
		xlab("") +
		xlim("all genes DOWN","cryptic exon genes DOWN", "other cryptic exon genes DOWN", "all genes UP", "cryptic exon genes UP","other cryptic exon genes UP") + 
		ylab("percentage of genes in group differentially expressed")

d.sum.graph <- subset(d.sum, group != "other cryptic exon genes")
p.de.no.other <- ggplot(d.sum.graph, aes(x = paste(group,direction), fill = dataset, label = N ) ) +
		geom_bar(aes(y = prop), stat = "identity",colour = "black") + 
		geom_errorbar(aes(ymin=CI.minus, ymax=CI.plus), width = 0.2 ) +
		#geom_text(aes(y = 1) ) +
		facet_wrap(~dataset,ncol=2) + 
		theme(legend.position="none",
		axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) + 
		scale_y_continuous(labels = scales::percent, limits = c(0,1)) + 
		xlab("") +
		xlim("all genes DOWN","cryptic exon genes DOWN", "all genes UP", "cryptic exon genes UP") + 
		ylab("percentage of genes in group differentially expressed")


pdf("Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Figure_2_differential_expression_enrichment.pdf")
print(p.de.no.other)
print(p.de.other)
dev.off()	

# check out the Kapeli dataset - is there evidence of cryptic exon mediated downregulation?
kapeli.rpkm <- "/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/TDP43_human/Kapeli/expression/deseq2/Control_TDP43/deseq_Kapeli_differential_expression.tab"
k.rpkm <- as.data.frame(fread(kapeli.rpkm))
union.cryptic.ids <- human.union$EnsemblID
union.cryptic.ids <- unique(union.cryptic.ids)
k.overlap <- filter(k.rpkm, grepl( paste(union.cryptic.ids, collapse = "|" ), EnsemblID) & padj < 0.1 )
k.sig <- filter(k.rpkm, padj < 0.1)


# rewrite tersely just what I want - the numbers of total and significant genes for each dataset, cryptic, all(expressed as much as cryptic) and other.
all.datasets <- list()
min.baseMean.list <- c()	
# for each dataset
for(i in 1:4){
	deseq <- as.data.frame(fread(deseq.res[i]))
	cryptic <- as.data.frame(fread(splicing.analysis.res[i]))	
	cryptic <- filter(cryptic, class == "SJ.SUPPORTED.UP")
	cryptic.names <- cryptic$fix.gene.names[cryptic$fix.gene.names != 0] # remove any empty gene names	
	cryptic.names <- paste0("^", cryptic.names, "$")
	cryptic.ids <- cryptic$EnsemblID[ cryptic$EnsemblID != 0 ] # remove any empty EnsemblIDs
	# keep only cryptic genes. match on both gene name (exact match) and EnsemblID (break up multiples with gsub)
	cryptic.genes <- filter(deseq, 
		grepl( paste(cryptic.names, collapse = "|" ), external_gene_id) |
		grepl( gsub("+", "|", paste(cryptic.ids, collapse = "|" ), fixed =T), EnsemblID ) )
	# remove any genes that do not have baseMean >= to the lowest cryptic exon
	min.baseMean <- min(cryptic.genes$baseMean)
	min.baseMean.list[i] <- min.baseMean
	all.genes <- filter(deseq, baseMean >= min.baseMean)
	
	all.genes.up <- filter(all.genes, log2FoldChange > 0)
	all.genes.down <- filter(all.genes, log2FoldChange < 0)
	cryptic.genes.up <- filter(cryptic.genes, log2FoldChange > 0 )
	cryptic.genes.down <- filter(cryptic.genes, log2FoldChange < 0)
	
	genes <- list(all.genes,cryptic.genes, all.genes, cryptic.genes)
	genes.directions <- list(all.genes.up, cryptic.genes.up, all.genes.down, cryptic.genes.down)
	N <- c()
	sig <- c()
	for(j in 1:4){
		N[j] <- nrow(genes[[j]])
		sig[j] <- nrow( filter( genes.directions[[j]], padj < 0.1) )
	}
	dataset <- data.frame(dataset = dataset.names[i], 
			      direction = c(rep("UP",2),rep("DOWN",2) ),
			      type = rep(c("all genes","cryptic exon genes"),2),
			      N = N,
			      sig = sig)		   
	all.datasets[[i]] <- dataset
}

all.datasets <- do.call(rbind, all.datasets)
all.datasets$prop <- all.datasets$sig / all.datasets$N
all.datasets$CI.plus <- 1
all.datasets$CI.minus <- 1
for(i in 1:nrow(all.datasets) ){
        all.datasets$CI.minus[i] <- binom.test(x = all.datasets$sig[i], n = all.datasets$N[i])$conf.int[1]
        all.datasets$CI.plus[i] <- binom.test(x = all.datasets$sig[i], n = all.datasets$N[i])$conf.int[2]
}
# conduct hypergeometric tests on each pair of values.
all.datasets$hypergeometric.res <- NA
for(i in seq(1, nrow(all.datasets), 2)){
	p <- phyper( q = all.datasets$sig[i+1] - 1, m = all.datasets$sig[i], n = (all.datasets$N[i] - all.datasets$sig[i] ), k = all.datasets$N[i + 1], lower.tail =F )
	all.datasets$hypergeometric.res[i] <- p
} 
all.datasets$sig.stars <- ifelse(test =all.datasets$hypergeometric.res >= 0.05, yes = NA, no = ifelse(test = all.datasets$hypergeometric.res < 0.05 & all.datasets$hypergeometric.res > 0.001, yes = "*", no = ifelse(all.datasets$hypergeometric.res < 0.001 & all.datasets$hypergeometric.res > 1e-16, yes = "**", no ="***")))


p.same.datasets <- ggplot(all.datasets, aes(x = paste(type,direction) ) ) +
                geom_bar(aes(y = prop), stat = "identity",colour = "black", fill = "deepskyblue2") +   
		geom_errorbar(aes(ymin=CI.minus, ymax=CI.plus), width = 0.2 ) +
                #geom_text(aes(y = 1) ) +
                facet_wrap(~dataset,ncol=2) +
                theme_bw() +
                theme(legend.position="none",
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) +
                scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
                xlab("") +
                xlim("all genes DOWN","cryptic exon genes DOWN", "all genes UP", "cryptic exon genes UP") +
                ylab("percentage of genes in group differentially expressed")
graph.out <- "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Figure_2_differential_expression_enrichment_min_baseMean.pdf"
ggsave(graph.out, plot = p.same.datasets)
same.table.out <- "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Figure_2_differential_expression_enrichment_min_baseMean_same.tab"
write.table(format(all.datasets, digits = 2),same.table.out, sep ="\t", row.names = F, quote=F)

# repeat but this time test the cryptic exon genes from the other dataset.
opposite.splicing.analysis.res <- splicing.analysis.res[c(2,1,4,3)]
opposite.datasets <- list()	
# for each dataset
for(i in 1:4){
	deseq <- as.data.frame(fread(deseq.res[i]))
	# get the gene names and EnsemblIDs of each dataset's cryptic exons
	cryptic <- as.data.frame(fread(splicing.analysis.res[i]))	
	cryptic <- filter(cryptic, class == "SJ.SUPPORTED.UP")
	cryptic.names <- cryptic$fix.gene.names[cryptic$fix.gene.names != 0] # remove any empty gene names	
	cryptic.names <- paste0("^", cryptic.names, "$")
	cryptic.ids <- cryptic$EnsemblID[ cryptic$EnsemblID != 0 ]

	# do the same for the opposite dataset's cryptic exons
	other.cryptic <- as.data.frame(fread(opposite.splicing.analysis.res[i]))	
	other.cryptic <- filter(other.cryptic, class == "SJ.SUPPORTED.UP")
	other.cryptic.names <- other.cryptic$fix.gene.names[other.cryptic$fix.gene.names != 0] # remove any empty gene names	
	other.cryptic.names <- paste0("^", other.cryptic.names, "$")
	other.cryptic.ids <- other.cryptic$EnsemblID[ other.cryptic$EnsemblID != 0 ] # remove any empty EnsemblIDs
	
	# retain only the names and IDs of genes unique to the other dataset
	unique.cryptic.names <- setdiff(other.cryptic.names, cryptic.names)
	unique.cryptic.ids <- setdiff(other.cryptic.ids, cryptic.ids)



	# keep only cryptic genes. match on both gene name (exact match) and EnsemblID (break up multiples with gsub)
	cryptic.genes <- filter(deseq, 
		grepl( paste(unique.cryptic.names, collapse = "|" ), external_gene_id) |
		grepl( gsub("+", "|", paste(unique.cryptic.ids, collapse = "|" ), fixed =T), EnsemblID ) )
	# remove any genes that do not have baseMean >= to the lowest cryptic exon
	min.baseMean <- min.baseMean.list[i]
	all.genes <- filter(deseq, baseMean >= min.baseMean)
	
	all.genes.up <- filter(all.genes, log2FoldChange > 0)
	all.genes.down <- filter(all.genes, log2FoldChange < 0)
	cryptic.genes.up <- filter(cryptic.genes, log2FoldChange > 0 )
	cryptic.genes.down <- filter(cryptic.genes, log2FoldChange < 0)
	
	genes <- list(all.genes,cryptic.genes, all.genes, cryptic.genes)
	genes.directions <- list(all.genes.up, cryptic.genes.up, all.genes.down, cryptic.genes.down)
	N <- c()
	sig <- c()
	for(j in 1:4){
		N[j] <- nrow(genes[[j]])
		sig[j] <- nrow( filter( genes.directions[[j]], padj < 0.1) )
	}
	dataset <- data.frame(dataset = dataset.names[i], 
			      direction = c(rep("UP",2),rep("DOWN",2) ),
			      type = rep(c("all genes","cryptic exon genes"),2),
			      N = N,
			      sig = sig)		   
	opposite.datasets[[i]] <- dataset
}

opposite.datasets <- do.call(rbind, opposite.datasets)
opposite.datasets$prop <- opposite.datasets$sig / opposite.datasets$N
opposite.datasets$CI.plus <- 1
opposite.datasets$CI.minus <- 1
for(i in 1:nrow(opposite.datasets) ){
        opposite.datasets$CI.minus[i] <- binom.test(x = opposite.datasets$sig[i], n = opposite.datasets$N[i])$conf.int[1]
        opposite.datasets$CI.plus[i] <- binom.test(x = opposite.datasets$sig[i], n = opposite.datasets$N[i])$conf.int[2]
}
# conduct hypergeometric tests on each pair of values.
opposite.datasets$hypergeometric.res <- NA
for(i in seq(1, nrow(opposite.datasets), 2)){
	p <- phyper( q = opposite.datasets$sig[i+1] - 1, m = opposite.datasets$sig[i], n = (opposite.datasets$N[i] - opposite.datasets$sig[i] ), k = opposite.datasets$N[i + 1], lower.tail =F )
	opposite.datasets$hypergeometric.res[i] <- p
} 
opposite.datasets$sig.stars <- ifelse(test =opposite.datasets$hypergeometric.res >= 0.05, yes = NA, no = ifelse(test = opposite.datasets$hypergeometric.res < 0.05 & opposite.datasets$hypergeometric.res > 0.001, yes = "*", no = ifelse(opposite.datasets$hypergeometric.res < 0.001 & opposite.datasets$hypergeometric.res > 1e-16, yes = "**", no ="***")))


p.opposite.datasets <- ggplot(opposite.datasets, aes(x = paste(type,direction) ) ) +
                geom_bar(aes(y = prop), stat = "identity",colour = "black", fill = "deepskyblue2") +   
		geom_errorbar(aes(ymin=CI.minus, ymax=CI.plus), width = 0.2 ) +
                #geom_text(aes(y = 1) ) +
                facet_wrap(~dataset,ncol=2) +
                theme_bw() +
                theme(legend.position="none",
                axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1) ) +
                scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
                xlab("") +
                xlim("all genes DOWN","cryptic exon genes DOWN", "all genes UP", "cryptic exon genes UP") +
                ylab("percentage of genes in group differentially expressed")
opposite.graph.out <- "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Figure_2_differential_expression_enrichment_min_baseMean_opposite.pdf"
opposite.table.out <- "Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Figure_2_differential_expression_enrichment_min_baseMean_opposite.tab"
write.table(format(opposite.datasets, digits = 2),opposite.table.out, sep ="\t", row.names = F, quote=F)
ggsave(opposite.graph.out, plot = p.opposite.datasets)

