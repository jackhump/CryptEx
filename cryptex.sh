# CRYPTEX

# Jack Humphrey 
# UCL

##Given a BAM file, extract the spliced reads and convert these to a BED file. Intersect with an intron bed file.
## Output should be a list of all spliced intronic reads. This will be the basis of my cryptic exon hunting
# this is to speed up debugging
set -euo pipefail

## these variables will be arguments passed to the pipeline script by the submission script
codeFolder=/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/CryptEx
Step2_master=${codeFolder}/Step2_master.sh

# R scripts for CryptEx
Rbin=/share/apps/R/bin/R
dexseqFinalProcessR=${codeFolder}/dexseq/forked_dexseq_pipeline_v2.R
countPrepareR=${codeFolder}/dexseq/forked_counts_prepare_pipeline.R
deseqFinalProcessR=${codeFolder}/dexseq/forked_deseq2_pipeline.R
R_support_chopper=${codeFolder}/support_frame_chopper.R

# R scripts for downstream analyses
R_splice_junction_analyzer=${codeFolder}/downstream_analyses/splice_junction_analyzer.R
R_functional_enrichment=${codeFolder}/downstream_analyses/Functional_Enrichment.R

# dexseq_count.py comes with HTSeq
pycount=${codeFolder}/dexseq_count.py
# the annotation files need to be flexible
# ideally the user should provide just one GFF file and then CryptEx extracts the exons and calculates the introns


#for testing only
species="mouse"
protein="FUS"
#support=${oFolder}/support/FUS_mouse_support.tab
submit="no"
splice_extractor=no
gff_creator=yes
read_counter=yes
DEXSeq=no
DESeq=no
cohort_merger=no


# Current strategy: merge all intronic splice features that are within 500bp of each other
strict_num=500


# Each step in the pipeline can be run independently of every other step but for efficiency purposes I'd like to be able to run each step in series automatically. 
# Each step's job script has unique ID and variable name

#case statement to match the argument variables coming in from the submission script
until [ -z $1 ];do
	case $1 in
	--species)
	shift
	species=$1;;
	--annotation_file)
	shift
	annotation_file=$1;;
	--protein)
	shift
	protein=$1;;
	--submit)
	shift
	submit=$1;;
	--code)
	shift
	code=$1;;
	--support)
	shift
	support=$1;;
	--gff)
	shift
	gff=$1;;
	--splice_extractor)
	shift
	splice_extractor=$1;;
	--gff_creator)
	shift
	gff_creator=$1;;
	--read_counter)
	shift
	read_counter=$1;;
	--DEXSeq)
	shift
	DEXSeq=$1;;
	--DESeq)
	shift
	DESeq=$1;;
	--intron_retainer)
	shift
	intron_retainer=$1;;
	--intron_DEXSeq)
	shift
	intron_DEXSeq=$1;;
	--cohort_merger)
	shift
	cohort_merger=$1;;
	--IGV)
	shift
	IGV=$1;;
	--splice_junction_analyzer)
	shift
	splice_junction_analyzer=$1;;
	--functional_enrichment)
	shift
	functional_enrichment=$1;;
	--hold_Step1)
	shift
	hold_Step1=$1;;
	--paired)
	shift
	paired=$1;;
	--stranded)
	shift
	stranded=yes
	libstrand=$1;;
	-* )
	echo "unrecognised argument: $1"
	exit 1;;
esac
shift
if [ "$#" = "0" ]; then break; fi
	echo $1 $2
done


oFolder=/cluster/project8/vyp/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/
results=${oFolder}/${protein}_${species}/${code}
reference=${oFolder}/reference
clusterFolder=${results}/cluster

for folder in ${oFolder} ${results} ${reference} ${clusterFolder} ${clusterFolder}/out ${clusterFolder}/error ${clusterFolder}/R  ${clusterFolder}/submission; do
    if [ ! -e $folder ]; then mkdir -p $folder; fi
done

# Users provide their own GFF file and their own GFF introns file. 
	# The CryptEx wiki will provide a link to Devon Ryan's intron script to create the GFF introns file.

#if [[ "$gff" == "hardcoded" ]]; then
	if [[ "$species" == "mouse" ]];then
	gff_base="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed"
	elif [[ "$species" == "human" ]];then
	gff_base="/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed"
	fi
#fi

#if [[ "$annotation_file" == "hardcoded" ]];then
	if [[ "$species" == "mouse" ]];then
	annotation_file=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/mouse/biomart/biomart_annotations_mouse.tab
	elif [[ "$species" == "human" ]];then
	 annotation_file=/cluster/project8/vyp/vincent/Software/RNASeq_pipeline/bundle/human_hg38/biomart/biomart_annotations_human.tab
	fi
#fi

intron_GFF=${gff_base}_introns_only.gff
exon_GFF=${gff_base}_exons_only.gff

# create intron_GFF and exon_GFF from the supplied GFF file. 




#The Intron GFF, created from the exon GFF using an R script written by Devon Ryan is to be used as a BED file.
#The BED file has information about strand,gene ID and crucially intron number. Each entry should in theory be unique.
intron_BED=${oFolder}/reference/${species}_full_introns.bed
intron_tweaked_GFF=${oFolder}/reference/${species}_introns_for_HTseq.gff
if [ ! -e $intron_BED ];then
	 awk 'BEGIN{OFS="\t"}{split($12,a,"\"");split($NF,b,"\"");print $1,$4,$5,$7,b[2]"_"a[2]}' $intron_GFF > $intron_BED
fi
# Devon Ryan's script labels each intron as "intronic sequence" in the GFF which HTSeq will throw out. So we have to alter it.
if [ ! -e $intron_tweaked_GFF ];then
	sed 's/intronic/exonic/g' ${intron_GFF} > ${intron_tweaked_GFF}
fi
intron_GFF=${intron_tweaked_GFF}


## Error Reporting
# I'd like each module of the pipeline to report to a central file so I can easily check the progress of a cohort.
report_file=${results}/report.txt
echo "CryptEx
Started at:	`date`
--Protein:	$protein
--Species:	$species
--Step1:	$splice_extractor
--hold_Step1:	$hold_Step1
--Step2:	$gff_creator
--Step3:	$read_counter
--Step4:	$DEXSeq
Support file:
" >> $report_file
cat $support >> $report_file

#################################
### STEP 1: SPLICE EXTRACTION ###
#################################

if [[ "$splice_extractor" == "yes" ]]

then 
echo "creating job scripts for spliced read extraction" 

#spliced intronic reads are extract for each dataset within a cohort - datasets no longer rely on each other!

#create an array job. This will contain a list of the individual jobs and execute them all together
## Step1_jobscript is the master job array.
# Each sample will have its own step1 script with the form splice_extract_${sample}
Step1_jobscript=${clusterFolder}/submission/Step1_${protein}_${species}.sh
# I want an array of jobs
sample_num_1=`wc -l ${support} | awk '{print $1 - 1}' `

Step1_jobID=Step1_${protein}_${species}_${code}
#strict doesn't need a different jobID right? If I want to the run the pipeline on regular or Strict mode at the same time then:
# Step 1 will run (called by regular pipeline and regular Step2 will wait for it
# Pipeline invoked with Strict parameter doesn't run Step1 but still waits for Step1 for finish so I've created a "hold_Step1" flag.
echo "
#$ -S /bin/bash
#$ -l h_vmem=5.90G
#$ -l tmem=5.9G
#$ -l h_rt=24:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -N ${Step1_jobID}
#$ -wd ${oFolder}
#$ -t 1-${sample_num_1}	
#$ -tc 20

if [[ \"\$SGE_TASK_ID\" == \"1\" ]];then
	echo \"
Step1 started at \`date +%H:%M:%S\`
\" >> $report_file
fi	
jobs=\"" > $Step1_jobscript


# load in information from the support file
awk 'NR >1{print $1,$2,$3}' $support | while read sample bam dataset; do
	#remember this is now for each bam file in the support file
	splicefolder=${results}/${dataset}/splice_extraction
        if [ ! -e $splicefolder ]; then mkdir -p $splicefolder; fi
	output=${splicefolder}/${dataset}_${sample}
	sample_jobscript=${clusterFolder}/submission/splice_extract_${protein}_${species}_${dataset}_${sample}.sh
	
	echo "
	(>&2 echo \$HOSTNAME)
	echo \$HOSTNAME
	mkdir /scratch0/CryptEx_tmp
	export TMPDIR=/scratch0/CryptEx_tmp
#extract spliced reads from bam, ignore header. -F 256 filters out any secondary aligning reads
#if [ ! -e ${output}.spliced.bam ]; then
samtools view -h -F 256 $bam | awk '\$1~/@/ || \$6~/N/' | samtools view -bh - > ${output}.spliced.bam 
#fi

#intersect spliced reads with the exon GFF file - only reads that intersect an exon are kept. 
bedtools intersect -a ${output}.spliced.bam -b ${exon_GFF} > ${output}.spliced.exons.bam

#if [ ! -e ${output}.spliced.bed ];then
#convert spliced bam to bed file
bedtools bamtobed -i ${output}.spliced.exons.bam -split | sort -k1,1 -k2,2n > ${output}.spliced.bed
#fi

#intersect with the exon gff file, outputting what falls into the introns as the intervals
bedtools intersect -a ${output}.spliced.bed -b ${exon_GFF} -v > ${output}.spliced.introns.bed

#remove intermediate files - this is put on hold while I develop the method
#rm ${output}.spliced.bed ${output}.spliced.bam ${output}.spliced.exons.bam

rm -rf /scratch0/CryptEx_tmp
	
echo \"Step 1 finished for $sample at \`date +%H:%M:%S\` \" >> $report_file 
" > $sample_jobscript


	echo $sample_jobscript >> $Step1_jobscript
done
echo "\"
script=\`echo \$jobs | cut -f\$SGE_TASK_ID -d \" \" \`
sh \$script
" >> $Step1_jobscript

fi
############################################
### STEP 2: BED MERGING AND GFF CREATION ###
############################################
# DEVELOPMENT VERSION
	# Merge now occurs within a dataset, not within a cohort.

if [ $gff_creator = "yes" ]
then
echo "creating job script for GFF creation"

Step2_jobscript=${clusterFolder}/submission/Step2_GFF_creator_${protein}_${species}.sh
#if Step2_jobscript already exists then remove it.
if [ -e $Step2_jobscript ]; then rm $Step2_jobscript;fi

Step2_jobID=Step2_${protein}_${species}_${code}
#extract number of datasets from the support file
sample_num_2=`cat $support | awk 'NR > 1{print $3}' | uniq | wc -l`
# Create the header for the master job script
echo "
#$ -S /bin/bash
#$ -l h_vmem=3.75G
#$ -l tmem=3.75G
#$ -l h_rt=24:00:00
#$ -pe smp 4
#$ -R y
#$ -o ${clusterFolder}/out
#$ -e ${report_file}
#$ -N ${Step2_jobID}
#$ -wd ${oFolder}
#$ -t 1-${sample_num_2} 
#$ -tc 20

if [[ \"\$SGE_TASK_ID\" == \"1\" ]];then
        echo \"
Step2 started at \`date +%H:%M:%S\`
\" >> $report_file
fi      
jobs=\"" > $Step2_jobscript

for dataset in `cat $support | awk 'NR > 1{print $3}' | uniq `;do
	if [ ! -e ${results}/${dataset}/GFF ];then 
		mkdir -p ${results}/${dataset}/GFF
	fi
	spliced_beds=${results}/${dataset}/splice_extraction/
	output=${results}/${dataset}/GFF/${protein}_${species}_${dataset}
	step2_dataset_script=${clusterFolder}/submission/GFF_creator_${protein}_${species}_${dataset}
# the "strict mode" should be hard coded  now.
	echo "
bash $Step2_master --dataset ${dataset} --output ${output} --intron_BED ${intron_BED} --exon_GFF ${exon_GFF} --spliced_beds ${spliced_beds}
	" > $step2_dataset_script

echo $step2_dataset_script >> $Step2_jobscript

done

echo "\"
script=\`echo \$jobs | cut -f\$SGE_TASK_ID -d \" \" \`
sh \$script
" >> $Step2_jobscript


fi


##############################################################
## STEP 3: READ COUNTING FOR EACH BAM WITH THE NEW GFF FILE###
##############################################################
# each subcohort within a cohort should be counted separately
#Submitted as an array job. Step3_master_jobscript is the master array job which is submitted.
# Step3_sample_jobscript are the individual jobs that are submitted by the master job.
# Is there an info.tab available with dataset-specific strand information?


if [[ "$read_counter" = "yes" ]];then 

		GFF=${results}/GFF/${protein}_${species}.total.cryptics.gff
Step3_master_jobscript=${clusterFolder}/submission/Step3_count_${protein}_${species}_${code}.sh

step3_num=`wc -l ${support} | awk '{print $1 -1}' `

Step3_jobID=Step3_${protein}_${species}_${code}

echo "
#$ -S /bin/bash
#$ -l h_vmem=3.6G
#$ -l tmem=3.6G
#$ -l h_rt=24:00:00
#$ -pe smp 2
#$ -R y
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -N $Step3_jobID
#$ -wd ${oFolder}
#$ -t 1-${step3_num}
#$ -tc 20

if [[ \"\$SGE_TASK_ID\" == \"1\" ]];then
        echo \"
Step3 started at \`date +%H:%M:%S\`
\" >> $report_file
fi
jobs=\"" > $Step3_master_jobscript

# separate cohort into subcohorts using sample table
## column 3 is 'dataset'. Each dataset should have its own count output
awk 'NR >1{print $1,$2,$3,$4}' $support | while read sample bam dataset condition;do
    GFF=${results}/${dataset}/GFF/${protein}_${species}_${dataset}.total.cryptics.gff

	info_table=${oFolder}/support/${protein}_${species}_info.tab
	echo "paired: $paired"
	echo "stranded: $stranded $libstrand"

	if [ -e ${info_table} ];then
		#echo "dataset:  $dataset
		#info_table: $info_table"
		#echo "Info table found for cohort"
		paired=`awk -v dataset=$dataset '$1==dataset {print $2}' $info_table`
		stranded=`awk -v dataset=$dataset '$1==dataset {print $3}' $info_table`
	else
	paired=$paired
	stranded=$stranded
	fi


# One jobscript per bam file
	countFolder=${results}/${dataset}/counts
	Step3_sample_jobscript=${clusterFolder}/submission/count_${dataset}_${sample}.sh

	output=${countFolder}/${sample}_dexseq_counts.txt
	if [ ! -e ${countFolder} ]; then mkdir -p ${countFolder};fi

# taken from the RNASeq pipelin 
	countStrand=no
	if [[ "$libstrand" == "fr-firststrand" ]]
	then
	  countStrand=yes
	  countStrandReverse=reverse
	elif [[ "$libstrand" == "fr-secondstrand" ]]
	then
	  countStrand=reverse
	  countStrandReverse=yes
	else
	    echo unknown libstrand $libstrand
	fi


	echo "
#BEDtools implementation
#bedtools multicov -bed $GFF -bams $bam -split | awk '{print \$14,\$12,\$NF}' | tr -d '\";' | awk '{print \$1 \":\" \$2 \"\t\" \$3}' > ${output}
#HTSeq implementation
#introduces paired and stranded arguments.
## 12/12/15 - hard-coding stranded=no as the counts were coming out near-empty. Try again with this.
## 14/10/16 - finally figured out why. Hah!

python $pycount --stranded=${countStrand} -p ${paired} -f bam -r pos $GFF $bam ${output}
echo \"Step 3 finished for $sample at \`date +%H:%M:%S\` \" >> $report_file 
" > $Step3_sample_jobscript
echo $Step3_sample_jobscript >> $Step3_master_jobscript
done

echo "\"
script=\`echo \$jobs | cut -f\$SGE_TASK_ID -d \" \"\`
sh \$script
" >> $Step3_master_jobscript
echo "creating job scripts for read counting"
fi

####################
## STEP 4: DEXSeq ##
####################

#### the R instructions are kept as a separate standalone script. The jobscript the pipeline outputs just invokes the R script with a bunch of preset variables
# the dexseq support frame needs to created from the cohort support frame for each subcohort within  
sanity_check=`awk 'NR == 1 {print $3}' $support `

if [ $DEXSeq = "yes" ] && [ $sanity_check = "dataset" ] ; then

Step4_master_jobscript=${clusterFolder}/submission/Step4_DEXSeq_${protein}_${species}.sh
dataset_num=`cat $support | awk 'NR > 1{print $3}' | uniq | wc -l`

Step4_jobID=Step4_${protein}_${species}_${code}

sample_num=`awk 'NR >1' $support | wc -l`
DEXSEQ_MEM=8
DEXSEQ_CORES=4

echo "sample number is: $sample_num "
if [ $sample_num -gt 8 ]; then
	DEXSEQ_CORES=12
	DEXSEQ_MEM=3.5
fi

echo "
#$ -S /bin/bash
#$ -l h_vmem=${DEXSEQ_MEM}G
#$ -l tmem=${DEXSEQ_MEM}G
#$ -l h_rt=36:00:00
#$ -pe smp ${DEXSEQ_CORES}
#$ -R y
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -N ${Step4_jobID}
#$ -wd ${oFolder}
#$ -l h_rt=24:00:00
#$ -t 1-${dataset_num}
#$ -tc 20

if [[ \"\$SGE_TASK_ID\" == \"1\" ]];then
        echo \"
Step4 started at \`date +%H:%M:%S\`
\" >> $report_file
fi
jobs=\"" > $Step4_master_jobscript

# do this on a per-cohort basis
for dataset in `cat $support | awk 'NR > 1{print $3}' | uniq `  
do

DEXSeqFolder=${results}/${dataset}/dexseq
if [ ! -e $DEXSeqFolder ]; then mkdir $DEXSeqFolder; fi
support_frame=${DEXSeqFolder}/${dataset}_dexseq_frame.tab
jobscript=${clusterFolder}/submission/DEXSeq_${dataset}.sh
GFF=${results}/${dataset}/GFF/${protein}_${species}_${dataset}.total.cryptics.gff
keepSex=TRUE
keepDups=FALSE
cryptic=TRUE
iFolder=${results}/${dataset}
error_file=${clusterFolder}/R/dexseq_${dataset}.out

## annotation_file is defined at the start

# slice up cohort support file into subcohort supports
awk -v dataset=$dataset 'NR == 1 {print $0} $3 == dataset {print $0}' $support > ${support_frame}

# for the downstream analysis of splice junctions. This uses the spliced.exons.bam files created in Step1

#bam_list=`cat $support_frame | awk 'BEGIN{ORS="\t"} NR>1 {print $2}' `
bam_list=`cat $support_frame | awk -F"/" -v results=$results -v dataset=$dataset 'BEGIN{ORS="\t"} NR>1 {split($NF,a,"_unique.bam");print results"/"dataset"/splice_extraction/"dataset"_" a[1]".spliced.exons.bam"}'`


sample_list=`cat $support_frame | awk 'BEGIN{ORS="\t"}NR>1{print $1}'`


echo "
# DEXSeq - Cryptic Exon Analysis

${Rbin}script ${dexseqFinalProcessR} --cryptic ${cryptic} --gff ${GFF} --keep.sex ${keepSex} --keep.dups ${keepDups} --support.frame ${support_frame} --code ${dataset} --annotation.file ${annotation_file} --iFolder ${iFolder}  > ${error_file} 2>&1

# #append the R output to the report
# cat ${error_file} >> ${report_file}


# #index the spliced bams
# for bam in $bam_list;do
#         samtools index \$bam
# done

# for i in `ls ${DEXSeqFolder}/`;do
# comparison=\`echo \$i | awk -F\"/\" '{print \$(NF-1)}'\`
#         if [ -e ${DEXSeqFolder}/\${comparison}/${dataset}_\${comparison}_CrypticExons.bed ]; then

#         bedtools multicov -bed ${DEXSeqFolder}/\${comparison}/${dataset}_\${comparison}_CrypticExons.bed -bams $bam_list -split > ${DEXSeqFolder}/\${comparison}/${dataset}_\${comparison}_SJ_analysis.bed 

#         echo \"chr start end gene.id strand intron.id log2FC FDR $sample_list\" | tr \" \" \"\t\" > ${DEXSeqFolder}/\${comparison}/header       
#         cat ${DEXSeqFolder}/\${comparison}/header ${DEXSeqFolder}/\${comparison}/${dataset}_\${comparison}_SJ_analysis.bed > tmp 
#         mv tmp  ${DEXSeqFolder}/\${comparison}/${dataset}_\${comparison}_SJ_analysis.bed  

#         fi
#done



" > $jobscript

echo $jobscript >> $Step4_master_jobscript

done

echo "\"
script=\`echo \$jobs | cut -f\$SGE_TASK_ID -d \" \"\`
sh \$script
" >> $Step4_master_jobscript
echo "creating job scripts for DEXSeq testing"
fi

############################
### STEP 4b: SJ ANALYZER ###
############################

if [[ "$splice_junction_analyzer" == "yes" ]]; then

Step4b_master_jobscript=${clusterFolder}/submission/Step4b_SJ_analyzer_${protein}_${species}.sh

dataset_num=`cat $support | awk 'NR > 1{print $3}' | uniq | wc -l`

Step4b_jobID=Step4b_${protein}_${species}

echo "
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=36:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -N ${Step4b_jobID}
#$ -wd ${oFolder}
#$ -l h_rt=24:00:00
#$ -t 1-${dataset_num}
#$ -tc 20

if [[ \"\$SGE_TASK_ID\" == \"1\" ]];then
        echo \"
Step4b started at \`date +%H:%M:%S\`
\" >> $report_file
fi
jobs=\"" > $Step4b_master_jobscript

# do this on a per-cohort basis and also by condition
for dataset in `cat $support | awk 'NR > 1{print $3}' | uniq `;do
	echo $dataset
	DEXSeqFolder=${results}/${dataset}/dexseq
	if [ ! -e $DEXSeqFolder ];then
		DEXSeqFolder=${results}/${dataset}/dexseq
	fi
	outFolder=${results}/${dataset}/splice_junction_analysis
		if [ ! -e $outFolder ]; then mkdir $outFolder; fi
	jobscript=${clusterFolder}/submission/SJ_analyzer_${dataset}.sh
			echo "
#SPLICE JUNCTION ANALYSIS FOR ${dataset}
" > $jobscript
#different datasets have different conditions (CTL vs HET, CTL vs HOM etc. Run SJ analysis for each.)
	for i in `ls ${results}/${dataset}/dexseq/`;do
		if [[ $i =~ .*[\.][a-z]+ ]];then
		echo "$i is not a valid set of conditions"
		continue;fi
		condition_names=`echo $i | awk -F"/" '{print $(NF-1)}'`
		echo $condition_names
        if [ -e ${DEXSeqFolder}/${condition_names}/${dataset}_${condition_names}_SignificantExons.csv ]; then
        	dexseq_res=${DEXSeqFolder}/${condition_names}/${dataset}_${condition_names}_SignificantExons.csv
        else 
        	echo "cannot find ${DEXSeqFolder}/${condition_names}/${dataset}_${condition_names}_CrypticExons.bed"
        	continue
        fi
		
		support_frame=${outFolder}/${dataset}_support_frame.tab
		error_file=${clusterFolder}/R/SJ_analyzer_${dataset}.out

		awk -v dataset=$dataset 'NR == 1 {print $0} $3 == dataset {print $0}' $support > ${support_frame}

		echo "
echo $condition_names
${Rbin}script ${R_splice_junction_analyzer} --support.frame ${support_frame} --code ${dataset} --species $species --condition.names ${condition_names} --dexseq.res ${dexseq_res} --outFolder ${outFolder}  > ${error_file} 2>&1
" >> $jobscript
			done

echo $jobscript >> $Step4b_master_jobscript

done

echo "\"
script=\`echo \$jobs | cut -f\$SGE_TASK_ID -d \" \"\`
sh \$script
" >> $Step4b_master_jobscript
echo "creating job scripts for splice junction analysis"
fi

#################################
## STEP 4c: FUNCTIONAL ENRICHMENT
##################################
# overlaps suspected cryptic exons with lists of repeat elements and peaks from iCLIP or eCLIP datasets. species specific

if [[ "$functional_enrichment" == "yes" ]]; then

Step4c_master_jobscript=${clusterFolder}/submission/Step4c_Functional_Enrichment_${protein}_${species}.sh

dataset_num=`cat $support | awk 'NR > 1{print $3}' | uniq | wc -l`

Step4c_jobID=Step4c_${protein}_${species}

echo "
#$ -S /bin/bash
#$ -l h_vmem=4G
#$ -l tmem=4G
#$ -l h_rt=36:00:00
#$ -pe smp 1
#$ -R y
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -N ${Step4c_jobID}
#$ -wd ${oFolder}
#$ -l h_rt=24:00:00
#$ -t 1-${dataset_num}
#$ -tc 20

if [[ \"\$SGE_TASK_ID\" == \"1\" ]];then
        echo \"
Step4c started at \`date +%H:%M:%S\`
\" >> $report_file
fi
jobs=\"" > $Step4c_master_jobscript

# do this on a per-cohort basis and also by condition
for dataset in `cat $support | awk 'NR > 1{print $3}' | uniq `;do
	echo $dataset
	SJAnalysisFolder=${results}/${dataset}/splice_junction_analysis
	if [ ! -e $SJAnalysisFolder ];then
		echo "cannot find splice junction analysis. Have you run Step 4b?"
	fi
	outFolder=${results}/${dataset}
	jobscript=${clusterFolder}/submission/functional_enrichment_${dataset}.sh
			echo "
#ENRICHMENT ANALYSIS FOR ${dataset}
" > $jobscript
#different datasets have different conditions (CTL vs HET, CTL vs HOM etc. Run enrichment analysis for each.)
	for i in `ls ${results}/${dataset}/strict_500/dexseq/`;do
		if [[ $i =~ .*[\.][a-z]+ ]];then
		echo "$i is not a valid set of conditions"
		continue;fi
		condition_names=`echo $i | awk -F"/" '{print $(NF-1)}'`
		echo $condition_names
		SJAnalysis_res=${SJAnalysisFolder}/${dataset}_${condition_names}_splicing_analysis.tab
        if [ ! -e ${SJAnalysis_res} ]; then
        	echo cannot find ${SJAnalysis_res}
        fi
		
		support_frame=${outFolder}/${dataset}_support_frame.tab
		error_file=${clusterFolder}/R/Functional_Enrichment_${dataset}.out

		awk -v dataset=$dataset 'NR == 1 {print $0} $3 == dataset {print $0}' $support > ${support_frame}

		echo "
echo $condition_names
${Rbin}script ${R_functional_enrichment} --support.frame ${support_frame} --code ${dataset} --species $species --condition.names ${condition_names} --outFolder ${outFolder}  > ${error_file} 2>&1
" >> $jobscript
			done

echo $jobscript >> $Step4c_master_jobscript

done

echo "\"
script=\`echo \$jobs | cut -f\$SGE_TASK_ID -d \" \"\`
sh \$script
" >> $Step4c_master_jobscript
echo "creating job scripts for splice junction analysis"
fi


#######################
### STEP 5: DESeq #####
#######################
# I want to see how gene expression (defined as numbers of reads falling into exons) correlates with cryptic exon presence.
# Step 3 creates counts that include the cryptic exons as well.
# Vincent's countsPrepareR script takes a folder full of counts and generates input files for DESeq2. I don't want to edit his script in particular.
# This will be carried out on each dataset in a cohort. There is no need for a "strict" parameter.

sanity_check=`awk 'NR == 1 {print $3}' $support `

if [ $DESeq = "yes" ] && [ $sanity_check = "dataset" ] ; then

Step5_master_jobscript=${clusterFolder}/submission/Step5_DESeq_${protein}_${species}.sh
dataset_num=`cat $support | awk 'NR > 1{print $3}' | uniq | wc -l`

Step5_jobID=Step5_${protein}_${species}

echo "
#$ -S /bin/bash
#$ -l h_vmem=8G
#$ -l tmem=8G
#$ -l h_rt=24:00:00
#$ -pe smp 4
#$ -R y
#$ -o ${clusterFolder}/out
#$ -e ${clusterFolder}/error
#$ -N ${Step5_jobID}
#$ -wd ${oFolder}
#$ -l h_rt=24:00:00
#$ -t 1-${dataset_num}
#$ -tc 20

if [[ \"\$SGE_TASK_ID\" == \"1\" ]];then
        echo \"
Step5 started at \`date +%H:%M:%S\`
\" >> $report_file
fi
jobs=\"" > $Step5_master_jobscript

# do this on a per-cohort basis
for dataset in `cat $support | awk 'NR > 1{print $3}' | uniq `
do

outFolder=${results}/${dataset}/expression
if [ ! -e $outFolder ]; then mkdir $outFolder; fi
support_frame=${outFolder}/${dataset}_deseq_frame.tab
jobscript=${clusterFolder}/submission/DESeq_${dataset}.sh
if [[ "$species" == "mouse" ]];then

GFF=/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Mouse/Mus_musculus.GRCm38.82_fixed.gff

elif [[ "$species" == "human" ]]; then
GFF=/cluster/scratch3/vyp-scratch2/reference_datasets/RNASeq/Human_hg38/Homo_sapiens.GRCh38.82_fixed.gff
fi
keepSex=TRUE
keepDups=FALSE
cryptic=TRUE

#Create the cryptic exon-excluded counts. 
# CountsPrepareR expects them to be like this: $iFolder/dexseq/sample_dexseq_counts.txt
# look for where the dexseq counts are
dexseq_counts=${results}/${dataset}/counts
if [ ! -e $dexseq_counts ]; then
dexseq_counts=${results}/${dataset}/strict_${strict_num}/counts
fi
#check this step hasn't already been attempted
new_countFolder=${outFolder}/dexseq/
if [ ! -e $new_countFolder ];then
mkdir $new_countFolder
fi

# slice up cohort support file into subcohort supports
awk -v dataset=$dataset 'BEGIN{OFS="\t"} NR == 1 {print $0} $3 == dataset {print $0}' $support > ${support_frame}
#removes dodgy columns full of NA values
Rscript ${R_support_chopper} ${support_frame}

echo "
#create counts for deseq from the cryptic exon counts - remove any lines with an "i"
for countfile in \`ls $dexseq_counts\`
        do awk '\$1 !~ /i/' ${dexseq_counts}/\${countfile} > $new_countFolder/\$countfile
done

#prepare counts
${Rbin}script ${countPrepareR} --gff ${GFF} --keep.dups ${keepDups} --support.frame ${support_frame} --code ${dataset} --annotation.file ${annotation_file} --iFolder ${outFolder}  > ${clusterFolder}/R/prepare_counts_${dataset}.out 2>&1

cat ${clusterFolder}/R/prepare_counts_${dataset}.out >> $report_file

#deseq
${Rbin}script ${deseqFinalProcessR} --keep.sex ${keepSex} --support.frame ${support_frame} --keep.dups ${keepDups} --code ${dataset} --annotation.file ${annotation_file} --iFolder ${outFolder} > ${clusterFolder}/R/deseq_${dataset}.out 2>&1

cat ${clusterFolder}/R/deseq_${dataset}.out >> $report_file
" > $jobscript

echo $jobscript >> $Step5_master_jobscript

done

echo "\"
script=\`echo \$jobs | cut -f\$SGE_TASK_ID -d \" \"\`
sh \$script
" >> $Step5_master_jobscript
echo "creating job scripts for DESeq testing"

fi

##############################
### FINAL STEP: SUBMISSION ###
##############################

# I want to run the pipeline script once and set all 4 steps running in series. This means that each job has to wait for the previous step to finish before executing.
#The 4 steps each have a variable assigned to them which is the -N option when submitting. I can use 'qsub -hold-jid $jobID' 
if [[ "$submit" == "yes" ]];then
	hold=""
	if [[ "$splice_extractor" == "yes" ]];then
		qsub $Step1_jobscript
		hold="-hold_jid $Step1_jobID"
	fi
	if [[ "$gff_creator" == "yes" ]];then
		if [[ "$hold_Step1" == "yes" ]];then
		Step1_jobID=Step1_${protein}_${species}			
		hold="-hold_jid $Step1_jobID"
		fi
		qsub $hold $Step2_jobscript
		hold="-hold_jid $Step2_jobID"
	fi
	if [[ "$read_counter" == "yes" ]]; then
		qsub $hold $Step3_master_jobscript
		hold="-hold_jid $Step3_jobID"
	fi
	if [[ "$DESeq" == "yes" ]]; then
            qsub $hold $Step5_master_jobscript
         	if [[ "$read_counter" == "yes" ]]; then
       		hold="-hold_jid $Step3_jobID";fi
	fi
	if [[ "$DEXSeq" == "yes" ]]; then
		qsub $hold $Step4_master_jobscript
		if [[ "$read_counter" == "yes" ]]; then
		hold="-hold_jid $Step3_jobID";fi
	fi
	if [[ "$splice_junction_analyzer" == "yes" ]];then
		if [[ "$DEXSeq" == "yes" ]]; then
			hold="-hold_jid $Step4_jobID"
		fi
		qsub $hold $Step4b_master_jobscript
		hold=""
	fi
	if [[ "$functional_enrichment" == "yes" ]];then
			if [[ "$splice_junction_analyzer" == "yes" ]];then
				hold="-hold_jid $Step4b_jobID"
			fi
		qsub $hold $Step4c_master_jobscript
		hold=""
	fi
	if [[ "$intron_retainer" == "yes" ]]; then
                qsub $Step6_master_jobscript
                hold="-hold_jid $Step6_jobID"
    fi
	if [[ "$intron_DEXSeq" == "yes" ]]; then
                qsub $hold $Step7_master_jobscript
        fi
fi
	


exit
