#Step 2 Master Script
# this controls the maximum distance between pairs of novel junctions to be merged into a cryptic tag
mergeNum=500
minSplice=5

#Case Statement to test arguments
until [ -z $1 ];do
        case $1 in
        --dataset)
	shift
	dataset=$1;;
	--resFolder)
	shift
	resFolder=$1;;
	--output)
        shift
        output=$1;;
	--spliced_beds)
	shift
	spliced_beds=$1;;
	--exon_GFF)
	shift
	exon_GFF=$1;;
	--intron_BED)
	shift
	intron_BED=$1;;
        -* )
        echo "unrecognised argument: $1"
        exit 1;;
esac
shift
if [ "$#" = "0" ]; then break; fi
        echo $1 $2
done

echo $output


## Concatenate all spliced intron bed files together and sort by start and coordinate
##
cat ${spliced_beds}/*spliced.introns.bed | sort -S 50% -k1,1 -k2,2n > ${output}.sorted.bed

#Merge cryptic tags if they are within 500bp of each other.
bedtools merge -i ${output}.sorted.bed -d ${mergeNum} -c 1 -o count > ${output}.overlap.merged.bed
# remove any tags with less than a set number of counts
awk -v MINSPLICE=$minSplice '$4 >= MINSPLICE {print } ' ${output}.overlap.merged.bed > ${output}.overlap.merged.reduced.bed
bedtools subtract -a ${output}.overlap.merged.reduced.bed -b ${exon_GFF} > ${output}.merged.bed

## Intersect again. this removes intergenic spliced intervals and adds information about the intron that it intersects.
bedtools intersect -a ${output}.merged.bed -b ${intron_BED} -wb | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$8,$9}' | sort -k1,1V -k5,5n > ${output}.cryptics.merged.bed

# add exon number
if [ -e ${output}.merged.annotated.bed ]; then rm ${output}.merged.annotated.bed;fi

## create unique list of gene/strand/intron_number and grep the merged list for each instance of that combination
## use awk to then append unique numbers on to each one
# for testing purposes just take the first 1000
cat ${output}.cryptics.merged.bed | awk '{print $6}' | sort -V | uniq > ${output}.unique_gene_introns.tab

# N is the number of allowed concurrently running forks
N=8

for entry in `cat ${output}.unique_gene_introns.tab`;do 
	((i=i%N)); ((i++==0)) && wait
	grep $entry ${output}.cryptics.merged.bed | awk 'BEGIN{s=1}{print $0"i"s;s+=1}' >> ${output}.merged.annotated.bed &
done

## Convert into a GFF file
cat ${output}.merged.annotated.bed | awk 'BEGIN{OFS="\t"}{split($6,a,"_");print $1, "mouse_iGenomes_GRCm38_with_ensembl.gtf", "exonic_part", $2, $3, ".", $5, ".", "transcripts \"cryptic_exon\"; exonic_part_number \""a[2]"\"; gene_id \""a[1]"\"" }' |
                sort -k1,1 -k2,2n > ${output}.cryptics.gff

outFile=`basename $output`

## Place cryptic exons within the total exon GFF
cat ${output}.cryptics.gff ${exon_GFF}| sort -k1,1V -k4,4n -k5,5n | awk '$14 ~ /ENS/' > ${resFolder}/${outFile}.total.cryptics.gff

# could we get rid of any genes that don't show any cryptic splicing?
# write list of genes that contain cryptic tags
# in R:

Rscript gff_reducer.R --output ${output} --resFolder ${resFolder} --exon_GFF ${exon_GFF} --outFile ${outFile}

# test whether the reduced GFF is as good as the the full!


exit
