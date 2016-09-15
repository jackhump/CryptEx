

outFolder="/Users/Jack/Documents/Cryptic_Exons/IGV_supplementary"

for species in human mouse;do

	if [[ "$species" == "human" ]];then
		data=(total mRNA)
	elif [[ "$species" == "mouse" ]];then
		data=(AB ES)
	fi

	if [ ! -e ${outFolder}/${species} ];then
		mkdir -p ${outFolder}/${species}
	fi

	for dataset in ${data[@]};do

		batchFile=${outFolder}/${species}_${dataset}_IGV_batch.txt

		BEDtable="/Users/Jack/project/Humphrey_RNASeq_brain/jack_git/Humphrey_RNASeq_brain/splice_junction_detection/extended_hunting/Figures/Union_Datasets/summary_tables/${species}_IGV_table.bed"

		echo snapshotDirectory ${outFolder}/${species} > $batchFile

		awk -v species_id=${species}_${dataset} '{
			print "goto "$1" "$2" "$3;
			print "snapshot  " NR"_"$4"_"$5"_"species_id".png";
			print "" }' $BEDtable >> $batchFile
	done

done


# run this after you've taken all the screenshots with IGV. 

# crop the top and bottom off of each screenshot
for species in human mouse; do
	mkdir -p ${outFolder}/${species}/cropped/
	for picture in `ls ${outFolder}/${species}`;do
		# crop top
		convert -crop +0+130 ${outFolder}/${species}/${picture} ${outFolder}/${species}/cropped/crop_${picture}
		# crop bottom
		convert -crop +0-410 ${outFolder}/${species}/cropped/crop_${picture} ${outFolder}/${species}/cropped/${picture} 
		# remove intermediate file
		rm ${outFolder}/${species}/cropped/crop_${picture} 
	done
done


# assemble all the pictures in a TeX file

for species in human mouse; do
	if [[ "$species" == "human" ]];then
		end_Num=95
		data=(total mRNA)
	elif [[ "$species" == "mouse" ]];then
		endNum=52
		data=(AB ES)
		data_long=("Mouse Adult Brain" "Mouse Embryonic Stem Cell")
	fi
	TeX=${outFolder}/${species}_all_images.tex

	echo "
\documentclass{article}

\usepackage{graphicx}
\usepackage{float}
\usepackage{caption}
\usepackage{fullpage}

\title{Mouse Cryptic Exons \\ from \\ \textit{Quantitative analysis of cryptic splicing associated with TDP-43 depletion}}
\author{Jack Humphrey, Warren Emmett, Pietro Fratta, Adrian M. Isaacs \& Vincent Plagnol}
\begin{document}
\maketitle
\newpage

" > $TeX
	endNum=2
	for i in `seq 1 $endNum`; do
		pic_one=`ls ${outFolder}/${species}/cropped/${i}_*_${data[0]}*png`
		pic_two=`ls ${outFolder}/${species}/cropped/${i}_*_${data[1]}*png`
		gene_name=`echo $pic_one | awk -F'/' '{split($NF,a,"_");print a[2]}' `

		echo "
\centering
\section{$gene_name}

\begin{figure}[H]
	\centering
	\caption*{${data_long[0]}}
	\includegraphics[width=\textwidth]{$pic_one}
	\newline
	\newline
	\caption*{${data_long[1]}}
	\includegraphics[width=\textwidth]{$pic_two}
\end{figure}

" >> $TeX
	 done

	 echo "
\end{document}
" >> $TeX
pdflatex $TeX
done