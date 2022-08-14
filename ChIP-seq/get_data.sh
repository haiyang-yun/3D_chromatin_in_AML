#!/bin/bash
######################################
# 2013-05-13
# David Ruau, Department of Haematology
# Gottgens lab
# CIMR, University of Cambridge
# All right reserved.
# Licence: GPL (>=2)
######################################
# Example usage of running this script: get_data.sh -g [GENOTYPE] -m [OUTPUT_FOLDER] -i [INPUT_FASTQ] -x mm10
# -g [GENOTYPE]: top folder named as genotype (e.g. WT, FLT3, NPM1, DM)for input fastq data 
# -m [OUTPUT_FOLDER]: sub folder (e.g. H3K4me1, ATAC-seq) where fastq data are stored and serve as output folder to store mapped reads and QC report
######################################
# GLOBAL VARIABLES
######################################
# This is the location where the reference genome for bowtie2 are stored
resource='/serenity/data/reference-genomes/'
# and this for bowtie1, only considered to process human samples
resource_bw1='/serenity/data/reference-genomes/'
# The GNU parallel utility is not working on all platforms.
# set to no if you do not want to use parallel
use_parallel='yes'
# Sensitivity level for trimGalore
# Choose either WARN or FAIL (default)
trimLevel="FAIL"

usage(){
	echo 
    echo -e '\033[1m ChIP-seq /ATAC  pipeline\033[0m'
	echo 
	echo -e '\033[1mNAME\033[0m'
	echo '     get_data'
	echo 
	echo -e '\033[1mSYNOPSIS\033[0m'
	echo '     usage: get_data  [-g top_folder] [-m sub_folder] [-x ref_genome] [-i original fastq] [--merged] [--paired]  [--fasta]'
	echo 
    echo -e '\033[1mDESCRIPTION\033[0m'
	echo '     This utility will process raw sequence data from ChIP-seq/ATAC experiments in fastQ format.'
	echo '     When processing data, it will run FastQC quality check on each fastQ file,'
	echo '     trim the adpaters if present and run the aligner Bowtie2.'
	echo '     Human samples processed with the hg19 reference genome will be aligned with Bowtie 0.12.9.'
	echo '     It can process fasta file with --fasta options'
	echo '     You have to specify if the data are paired-end using the --paired option.'
	echo 
	echo '     The options are as follows:'
	echo 
    echo -e '     \033[1m-h, --help\033[0m'
    echo '             Show this message'
	echo 
    echo -e '     \033[1m-g \033[4mcommand\033[0m'
	echo '             The top folder to fastQ files'
	echo 
    echo -e '     \033[1m-m \033[4mcommand\033[0m'
    echo '             The sub folder to fastQ files'
	echo 
	echo -e '     \033[1m-i \033[4mcommand\033[0m'
	echo '             Sample Name'
	echo 
	echo -e '     \033[1m-x \033[4mcommand\033[0m'
	echo '             [optional] Reference genome: mm10, mm9, hg19, rn4 or rn5. If empty program will quit after fastQC'
	echo 
	echo -e '     \033[1m--paired\033[0m'
	echo '             [optional] Flag indicating that we are dealing with paired sequence reads.'
	echo
	echo -e '     \033[1m--fasta\033[0m'
	echo '             [optional] Flag indicating that we are processing fasta files. Fasta files need to be merged before hand using cat and placed in the GSE/GSM/ folder. The --fasta flag imply --local. It is used to specify that the "local" files are in fasta format and not in fastQ format.'
	echo 
	echo '             # To process local reads in fastQ or fasta formats.'
	echo '             # The files need to be in: top_folder/sub_folder/file.fastq or top_folder/sub_folder/file.fa'
	echo '             get_data  -g top_folder -m sub_folder -i original_fastq -x mm10'
	echo '             get_data  -g top_folder -m sub_folder -i original_fastq   -x mm10 --paired'		
	echo '             get_data --fasta -g top_folder -m sub_folder -i original_fastq  -x hg19'
	echo 
	echo  -e '\033[1mAUTHOR\033[0m'
	echo '     David Ruau <davidruau@gmail.com>'
	echo '     Department of Haematology, Gottgens lab'
	echo '     CIMR, University of Cambridge'
	echo '     Licence: GPL (>=2)'
	echo
	echo '02 Oct, 2013'
	exit 1
}

do_bowtie(){
	case "$GENOME" in
		mm10|mm9|hg19|rn4|rn5 )
			
			if [[ $GENOME == "hg19" ]] ; then
				echo '==> Running Bowtie 0.12.9...'
			else
				echo '==> Running Bowtie2...'
			fi
			
			REPORT="report_bowtie_"$GSM".txt"
			samFile=$GSM'_full.sam'
			
			if [[ $PAIRED == 0 ]]; then
				## ALIGNMENT TO REFERENCE GENOME STEP
				# For mouse dataset we use bowtie2 with the options described at
				# http://haemcode.stemcells.cam.ac.uk/doc_pipeline.php
				# -k 2 tell bowtie2 to keep reads with at most 2 match
				# -N 1 allow for one mismatch in the seed (20bp default).
				# the rest of the options is just for speed
				if [[ $FASTA == 1 ]];then
					echo '   - Processing FASTA sequence reads'
					if [[ $GENOME == "hg19" ]] ; then
						echo 'bowtie -m 2 -v 1 --best --strata --seed 0 '$resource_bw1$GENOME '-f *.fa -S '$samFile '2> '$REPORT
						bowtie -m 2 -v 1 --best --strata --seed 0 $resource_bw1$GENOME -f *.fa -S $samFile 2> $REPORT
					else
						echo 'bowtie2 -k 2 -N 1 --mm -x '$resource$GENOME '-f -U *.fa -S' $samFile '-p 20 2>' $REPORT
						bowtie2 -k 2 -N 1 --mm -x $resource$GENOME -f -U *.fa -S $samFile -p 20 2> $REPORT
					fi
										
				else
					echo '   - Processing '$GSM'.fastq'
					if [[ $GENOME == "hg19" ]] ; then
						echo 'bowtie -m 2 -v 1 --best --strata --seed 0'  $resource_bw1$GENOME '-q' $GSM'.fastq -S' $samFile '2>' $REPORT
						bowtie -m 2 -v 1 --best --strata --seed 0  $resource_bw1$GENOME -q $GSM.fastq -S $samFile 2> $REPORT
					else
						echo 'bowtie2 -k 2 -N 1 --mm -x' $resource$GENOME '-U' $GSM'.fastq' '-S' $samFile '-p 20 2> '$REPORT
						bowtie2 -k 2 -N 1 --mm -x $resource$GENOME -U $GSM'.fastq' -S $samFile -p 20  2> $REPORT
					fi
				fi
				## We remove the lines of the SAM file indicating the read matched multiple location. Keeping 
				# only uniquely mappable reads. The awk command filter for reads that actually map the forward 
				# strand (0) and the reverse strand (16) and also keep the header (@).
				grep -v XS:i: $samFile | awk '{if($2==0 || $2==16 || $1~/^@/) print $0}' > $GSM'.sam'
				
			else
				# IDENTIFY the fastQ files
				fileArray=(`ls *.fastq`)
				echo "   - Processing paired data ${fileArray[0]} and ${fileArray[1]}"
				if [[ $GENOME == "hg19" ]] ; then
					bowtie -m 2 -v 1 --best --strata --seed 0 $resource_bw1$GENOME -q -1 ${fileArray[0]} -2 ${fileArray[1]} -S $samFile 2> $REPORT
				else
					bowtie2 -k 2 -N 1 --mm --no-mixed --no-discordant  -x $resource$GENOME -1 ${fileArray[0]} -2 ${fileArray[1]} -S $samFile -p 20  2> $REPORT
				fi			
				## We remove the lines of the SAM file indicating the read matched multiple location. Keeping 
				# only uniquely mappable reads. The awk command filter for reads mate mapping uniquely the 
				# forward and reverse strand and also keep the headers.
				grep chr $samFile | grep -v "XS:i:" | awk '{if($2==147 || $2==83 || $2==99 || $2==163 || $2==81 || $2==97 || $2==145 || $2==161 || $1~/^@/) print $0}' > $GSM.sam
			fi
			rm $samFile
			if [[ $GENOME == "hg19" ]] ; then
				echo 
				echo -e "\033[1m\e[00;31m***\e[00m Human samples are processed with Bowtie 0.12.9!!! \e[00;31m***\e[00m\033[0m"
				echo 
				echo -e "\033[1m\e[00;31m***\e[00m Don't forget to copy the Bowtie report into the staging DB !!! \e[00;31m***\e[00m\033[0m"
				echo
				echo -e "   - Bowtie report in \e[00;31m$REPORT\e[00m"
			else
				echo 
				echo -e "\033[1m\e[00;31m***\e[00m Don't forget to copy the Bowtie2 report into the staging DB !!! \e[00;31m***\e[00m\033[0m"
				echo
				echo -e "   - Bowtie2 report in \e[00;31m$REPORT\e[00m"
			fi		
			;;
		* )
			echo "==> Unknown or No reference genome given. Exiting... Bowtie2 not run."
			;;
	esac
}

doQC(){
	FASTQ=$1
	PAIRED=$2
	GENOME=$3
	
	if [[ $PAIRED == 0 ]]
	then
		fastqc --quiet -f fastq $GSM'.fastq'
	else
		if [[ $use_parallel == "yes" ]]; then
			parallel --gnu 'fastqc --quiet -f fastq {}' ::: *.fastq
		else
			for filename in ./*.fastq; do 
				fastqc --quiet -f fastq "${filename}" 
			done
		fi
	fi
	
	# In case of multiple fastQ files for paired data you have to process the fastQC/trim_galore step differently.
	# In case of merged files you already ran fastQC and rtimgalore so no need for this step.
	
	if [[ $trimLevel == "WARN" ]]; then
		ERROR="$(grep -h -E $trimLevel$'\tOverrepresented sequences' *_fastqc/summary.txt | wc -l)"
		if [[ $ERROR == 0 ]]; then
			ERROR="$(grep -h -E FAIL$'\tOverrepresented sequences' *_fastqc/summary.txt | wc -l)"
		fi
	else
		ERROR="$(grep -h -E $trimLevel$'\tOverrepresented sequences' *_fastqc/summary.txt | wc -l)"
	fi
	
	if [[ $PAIRED == 0 ]]; then
		
		if [[ $ERROR != 0 ]]; then
			echo -e '\033[1m\e[00;31m   - fastQC FAIL   Overrepresented sequences\e[00m (adapters)\033[0m'
			echo  
			echo '==> Running trimGalore'
			
			# cleaning if get_data was interupted.
			FILE_TEST=`test -f *_trimmed.fq && echo 1`
			if [[ ! -z $FILE_TEST ]]; then
				rm -f *_trimmed.fq
			fi
			
			adapterArray=( $(awk '/#Sequence\tCount\tPercentage\tPossible Source/{flag=1;next} />>END_MODULE/{flag=0} flag {print}' */fastqc_data.txt | cut -f 1))
			
			## remove the unspecific sequence from array
			## those pattern remove all sequences from fastQ files
			# declare -a adapterArray=( ${adapterArray[@]/*NNNNNNNNNNNNNNNNNNNNNNNNN/} )			
			
			## Copy the original fastQ files for backup
			echo "==> Making copy of original fastQ file before trimming..."
			mkdir originals
			cp $GSM.fastq originals/"$GSM"_original.fq
			cp "$GSM"_fastqc.zip originals/"$GSM"_original-fastqc.zip
			
			for (( i=0; i<=$(( ${#adapterArray[@]} -1 )); i++ ))
			do
				echo -e "   - Removing adapter $(($i + 1)) out of ${#adapterArray[@]}"
				echo -e "        removing: ${adapterArray[$i]}"
				# trim_galore remove the adpater if the full lengh of the adpator is found (aka: stringency)
				trim_galore --no_report_file --quality 20 -a ${adapterArray[$i]} --stringency ${#adapterArray[$i]} $GSM.fastq &>> $GSM"_trimming_report.logFile"
				mv "$GSM"_trimmed.fq "$GSM".fastq
			done
			echo '   - trim_galore done'
			# erase old fastQC report
			rm -Rf "$GSM"_fastqc*
			# run fastQC on trimmed file
			echo '   - fastQC on trimmed fastQ file...'
			fastqc --quiet -f fastq $GSM'.fastq'
			
			# If reads are still too long
			# seqtk trimfq -b 1 -e 50 GSM1074872_trimmed.fq > GSM1074872_trimmed_1-50.fq
		fi
		echo '   - fastQC reported NO overrepresented sequences'
	else
		# We are dealing with non merged data but multiple fastQ file because paired reads.
		# Loop through the fastQC report and apply trim galore on both fastQ files
		# fastq file are named SRRxxx.fastq
		if [[ $ERROR != 0 ]]; then
			echo -e '\033[1m\e[00;31m   - fastQC FAIL   Overrepresented sequences\e[00m (adapters)\033[0m'
			echo 
			
			FAULTY_FASTQ=($(grep -h -E $trimLevel$'\tOverrepresented sequences' */summary.txt | perl -pe 's/.*\t(.*)\.fastq/\1/g'))
			
			for i in $( seq 0 $((${#FAULTY_FASTQ[@]} - 1)) ); do
				echo "==> Running trimGalore on ${FAULTY_FASTQ[$i]}"
				
				# cleaning if get_data was interupted.
				FILE_TEST=`test -f ${FAULTY_FASTQ[$i]}_trimmed.fq && echo 1`
				if [[ ! -z $FILE_TEST ]]; then
					rm -f *_trimmed.fq
				fi
				
				## get fastq file second reads
				# this is convoluted... trim_galore process paired end file together and one need to give
				# the name of the other file. Basically here we loop over all the adapter sequences for the
				# first file and then over the adpater seq for the second file. I do not use the -a2 option yet.
				index=$(echo ${FAULTY_FASTQ[$i]} | perl -pe 's/.*_(.*)/\1/g')
				fileroot=$(echo ${FAULTY_FASTQ[$i]} | perl -pe 's/(.*)_.*/\1/g')
				if [[ $index == 1 ]]; then
					OTHER_SRA_FILE=$fileroot"_2.fastq"
					theFile=$(echo ${FAULTY_FASTQ[$i]} | perl -pe 's/(.*)\.fastq/\1/g')
					otherFile=$(echo $OTHER_SRA_FILE | perl -pe 's/(.*)\.fastq/\1/g')
				else
					OTHER_SRA_FILE=$fileroot"_1.fastq"
					theFile=$(echo ${FAULTY_FASTQ[$i]} | perl -pe 's/(.*)\.fastq/\1/g')
					otherFile=$(echo $OTHER_SRA_FILE | perl -pe 's/(.*)\.fastq/\1/g')
				fi
				
				## Copy the original fastQ files for backup
				# Even if only one of the paired-end file contain adapter we treat both file at the same time.
				# We follow the recommendation of trimGalore. Because we trim both file_1 and file_2  we backup both.
				echo "==> Making copy of original fastQ file before trimming..."
				mkdir originals
				cp ${FAULTY_FASTQ[$i]}.fastq originals/${FAULTY_FASTQ[$i]}_original.fq
				cp ${FAULTY_FASTQ[$i]}_fastqc.zip originals/${FAULTY_FASTQ[$i]}_original-fastqc.zip
				cp $OTHER_SRA_FILE originals/$otherFile"_original.fq"
				cp $otherFile"_fastqc.zip" originals/$otherFile"_original-fastqc.zip"
									
				adapterArray=( $(awk '/#Sequence\tCount\tPercentage\tPossible Source/{flag=1;next} />>END_MODULE/{flag=0} flag {print}' ${FAULTY_FASTQ[$i]}"_fastqc"/fastqc_data.txt | cut -f 1))
				
				## remove the unspecific sequence from array
				# those pattern remove al sequences from fastQ files
				# declare -a adapterArray=( ${adapterArray[@]/*NNNNNNNNNNNNNNNNNNNNNNNNN/} )			
				
				for (( j=0; j<=$(( ${#adapterArray[@]} -1 )); j++ ))
				do
					echo -e "   - Removing adapter $(($j + 1)) out of ${#adapterArray[@]}"
					echo -e "        removing: ${adapterArray[$i]}"
					# trim_galore remove the adpater if the full lengh of the adpator is found (aka: stringency)			
					trim_galore --no_report_file --paired --quality 20 -a ${adapterArray[$j]} --stringency ${#adapterArray[$j]} ${FAULTY_FASTQ[$i]}.fastq $OTHER_SRA_FILE &>> $theFile"_trimming_report.logFile"
					
					## here we need the suffix of the *_val_* file
					index=$(echo ${FAULTY_FASTQ[$i]}_val_*.fq | perl -pe 's/.*_val_(.*).fq/\1/g')
					if [[ $index == 1 ]]; then
						mv "${FAULTY_FASTQ[$i]}"_val_1.fq "${FAULTY_FASTQ[$i]}".fastq
						mv $otherFile"_val_2.fq" $OTHER_SRA_FILE
					else
						mv "${FAULTY_FASTQ[$i]}"_val_2.fq "${FAULTY_FASTQ[$i]}".fastq
						mv $otherFile"_val_1.fq" $OTHER_SRA_FILE
					fi
				done
				echo '   - trim_galore done'
				
				# erase old fastQC report. Original fastqc report is left intact.
				rm -Rf *_fastqc*
				
			done # for faulty
			
			echo '   - fastQC on trimmed fastQ file...'
			if [[ $use_parallel == "yes" ]]; then
				parallel --gnu 'fastqc --quiet -f fastq {}' ::: *.fastq
			else
				for filename in $PWD/*.fastq; do 
					fastqc --quiet -f fastq $filename
				done
			fi
		else
			echo '   - fastQC reported NO overrepresented sequences'
		fi # if ERROR
		
	fi # IF PAIRED
	
	## get the PHRED info from fastQC report
	encoding=($(grep -h -E "^Encoding(.*)" */fastqc_data.txt | sort | uniq | perl -pe 's/Encoding\s+(.*)/\1/g'))
	# transform array to string
	encoding=${encoding[@]}
	echo "   - fastQ format: "$encoding"<<"
	if [[ $encoding == "Sanger / Illumina 1.9" ]] || [[ $encoding == "Sanger" ]]; then
		echo '   - Quality encoding: phred33'
		QUALITY="--phred33"
	elif [[ $encoding == "Illumina 1.5" ]] || [[ $encoding == "Illumina <1.3" ]] || [[ $encoding == "Illumina 1.3" ]]; then
		echo '   - Quality encoding: phred64'
		QUALITY="--phred64"
	else
		echo "\033[1m\e[00;31m***\e[00m fastQ format encoding not recognised! Report fastQC report to developer. Assuming default --phred33 \e[00;31m***\e[00m\033[0m"
		QUALITY="--phred33"
		# QUALITY="--solexa-quals"
	fi
}

GSE=
GSM=
GENOME=
FILE=
FASTQ=
PAIRED=0
FASTA=0

# Note that the : after an option flag means that it should have a value instead of
# just being the boolean flag that a is.
# OPTS=`getopt -o hg:m:s:x: --long help,trimmed -- "$@"`
OPTS=`getopt -o hg:m:x:i: --long help,paired,fasta -- "$@"`
if [ $? != 0 ]
then
	# something went wrong, getopt will put out an error message for us
    exit 1
fi

eval set -- "$OPTS"

while true
do
	case "$1" in
        -h | --help)
            usage
			;;
		# for options with required arguments, an additional shift is required
        -g)
            GSE=$2
			shift 2;;
        -m)
            GSM=$2
			shift 2;;
		-x)
			GENOME=$2
			shift 2;;
		-i) 
			FILE=$2
			shift 2;;
		--paired)
			PAIRED=1
			shift;;
		--fasta)
			FASTA=1
			shift;;
		--) break;;
		--*) break;;
        -?)
            usage
            ;;
    esac
done


if [[ -z $GSE ]] || [[ -z $GSM ]] || [[ -z $FILE ]]|| [[ -z $GENOME ]]; then
	   usage
else
	if [[ $PAIRED == 0 ]]; then

		mkdir -p $GSE/$GSM
		mv $FILE $GSE/$GSM/$GSM".fastq"	
		cd $GSE/$GSM
	else
		mkdir -p $GSE/$GSM
		mv $FILE".r_1.fastq" $GSE/$GSM/$GSM".r_1.fastq" 
		mv $FILE".r_2.fastq" $GSE/$GSM/$GSM".r_2.fastq"
		cd $GSE/$GSM
	fi
	##############
	## fastQC
	##############
	echo 
	echo '==> fastQC on fastQ file...'
	doQC $GSM $PAIRED $GENOME

	##############
	## BOWTIE2
	##############
	do_bowtie
fi	


echo 
echo '==> Summary:'
echo 
cd ../..
tree -shFL 1 $GSE/$GSM

echo 
echo -e "\033[1m\e[00;31mdone.\e[00m\033[0m"
