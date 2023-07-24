#!bin/bash

#### All the script was written take into account we are inside my RNAmapping conda env
### First all indexes were computed

### Variables: tune in a species-specific fashion

PLANTINA=/home/rocesv/
UTILS=/mnt/h/Taiotactical/Taiotactical_RNAseq/Utils/
HOME=/mnt/h/Taiotactical/Taiotactical_RNAseq/0.Taiotactical/0.At/
GENOME=/mnt/h/Taiotactical/Taiotactical_RNAseq/Utils/Genomes/0.At/At_genome.fa
GSIZE=91254070
INDEX=/mnt/h/Taiotactical/Taiotactical_RNAseq/Utils/Index/0.At/STAR/
INTRONM=11602
INTRONm=4
ITER=$(ls $HOME/0.SRA/ | cut -d'_' -f7 | sort -u -n)
echo $ITER

# loop through group of files that should be processed together

for i in $ITER
do

FILES=$(find $HOME/0.SRA/* -type f -iname "*_$i")
FILESN=$(ls $FILES | wc -l)

# Main code

# fastq extracion

cd $HOME/1.fastq
echo "STEP0: All iter ${i} files fastq-dumping ... | $(date)"
for f in $FILES
do
$UTILS/software_scripts/sratoolkit.2.11.3-ubuntu64/bin/fasterq-dump-orig.2.11.2 -p -e 6 -m 30G --split-3 $f
done

# If paired or single conditional
TYPE=$(find $HOME/1.fastq/* -type f -iname "*_${i}_1.fastq")
	if [ -z "$TYPE" ]; then
	echo "All iter ${i} files are SINGLE 100%"
	FASTQs=$(find $HOME/1.fastq/* -type f -iname "*_${i}.fastq")
	echo "STEP1: All iter ${i} files pre-processing ... | $(date)"
	for q in $FASTQs
	do
	# pre-process
	$UTILS/software_scripts/TrimGalore-0.6.6/trim_galore $q --fastqc -j 8 -o $HOME/2.trimmed/
	done
	# extract read length: careful this changes depending on SE or PE
	REPORTS=$(echo $(find $HOME/2.trimmed/* -type f -iname "*_${i}_trimmed_fastqc.zip") $(find $HOME/2.trimmed/* -type f -iname "*_${i}.fastq_trimming_report.txt"))
	echo "STEP2: All iter ${i} files multiqc ... | $(date)"
	multiqc $REPORTS -o $HOME/multiqc/ -n "All_${i}_files_multiqc"
	RLEN=$(awk 'FNR == 2 {print $12}' $HOME/multiqc/"All_${i}_files_multiqc_data/multiqc_fastqc.txt" | cut -d'-' -f2)
	RLEN=$(expr $RLEN - 1)
	# build BASENAME and mapping: careful this changes depending on SE or PE
	BASENAME=$(echo $REPORTS | cut -d'_' -f1 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f2 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f3 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f4 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f5 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f6 | sort -u)"_${i}"
	TRIMMED=$(find $HOME/2.trimmed/* -type f -iname "*_${i}_trimmed.fq")
	echo $TRIMMED > $HOME/dumpy_file_single.txt
	sed -e 's/\s\+/,/g' $HOME/dumpy_file_single.txt > $HOME/dumpy_file_single_comma.txt
	ReadFilesIn=$(cat $HOME/dumpy_file_single_comma.txt)
	echo "STEP3: All iter ${i} files mapping ... | $(date)"
	cd /home/rocesv/ # issue with fifo files while mapping in non Linux partitions and multiple fastq files
	STAR --genomeDir $INDEX --runThreadN 16 --genomeLoad LoadAndKeep --outFilterMismatchNmax 2 --alignIntronMin $INTRONm --alignIntronMax $INTRONM --outFileNamePrefix $HOME/3.mapping/$(basename $BASENAME) --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 16 --limitBAMsortRAM 60000000000 --sjdbOverhang $RLEN --readFilesIn $ReadFilesIn --outTmpDir /tmp/tmp/star
	# samtools idx .bam and .bam to bigwig: careful this changes depending on SE or PE
	echo "STEP4: All iter ${i} files formatting ... | $(date)"
	samtools index -b $HOME/3.mapping/$(basename $BASENAME)"Aligned.sortedByCoord.out.bam" -@ 16
	$PLANTINA/deepTools/bin/bamCoverage -b $HOME/3.mapping/$(basename $BASENAME)"Aligned.sortedByCoord.out.bam" -o $HOME/4.outputs/$(basename $BASENAME)".bw" --skipNonCoveredRegions -p 16 -v --normalizeUsing RPGC --effectiveGenomeSize $GSIZE
	# remove and clean up
	echo "STEP5: All iter ${i} files cleanning ... | $(date)"
	rm $HOME/1.fastq/*.fastq
	rm $HOME/2.trimmed/*.fq
	else
	TYPEN=$(ls $TYPE | wc -l) # number of paired versus total number to check mix single paired in i iter
		if [ $FILESN -eq $TYPEN ]; then
		echo "All iter ${i} files are PAIRED 100%"
		FASTQs=$(find $HOME/1.fastq/* -type f -iname "*_${i}_1.fastq")
		echo "STEP1: All iter ${i} files pre-processing ... | $(date)"
		for q in $FASTQs
		do
		# pre-process
		BASE=${q%_1.fastq}
		FORWARD="${BASE}_1.fastq"
		REVERSE="${BASE}_2.fastq"
		$UTILS/software_scripts/TrimGalore-0.6.6/trim_galore $FORWARD $REVERSE --paired --fastqc -j 8 -o $HOME/2.trimmed/
		done
		# extract read length: careful this changes depending on SE or PE
		REPORTS=$(echo $(find $HOME/2.trimmed/* -type f -iname "*_${i}_1_val_1_fastqc.zip") $(find $HOME/2.trimmed/* -type f -iname "*_${i}_2_val_2_fastqc.zip") $(find $HOME/2.trimmed/* -type f -iname "*_${i}_1.fastq_trimming_report.txt") $(find $HOME/2.trimmed/* -type f -iname "*_${i}_2.fastq_trimming_report.txt"))
		multiqc $REPORTS -o $HOME/multiqc/ -n "All_${i}_files_multiqc"
		RLEN=$(awk 'FNR == 2 {print $12}' $HOME/multiqc/"All_${i}_files_multiqc_data/multiqc_fastqc.txt" | cut -d'-' -f2)
		RLEN=$(expr $RLEN - 1)
		# build BASENAME and mapping: careful this changes depending on SE or PE
		BASENAME=$(echo $REPORTS | cut -d'_' -f1 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f2 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f3 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f4 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f5 | sort -u)"_"$(echo $REPORTS | cut -d'_' -f6 | sort -u)"_${i}"
		TRIMMEDf=$(find $HOME/2.trimmed/* -type f -iname "*_${i}_1_val_1.fq")
		TRIMMEDr=$(find $HOME/2.trimmed/* -type f -iname "*_${i}_2_val_2.fq")
		echo $TRIMMEDf > $HOME/dumpy_file_forward.txt
		echo $TRIMMEDr > $HOME/dumpy_file_reverse.txt
		sed -e 's/\s\+/,/g' $HOME/dumpy_file_forward.txt > $HOME/dumpy_file_forward_comma.txt
		sed -e 's/\s\+/,/g' $HOME/dumpy_file_reverse.txt > $HOME/dumpy_file_reverse_comma.txt
		ReadFilesInF=$(cat $HOME/dumpy_file_forward_comma.txt)
		ReadFilesInR=$(cat $HOME/dumpy_file_reverse_comma.txt)
		echo "STEP3: All iter ${i} files mapping ... | $(date)"
		cd /home/rocesv/ # issue with fifo files while mapping in non Linux partitions and multiple fastq files
		STAR --genomeDir $INDEX --runThreadN 16 --genomeLoad LoadAndKeep --outFilterMismatchNmax 2 --alignIntronMin $INTRONm --alignIntronMax $INTRONM --outFileNamePrefix $HOME/3.mapping/$(basename $BASENAME) --outSAMtype BAM SortedByCoordinate --outBAMsortingThreadN 16 --limitBAMsortRAM 60000000000 --sjdbOverhang $RLEN --readFilesIn $ReadFilesInF $ReadFilesInR --outTmpDir /tmp/tmp/star
		# samtools idx .bam and .bam to bigwig: careful this changes depending on SE or PE
		echo "STEP4: All iter ${i} files formatting ... | $(date)"
		samtools index -b $HOME/3.mapping/$(basename $BASENAME)"Aligned.sortedByCoord.out.bam" -@ 16
		$PLANTINA/deepTools/bin/bamCoverage -b $HOME/3.mapping/$(basename $BASENAME)"Aligned.sortedByCoord.out.bam"  -o $HOME/4.outputs/$(basename $BASENAME)".bw" --skipNonCoveredRegions -p 16 -v --normalizeUsing RPGC --effectiveGenomeSize $GSIZE --samFlagExclude 16
		# remove and clean up
		echo "STEP5: All iter ${i} files cleanning ... | $(date)"
		rm $HOME/1.fastq/*.fastq
		rm $HOME/2.trimmed/*.fq
		else 
		echo "Iter ${i} files are a MIX of paired and single: This should not be happening CHECK"
		exit 0
		fi
	fi

done
