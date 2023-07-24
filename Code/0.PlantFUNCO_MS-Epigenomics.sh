#!bin/bash

#### All the script was written take into account we are inside plantina-chiplike docker container
### First activate conda and change to desire conda env
### Next modify vdb-config -i to uncompress sra files in current directory nad paste .ncbi in the desired folder
### Set sra-toolkit into PATH

### Variables

PLANTINA=/home/Plantina_Utils/
UTILS=/home/Transfer/Utils/
WD=/home/Transfer/0.Taiotactical/
log=/home/Transfer/0.Taiotactical/Taiotactical_log.txt
SPECIES=$(ls $WD)

### Master Script

for sp in $SPECIES # general loop
	do
	echo "Processing $sp directoy ... | $(date)"

	## conditional Arabidopsis thaliana

	if [ "${sp}" = "0.Zm" ]; then

	# Species-specific variables needed

	GENOME=$UTILS/Genomes/$sp/*_genome.fa
	FSIZE=$UTILS/Genomes/$sp/*_genome.sizes
	GSIZE=$(cat $UTILS/gEffectiveSize/$sp/gsize.txt)
	INDEX=$UTILS/Index/$sp/bowtie2/$(echo $sp | cut -d'.' -f2)"_genome"
	HOME=$WD/$sp/hiHMM_marks/
	ITER=$(ls $HOME/0.SRA/ | cut -d'_' -f7 | sort -u -n)
	echo $ITER
	# loop thorugh group of files that should be processed together

	for i in $ITER
	do
        # Variables needed for the commands

	FILES=$(find $HOME/0.SRA/* -type f -iname "*_$i")
	FILESN=$(ls $FILES | wc -l)

	# Main code

	# fastq extraction
	cd $HOME/1.fastq/
	for f in $FILES
	do
	fasterq-dump -p -e 6 -m 10G --split-3 $f
	done

	# If paired or single conditional
	TYPE=$(find $HOME/1.fastq/* -type f -iname "*_${i}_1.fastq")
	if [ -z "$TYPE" ]; then
	echo "All iter ${i} files are SINGLE 100%"
	INPUT=$(find $HOME/1.fastq/* -type f -iname "*_Input_${i}.fastq")
	IP=$(find $HOME/1.fastq/* -type f -iname "*_IP_${i}.fastq")
	OPEN=$(find $HOME/1.fastq/* -type f -iname "*__${i}.fastq")
	FASTQs=$(echo $INPUT $IP $OPEN)
	# pre-process: for start
	for q in $FASTQs # pre-process for start
	do
	# pre-process
	$PLANTINA/TrimGalore-0.6.6/trim_galore $q --fastqc -j 8 -o $HOME/2.trimmed/
	# mapping
	BASE=${q%.fastq}
	BASENAME=$(basename $BASE)
	TRIMMED=$(find $HOME/2.trimmed/* -type f -iname "${BASENAME}_trimmed.fq")
	bowtie2 --very-sensitive -N 1 -p 16 -x $INDEX -U $TRIMMED -S $HOME/3.sbam/"${BASENAME}.sam" 
	# convert to bam and filter bam
	samtools view --threads 16 -bS $HOME/3.sbam/"${BASENAME}.sam" -o $HOME/3.sbam/"${BASENAME}.bam"
	samtools sort $HOME/3.sbam/"${BASENAME}.bam" --threads 16 -o $HOME/3.sbam/"${BASENAME}.bam"
	java -jar $PLANTINA/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$HOME/3.sbam/"${BASENAME}.bam" O=$HOME/3.sbam/"${BASENAME}.rmdup.bam" M=$HOME/3.sbam/"${BASENAME}.rmdup.bam.txt" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
	samtools view -b $HOME/3.sbam/"${BASENAME}.rmdup.bam" -q 30 -o $HOME/3.sbam/"${BASENAME}.unique.bam"
	# filtered bigwig and bed from filtered bam - raw signals: for end
	samtools index -b $HOME/3.sbam/"${BASENAME}.unique.bam"
	$PLANTINA/deepTools/bin/bamCoverage -b $HOME/3.sbam/"${BASENAME}.unique.bam" -o $HOME/4.outputs/"${BASENAME}.unique.bw" -of bigwig -p 12
	$PLANTINA/deepTools/bin/bamCoverage -b $HOME/3.sbam/"${BASENAME}.unique.bam" -o $HOME/4.outputs/"${BASENAME}.unique.bg" -of bedgraph -p 12
	done
	# If epigenomic mark
	MARKS=$(echo $FASTQs | cut -d'_' -f4)
	if [ "${MARKS}" = "5mC" ]; then
		# If input availability and several files
		IPp=$(find $HOME/3.sbam/* -type f -iname "*_IP_${i}.unique.bam")
		INPUTp=$(find $HOME/3.sbam/* -type f -iname "*_Input_${i}.unique.bam")
		if [ -z "$INPUTp" ]; then
		echo "5mC INPUT not available"
		# peak calling
		NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
		macs2 callpeak -t $IPp -f BAM -n "${NAMEp}_Single_noInput_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		macs2 callpeak -t $IPp -f BAM -n "${NAMEp}_Single_noInput_broad" --SPMR -B --broad --nolambda --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		# signalt tracks: FE and log10FE + 0.1
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_FE.bdg"
		# always pref narrow if not broad
		# pass to bg-bw
		cd $PLANTINA
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_FE.bdg" $FSIZE
		mv ./*.bw $HOME/6.SignalTracks/
		else
		echo "5mC INPUT available"
                # peak calling
		NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
		macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_broad" --SPMR -B --broad --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                # signal tracks: FE and log10FE + 0.1
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_FE.bdg"
		# always pref narrow if not broad
		# pass to bg-bw
		cd $PLANTINA
                ./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_logFE05.bdg" $FSIZES
                ./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_logFE05.bdg" $FSIZES
                ./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_FE.bdg" $FSIZES
                ./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_FE.bdg" $FSIZES
                mv ./*.bw $HOME/6.SignalTracks/
		fi
	elif [ "${MARKS}" = "ATAC-seq" ] || [ "${MARKS}" = "DNase-seq" ]; then
		# If IP/Input availability or just IP and several files (no Input ones has __ instead IP_ or Input_)
		IPp=$(find $HOME/3.sbam/* -type f -iname "*_IP_${i}.unique.bam")
		INPUTp=$(find $HOME/3.sbam/* -type f -iname "*_Input_${i}.unique.bam")
		OPENp=$(find $HOME/3.sbam/* -type f -iname "*__${i}.unique.bam")
		if [ -z "$IPp" ]; then
		echo "Chromatin acc. IP/INPUT not available, just IPs"
		# peak calling
		NAMEp=$(basename $(echo $OPENp | cut -d'_' -f1,2,3,4,5,8))
		macs2 callpeak -t $OPENp -f BAM -n "${NAMEp}_Single_noInput_narrow" --SPMR -B --keep-dup all --nomodel --shift 75 --extsize 150 -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		macs2 callpeak -t $OPENp -f BAM -n "${NAMEp}_Single_noInput_broad" --SPMR -B --broad --nolambda --keep-dup all --nomodel --shift 75 --extsize 150 -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		# signal tracks: FE and log10FE + 0.1
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_FE.bdg"
		# always pref narrow if not broad
		# pass to bg-bw
		cd $PLANTINA
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_narrow_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_broad_signaltrack_FE.bdg" $FSIZE
		mv ./*.bw $HOME/6.SignalTracks/
		else
		echo "Chromatin acc. IP/INPUT available"
		# peak calling
		NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
		macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_narrow" --SPMR -B --keep-dup all --nomodel --shift 75 --extsize 150 -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_broad" --SPMR -B --broad --keep-dup all --nomodel --shift 75 --extsize 150 -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		# signal tracks: FE and log10FE + 0.1
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_FE.bdg"
		# always pref narrow if not broad
		# pass to bg-bw
		cd $PLANTINA
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_narrow_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_broad_signaltrack_FE.bdg" $FSIZE
		mv ./*.bw $HOME/6.SignalTracks/
		fi
	else
                # If input availability and several files
                IPp=$(find $HOME/3.sbam/* -type f -iname "*_IP_${i}.unique.bam")
                INPUTp=$(find $HOME/3.sbam/* -type f -iname "*_Input_${i}.unique.bam")
                if [ -z "$INPUTp" ]; then
                echo "Histone mod. INPUT not available"
                # peak calling
		NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
		macs2 callpeak -t $IPp -f BAM -n "${NAMEp}_Single_noInput_nomodel_narrow" --SPMR -B --keep-dup all --nomodel -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		macs2 callpeak -t $IPp -f BAM -n "${NAMEp}_Single_noInput_nomodel_broad" --SPMR -B --broad --nolambda --keep-dup all --nomodel -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		macs2 callpeak -t $IPp -f BAM -n "${NAMEp}_Single_noInput_model_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
		macs2 callpeak -t $IPp -f BAM -n "${NAMEp}_Single_noInput_model_broad" --SPMR -B --broad --nolambda --keep-dup all  -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                # signalt tracks: FE and log10FE + 0.1
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_nomodel_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_broad_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_model_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_narrow_control_lambda.bdg" -m FE -o $HOME//6.SignalTracks/"${NAMEp}_Single_noInput_model_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_model_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_noInput_model_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_model_broad_signaltrack_FE.bdg"
		# pass to bg-bw
		cd $PLANTINA
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_narrow_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_narrow_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_broad_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_nomodel_broad_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_model_narrow_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_model_narrow_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_model_broad_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_noInput_model_broad_signaltrack_FE.bdg" $FSIZE
		mv ./*.bw $HOME/6.SignalTracks/
                else
                echo "Histone mod. INPUT available"
                # peak calling
		NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
                macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_model_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_model_broad" --SPMR -B --broad --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_nomodel_narrow" --SPMR -B --keep-dup all --nomodel -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                macs2 callpeak -t $IPp -c $INPUTp -f BAM -n "${NAMEp}_Single_yesInput_nomodel_broad" --SPMR -B --broad --keep-dup all --nomodel -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                # signal tracks: FE and log10FE + 0.1
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_nomodel_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_broad_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_narrow_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_narrow_signaltrack_FE.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_broad_signaltrack_logFE05.bdg"
		macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Single_yesInput_model_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_broad_signaltrack_FE.bdg"
		# pass to bg-bw
		cd $PLANTINA
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_narrow_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_narrow_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_broad_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_nomodel_broad_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_narrow_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_narrow_signaltrack_FE.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_broad_signaltrack_logFE05.bdg" $FSIZE
		./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Single_yesInput_model_broad_signaltrack_FE.bdg" $FSIZE
		mv ./*.bw $HOME/6.SignalTracks/
                fi
	fi
	else
	TYPEN=$(ls $TYPE | wc -l) # number of paired versus total number to check mix single paired in i iter
		if [ $FILESN -eq $TYPEN ]; then
		echo "All iter ${i} files are PAIRED 100%"
		FASTQs=$(find $HOME/1.fastq/* -type f -iname "*_${i}_1.fastq")
		# pre-process: paired for start
		for q in $FASTQs # pre-process for start
		do
		BASE=${q%_1.fastq}
		FORWARD="${BASE}_1.fastq"
		REVERSE="${BASE}_2.fastq"
		$PLANTINA/TrimGalore-0.6.6/trim_galore $FORWARD $REVERSE --paired --fastqc -j 8 -o $HOME/2.trimmed/
		# cat and mapping
		BASENAME=$(basename $BASE)
		FORWARDt=$(find $HOME/2.trimmed/* -type f -iname "${BASENAME}_1_val_1.fq")
		REVERSEt=$(find $HOME/2.trimmed/* -type f -iname "${BASENAME}_2_val_2.fq")
		bowtie2 --very-sensitive -N 1 -p 16 -x $INDEX -1 $FORWARDt -2 $REVERSEt -S $HOME/3.sbam/"${BASENAME}.sam"
		# convert to bam and filter bam
		samtools view --threads 16 -bS $HOME/3.sbam/"${BASENAME}.sam" -o $HOME/3.sbam/"${BASENAME}.bam"
       		samtools sort $HOME/3.sbam/"${BASENAME}.bam" --threads 16 -o $HOME/3.sbam/"${BASENAME}.bam"
       		java -jar $PLANTINA/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$HOME/3.sbam/"${BASENAME}.bam" O=$HOME/3.sbam/"${BASENAME}.rmdup.bam" M=$HOME/3.sbam/"${BASENAME}.rmdup.bam.txt" USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
		samtools view -b $HOME/3.sbam/"${BASENAME}.rmdup.bam" -q 30 -o $HOME/3.sbam/"${BASENAME}.unique.bam"
		# filtered bigwig and bed from filtered bam - raw signals: paired for end
		samtools index -b $HOME/3.sbam/"${BASENAME}.unique.bam"
		$PLANTINA/deepTools/bin/bamCoverage -b $HOME/3.sbam/"${BASENAME}.unique.bam" -o $HOME/4.outputs/"${BASENAME}.unique.bw" -of bigwig -p 12
		$PLANTINA/deepTools/bin/bamCoverage -b $HOME/3.sbam/"${BASENAME}.unique.bam" -o $HOME/4.outputs/"${BASENAME}.unique.bg" -of bedgraph -p 12
		done
		# If epigenomic mark
		MARKS=$(echo $FASTQs | cut -d'_' -f4)
	        if [ "${MARKS}" = "5mC" ]; then
			IPp=$(find $HOME/3.sbam/* -type f -iname "*_IP_${i}.unique.bam")
                	INPUTp=$(find $HOME/3.sbam/* -type f -iname "*_Input_${i}.unique.bam")
                	if [ -z "$INPUTp" ]; then
                	echo "5mC INPUT not available"
                	# peak calling
                	NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
                	macs2 callpeak -t $IPp -f BAMPE -n "${NAMEp}_Paired_noInput_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
			macs2 callpeak -t $IPp -f BAMPE -n "${NAMEp}_Paired_noInput_broad" --SPMR -B --broad --nolambda --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
			# signalt tracks: FE and log10FE + 0.1
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_FE.bdg"
                	macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_FE.bdg"
			# always pref narrow if not broad
			# pass to bg-bw
			cd $PLANTINA
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_FE.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_FE.bdg" $FSIZE
			mv ./*.bw $HOME/6.SignalTracks/
                	else
                	echo "5mC INPUT available"
                	# peak calling
			NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
			macs2 callpeak -t $IPp -c $INPUTp -f BAMPE -n "${NAMEp}_Paired_yesInput_narrow" --SPMR -B --keep-dup all  -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
			macs2 callpeak -t $IPp -c $INPUTp -f BAMPE -n "${NAMEp}_Paired_yesInput_broad" --SPMR -B --broad --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/ 
                	# signal tracks: FE and log10FE + 0.1
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_FE.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_FE.bdg"
			# always pref narrow if not broad
			# pass to bg-bw
                        cd $PLANTINA
                        ./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_logFE05.bdg" $FSIZE
                        ./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_FE.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_FE.bdg" $FSIZE
                        mv ./*.bw $HOME/6.SignalTracks/
                	fi
        	elif [ "${MARKS}" = "ATAC-seq" ] || [ "${MARKS}" = "DNase-seq" ]; then
			# If IP/Input availability or just IP and several files (no Input ones has __ instead IP_ or Input_)
                	IPp=$(find $HOME/3.sbam/* -type f -iname "*_IP_${i}.unique.bam")
                	INPUTp=$(find $HOME/3.sbam/* -type f -iname "*_Input_${i}.unique.bam")
			OPENp=$(find $HOME/3.sbam/* -type f -iname "*__${i}.unique.bam")
                	if [ -z "$IPp" ]; then
                	echo "Chromatin acc. IP/INPUT not available, just IPs"
                	# peak calling
			NAMEp=$(basename $(echo $OPENp | cut -d'_' -f1,2,3,4,5,8))
                	macs2 callpeak -t $OPENp -f BAMPE -n "${NAMEp}_Paired_noInput_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
			macs2 callpeak -t $OPENp -f BAMPE -n "${NAMEp}_Paired_noInput_broad" --SPMR -B --broad --nolambda --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                	# signal tracks: FE and log10FE + 0.1
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_FE.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_FE.bdg"
			# always pref narrow if not broad
			# pass to bg-bw
			cd $PLANTINA
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_FE.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_FE.bdg" $FSIZE
			mv ./*.bw $HOME/6.SignalTracks/
                	else
                	echo "Chromatin acc. IP/INPUT available"
                	# peak calling
			NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
			macs2 callpeak -t $IPp -c $INPUTp -f BAMPE -n "${NAMEp}_Paired_yesInput_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
			macs2 callpeak -t $IPp -c $INPUTp -f BAMPE -n "${NAMEp}_Paired_yesInput_broad" --SPMR -B --broad --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                	# signal tracks: FE and log10FE + 0.1
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_FE.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_FE.bdg"
			# always pref narrow if not broad
			# pass to bg-bw
			cd $PLANTINA
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_FE.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_FE.bdg" $FSIZE
			mv ./*.bw $HOME/6.SignalTracks/
                	fi
        	else
			# If input availability and several files
                	IPp=$(find $HOME/3.sbam/* -type f -iname "*_IP_${i}.unique.bam")
                	INPUTp=$(find $HOME/3.sbam/* -type f -iname "*_Input_${i}.unique.bam")
                	if [ -z "$INPUTp" ]; then
                	echo "Histone mod. INPUT not available"
                	# peak calling
			NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
                	macs2 callpeak -t $IPp -f BAMPE -n "${NAMEp}_Paired_noInput_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
			macs2 callpeak -t $IPp -f BAMPE -n "${NAMEp}_Paired_noInput_broad" --SPMR -B --broad --nolambda --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                	# signalt tracks: FE and log10FE + 0.1
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_FE.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_noInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_FE.bdg"
			# pass to bg-bw
			cd $PLANTINA
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_narrow_signaltrack_FE.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_noInput_broad_signaltrack_FE.bdg" $FSIZE
			mv ./*.bw $HOME/6.SignalTracks/
                	else
                	echo "Histone mod. INPUT available"
                	# peak calling
			NAMEp=$(basename $(echo $IPp | cut -d'_' -f1,2,3,4,5,8))
			macs2 callpeak -t $IPp -c $INPUTp -f BAMPE -n "${NAMEp}_Paired_yesInput_narrow" --SPMR -B --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
			macs2 callpeak -t $IPp -c $INPUTp -f BAMPE -n "${NAMEp}_Paired_yesInput_broad" --SPMR -B --broad --keep-dup all -g $GSIZE --outdir $HOME/5.Peaks/$NAMEp/
                	# signal tracks: FE and log10FE + 0.1
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_narrow_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_FE.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_control_lambda.bdg" -p 0.1 -m logFE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_logFE05.bdg"
			macs2 bdgcmp -t $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_treat_pileup.bdg" -c $HOME/5.Peaks/$NAMEp/"${NAMEp}_Paired_yesInput_broad_control_lambda.bdg" -m FE -o $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_FE.bdg"
			# pass to bg-bw
			cd $PLANTINA
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_narrow_signaltrack_FE.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_logFE05.bdg" $FSIZE
			./bdg2bw_2 $HOME/6.SignalTracks/"${NAMEp}_Paired_yesInput_broad_signaltrack_FE.bdg" $FSIZE
			mv ./*.bw $HOME/6.SignalTracks/
                	fi
        	fi
		else
		echo "MIX: This not should be happening CHECK"
		exit 0
		fi
	fi

	# fill log with sp, ITER, mark, SRR, Paired or Single, Reps, Input Avail

	if [ -z "$TYPE" ]; then
		if [ "${MARKS}" = "5mC" ]; then
			if [ -z "$INPUTp" ]; then
			echo "$sp $i $MARK Single noInput" >> $log
			else
			echo "$sp $i $MARK Single yesInput" >> $log
			fi
		elif [ "${MARKS}" = "ATAC-seq" ] || [ "${MARKS}" = "DNase-seq" ]; then
			if [ -z "$IPp" ]; then
			echo "$sp $i $MARK Single noInput" >> $log
			else
			echo "$sp $i $MARK Single yesInput" >> $log
			fi
		else
			if [ -z "$INPUTp" ]; then
			echo "$sp $i $MARK Single noInput" >> $log
			else
			echo "$sp $i $MARK Single yesInput" >> $log
			fi
		fi
	else
		if [ "${MARKS}" = "5mC" ]; then
			if [ -z "$INPUTp" ]; then
			echo "$sp $i $MARK Paired noInput" >> $log
			else
			echo "$sp $i $MARK Paired yesInput" >> $log
			fi
		elif [ "${MARKS}" = "ATAC-seq" ] || [ "${MARKS}" = "DNase-seq" ]; then
			if [ -z "$IPp" ]; then
			echo "$sp $i $MARK Paired noInput" >> $log
			else
			echo "$sp $i $MARK Paired yesInput" >> $log
			fi
		else
			if [ -z "$INPUTp" ]; then
			echo "$sp $i $MARK Paired noInput" >> $log
			else
			echo "$sp $i $MARK Paired yesInput" >> $log
			fi
		fi
	fi

	# clear and remove: fastqs, trimmed fastqs, sams, non filtered bams, empty files

	mv $HOME/3.sbam/*.unique.bam $HOME/4.outputs/
	mv $HOME/3.sbam/*.rmdup.bam.txt $HOME/4.outputs/
	mv $HOME/3.sbam/*.unique.bam.bai $HOME/4.outputs/
	rm $HOME/1.fastq/*
	rm $HOME/2.trimmed/*.fq
	rm $HOME/3.sbam/*.sam
	rm $HOME/3.sbam/*.rmdup.bam
	rm $HOME/3.sbam/*.bam
	rm $(find /home/Transfer/0.Taiotactical/0.At/hiHMM_marks/6.SignalTracks/* -size 1)

	done # iter loop

	# multiqc all i iter files

        multiqc $HOME/2.trimmed/* -o $HOME/2.trimmed/ -n "Alliter_${sp}_prebamfilter_multiqc"

	## hiHMM input: signal tracks FE binned by 200 bp: all
	# For each epigenomic mark merge logFoldChange BedGraphs: log2FE + pseudocount 0.5 binned by each 200 bp and averaged | from directly bedgraph outpot and signal tracks log10+0.1
	MODELfbdg10=$(find $HOME/6.SignalTracks/* -type f -iname "*_model_*_signaltrack_logFE05.bdg")
	MODELfbw10=$(find $HOME/6.SignalTracks/* -type f -iname "*_model_*_signaltrack_logFE05.bw")
	MODELfbdg=$(find $HOME/6.SignalTracks/* -type f -iname "*_model_*_signaltrack_FE.bdg")
        MODELfbw=$(find $HOME/6.SignalTracks/* -type f -iname "*_model_*_signaltrack_FE.bw")
	mkdir $HOME/6.SignalTracks_model/
	mv $MODELfbdg10 $HOME/6.SignalTracks_model/
	mv $MODELfbw10 $HOME/6.SignalTracks_model/
        mv $MODELfbdg $HOME/6.SignalTracks_model/
        mv $MODELfbw $HOME/6.SignalTracks_model/
	MARKSf=$(ls $HOME/6.SignalTracks/ | cut -d'_' -f3 | sort -u)
	bedtools makewindows -g $FSIZE -w 200 > $HOME/7.hiHMMinputs/"${sp}_Genome_windows_unsorted.bed"
	LC_COLLATE=C sort -k1,1 -k2,2n $HOME/7.hiHMMinputs/"${sp}_Genome_windows_unsorted.bed" > $HOME/7.hiHMMinputs/"${sp}_Genome_windows.bed"
	echo $MARKSf
	for k in $MARKSf
		do
                FILESf=$(find $HOME/6.SignalTracks* -type f -iname "*_${k}_*")
		if [ "${k}" = "5mC" ]; then
		# 5mC: narrow (signal track)
		fivemClog=$(find $FILESf -type f -iname "*_narrow_signaltrack_logFE05.bdg")
		fivemCFE=$(find $FILESf -type f -iname "*_narrow_signaltrack_FE.bdg")
		# 0. Merge bed graphs from same epigenomic mark | bedtools unionbedg and awk mean | add conditional for single files no unionbedg
		Llog=$(ls $fivemClog | wc -l)
		Lfe=$(ls $fivemCFE | wc -l)
			if [ "$Llog" -eq "1" ]; then
			cp $fivemClog $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg"
			else
			echo "Several Signal Tracks Available computing mean ..."
			bedtools unionbedg -i $fivemClog | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg"
			fi
			if [ "$Lfe" -eq "1" ]; then
                        cp $fivemCFE $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg"
                        else
			echo "Several Signal Tracks Available computing mean ..."
			bedtools unionbedg -i $fivemCFE | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg"
			fi
		# log2(x + 0.5) FE
		awk '{print $1"\t"$2"\t"$3"\t"log($4 +0.5)/log(2)}' $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean.bdg"
		# 1. Map bedgrpah scores to 200 bp bin
		LC_COLLATE=C sort -k1,1 -k2,2n $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_sorted.bdg"
		LC_COLLATE=C sort -k1,1 -k2,2n $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_sorted.bdg"
		bedtools map -a $HOME/7.hiHMMinputs/"${sp}_Genome_windows.bed" -b $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_sorted.bdg" -c 4 -o mean -null -1 > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_binned.bdg"
		bedtools map -a $HOME/7.hiHMMinputs/"${sp}_Genome_windows.bed" -b $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_sorted.bdg" -c 4 -o mean -null -1 > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_binned.bdg"
		elif [ "${MARKS}" = "ATAC-seq" ] || [ "${MARKS}" = "DNase-seq" ]; then
		# Open: narrow and broad
                Openlog=$(find $FILESf -type f -iname "*_signaltrack_logFE05.bdg")
                OpenFE=$(find $FILESf -type f -iname "*_signaltrack_FE.bdg")
                # 0. Merge bed graphs from same epigenomic mark | bedtools unionbedg and awk mean | single file conditional
                Llog=$(ls $Openlog | wc -l)
                Lfe=$(ls $OpenFE | wc -l)
			if [ "$Llog" -eq "1" ]; then
			cp $Openlog $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg"
			else
			echo "Several Signal Tracks Available computing mean ..."
			bedtools unionbedg -i $Openlog | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg"
			fi
			if [ "$Lfe" -eq "1" ]; then
			cp $OpenFE $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg"
			else
			echo "Several Signal Tracks Available computing mean ..."
			bedtools unionbedg -i $OpenFE | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg"
			fi
                # log2(x + 0.5) FE
		awk '{print $1"\t"$2"\t"$3"\t"log($4 +0.5)/log(2)}' $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean.bdg"
                # 1. Map bedgrpah scores to 200 bp bin
		LC_COLLATE=C sort -k1,1 -k2,2n $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_sorted.bdg"
		LC_COLLATE=C sort -k1,1 -k2,2n $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_sorted.bdg"
		bedtools map -a $HOME/7.hiHMMinputs/"${sp}_Genome_windows.bed" -b $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_sorted.bdg" -c 4 -o mean -null -1 > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_binned.bdg"
		bedtools map -a $HOME/7.hiHMMinputs/"${sp}_Genome_windows.bed" -b $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_sorted.bdg" -c 4 -o mean -null -1 > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_binned.bdg"
		else
		# Histones: narrow and broad (or from .*Peaks narrow broad and gapped nomodel)
                Histoneslog=$(find $FILESf -type f -iname "*_signaltrack_logFE05.bdg")
                HistonesFE=$(find $FILESf -type f -iname "*_signaltrack_FE.bdg")
                # 0. Merge bed graphs from same epigenomic mark | bedtools unionbedg and awk mean | single file coniditional
                Llog=$(ls $Histoneslog | wc -l)
                Lfe=$(ls $HistonesFE | wc -l)
			if [ "$Llog" -eq "1" ]; then
			cp $Histoneslog $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg"
			else
			echo "Several Signal Tracks Available computing mean ..."
			bedtools unionbedg -i $Histoneslog | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg"
			fi
			if [ "$Lfe" -eq "1" ]; then
			cp $HistonesFE $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg"
			else
			echo "Several Signal Tracks Available computing mean ..."
			bedtools unionbedg -i $HistonesFE | awk '{sum=0; for (col=4; col<=NF; col++) sum += $col; print $1"\t"$2"\t"$3"\t"sum/(NF-4+1); }' > $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg"
			fi
                # log2(x + 0.5) FE
		awk '{print $1"\t"$2"\t"$3"\t"log($4 +0.5)/log(2)}' $HOME/7.hiHMMinputs/"${sp}_${k}_all_FE_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean.bdg"
                # 1. Map bedgrpah scores to 200 bp bin
		LC_COLLATE=C sort -k1,1 -k2,2n $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_sorted.bdg"
		LC_COLLATE=C sort -k1,1 -k2,2n $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean.bdg" > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_sorted.bdg"
		bedtools map -a $HOME/7.hiHMMinputs/"${sp}_Genome_windows.bed" -b $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_sorted.bdg" -c 4 -o mean -null -1 > $HOME/7.hiHMMinputs/"${sp}_${k}_all_log2FE05_mean_binned.bdg"
		bedtools map -a $HOME/7.hiHMMinputs/"${sp}_Genome_windows.bed" -b $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_sorted.bdg" -c 4 -o mean -null -1 > $HOME/7.hiHMMinputs/"${sp}_${k}_all_logFE05_mean_binned.bdg"
		fi
		done

	# 2. Concat all marks with unionbdg
	log2FE05=$(find $HOME/7.hiHMMinputs/* -type f -iname "*_log2FE05_mean_binned.bdg")
	logFE05=$(find $HOME/7.hiHMMinputs/* -type f -iname "*_logFE05_mean_binned.bdg")
	bedtools unionbedg -i $log2FE05 -header -names $MARKSf > $HOME/7.hiHMMinputs/"${sp}_all_log2FE05_hiHMMinput.bdg"
	bedtools unionbedg -i $logFE05 -header -names $MARKSf > $HOME/7.hiHMMinputs/"${sp}_all_logFE05_hiHMMinput.bdg"

	## hiHMM input: empirical mappability genome: all
	# merge all Input processed bams and with bamCoverage to bedGraph extending it
	samtools merge -O BAM --threads 16 $HOME/7.hiHMMinputs/Mappability_files.bam $(find $HOME/4.outputs/* -type f -iname "*.unique.bam")
	samtools index -b $HOME/7.hiHMMinputs/Mappability_files.bam
	$PLANTINA/deepTools/bin/bamCoverage -b $HOME/7.hiHMMinputs/Mappability_files.bam -o $HOME/7.hiHMMinputs/Merged_inputs.bdg -of bedgraph -p 16 --extendReads --skipNonCoveredRegions
	bedtools merge -i $HOME/7.hiHMMinputs/Merged_inputs.bdg > $HOME/7.hiHMMinputs/"mappable_${sp}.bed"
	rm $HOME/7.hiHMMinputs/Mappability_files.bam
	rm $HOME/7.hiHMMinputs/Merged_inputs.bdg

	## hiHMM input: finish signal tracks
	# Mean centering and scaling to 1 standard deviation in R with scale() and split by structural unit

	fi # sp conditional

	# conditional Oryza sativa 
	if [ "${sp}" = "0.At" ]; then
	echo $sp
	fi

	# condtional Zea mays
	if [ "${sp}" = "0.Zm" ]; then
	echo $sp
	fi

	done # sp loop

