#!/bin/bash

pushd `dirname $0` > /dev/null
MEA_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

# MANUAL INSTALLATION
#source $MEA_DIR_TOOLS/mea.config
# DOCKER INSTALLATION
source /mea-data/mea.config

##############################################################################
#############   Module 3: allele-specific alignment
##############################################################################

PARAM_VALID=1
PARAM_SINGLE_READS=1
if [ $MEA_USE_CONCATENATED_GENOME = 1 ]; then
    if [ "$1" = "-s" ]; then
        if [ $# -eq 6 ]; then
            PARAM_FASTQ_FILE=$2
            PARAM_GENOME=$3
            PARAM_REFERENCE_GENOME=$MEA_DIR_REFERENCES
            PARAM_STRAIN1=$4
            PARAM_STRAIN2=$5
            PARAM_BAM_PREFIX=$6
        else
            PARAM_VALID=0
        fi
    elif [ "$1" = "-p" ]; then
        if [ $# -eq 7 ]; then
            PARAM_SINGLE_READS=0
            PARAM_FASTQ_FILE1=$2
            PARAM_FASTQ_FILE2=$3
            PARAM_GENOME=$4
            PARAM_REFERENCE_GENOME=$MEA_DIR_REFERENCES
            PARAM_STRAIN1=$5
            PARAM_STRAIN2=$6
            PARAM_BAM_PREFIX=$7
        else
            PARAM_VALID=0
        fi
    else
        PARAM_VALID=0
    fi
fi

if [ $PARAM_VALID = 0 ]; then
    echo "
Usage:
        $MEA alignReads <-s/-p> input_reads_1 [input_reads_2] genome_concat strain1 strain2 outputPrefix

Options:
    -s              to align single-end reads (requires one input file)
    -p              to align paired-end reads (requires two input files)
                    
    input_reads_1   the 1st input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)
                    
    input_reads_2   (paired end) the 2nd input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)
                    
    genome_concat   path to the indexed reference for concatenated insilico genome.
                    for BWA, specifiy path to the fasta.
                    for Bowtie2 and Tophat2, specify path and basename of index files
                    for Bismark, specify genome folder, excluding <Bisulphite_Genome>
                    
    strain1         name of strain1
                    (e.g. hap1 or CASTEiJ)
                    
    strain2         name of strain2
                    (e.g. hap2 or C57BL6J)
                    
    outputPrefix    prefix for output files, including the full path, without an extension
                    (e.g. ./TSC_H3K36me3 )

Output:
    
    outputPrefix_strain1_strain2.sam   (when MEA_USE_CONCATENATED_GENOME=1)
                                       all reads aligned to the concatenated insilico genome
                                       
    outputPrefix_strain1_all.sam       (when MEA_USE_CONCATENATED_GENOME=0)
                                       all reads aligned to the first insilico genome
                                       
    outputPrefix_strain2_all.sam       (when MEA_USE_CONCATENATED_GENOME=0)
                                       all reads aligned to the second insilico genome
                                       
    outputPrefix_strain1.bam           allelic reads for strain1 genome (sorted bam)
                                       
    outputPrefix_strain2.bam           allelic reads for strain2 genome (sorted bam)

Examples:
    (MEA_USE_BWA=1)
    $MEA alignReads -s H3K36me3.fastq createGenomeDir/hap1_hap2.fasta hap1 hap2 ./H3K36me3
    $MEA alignReads -p H3K36me3_1.fastq H3K36me3_2.fastq createGenomeDir/hap1_hap2.fasta hap1 hap2 ./hap1_hap2_H3K36me3
    
    (MEA_USE_STAR=1)
    $MEA alignReads -s H3K36me3.fastq createGenomeDir/STAR-index/hap1_hap2/ hap1 hap2 ./hap1_hap2_H3K36me3
    
    (MEA_USE_BOWTIE2=1 or MEA_USE_TOPHAT2=1)
    $MEA alignReads -s H3K36me3.fastq createGenomeDir/bowtie2-index/hap1_hap2 hap1 hap2 ./hap1_hap2_H3K36me3

    (MEA_USE_BISMARK=1)
    $MEA alignreads -p WGBS_1.fastq WGBS_2.fastq createGenomeDir/hap1_hap2/ hap1 hap2 ./hap1_hap2_WGBS
"
exit 1
fi

#------------------------------------------------------------------------------------
# Aligns single-end reads using bwa and creates a bam file sorted by position.
#------------------------------------------------------------------------------------
function alignReadsBWA_S {
    local PARAM_FASTQ_FILE=$1
    local PARAM_FASTA_FILE=$2
    local PARAM_BAM_PREFIX=$3
    printProgress "[alignReadsBWA_S] Started"
    
    local PARAM_IS_BAM=""
    if [[ "$PARAM_FASTQ_FILE" == *.bam ]]
    then
        PARAM_IS_BAM="-b"
    fi
    
    # align to reference genome
    $MEA_BIN_BWA aln $MEA_BWA_ALN_PARAMS $PARAM_IS_BAM $PARAM_FASTA_FILE $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX".sai
    
    # generate alignment (single ended)
    $MEA_BIN_BWA samse $PARAM_FASTA_FILE "$PARAM_BAM_PREFIX".sai $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX".sam
    
    # convert SAM to BAM
    #$MEA_BIN_SAMTOOLS view -b -S -o "$PARAM_BAM_PREFIX".bam "$PARAM_BAM_PREFIX".sam
    
    if [ $MEA_DEBUG = 0 ]; then
        rm -f "$PARAM_BAM_PREFIX".sai
    fi
    
    # remove temporary files
    if [ -f  "$PARAM_BAM_PREFIX".sam ]
    then
        printProgress "[alignReadsBWA] SAM file created successfully."
    else
        printProgress "[alignReadsBWA] ERROR: SAM file was not created."
    fi
    
    printProgress "[alignReadsBWA_S] Done"
}

#------------------------------------------------------------------------------------
# Aligns paired-end reads using bwa and creates a bam file sorted by position.
#------------------------------------------------------------------------------------
function alignReadsBWA_P {
    local PARAM_FASTQ_FILE1=$1
    local PARAM_FASTQ_FILE2=$2
    local PARAM_FASTA_FILE=$3
    local PARAM_BAM_PREFIX=$4
    
    printProgress "[alignReadsBWA_P] Started"
    
    local PARAM_IS_BAM1=""
    if [[ "$PARAM_FASTQ_FILE1" == *.bam ]]
    then
        PARAM_IS_BAM1="-b1"
    fi
    
    local PARAM_IS_BAM2=""
    if [[ "$PARAM_FASTQ_FILE2" == *.bam ]]
    then
        PARAM_IS_BAM2="-b2"
    fi
    
    # align to reference genome
    $MEA_BIN_BWA aln $MEA_BWA_ALN_PARAMS $PARAM_IS_BAM1 $PARAM_FASTA_FILE $PARAM_FASTQ_FILE1 > "$PARAM_BAM_PREFIX"_1.sai
    $MEA_BIN_BWA aln $MEA_BWA_ALN_PARAMS $PARAM_IS_BAM2 $PARAM_FASTA_FILE $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_2.sai
    
    # generate alignment (paired-ended)
    $MEA_BIN_BWA sampe $PARAM_FASTA_FILE "$PARAM_BAM_PREFIX"_1.sai "$PARAM_BAM_PREFIX"_2.sai $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX".sam
    
    # convert SAM to BAM and sort
    #$MEA_BIN_SAMTOOLS view -b -S -o "$PARAM_BAM_PREFIX".bam "$PARAM_BAM_PREFIX".sam
    
    if [ $MEA_DEBUG = 0 ]; then
        rm -f "$PARAM_BAM_PREFIX"_1.sai
        rm -f "$PARAM_BAM_PREFIX"_2.sai
    fi
    
    if [ -f  "$PARAM_BAM_PREFIX".sam ]
    then
        printProgress "[alignReadsBWA] SAM file created successfully."
    else
        printProgress "[alignReadsBWA] ERROR: SAM file was not created."
    fi
    
    printProgress "[alignReadsBWA_P] Done"
}

#------------------------------------------------------------------------------------
# detect the allelic reads from the reads aligned to the concatenated genome
#------------------------------------------------------------------------------------
function detectAllelicConcatenated {
    
    printProgress "[detectAllelicConcatenated] Started"
    
    local PARAM_INPUT_SAM=$1
    local PARAM_STRAIN=$2
    local PARAM_REFFASTA=$3
    local PARAM_QUALITY=$4
    local PARAM_OUT_PREFIX=$5
    
    # output header first
    $MEA_BIN_SAMTOOLS view -SH "$PARAM_INPUT_SAM" \
        | awk -v ref="$PARAM_STRAIN" '($0 ~ ref) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_//g' \
        > "$PARAM_OUT_PREFIX".sam
    
    # append reads
    $MEA_BIN_SAMTOOLS view -S $PARAM_INPUT_SAM \
        | awk -v ref="$PARAM_STRAIN" '(($3 ~ ref)&&($5'"$PARAM_QUALITY"')) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_//g' \
        >> "$PARAM_OUT_PREFIX".sam
    
    # convert to bam
    $MEA_BIN_SAMTOOLS view -bt $PARAM_REFFASTA "$PARAM_OUT_PREFIX".sam > "$PARAM_OUT_PREFIX".unsorted.bam
    
    # sort by coordinates
    $MEA_BIN_SAMTOOLS sort "$PARAM_OUT_PREFIX".unsorted.bam "$PARAM_OUT_PREFIX" 
    
    if [ -f "$PARAM_OUT_PREFIX".bam ]
    then
        printProgress "[detectAllelicConcatenated] Filtered BAM file created."
        # remove temp files
        if [ $MEA_DEBUG = 0 ]; then
            printProgress "[detectAllelicConcatenated] Removing temporary SAM files."
            rm -f "$PARAM_OUT_PREFIX".sam
            rm -f "$PARAM_OUT_PREFIX".unsorted.bam
        fi
    else
        printProgress "[detectAllelicConcatenated] ERROR: Filtered BAM file was not created."
    fi
    printProgress "[detectAllelicConcatenated] Done"
}

#------------------------------------------------------------------------------------
# detect the allelic cytosine information in the concatenated genome method
#------------------------------------------------------------------------------------
function detectAllelicCytoConcatenated {
    
    printProgress "Started detectAllelicCytoConcatenated"
    
    local PARAM_INPUT_METHYL=$1
    local PARAM_STRAIN=$2
    local PARAM_OUT_PREFIX=$3
    
    meaCheckFileExists "$PARAM_INPUT_METHYL"
    
    cat "$PARAM_INPUT_METHYL" \
        | awk -v ref="$PARAM_STRAIN" '($0 ~ ref) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_//g' \
        >> "$PARAM_OUT_PREFIX".CpG_report.txt
    
    printProgress "Finished detectAllelicCytoConcatenated"
}

#------------------------------------------------------------------------------------
# alignReads
#------------------------------------------------------------------------------------

if [ $MEA_USE_CONCATENATED_GENOME = 1 ]; then
    
    if [ $PARAM_SINGLE_READS = 1 ]; then
        #align single-end reads to concatenated insilico genome
        meaCheckFileExists "$PARAM_FASTQ_FILE"
        
        if [ $MEA_USE_BWA = 1 ]; then
            meaCheckFileExists "$PARAM_GENOME"
            printProgress "aligning to concatenated insilico genome"

#			replace this
#			alignReadsBWA_S "$PARAM_FASTQ_FILE" "$PARAM_GENOME" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
#			with this
			$MEA_BIN_BWA aln $MEA_BWA_ALN_PARAMS "$PARAM_GENOME" $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX".sai
			$MEA_BIN_BWA samse "$PARAM_GENOME" "$PARAM_BAM_PREFIX".sai $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
#			rm "$PARAM_BAM_PREFIX".sai
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".bam

            meaCheckFileExists "$MEA_REFERENCE_GENOME"
            printProgress "aligning to reference genome"
            $MEA_BIN_BWA aln "$MEA_BWA_ALN_TOTAL_PARAMS" $MEA_REFERENCE_GENOME "$PARAM_FASTQ_FILE" > "$PARAM_BAM_PREFIX"_total.sai
            $MEA_BIN_BWA samse $MEA_REFERENCE_GENOME "$PARAM_BAM_PREFIX"_total.sai $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_total.sam
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
        
        elif [ $MEA_USE_BOWTIE2 = 1 ]; then
            meaCheckFileExists "$PARAM_GENOME".1.bt2*
            printProgress "align to insilico concatenated genome"
            $MEA_BIN_BOWTIE2 "$MEA_BOWTIE2_ALN_PARAMS" -x "$PARAM_GENOME" $PARAM_FASTQ_FILE -S "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
            printProgress "align to reference genome"
            meaCheckFileExists "$MEA_DIR_REFERENCES"/bowtie2-index/"$MEA_BUILD".1.bt2*
            $MEA_BIN_BOWTIE2 "$MEA_BOWTIE2_ALN_TOTAL_PARAMS" -x "$MEA_DIR_REFERENCES"/bowtie2-index/"$MEA_BUILD" $PARAM_FASTQ_FILE -S "$PARAM_BAM_PREFIX"_total.sam
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
        
        elif [ $MEA_USE_BISMARK = 1 ]; then
            VAR_BASENAME=`basename $PARAM_BAM_PREFIX`
            VAR_DIRNAME=`dirname $PARAM_BAM_PREFIX`
            VAR_q=$MEA_BAM2WIG_PARAM_MIN_QUALITY
            VAR_F=$MEA_BAM2WIG_PARAM_FILTERING_FLAG
            printProgress "align to the insilico concatenated genome"
            meaCheckDirExists "$PARAM_GENOME"/Bisulfite_Genome
            $MEA_BIN_BISMARK $MEA_BISMARK_ALN_PARAMS --basename "$VAR_BASENAME"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2" -o $VAR_DIRNAME $PARAM_GENOME $PARAM_FASTQ_FILE
            printProgress "aligning to reference genome"
            meaCheckDirExists "$MEA_DIR_REFERENCES"/"$MEA_BUILD"/Bisulfite_Genome
            $MEA_BIN_BISMARK $MEA_BISMARK_ALN_TOTAL_PARAMS --basename "$VAR_BASENAME"_total -o $VAR_DIRNAME "$MEA_DIR_REFERENCES"/"$MEA_BUILD" $PARAM_FASTQ_FILE
            $MEA_BIN_SAMTOOLS view -bShu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" ">= 0" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" ">= 0" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            printProgress "extract methylation of the insilico concatenated genome"
            $MEA_BIN_SAMTOOLS view -Sbh -F $VAR_F -q $VAR_q "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam \
                > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".bam
            $MEA_BIN_BISMARK_EXTRACT -s --comprehensive --cytosine_report --samtools_path $MEA_DIR_TOOLS -o $VAR_DIRNAME --genome_folder $PARAM_GENOME \
                "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".bam
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt \
                "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            sort -k1,1 -k2,2n "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt \
                > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt
            if [ $MEA_DEBUG = 0 ]; then
                rm "$VAR_DIRNAME"/C??_context_*.txt
                rm "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            fi
            printProgress "extract methylation of reference genome"
            $MEA_BIN_SAMTOOLS view -Sbh -F $VAR_F -q $VAR_q "$PARAM_BAM_PREFIX"_total.sam > "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".bam
            $MEA_BIN_BISMARK_EXTRACT -s --comprehensive --cytosine_report --samtools_path $MEA_DIR_TOOLS -o $VAR_DIRNAME --genome_folder $MEA_DIR_REFERENCES/$MEA_BUILD \
                "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".bam
            mv "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.txt "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            sort -k1,1 -k2,2n "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt > "$PARAM_BAM_PREFIX"_total.CpG_report.txt
            if [ $MEA_DEBUG = 0 ]; then
                rm "$VAR_DIRNAME"/C??_context_*.txt
                rm "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            fi
            printProgress "detecting allelic cytosines information"
            detectAllelicCytoConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt \
                $PARAM_STRAIN1 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicCytoConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt \
                $PARAM_STRAIN2 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1".CpG_report.txt "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_preProject.CpG_report.txt
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2".CpG_report.txt "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_preProject.CpG_report.txt
        
        elif [ $MEA_USE_STAR = 1 ]; then
            printProgress "align to concatenated insilico genome"
            meaCheckDirExists "$PARAM_GENOME"
            if [ "$PARAM_FASTQ_FILE" == "*.gz" ]
            then
              $MEA_BIN_STAR --runMode alignReads --genomeDir "$PARAM_GENOME" "$MEA_STAR_ALN_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE" --readFilesCommand gunzip -c
            else 
              $MEA_BIN_STAR --runMode alignReads --genomeDir "$PARAM_GENOME" "$MEA_STAR_ALN_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE"
            fi
            mv Aligned.out.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
            mv Log.out "$PARAM_BAM_PREFIX"_STAR_RunParameters.tsv
            mv Log.final.out "$PARAM_BAM_PREFIX"_STAR_AlignmentSummary.tsv
            printProgress "align to reference genome"
            if [ "$PARAM_FASTQ_FILE" == "*.gz" ]
            then
              $MEA_BIN_STAR --runMode alignReads --genomeDir "$MEA_DIR_REFERENCES"/STAR-index/"$MEA_BUILD"/ "$MEA_STAR_ALN_TOTAL_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE" --readFilesCommand gunzip -c
            else 
              $MEA_BIN_STAR --runMode alignReads --genomeDir "$MEA_DIR_REFERENCES"/STAR-index/"$MEA_BUILD"/ "$MEA_STAR_ALN_TOTAL_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE"
            fi
            meaCheckDirExists "$MEA_DIR_REFERENCES"/STAR-index/"$MEA_BUILD"
            mv Aligned.out.sam "$PARAM_BAM_PREFIX"_total.sam
            $MEA_BIN_SAMTOOLS view -bShu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            mv Log.out "$PARAM_BAM_PREFIX"_total_STAR_referenceRunParameters.tsv
            mv Log.final.out "$PARAM_BAM_PREFIX"_total_STAR_referenceAlignmentSummary.tsv
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            if [ $MEA_DEBUG = 0 ]; then
                rm SJ.out.tab Log.progress.out
            fi
        
        elif [ $MEA_USE_TOPHAT2 = 1 ]; then
            printProgress "align to insilico concatenated genome"
            meaCheckFileExists "$PARAM_GENOME".1.bt2*
            $MEA_BIN_TOPHAT2 $MEA_TOPHAT2_ALN_PARAMS --output-dir ./tophat_out_concat "$PARAM_GENOME" $PARAM_FASTQ_FILE
            mv ./tophat_out_concat/accepted_hits.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
            mv ./tophat_out_concat/align_summary.txt ./"$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_alignSummary.txt
            printProgress "align to reference genome"
            meaCheckFileExists "$MEA_DIR_REFERENCES"/bowtie2-index/"$MEA_BUILD".1.bt2*
            $MEA_BIN_TOPHAT2 $MEA_TOPHAT2_ALN_TOTAL_PARAMS --output-dir ./tophat_out_total "$MEA_DIR_REFERENCES"/bowtie2-index/"$MEA_BUILD" $PARAM_FASTQ_FILE
            mv ./tophat_out_total/accepted_hits.sam "$PARAM_BAM_PREFIX"_total.sam
            mv ./tophat_out_total/align_summary.txt ./"$PARAM_BAM_PREFIX"_total_alignSummary.txt
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            if [ $MEA_DEBUG = 0 ]; then
                rm -r tophat_out_*
            fi
        
        fi
    
    else #[ $PARAM_SINGLE_READS = 0 ]
        #align paired-end reads to concatenated insilico genome
        meaCheckFileExists $PARAM_FASTQ_FILE1
        meaCheckFileExists $PARAM_FASTQ_FILE2
        
        if [ $MEA_USE_BWA = 1 ]; then
            meaCheckFileExists $PARAM_GENOME
            printProgress "aligning reads to concatenated insilico genome"
			$MEA_BIN_BWA aln $MEA_BWA_ALN_PARAMS "$PARAM_GENOME" $PARAM_FASTQ_FILE1 > "$PARAM_BAM_PREFIX"_1.sai
			$MEA_BIN_BWA aln $MEA_BWA_ALN_PARAMS "$PARAM_GENOME" $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_2.sai
			$MEA_BIN_BWA sampe "$PARAM_GENOME" "$PARAM_BAM_PREFIX"_1.sai "$PARAM_BAM_PREFIX"_2.sai $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
			rm "$PARAM_BAM_PREFIX"_1.sai "$PARAM_BAM_PREFIX"_2.sai
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".bam
#            alignReadsBWA_P $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 $PARAM_GENOME "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
            


			meaCheckFileExists "$MEA_REFERENCE_GENOME"
            printProgress "aligning reads to reference genome"
            $MEA_BIN_BWA aln "$MEA_BWA_ALN_TOTAL_PARAMS" $MEA_REFERENCE_GENOME "$PARAM_FASTQ_FILE1" > "$PARAM_BAM_PREFIX"_1.sai
            $MEA_BIN_BWA aln "$MEA_BWA_ALN_TOTAL_PARAMS" $MEA_REFERENCE_GENOME "$PARAM_FASTQ_FILE2" > "$PARAM_BAM_PREFIX"_2.sai
            $MEA_BIN_BWA sampe $MEA_REFERENCE_GENOME "$PARAM_BAM_PREFIX"_1.sai "$PARAM_BAM_PREFIX"_2.sai "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2" > "$PARAM_BAM_PREFIX"_total.sam
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
        

       	elif [ $MEA_USE_BOWTIE2 = 1 ]; then
            meaCheckFileExists "$PARAM_GENOME".1.bt2*
            printProgress "align to insilico concatenated genome"
            $MEA_BIN_BOWTIE2 $MEA_BOWTIE2_ALN_PARAMS -x "$PARAM_GENOME" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 -S "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
            printProgress "align to the reference genome"
            meaCheckFileExists "$MEA_DIR_REFERENCES"/"bowtie2-index"/"$MEA_BUILD".1.bt2*
            $MEA_BIN_BOWTIE2 "$MEA_BOWTIE2_ALN_TOTAL_PARAMS" -x "$MEA_DIR_REFERENCES"/"bowtie2-index"/"$MEA_BUILD" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 -S "$PARAM_BAM_PREFIX"_total.sam
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
        
        elif [ $MEA_USE_BISMARK = 1 ]; then
            VAR_BASENAME=`basename $PARAM_BAM_PREFIX`
            VAR_DIRNAME=`dirname $PARAM_BAM_PREFIX`
            VAR_q=$MEA_BAM2WIG_PARAM_MIN_QUALITY
            VAR_F=$MEA_BAM2WIG_PARAM_FILTERING_FLAG
            printProgress "align to the insilico concatenated genome"
            meaCheckDirExists "$PARAM_GENOME"/Bisulfite_Genome
            $MEA_BIN_BISMARK $MEA_BISMARK_ALN_PARAMS --basename "$VAR_BASENAME"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2" -o $VAR_DIRNAME $PARAM_GENOME -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_pe.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
            printProgress "aligning to reference genome"
            meaCheckDirExists "$MEA_DIR_REFERENCES"/"$MEA_BUILD"/Bisulfite_Genome
            $MEA_BIN_BISMARK $MEA_BISMARK_ALN_TOTAL_PARAMS --basename "$VAR_BASENAME"_total -o $VAR_DIRNAME "$MEA_DIR_REFERENCES"/"$MEA_BUILD" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2
            mv "$PARAM_BAM_PREFIX"_total_pe.sam "$PARAM_BAM_PREFIX"_total.sam
            $MEA_BIN_SAMTOOLS view -bShu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" ">= 0" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" ">= 0" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            printProgress "extract methylation of the insilico concatenated genome"
            $MEA_BIN_SAMTOOLS view -Sbh -F $VAR_F -q $VAR_q "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam \
                > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".bam
            $MEA_BIN_BISMARK_EXTRACT -p --comprehensive --cytosine_report --samtools_path $MEA_DIR_TOOLS -o $VAR_DIRNAME --genome_folder $PARAM_GENOME \
                "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".bam
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt \
                "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            sort -k1,1 -k2,2n "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt \
                > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt
            if [ $MEA_DEBUG = 0 ]; then
                rm "$VAR_DIRNAME"/C??_context_*.txt
                rm "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            fi
            printProgress "extract methylation of reference genome"
            $MEA_BIN_SAMTOOLS view -Sbh -F $VAR_F -q $VAR_q "$PARAM_BAM_PREFIX"_total.sam > "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".bam
            $MEA_BIN_BISMARK_EXTRACT -p --comprehensive --cytosine_report --samtools_path $MEA_DIR_TOOLS -o $VAR_DIRNAME --genome_folder $MEA_DIR_REFERENCES/$MEA_BUILD \
                "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".bam
            mv "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.txt "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            sort -k1,1 -k2,2n "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt > "$PARAM_BAM_PREFIX"_total.CpG_report.txt
            if [ $MEA_DEBUG = 0 ]; then
                rm "$VAR_DIRNAME"/C??_context_*.txt
                rm "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".CpG_report.unsorted.txt
            fi
            printProgress "detecting allelic cytosines information"
            detectAllelicCytoConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt \
                $PARAM_STRAIN1 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicCytoConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_F"$VAR_F"_q"$VAR_q".CpG_report.txt \
                $PARAM_STRAIN2 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1".CpG_report.txt "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_preProject.CpG_report.txt
            mv "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2".CpG_report.txt "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_preProject.CpG_report.txt
        
        elif [ $MEA_USE_STAR = 1 ]; then
            printProgress "align to the insilico concatenated genome"
            meaCheckDirExists "$PARAM_GENOME"
            $MEA_BIN_STAR --runMode alignReads --genomeDir "$PARAM_GENOME" "$MEA_STAR_ALN_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2"
            mv Aligned.out.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
            mv Log.out "$PARAM_BAM_PREFIX"_STAR_RunParameters.tsv
            mv Log.final.out "$PARAM_BAM_PREFIX"_STAR_AlignmentSummary.tsv
            printProgress "align to the reference genome"
            meaCheckDirExists "$MEA_DIR_REFERENCES"/STAR-index/"$MEA_BUILD"
            $MEA_BIN_STAR --runMode alignReads --genomeDir "$MEA_DIR_REFERENCES"/STAR-index/"$MEA_BUILD"/ "MEA_STAR_ALN_TOTAL_PARAMS" --readFilesIn "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2"
            mv Aligned.out.sam "$PARAM_BAM_PREFIX"_total.sam
            $MEA_BIN_SAMTOOLS view -bShu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            mv Log.out "$PARAM_BAM_PREFIX"_total_STAR_referenceRunParameters.tsv
            mv Log.final.out "$PARAM_BAM_PREFIX"_total_STAR_referenceAlignmentSummary.tsv
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" "== 255" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            if [ $MEA_DEBUG = 0 ]; then
                rm SJ.out.tab Log.progress.out
            fi
        
        elif [ $MEA_USE_TOPHAT2 = 1 ]; then
            printProgress "align to the insilico concatenated genome"
            meaCheckFileExists "$PARAM_GENOME".1.bt2*
            $MEA_BIN_TOPHAT2 $MEA_TOPHAT2_ALN_PARAMS --output-dir ./tophat_out_concat "$PARAM_GENOME" "$PARAM_FASTQ_FILE1" "$PARAM_FASTQ_FILE2"
            mv ./tophat_out_concat/accepted_hits.sam ./"$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
            mv ./tophat_out_concat/align_summary.txt ./"$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"_alignSummary.txt
            printProgress "align to the reference genome"
            meaCheckFileExists "$MEA_DIR_REFERENCES"/"bowtie2-index"/"$MEA_BUILD".1.bt2*
            $MEA_BIN_TOPHAT2 $MEA_TOPHAT2_ALN_TOTAL_PARAMS --output-dir ./tophat_out_total "$MEA_DIR_REFERENCES"/bowtie2-index/"$MEA_BUILD" $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2
            mv ./tophat_out_total/accepted_hits.sam "$PARAM_BAM_PREFIX"_total.sam
            mv ./tophat_out_total/align_summary.txt ./"$PARAM_BAM_PREFIX"_total_alignSummary.txt
            $MEA_BIN_SAMTOOLS view -bSu "$PARAM_BAM_PREFIX"_total.sam | $MEA_BIN_SAMTOOLS sort - "$PARAM_BAM_PREFIX"_total
            $MEA_BIN_SAMTOOLS index "$PARAM_BAM_PREFIX"_total.bam
            printProgress "detecting allelic reads"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
            detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" ">= 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"
            if [ $MEA_DEBUG = 0 ]; then
                rm -r tophat_out_*
            fi
        fi
    fi
fi
