#!/bin/bash

pushd `dirname $0` > /dev/null
AL_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

source $AL_DIR_TOOLS/alea.config

##############################################################################
#############   Module 3: allele-specific alignment
##############################################################################

PARAM_VALID=1
PARAM_SINGLE_READS=1
if [ $AL_USE_CONCATENATED_GENOME = 1 ]; then
    if [ "$1" = "-s" ]; then
        if [ $# -eq 6 ]; then
            PARAM_FASTQ_FILE=$2
            PARAM_GENOME=$3
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
            PARAM_STRAIN1=$5
            PARAM_STRAIN2=$6
            PARAM_BAM_PREFIX=$7
        else
            PARAM_VALID=0
        fi
    else
        PARAM_VALID=0
    fi
else #[ $AL_USE_CONCATENATED_GENOME != 1 ]
    if [ "$1" = "-s" ]; then
        if [ $# -eq 7 ]; then
            PARAM_FASTQ_FILE=$2
            PARAM_GENOME1=$3
            PARAM_GENOME2=$4
            PARAM_STRAIN1=$5
            PARAM_STRAIN2=$6
            PARAM_BAM_PREFIX=$7
        else
            PARAM_VALID=0
        fi
    elif [ "$1" = "-p" ]; then
        if [ $# -eq 8 ]; then
            PARAM_SINGLE_READS=0
            PARAM_FASTQ_FILE1=$2
            PARAM_FASTQ_FILE2=$3
            PARAM_GENOME1=$4
            PARAM_GENOME2=$5
            PARAM_STRAIN1=$6
            PARAM_STRAIN2=$7
            PARAM_BAM_PREFIX=$8
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
    using concatenated genome method (AL_USE_CONCATENATED_GENOME=1):
        alea alignReads <-s/-p> <input_reads_1 [input_reads_2]> <genome_concat> <strain1 strain2> <outputPrefix>

    using separate insilico genomes method (AL_USE_CONCATENATED_GENOME=0):
        alea alignReads <-s/-p> <input_reads_1 [input_reads_2]> <genome1  genome2> <strain1 strain2> <outputPrefix>

Options:
    -s              to align single-end reads (requires one input file)
    -p              to align paired-end reads (requires two input files)
         
    input_reads_1   the 1st input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)
                    
    input_reads_2   (paired end) the 2nd input reads file in fastq.
                    (fastq.gz or bam is supported when using BWA)
                    
    genome_concat   (when AL_USE_CONCATENATED_GENOME=1)
                    path to the indexed reference for concatenated insilico genome.
                    for BWA, specifiy path to the fasta.
                    for Bowtie, specify basename of index file
                    
    genome1         (when AL_USE_CONCATENATED_GENOME=0)
                    path to the indexed reference for 1st insilico genome (of strain1).
                    for BWA, specifiy the fasta file.
                    for Bowtie, specify index filename prefix (minus trailing .X.ebwt or .X.bt2)
                    
    genome1         (when AL_USE_CONCATENATED_GENOME=0)
                    path to the indexed reference for 2nd insilico genome (of strain2).
                    for BWA, specifiy the fasta file.
                    for Bowtie, specify index filename prefix (minus trailing .X.ebwt or .X.bt2)
                    
    strain1         name of strain1
                    (e.g. hap1 or CASTEiJ)
                    
    strain2         name of strain2
                    (e.g. hap2 or C57BL6J)
                    
    outputPrefix    prefix for output files, including the full path, without an extension
                    (e.g. ./TSC_H3K36me3 )

Output:
        
    outputPrefix_strain1_starin2.sam   (when AL_USE_CONCATENATED_GENOME=1)
                                       all reads aligned to the concatenated insilico genome
                                       
    outputPrefix_strain1_all.sam       (when AL_USE_CONCATENATED_GENOME=0)
                                       all reads aligned to the first insilico genome
                                       
    outputPrefix_strain2_all.sam       (when AL_USE_CONCATENATED_GENOME=0)
                                       all reads aligned to the second insilico genome
                                       
    outputPrefix_strain1.bam           allelic reads for strain1 genome (sorted bam)
    
    outputPrefix_strain2.bam           allelic reads for strain2 genome (sorted bam)

Examples:
    (AL_USE_CONCATENATED_GENOME=1, AL_USE_BWA=1)
    alea alignReads -s H3K36me3.fastq CASTEiJ_C57BL6J.fasta CASTEiJ C57BL6J ./H3K36me3
    alea alignReads -p H3K36me3_1.fastq H3K36me3_2.fastq CASTEiJ_C57BL6J.fasta CASTEiJ C57BL6J ./H3K36me3

    (AL_USE_CONCATENATED_GENOME=0, AL_USE_BWA=1)
    alea alignReads -s H3K36me3.fastq CASTEiJ.fasta C57BL6J.fasta CASTEiJ C57BL6J ./H3K36me3
    
    (AL_USE_CONCATENATED_GENOME=0, AL_USE_BOWTIE1=1)
    alea alignReads -s H3K36me3.fastq bowtie1-index/CASTEiJ bowtie1-index/C57BL6J CASTEiJ C57BL6J ./H3K36me3  
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
    $AL_BIN_BWA aln $AL_BWA_ALN_PARAMS $PARAM_IS_BAM $PARAM_FASTA_FILE $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX".sai
  
    # generate alignment (single ended)
    $AL_BIN_BWA samse $PARAM_FASTA_FILE "$PARAM_BAM_PREFIX".sai $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX".sam
  
    # convert SAM to BAM
    #$AL_BIN_SAMTOOLS view -b -S -o "$PARAM_BAM_PREFIX".bam "$PARAM_BAM_PREFIX".sam
  
    if [ $AL_DEBUG = 0 ]; then
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
    $AL_BIN_BWA aln $AL_BWA_ALN_PARAMS $PARAM_IS_BAM1 $PARAM_FASTA_FILE $PARAM_FASTQ_FILE1 > "$PARAM_BAM_PREFIX"_1.sai
    $AL_BIN_BWA aln $AL_BWA_ALN_PARAMS $PARAM_IS_BAM2 $PARAM_FASTA_FILE $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_2.sai
  
    # generate alignment (paired-ended)
    $AL_BIN_BWA sampe $PARAM_FASTA_FILE "$PARAM_BAM_PREFIX"_1.sai "$PARAM_BAM_PREFIX"_2.sai $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX".sam
  
    # convert SAM to BAM and sort
    #$AL_BIN_SAMTOOLS view -b -S -o "$PARAM_BAM_PREFIX".bam "$PARAM_BAM_PREFIX".sam
  
    if [ $AL_DEBUG = 0 ]; then
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
    $AL_BIN_SAMTOOLS view -SH "$PARAM_INPUT_SAM" \
        | awk -v ref="$PARAM_STRAIN" '($0 ~ ref) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_chr//g' \
        > "$PARAM_OUT_PREFIX".sam
    
    # append reads
    $AL_BIN_SAMTOOLS view -S $PARAM_INPUT_SAM \
        | awk -v ref="$PARAM_STRAIN" '(($3 ~ ref)&&($5'"$PARAM_QUALITY"')) {print $0}' \
        | sed 's/'"$PARAM_STRAIN"'_chr//g' \
        >> "$PARAM_OUT_PREFIX".sam

    # convert to bam
    $AL_BIN_SAMTOOLS view -bt $PARAM_REFFASTA "$PARAM_OUT_PREFIX".sam > "$PARAM_OUT_PREFIX".unsorted.bam
    
    # sort by coordinates
    $AL_BIN_SAMTOOLS sort "$PARAM_OUT_PREFIX".unsorted.bam "$PARAM_OUT_PREFIX" 
    
    if [ -f "$PARAM_OUT_PREFIX".bam ]
    then
        printProgress "[detectAllelicConcatenated] Filtered BAM file created."
        # remove temp files
        if [ $AL_DEBUG = 0 ]; then
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
# detect the allelic reads from the reads aligned to separate insilico genomes
#------------------------------------------------------------------------------------
function detectAllelicSeparate {

    local PARAM_INPUT_SAM1=$1
    local PARAM_INPUT_SAM2=$2
    local PARAM_OUT_PREFIX1=$3
    local PARAM_OUT_PREFIX2=$4

    # extract MAPQ and CIGAR columns to one file
    awk '{print $5"\t"$6}' "$PARAM_INPUT_SAM1" > "$PARAM_INPUT_SAM1".c56
    awk '{print $5"\t"$6}' "$PARAM_INPUT_SAM2" > "$PARAM_INPUT_SAM2".c56
    paste "$PARAM_INPUT_SAM1".c56 "$PARAM_INPUT_SAM2".c56 > "$PARAM_INPUT_SAM1".c56c56
    
    # create output sam headers
    samtools view -SH "$PARAM_INPUT_SAM1" > "$PARAM_OUT_PREFIX1".unsorted.sam
    samtools view -SH "$PARAM_INPUT_SAM2" > "$PARAM_OUT_PREFIX2".unsorted.sam
    
    # detect allelic reads
    awk '{if (NR==FNR) {if ($1 > $3 && $4 == "*") a[NR]=1;} else {if (FNR in a) print $0;}}' "$PARAM_INPUT_SAM1".c56c56 "$PARAM_INPUT_SAM1" >> "$PARAM_OUT_PREFIX1".unsorted.sam
    awk '{if (NR==FNR) {if ($3 > $1 && $2 == "*") a[NR]=1;} else {if (FNR in a) print $0;}}' "$PARAM_INPUT_SAM1".c56c56 "$PARAM_INPUT_SAM2" >> "$PARAM_OUT_PREFIX2".unsorted.sam
    
    # create .bam
    samtools view -bS "$PARAM_OUT_PREFIX1".unsorted.sam > "$PARAM_OUT_PREFIX1".unsorted.bam
    samtools view -bS "$PARAM_OUT_PREFIX2".unsorted.sam > "$PARAM_OUT_PREFIX2".unsorted.bam
    
    # sort .bam
    samtools sort "$PARAM_OUT_PREFIX1".unsorted.bam "$PARAM_OUT_PREFIX1"
    samtools sort "$PARAM_OUT_PREFIX2".unsorted.bam "$PARAM_OUT_PREFIX2"
    
    # remove temp files
    if [ $AL_DEBUG = 0 ]; then
        rm -f "$PARAM_INPUT_SAM1".c56
        rm -f "$PARAM_INPUT_SAM2".c56
        rm -f "$PARAM_INPUT_SAM1".c56c56
        rm -f "$PARAM_OUT_PREFIX1".unsorted.*
        rm -f "$PARAM_OUT_PREFIX2".unsorted.*
    fi
    
    
    #TODO: check for ($0 ~ "XM:i:0" && $0 ~ "XG:i:0" && $0 ~ "XO:i:0" ) ?
}

#------------------------------------------------------------------------------------
# alignReads
#------------------------------------------------------------------------------------

if [ $AL_USE_CONCATENATED_GENOME = 1 ]; then
    
    if [ $PARAM_SINGLE_READS = 1 ]; then
        #align single-end reads to concatenated insilico genome
        aleaCheckFileExists "$PARAM_FASTQ_FILE"

        if [ $AL_USE_BWA = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME"
            alignReadsBWA_S "$PARAM_FASTQ_FILE" "$PARAM_GENOME" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
        elif [ $AL_USE_BOWTIE1 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.ebwt
            $AL_BIN_BOWTIE1 $AL_BOWTIE1_ALN_PARAMS "$PARAM_GENOME" $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
        elif [ $AL_USE_BOWTIE2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.bt2
            $AL_BIN_BOWTIE2 $AL_BOWTIE2_ALN_PARAMS -x "$PARAM_GENOME" $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
        fi
    else #[ $PARAM_SINGLE_READS = 0 ]
        #align paired-end reads to concatenated insilico genome
        
        aleaCheckFileExists $PARAM_FASTQ_FILE1
        aleaCheckFileExists $PARAM_FASTQ_FILE2

        if [ $AL_USE_BWA = 1 ]; then
            aleaCheckFileExists $PARAM_GENOME
            alignReadsBWA_P $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 $PARAM_GENOME "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
        elif [ $AL_USE_BOWTIE1 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.ebwt
            $AL_BIN_BOWTIE1 $AL_BOWTIE1_ALN_PARAMS "$PARAM_GENOME" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
        elif [ $AL_USE_BOWTIE2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME".1.bt2
            $AL_BIN_BOWTIE2 $AL_BOWTIE2_ALN_PARAMS -x "$PARAM_GENOME" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam
        fi
    fi
    
    # extract allelic reads
    detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN1" "$PARAM_GENOME" "> 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"
    detectAllelicConcatenated "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2".sam "$PARAM_STRAIN2" "$PARAM_GENOME" "> 1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"

else # [ $AL_USE_CONCATENATED_GENOME != 1 ] 

    if [ $PARAM_SINGLE_READS = 1 ]; then
        #align single-end reads to both insilico genome

        aleaCheckFileExists $PARAM_FASTQ_FILE
        
        if [ $AL_USE_BWA = 1 ]; then
            
            aleaCheckFileExists $PARAM_GENOME1
            aleaCheckFileExists $PARAM_GENOME2
            alignReadsBWA_S $PARAM_FASTQ_FILE $PARAM_GENOME1 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_all
            alignReadsBWA_S $PARAM_FASTQ_FILE $PARAM_GENOME2 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_all
            
        elif [ $AL_USE_BOWTIE1 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME1".1.ebwt
            aleaCheckFileExists "$PARAM_GENOME2".1.ebwt
            $AL_BIN_BOWTIE1 $AL_BOWTIE1_ALN_PARAMS "$PARAM_GENOME1" $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_all.sam
            $AL_BIN_BOWTIE1 $AL_BOWTIE1_ALN_PARAMS "$PARAM_GENOME2" $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_all.sam
        elif [ $AL_USE_BOWTIE2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME1".1.bt2
            aleaCheckFileExists "$PARAM_GENOME2".1.bt2
            $AL_BIN_BOWTIE2 $AL_BOWTIE2_ALN_PARAMS -x "$PARAM_GENOME1" $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_all.sam
            $AL_BIN_BOWTIE2 $AL_BOWTIE2_ALN_PARAMS -x "$PARAM_GENOME2" $PARAM_FASTQ_FILE > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_all.sam
        fi
        
    else #[ $PARAM_SINGLE_READS = 0 ]
        #align paired-end reads to both insilico genome
        
        aleaCheckFileExists $PARAM_FASTQ_FILE1
        aleaCheckFileExists $PARAM_FASTQ_FILE2

        if [ $AL_USE_BWA = 1 ]; then
            aleaCheckFileExists $PARAM_GENOME1
            aleaCheckFileExists $PARAM_GENOME2
            alignReadsBWA_P $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 $PARAM_GENOME1 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_all
            alignReadsBWA_P $PARAM_FASTQ_FILE1 $PARAM_FASTQ_FILE2 $PARAM_GENOME2 "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_all
        elif [ $AL_USE_BOWTIE1 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME1".1.ebwt
            aleaCheckFileExists "$PARAM_GENOME2".1.ebwt
            $AL_BIN_BOWTIE1 $AL_BOWTIE1_ALN_PARAMS "$PARAM_GENOME1" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_all.sam
            $AL_BIN_BOWTIE1 $AL_BOWTIE1_ALN_PARAMS "$PARAM_GENOME2" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_all.sam
        elif [ $AL_USE_BOWTIE2 = 1 ]; then
            aleaCheckFileExists "$PARAM_GENOME1".1.bt2
            aleaCheckFileExists "$PARAM_GENOME2".1.bt2
            $AL_BIN_BOWTIE2 $AL_BOWTIE2_ALN_PARAMS -x "$PARAM_GENOME1" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_all.sam
            $AL_BIN_BOWTIE2 $AL_BOWTIE2_ALN_PARAMS -x "$PARAM_GENOME2" -1 $PARAM_FASTQ_FILE1 -2 $PARAM_FASTQ_FILE2 > "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_all.sam
        fi
    fi
    
    detectAllelicSeparate "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_all.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_all.sam "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1" "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"

fi    
