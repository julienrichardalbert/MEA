#!/bin/bash

pushd `dirname $0` > /dev/null
AL_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

source $AL_DIR_TOOLS/alea.config


if test $# -ne 8
then
    echo "
Usage:   
         alea createTracks <-s/-p> bamPrefix strain1 strain2 genome1.refmap genome2.refmap chrom.sizes outputDIR
         
Options:
         -s to create tracks for the single-end aligned reads
         -p to create tracks for the paired-end aligned reads

         bamPrefix      prefix used for the output of alignReads command
         strain1        name of strain1 (e.g. hap1)
         strain2        name of strain2 (e.g. hap2)
         genome1.refmap path to the refmap file created for insilico genome 1
         genome1.refmap path to the refmap file created for insilico genome 2
         chrom.sizes    path to the chromosome size file (required for creating .bw)
         outputDIR      output directory (where to create track files)
         
Output:
         outputDIR/outputPrefix_strain1.bedGraph
         outputDIR/outputPrefix_strain1.bw        read profiles for strain1 projected to reference genome
         
         outputDIR/outputPrefix_strain2.bedGraph 
         outputDIR/outputPrefix_strain2.bw        read profiles for strain2 projected to reference genome
         
         outputDIR/outputPrefix_strain1.wig.gz
         outputDIR/outputPrefix_strain2.wig.gz    unprojected read profiles for strain1 and strain2
"
exit 1
fi

##############################################################################
#############   Module 4: projection to reference genome
##############################################################################

### Converts filtered bam files to wig
function convertBam2WigSE {
    local PARAM_INPUT_PREFIX=$1
    local PARAM_OUTPUT_DIR=$2

    local VAR_q=$AL_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
    local VAR_F=$AL_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
    local VAR_x=$AL_BAM2WIG_PARAM_SE_EXTENSION    # average fragment length used for fixed length of the read extension [0]. used for ChIP-seq (SET) only
    local VAR_INPUT_BASENAME=`basename $PARAM_INPUT_PREFIX`
    
    aleaCheckFileExists "$PARAM_INPUT_PREFIX".bam

    # Create a wig profile from the bam file
    $AL_BIN_BAM2WIG \
        -samtools $AL_BIN_SAMTOOLS \
        -bamFile "$PARAM_INPUT_PREFIX".bam \
        -out $PARAM_OUTPUT_DIR/ \
        -q $VAR_q \
        -F $VAR_F \
        -cs \
        -x $VAR_x
    
    mv $PARAM_OUTPUT_DIR/$VAR_INPUT_BASENAME.q"$VAR_q".F"$VAR_F".SET_"$VAR_x".wig.gz $PARAM_OUTPUT_DIR/$VAR_INPUT_BASENAME.wig.gz
}

function convertBam2WigPE {
    local PARAM_INPUT_PREFIX=$1
    local PARAM_OUTPUT_DIR=$2

    local VAR_q=$AL_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
    local VAR_F=$AL_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
    local VAR_x=$AL_BAM2WIG_PARAM_SE_EXTENSION    # average fragment length used for fixed length of the read extension [0]. used for ChIP-seq (SET) only
    local VAR_INPUT_BASENAME=`basename $PARAM_INPUT_PREFIX`

    aleaCheckFileExists "$PARAM_INPUT_PREFIX".bam
    
    # Create a wig profile from the bam file
    $AL_BIN_BAM2WIG \
        -samtools $AL_BIN_SAMTOOLS \
        -bamFile "$PARAM_INPUT_PREFIX".bam \
        -out $PARAM_OUTPUT_DIR/ \
        -q $VAR_q \
        -F $VAR_F \
        -cp \
        -x $VAR_x

    mv "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".q"$VAR_q".F"$VAR_F".PET.wig.gz "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".wig.gz
}


### projects a wig profile to reference genome
function projectToReferenceGenome {
    local PARAM_WIG_FILE=$1
    local PARAM_REFMAP_FILE=$2
    local PARAM_BEDGRAPH_FILE=$3
    
    aleaCheckFileExists "$PARAM_WIG_FILE"
    aleaCheckFileExists "$PARAM_REFMAP_FILE"

    printProgress "Started projectToReferenceGenome"

    $AL_BIN_ALEA project\
        --input-wig=$PARAM_WIG_FILE\
        --input-refmap=$PARAM_REFMAP_FILE\
        --output-bedgraph=$PARAM_BEDGRAPH_FILE
    
    printProgress "Finished projectToReferenceGenome"
}

VAR_OPTION=$1
shift

#function generateAllelicTracks {
    PARAM_BAM_PREFIX=$1
    PARAM_STRAIN1=$2
    PARAM_STRAIN2=$3
    PARAM_REFMAP_FILE1=$4
    PARAM_REFMAP_FILE2=$5
    PARAM_CHROM_SIZES=$6
    PARAM_OUTPUT_DIR=$7
    
    aleaCheckFileExists "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1".bam
    aleaCheckFileExists "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2".bam
    aleaCheckFileExists "$PARAM_REFMAP_FILE1"
    aleaCheckFileExists "$PARAM_REFMAP_FILE2"
    aleaCheckFileExists "$PARAM_CHROM_SIZES"
    aleaCreateDir "$PARAM_OUTPUT_DIR"
    
    
    
    if [ "$VAR_OPTION" = "-s" ]; then
        convertBam2WigSE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1" "$PARAM_OUTPUT_DIR"
        convertBam2WigSE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR"
    elif [ "$VAR_OPTION" = "-p" ]; then
        convertBam2WigPE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1" "$PARAM_OUTPUT_DIR"
        convertBam2WigPE "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR"
    else
        echo "Invalid option $VAR_OPTION"
        exit 1
    fi
    
    VAR_OUTPUT_BASENAME=`basename $PARAM_BAM_PREFIX`
    VAR_OUTPUT_PREFIX1="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"
    VAR_OUTPUT_PREFIX2="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"

    projectToReferenceGenome "$VAR_OUTPUT_PREFIX1".wig.gz "$PARAM_REFMAP_FILE1" "$VAR_OUTPUT_PREFIX1".bedGraph
    $AL_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX1".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX1".bw

    projectToReferenceGenome "$VAR_OUTPUT_PREFIX2".wig.gz "$PARAM_REFMAP_FILE2" "$VAR_OUTPUT_PREFIX2".bedGraph
    $AL_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX2".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX2".bw

#}
