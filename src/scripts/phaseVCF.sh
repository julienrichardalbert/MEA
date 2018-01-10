#!/bin/bash

pushd `dirname $0` > /dev/null
MEA_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

# MANUAL INSTALLATION
#source $MEA_DIR_TOOLS/mea.config
# DOCKER INSTALLATION
source /mea-data/mea.config

##############################################################################
#############   Module 1: creating phased genotypes
##############################################################################


if test $# -ne 3
then
    echo "
Usage:
         mea phaseVCF hapsDIR unphased.vcf outputPrefix

Options:
         hapsDIR        path to the directory containing the .haps files
         unphased.vcf   path to the vcf file containing unphased SNPs and Indels
         outputPrefix   output file prefix including the path but not the extension

Output:
         creates the file outputPrefix.vcf.gz 
"
exit 1
fi


# merges the SHAPEIT2 haplotype outputs (.haps). Adds the chromosome name to the records in the file
function mergeHaps {
    printProgress "[mergeHaps] Started"

    local PARAM_INPUT_DIR=$1
    local PARAM_OUTPUT_FILE=$2
    
    
    for hapfile in $PARAM_INPUT_DIR/*.haps
    do
        filename=`basename "$hapfile"`
        # get chromosome name, remove 'chr':  chrX.phased.haps --> X
        chr=`echo "${filename%.*.*}" | sed 's/chr//g'`
        echo $chr
        cat $hapfile | awk -v chr="$chr" '{ s = chr; for (i = 2; i <= NF; i++) s = s " " $i; print s }' >> $PARAM_OUTPUT_FILE
    done
    
    printProgress "[mergeHaps] Done"
}

# creates a phased variant (phased.vcf.gz) file given the SHAPEIT2 output directory and input variant (.vcf) file
printProgress "[phaseVCF] Started"

PARAM_HAPS_DIR=$1
PARAM_INPUT_VCF=$2
PARAM_OUT_PREFIX=$3

meaCheckDirExists $PARAM_HAPS_DIR
meaCheckFileExists $PARAM_INPUT_VCF


## merge the .haps files and fix the chromosome name field
mergeHaps "$PARAM_HAPS_DIR" "$PARAM_OUT_PREFIX".haps

## create phased vcf files
$MEA_BIN_PHASEVCF "$PARAM_OUT_PREFIX".haps "$PARAM_INPUT_VCF" "$PARAM_OUT_PREFIX".vcf
$MEA_BIN_BGZIP "$PARAM_OUT_PREFIX".vcf
$MEA_BIN_TABIX -p vcf "$PARAM_OUT_PREFIX".vcf.gz

printProgress "[phaseVCF] Done"
