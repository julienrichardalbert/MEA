#!/bin/bash

pushd `dirname $0` > /dev/null
AL_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

source $AL_DIR_TOOLS/alea.config

##############################################################################
#############   Module 2: creating insilico genome
##############################################################################

if test $# -lt 5
then
    echo "
Usage:
         alea createGenome reference.fasta phased.vcf.gz strain1 strain2 outputDir
         
         alea createGenome -snps-indels-separately reference.fasta phased_snps.vcf.gz phased_indels.vcf.gz strain1 strain2 outputDir

Options:
         reference.fasta        the reference genome fasta file
         phased.vcf.gz          the phased variants vcf file (including SNPs and Indels)
         strain1                name of strain1 exactly as specified in the vcf file (e.g. hap1)
         strain2                name of strain2 exactly as specified in the vcf file (e.g. hap2)
         outputDir              location of the output directory
         
         -snps-indels-separately    use if SNPs and Indels are in two separate vcf files
         phased-snps.vcf.gz         the phased SNPs (should be specified first)
         phased-indels.vcf.gz       the phased Indels  (should be specified second)

Output:
         Creates two parental insilico genomes strain1.fasta and strain2.fasta as well
         as alignment indeces.
         A concatenated genome strain1_strain2.fasta will be created if
         AL_USE_CONCATENATED_GENOME=1 is set in alea.config
         
Note:
         It is possible to have SNPs and Indels in two separate vcf files. In that case
         use -snps-indels-separately option, and make sure you specify SNPs before Indels.         
"
exit 1
fi


#concatenate two insilico genomes to make the 
function concatFasta {
    printProgress "[concatFasta] Started"

    local PARAM_FASTA1=$1
    local PARAM_FASTA2=$2
    local PARAM_STRAIN1=$3
    local PARAM_STRAIN2=$4
    local PARAM_FASTA_CONCAT=$5
    
    #first make the chromosome names unique
    cat "$PARAM_FASTA1" | sed 's/>/>'"$PARAM_STRAIN1"'_chr/g' > "$PARAM_FASTA_CONCAT".1
    cat "$PARAM_FASTA2" | sed 's/>/>'"$PARAM_STRAIN2"'_chr/g' > "$PARAM_FASTA_CONCAT".2
    
    # concatenate the two fasta files. (add a newline inbetween)
    echo >> "$PARAM_FASTA_CONCAT".1 
    cat "$PARAM_FASTA_CONCAT".2 >> "$PARAM_FASTA_CONCAT".1 
    mv "$PARAM_FASTA_CONCAT".1 "$PARAM_FASTA_CONCAT"
    rm -f "$PARAM_FASTA_CONCAT".2
    
    printProgress "[concatFasta] Done"
}

function createFastaIndex {
    local PARAM_FASTA=$1
    local PARAM_STRAIN=$2
    local PARAM_OUTPUT_DIR=$3

    printProgress "[createFastaIndex] Started"

    if [ $AL_USE_BWA = 1 ]; then
        $AL_BIN_BWA_INDEX $PARAM_FASTA
    elif [ $AL_USE_BOWTIE1 = 1 ]; then
        aleaCreateDir "$PARAM_OUTPUT_DIR"/bowtie1-index
        $AL_BIN_BOWTIE1_INDEX "$PARAM_FASTA" "$PARAM_OUTPUT_DIR"/bowtie1-index/"$PARAM_STRAIN"
    elif [ $AL_USE_BOWTIE2 = 1 ]; then
        aleaCreateDir "$PARAM_OUTPUT_DIR"/bowtie2-index
        $AL_BIN_BOWTIE2_INDEX "$PARAM_FASTA" "$PARAM_OUTPUT_DIR"/bowtie2-index/"$PARAM_STRAIN"
    fi

    printProgress "[createFastaIndex] Done"
}

## creates the insilico genome for each haplotype
#function createInsilicoGenome {
    printProgress "[createGenome] Started"
    
    if [ "$1" = "-snps-indels-separately" ]; then
        
        PARAM_INPUT_FASTA=$2
        PARAM_INPUT_VCF_SNPS=$3
        PARAM_INPUT_VCF_INDELS=$4
        PARAM_STRAIN1=$5
        PARAM_STRAIN2=$6
        PARAM_OUTPUT_DIR=$7
        
        aleaCheckFileExists $PARAM_INPUT_VCF_SNPS
        aleaCheckFileExists $PARAM_INPUT_VCF_INDELS
        aleaCheckFileExists $PARAM_INPUT_FASTA
        aleaCreateDir $PARAM_OUTPUT_DIR
        
        if [ $PARAM_STRAIN1 = $VAR_REFERENCE_STRAIN ]; then
            # nothing to do. it is the reference genome
            VAR_FASTA1="$PARAM_INPUT_FASTA"
        else
            VAR_GENOME1_SNPS="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".snps.fasta
            # create insilico genome for strain 1 snps using reference
            $AL_BIN_ALEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF_SNPS" \
                --strain="$PARAM_STRAIN1" \
                --output-fasta="$VAR_GENOME1_SNPS"
            
            VAR_FASTA1="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".fasta
            # create insilico genome for strain 1 indels using ref+snps (previously created)
            $AL_BIN_ALEA insilico \
                --input-fasta="$VAR_GENOME1_SNPS" \
                --input-vcf="$PARAM_INPUT_VCF_INDELS" \
                --strain="$PARAM_STRAIN1" \
                --output-fasta="$VAR_FASTA1"
            
            if [ $AL_DEBUG = 0 ]; then
                #remove partial files
                rm -f "$VAR_GENOME1_SNPS"*
            fi
        fi
        
        
        if [ $PARAM_STRAIN2 = $VAR_REFERENCE_STRAIN ]; then
            # nothing to do. it is the reference genome
            VAR_FASTA2="$PARAM_INPUT_FASTA"
        else
            VAR_GENOME2_SNPS="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".snps.fasta
            # create insilico genome for strain 2 snps using reference
            $AL_BIN_ALEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF_SNPS" \
                --strain="$PARAM_STRAIN2" \
                --output-fasta="$VAR_GENOME1_SNPS"
            
            VAR_FASTA2="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".fasta
            
            # create insilico genome for strain 2 indels using ref+snps (previously created)
            $AL_BIN_ALEA insilico \
                --input-fasta="$VAR_GENOME2_SNPS" \
                --input-vcf="$PARAM_INPUT_VCF_INDELS" \
                --strain="$PARAM_STRAIN2" \
                --output-fasta="$VAR_FASTA2"
            
            if [ $AL_DEBUG = 0 ]; then
                #remove partial files
                rm -f "$VAR_GENOME2_SNPS"*
            fi
        fi
        
    else
        #all varients (snps and indels) are in a single vcf file
        
        PARAM_INPUT_FASTA=$1
        PARAM_INPUT_VCF=$2
        PARAM_STRAIN1=$3
        PARAM_STRAIN2=$4
        PARAM_OUTPUT_DIR=$5
        
        aleaCheckFileExists $PARAM_INPUT_VCF
        aleaCheckFileExists $PARAM_INPUT_FASTA
        aleaCreateDir $PARAM_OUTPUT_DIR
        
        if [ $PARAM_STRAIN1 = $VAR_REFERENCE_STRAIN ]; then
            # nothing to do. it is the reference genome
            VAR_FASTA1="$PARAM_INPUT_FASTA"
        else
            VAR_FASTA1="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".fasta
            # create insilico genome for strain 1    
            $AL_BIN_ALEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF" \
                --strain="$PARAM_STRAIN1" \
                --output-fasta="$VAR_FASTA1"
        fi
        
        if [ $PARAM_STRAIN2 = $VAR_REFERENCE_STRAIN ]; then
            # nothing to do. it is the reference genome
            VAR_FASTA2="$PARAM_INPUT_FASTA"
        else
            VAR_FASTA2="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".fasta
            # create insilico genome for strain 2
            $AL_BIN_ALEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF" \
                --strain="$PARAM_STRAIN2" \
                --output-fasta="$VAR_FASTA2"
        fi
    fi
    
    createFastaIndex "$VAR_FASTA1" "$PARAM_STRAIN1" "$PARAM_OUTPUT_DIR"
    createFastaIndex "$VAR_FASTA2" "$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR"

    if [ $AL_USE_CONCATENATED_GENOME = 1 ]; then
        VAR_FASTA_CONCAT="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1"_"$PARAM_STRAIN2".fasta
        concatFasta "$VAR_FASTA1" "$VAR_FASTA2" "$PARAM_STRAIN1" "$PARAM_STRAIN2" "$VAR_FASTA_CONCAT"
    
        $AL_BIN_SAMTOOLS faidx "$VAR_FASTA_CONCAT"
        createFastaIndex "$VAR_FASTA_CONCAT" "$PARAM_STRAIN1"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR"
    fi
    
    printProgress "[createGenome] Done"
#}
