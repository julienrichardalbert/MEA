#!/bin/bash

pushd `dirname $0` > /dev/null
MEA_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

# MANUAL INSTALLATION
#source $MEA_DIR_TOOLS/mea.config
# DOCKER INSTALLATION
source /mea-data/mea.config

##############################################################################
#############   Module 2: creating insilico genome
##############################################################################

if test $# -lt 4
then
    echo "
Usage:
         $MEA createGenome phased.vcf.gz strain1 strain2 outputDir
         $MEA createGenome -snps-indels-separately phased_snps.vcf.gz phased_indels.vcf.gz strain1 strain2 outputDir

Options:
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
    cat "$PARAM_FASTA1" | sed 's/>/>'"$PARAM_STRAIN1"'_/g' > "$PARAM_FASTA_CONCAT".1
    cat "$PARAM_FASTA2" | sed 's/>/>'"$PARAM_STRAIN2"'_/g' > "$PARAM_FASTA_CONCAT".2
    
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
    local cpus_tmp=$(grep -c ^processor /proc/cpuinfo)
	local cpus=$(echo "scale=0;$cpus_tmp-($cpus_tmp/15)" | bc)
    
    printProgress "[createFastaIndex] Started"

    if [ $MEA_USE_BWA = 1 ]; then
        $MEA_BIN_BWA_INDEX $PARAM_FASTA
    elif [ $MEA_USE_BOWTIE2 = 1 ]; then
        meaCreateDir "$PARAM_OUTPUT_DIR"/bowtie2-index
	ln -s "$PARAM_FASTA" "$PARAM_OUTPUT_DIR"/bowtie2-index/"${PARAM_FASTA%.*}.fa"
        $MEA_BIN_BOWTIE2_INDEX  --threads $cpus --large-index "$PARAM_FASTA" "$PARAM_OUTPUT_DIR"/bowtie2-index/"$PARAM_STRAIN"
    elif [ $MEA_USE_BISMARK = 1 ]; then
        meaCreateDir "$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN"
        cp "$PARAM_FASTA" "$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN"
        $MEA_BIN_BISMARK_INDEX --bowtie2 --path_to_bowtie $MEA_DIR_TOOLS "$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN"
    elif [ $MEA_USE_STAR = 1 ]; then
        meaCreateDir "$PARAM_OUTPUT_DIR"/STAR-index/"$PARAM_STRAIN"
		local PATH=$PATH:/opt/mea/bin/STAR
        local totalk=$[600 * $(awk '/^MemTotal:/{print $2}' /proc/meminfo)]
 #       $MEA_BIN_STAR --runMode genomeGenerate --runThreadN $MEA_THREADS --genomeDir "$PARAM_OUTPUT_DIR"/STAR-index/"$PARAM_STRAIN" --genomeFastaFiles "$PARAM_FASTA" --genomeSAsparseD 2 --limitGenomeGenerateRAM $totalk --genomeChrBinNbits 14
		$MEA_BIN_STAR --runMode genomeGenerate --runThreadN $cpus --genomeDir "$PARAM_OUTPUT_DIR"/STAR-index/"$PARAM_STRAIN" --genomeFastaFiles "$PARAM_FASTA" --genomeSAsparseD 2
        cp /*.out "$PARAM_OUTPUT_DIR"/. 
    elif [ $MEA_USE_TOPHAT2 = 1 ]; then
        meaCreateDir "$PARAM_OUTPUT_DIR"/bowtie2-index
        ln -s "$PARAM_FASTA" "$PARAM_OUTPUT_DIR"/bowtie2-index/"${PARAM_FASTA%.*}.fa"
        $MEA_BIN_BOWTIE2_INDEX --threads $cpus --large-index "$PARAM_FASTA" "$PARAM_OUTPUT_DIR"/bowtie2-index/"$PARAM_STRAIN" 
    fi
    
    printProgress "[createFastaIndex] Done"
}

function createRefStrRefmap {
    local PARAM_FASTA_IDX=$1
    local PARAM_OUTPUT_REFMAP=$2
    
    printProgress "[createRefStrRefmap] Started"
    
    awk '
        BEGIN{
            FS = "\t"
            OFS = "\t"
        }
        {
            print ">" $1
            print "0", "0", $2
        }
    ' "$PARAM_FASTA_IDX" > "$PARAM_OUTPUT_REFMAP"
    
    printProgress "[createRefStrRefmap] Done"
}

## creates the insilico genome for each haplotype
#function createInsilicoGenome {
    printProgress "[createGenome] Started"
    
    if [ "$1" = "-snps-indels-separately" ]; then
        
        PARAM_INPUT_FASTA=$MEA_REFERENCE_GENOME
        PARAM_INPUT_VCF_SNPS=$2
        PARAM_INPUT_VCF_INDELS=$3
        PARAM_STRAIN1=$4
        PARAM_STRAIN2=$5
        PARAM_OUTPUT_DIR=$6
        
        meaCheckFileExists $PARAM_INPUT_VCF_SNPS
        meaCheckFileExists $PARAM_INPUT_VCF_INDELS
        meaCheckFileExists $PARAM_INPUT_FASTA
        meaCreateDir $PARAM_OUTPUT_DIR
        
        if [ ! -f "$PARAM_INPUT_FASTA".fai ]; then
            printProgress "[createGenome] Indexing Reference"
            $MEA_BIN_SAMTOOLS faidx "$PARAM_INPUT_FASTA"
        fi
        
        if [ $PARAM_STRAIN1 = $VAR_REFERENCE_STRAIN ]; then
            # copy reference genome as reference stain name
            VAR_FASTA1="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".fasta
            cp "$PARAM_INPUT_FASTA" "$VAR_FASTA1"
            cp "$PARAM_INPUT_FASTA".fai "$VAR_FASTA1".fai
            # create a refmap file for reference strain
            createRefStrRefmap "$VAR_FASTA1".fai "$VAR_FASTA1".refmap
        else
            VAR_GENOME1_SNPS="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".snps.fasta
            # create insilico genome for strain 1 snps using reference
            $MEA_BIN_MEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF_SNPS" \
                --strain="$PARAM_STRAIN1" \
                --output-fasta="$VAR_GENOME1_SNPS"
            
            VAR_FASTA1="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".fasta
            # create insilico genome for strain 1 indels using ref+snps (previously created)
            $MEA_BIN_MEA insilico \
                --input-fasta="$VAR_GENOME1_SNPS" \
                --input-vcf="$PARAM_INPUT_VCF_INDELS" \
                --strain="$PARAM_STRAIN1" \
                --output-fasta="$VAR_FASTA1"
            
            if [ $MEA_DEBUG = 0 ]; then
                #remove partial files
                rm -f "$VAR_GENOME1_SNPS"*
            fi
        fi
        
        if [ $PARAM_STRAIN2 = $VAR_REFERENCE_STRAIN ]; then
            # copy reference genome as reference stain name
            VAR_FASTA2="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".fasta
            cp "$PARAM_INPUT_FASTA" "$VAR_FASTA2"
            cp "$PARAM_INPUT_FASTA".fai "$VAR_FASTA2".fai
            # create a refmap file for reference strain
            createRefStrRefmap "$VAR_FASTA2".fai "$VAR_FASTA2".refmap
        else
            VAR_GENOME2_SNPS="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".snps.fasta
            # create insilico genome for strain 2 snps using reference
            $MEA_BIN_MEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF_SNPS" \
                --strain="$PARAM_STRAIN2" \
                --output-fasta="$VAR_GENOME2_SNPS"
            
            VAR_FASTA2="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".fasta
            
            # create insilico genome for strain 2 indels using ref+snps (previously created)
            $MEA_BIN_MEA insilico \
                --input-fasta="$VAR_GENOME2_SNPS" \
                --input-vcf="$PARAM_INPUT_VCF_INDELS" \
                --strain="$PARAM_STRAIN2" \
                --output-fasta="$VAR_FASTA2"
            
            if [ $MEA_DEBUG = 0 ]; then
                #remove partial files
                rm -f "$VAR_GENOME2_SNPS"*
            fi
        fi
        
    else
        #all varients (snps and indels) are in a single vcf file
        
        PARAM_INPUT_FASTA=$MEA_REFERENCE_GENOME
        PARAM_INPUT_VCF=$1
        PARAM_STRAIN1=$2
        PARAM_STRAIN2=$3
        PARAM_OUTPUT_DIR=$4
        
        meaCheckFileExists $PARAM_INPUT_VCF
        meaCheckFileExists $PARAM_INPUT_FASTA
        meaCreateDir $PARAM_OUTPUT_DIR
        
        if [ ! -f "$PARAM_INPUT_FASTA".fai ]; then
            printProgress "[createGenome] Indexing Reference"
            $MEA_BIN_SAMTOOLS faidx "$PARAM_INPUT_FASTA"
        fi
        
        if [ $PARAM_STRAIN1 = $VAR_REFERENCE_STRAIN ]; then
            # copy reference genome as reference stain name
            VAR_FASTA1="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".fasta
            cp "$PARAM_INPUT_FASTA" "$VAR_FASTA1"
            cp "$PARAM_INPUT_FASTA".fai "$VAR_FASTA1".fai
            # create a refmap file for reference strain
            createRefStrRefmap "$VAR_FASTA1".fai "$VAR_FASTA1".refmap
        else
            VAR_FASTA1="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1".fasta
            # create insilico genome for strain 1    
            $MEA_BIN_MEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF" \
                --strain="$PARAM_STRAIN1" \
                --output-fasta="$VAR_FASTA1"
        fi
        
        if [ $PARAM_STRAIN2 = $VAR_REFERENCE_STRAIN ]; then
            # copy reference genome as reference stain name
            VAR_FASTA2="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".fasta
            cp "$PARAM_INPUT_FASTA" "$VAR_FASTA2"
            cp "$PARAM_INPUT_FASTA".fai "$VAR_FASTA2".fai
            # create a refmap file for reference strain
            createRefStrRefmap "$VAR_FASTA2".fai "$VAR_FASTA2".refmap
        else
            VAR_FASTA2="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN2".fasta
            # create insilico genome for strain 2
            $MEA_BIN_MEA insilico \
                --input-fasta="$PARAM_INPUT_FASTA" \
                --input-vcf="$PARAM_INPUT_VCF" \
                --strain="$PARAM_STRAIN2" \
                --output-fasta="$VAR_FASTA2"
        fi
    fi
    
VAR_FASTA_CONCAT="$PARAM_OUTPUT_DIR"/"$PARAM_STRAIN1"_"$PARAM_STRAIN2".fasta
concatFasta "$VAR_FASTA1" "$VAR_FASTA2" "$PARAM_STRAIN1" "$PARAM_STRAIN2" "$VAR_FASTA_CONCAT"

if [ $MEA_USE_STAR = 1 ]; then
	FOLDER="$PARAM_OUTPUT_DIR/reorder_fasta"
	meaCreateDir "$FOLDER"
	awk -v FOLDER=$FOLDER 'BEGIN {OFS="\t"} {
		do {
			if ($1 ~ ">" ) {
				filename = substr($1,2);
				print $0 > FOLDER"/"filename;
			} else {
				print $0 > FOLDER"/"filename;
			}
		} while (getline)
	}' "$VAR_FASTA_CONCAT" 


	for FILE2 in $FOLDER/*; do
		FILE=$(basename $FILE2)
		name=${FILE##*_}
		name2=${FILE%_*}
		mv $FILE2 $FOLDER/"$name"_"$name2"
	done

	cat `ls $FOLDER/* | sort -V` > "$VAR_FASTA_CONCAT"
	rm -r $FOLDER
fi

$MEA_BIN_SAMTOOLS faidx "$VAR_FASTA_CONCAT"

createFastaIndex "$VAR_FASTA_CONCAT" "$PARAM_STRAIN1"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR"
createFastaIndex "$MEA_REFERENCE_GENOME" "$MEA_BUILD" "$MEA_DIR_REFERENCES"
printProgress "[createGenome] Done"
#}
