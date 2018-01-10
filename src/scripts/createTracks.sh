#!/bin/bash

pushd `dirname $0` > /dev/null
MEA_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null

# MANUAL INSTALLATION
#source $MEA_DIR_TOOLS/mea.config
# DOCKER INSTALLATION
source /mea-data/mea.config


if test $# -ne 7
then
    echo "
Usage:   
         $MEA createTracks <-s/-p> bamPrefix strain1 strain2 strain1.refmap strain2.refmap outputDIR
         
Options:
         -s to create tracks for the single-end aligned reads
         -p to create tracks for the paired-end aligned reads

         bamPrefix      prefix used for the output of alignReads command
         strain1        name of strain1 (e.g. hap1)
         strain2        name of strain2 (e.g. hap2)
         genome1.refmap path to the refmap file created for insilico genome 1
         genome1.refmap path to the refmap file created for insilico genome 2
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

###converts bam to wig using bedtools
function BAM2WIGbedtools {
    local PARAM_INPUT_PREFIX=$1
    local PARAM_OUTPUT_DIR=$2
    local PARAM_CHROM_SIZES=$MEA_REFERENCE_CHROM_SIZES
    
    local VAR_q=$MEA_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
    local VAR_F=$MEA_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
    local VAR_x=$MEA_BAM2WIG_PARAM_SE_EXTENSION    # average fragment length used for fixed length of the read extension [0]. used for ChIP-seq (SET) only
    local VAR_INPUT_BASENAME=`basename $PARAM_INPUT_PREFIX`
    
    meaCheckFileExists "$PARAM_INPUT_PREFIX".bam
    printProgress "[bam2wig of $VAR_INPUT_BASENAME] started"
    $MEA_BIN_SAMTOOLS view -bh -F "$VAR_F" -q "$VAR_q" "$PARAM_INPUT_PREFIX".bam > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam
    
    if [ "$VAR_OPTION" = "-s" ]; then
        if [ $MEA_USE_BWA = 1 ] || [ $MEA_USE_BOWTIE1 = 1 ] || [ $MEA_USE_BOWTIE2 = 1 ]; then
            $MEA_BIN_BEDTOOLS genomecov -bg -fs "$VAR_x" -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
                > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
        elif [ $MEA_USE_BISMARK = 1 ]; then
            $MEA_BIN_BEDTOOLS genomecov -bg -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
                > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
        elif [ $MEA_USE_STAR = 1 ] || [ $MEA_USE_TOPHAT2 = 1 ]; then
            $MEA_BIN_BEDTOOLS genomecov -bg -split -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
                > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
        fi
    elif [ "$VAR_OPTION" = "-p" ]; then
        if [ $MEA_USE_BWA = 1 ] || [ $MEA_USE_BOWTIE1 = 1 ] || [ $MEA_USE_BOWTIE2 = 1 ]; then
	#removed -pc flag as it produced an empty bedgraph (bedtools v2.26.0, JRA)
            $MEA_BIN_BEDTOOLS genomecov -bg -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
                > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
        elif [ $MEA_USE_BISMARK = 1 ]; then
            $MEA_BIN_BEDTOOLS genomecov -bg -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
                > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
        elif [ $MEA_USE_STAR = 1 ] || [ $MEA_USE_TOPHAT2 = 1 ]; then
            $MEA_BIN_BEDTOOLS genomecov -bg -split -ibam "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bam -g "$PARAM_CHROM_SIZES" \
                > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph
        fi
    else
        echo "Invalid option $VAR_OPTION"
        exit 1
    fi
    sort -k1,1 -k2,2n "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q"_unsorted.bedGraph > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph
    $MEA_BIN_BEDGRAPH_TO_BW "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph "$PARAM_CHROM_SIZES" "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bw #output this file for viz
    
    awk 'BEGIN {
        print "track type=wiggle_0"
    }
    NF == 4 {
        print "fixedStep chrom="$1" start="$2+1" step=1 span=1"
        for(i = 0; i < $3-$2; i++) {
            print $4
        }
    }' "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph > "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".wig
    $MEA_BIN_BGZIP -c "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".wig > "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".wig.gz
    mv "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bedGraph "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".bedGraph
    mv "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".bw "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".bw
    if [ $MEA_DEBUG = 0 ]; then
        rm "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".wig
    else
        mv "$PARAM_INPUT_PREFIX"_F"$VAR_F"_q"$VAR_q".wig "$PARAM_OUTPUT_DIR"/"$VAR_INPUT_BASENAME".wig
    fi
}

### projects a wig profile to reference genome
function projectToReferenceGenome {
    local PARAM_WIG_FILE=$1
    local PARAM_REFMAP_FILE=$2
    local PARAM_BEDGRAPH_FILE=$3
    
    meaCheckFileExists "$PARAM_WIG_FILE"
    meaCheckFileExists "$PARAM_REFMAP_FILE"
    
    printProgress "Started projectToReferenceGenome"
    
    $MEA_BIN_MEA project\
        --input-wig=$PARAM_WIG_FILE\
        --input-refmap=$PARAM_REFMAP_FILE\
        --output-bedgraph=$PARAM_BEDGRAPH_FILE
    
    printProgress "Finished projectToReferenceGenome"
}

### merge a cytosine report into a CpG site report
function mergeTwoStrandMethylation {
    local PARAM_CYTESINE_REPORT_FILE=$1
    local PARAM_SITE_REPORT_FILE=$2
    
    meaCheckFileExists "$PARAM_CYTESINE_REPORT_FILE"
    
    printProgress "Started mergeTwoStrandMethylation"
    
    awk '
        BEGIN{
            FS = "\t"
            OFS = "\t"
            
            FIRST_POS = 0
            FIRST_METHYL = 0
            FIRST_UNMETHYL = 0
            METHYL = 0
            UNMETHYL = 0
            FIRST_TRI = ""
        }
        $3 == "+"{
            FIRST_POS = $2
            FIRST_METHYL = $4
            FIRST_UNMETHYL = $5
            FIRST_TRI = $7
        }
        $3 == "-"{
            if ($2 == FIRST_POS + 1) {
                METHYL = FIRST_METHYL + $4
                UNMETHYL = FIRST_UNMETHYL + $5
                
                if (METHYL + UNMETHYL > 0) {
                    printf $1 "\t" FIRST_POS "\t" $2 "\t"
                    printf "%6f\t", METHYL / (METHYL + UNMETHYL) * 100.0
                    print METHYL, UNMETHYL, $6, FIRST_TRI
                }
                else {
                    print $1, FIRST_POS, $2, "NA", METHYL, UNMETHYL, $6, FIRST_TRI
                }
            }
            FIRST_POS = 0
            FIRST_METHYL = 0
            FIRST_UNMETHYL = 0
            FIRST_TRI = ""
            METHYL = 0
            UNMETHYL = 0
        }
    ' "$PARAM_CYTESINE_REPORT_FILE" > "$PARAM_SITE_REPORT_FILE"
    
    printProgress "Finished mergeTwoStrandMethylation"
}

function methyl2LOC {
    local INPUT=$1
    grep -v "track" $INPUT | awk '{OFS="\t";FS="\t"} {print $1, $2, $3, "1"}' > ${INPUT//.bedGraph/_LOC.bedGraph}
    $MEA_BIN_BEDGRAPH_TO_BW ${INPUT//.bedGraph/_LOC.bedGraph} "$PARAM_CHROM_SIZES" ${INPUT//.bedGraph/_LOC.bw}
    rm ${INPUT//.bedGraph/_LOC.bedGraph}
}

### convert a CpG site report into a Wig
function convertMethylationToWig {
    local PARAM_SITE_REPORT_FILE=$1
    local PARAM_WIG_FILE=$2
    
    local VAR_MIN_DEPTH=$MEA_METH2WIG_PARAM_MIN_DEPTH
    
    meaCheckFileExists "$PARAM_SITE_REPORT_FILE"
    
    printProgress "Started convertMethylationToWig"
    
    awk -v MIN_DEPTH=$VAR_MIN_DEPTH '
        BEGIN{
            FS = "\t"
            OFS = "\t"
            
            if (MIN_DEPTH < 1)
                MIN_DEPTH = 1
            
            print "track type=wiggle_0"
        }
        $5 + $6 >= MIN_DEPTH{
            print "fixedStep chrom=" $1 " start=" $2 " step=1 span=1"
            for(i = 0; i < $3-$2+1; i++)
                print $4
        }
    ' "$PARAM_SITE_REPORT_FILE" > "$PARAM_WIG_FILE"
    
    printProgress "Finished convertMethylationToWig"
}

### convert a CpG site report into a BedGraph
function convertMethylationToBedgraph {
    local PARAM_SITE_REPORT_FILE=$1
    local PARAM_BEDGRAPH_FILE=$2
    
    local VAR_MIN_DEPTH=$MEA_METH2WIG_PARAM_MIN_DEPTH
    
    meaCheckFileExists "$PARAM_SITE_REPORT_FILE"
    
    printProgress "Started convertMethylationToBedgraph"
    
    awk -v MIN_DEPTH=$VAR_MIN_DEPTH '
        BEGIN{
            FS = "\t"
            OFS = "\t"
            
            if (MIN_DEPTH < 1)
                MIN_DEPTH = 1
            
            print "track type=bedGraph"
        }
        $5 + $6 >= MIN_DEPTH{
            print $1, $2 - 1, $3, $4
        }
    ' "$PARAM_SITE_REPORT_FILE" > "$PARAM_BEDGRAPH_FILE"
    
    printProgress "Finished convertMethylationToBedgraph"
}


VAR_OPTION=$1
shift

#function generateAllelicTracks {
    PARAM_BAM_PREFIX=$1
    PARAM_STRAIN1=$2
    PARAM_STRAIN2=$3
    PARAM_REFMAP_FILE1=$4
    PARAM_REFMAP_FILE2=$5
    PARAM_CHROM_SIZES=$MEA_REFERENCE_CHROM_SIZES
    PARAM_OUTPUT_DIR=$6
    
    meaCheckFileExists "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1".bam
    meaCheckFileExists "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2".bam
    meaCheckFileExists "$PARAM_REFMAP_FILE1"
    meaCheckFileExists "$PARAM_REFMAP_FILE2"
    meaCheckFileExists "$PARAM_CHROM_SIZES"
    meaCreateDir "$PARAM_OUTPUT_DIR"
    
    VAR_OUTPUT_BASENAME=`basename $PARAM_BAM_PREFIX`
    VAR_OUTPUT_PREFIX1="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"
    VAR_OUTPUT_PREFIX2="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"
    
    BAM2WIGbedtools "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1" "$PARAM_OUTPUT_DIR" "$PARAM_CHROM_SIZES"
    BAM2WIGbedtools "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2" "$PARAM_OUTPUT_DIR" "$PARAM_CHROM_SIZES"
    mv "$VAR_OUTPUT_PREFIX1".wig.gz "$VAR_OUTPUT_PREFIX1"_preProject.wig.gz
    mv "$VAR_OUTPUT_PREFIX2".wig.gz "$VAR_OUTPUT_PREFIX2"_preProject.wig.gz
    mv "$VAR_OUTPUT_PREFIX1".bedGraph "$VAR_OUTPUT_PREFIX1"_preProject.bedGraph
    mv "$VAR_OUTPUT_PREFIX2".bedGraph "$VAR_OUTPUT_PREFIX2"_preProject.bedGraph
    mv "$VAR_OUTPUT_PREFIX1".bw "$VAR_OUTPUT_PREFIX1"_preProject.bw
    mv "$VAR_OUTPUT_PREFIX2".bw "$VAR_OUTPUT_PREFIX2"_preProject.bw
    if [ $MEA_DEBUG != 0 ]; then
        mv "$VAR_OUTPUT_PREFIX1".wig "$VAR_OUTPUT_PREFIX1"_preProject.wig
        mv "$VAR_OUTPUT_PREFIX2".wig "$VAR_OUTPUT_PREFIX2"_preProject.wig
    fi
    
    projectToReferenceGenome "$VAR_OUTPUT_PREFIX1"_preProject.wig.gz "$PARAM_REFMAP_FILE1" "$VAR_OUTPUT_PREFIX1".bedGraph
    $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX1".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX1".bw
    
    projectToReferenceGenome "$VAR_OUTPUT_PREFIX2"_preProject.wig.gz "$PARAM_REFMAP_FILE2" "$VAR_OUTPUT_PREFIX2".bedGraph
    $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX2".bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX2".bw
    
    BAM2WIGbedtools "$PARAM_BAM_PREFIX"_total "$PARAM_OUTPUT_DIR" "$PARAM_CHROM_SIZES"
    printProgress "[bam2wig] finished successfully"
#}


    if [ $MEA_USE_BISMARK = 1 ]; then
        :
    else

### generate normalized bigwigs
### and place them in UCSC track hub
#function TrackHubGenerate {
#    VAR_OUTPUT_BASENAME=`basename $PARAM_BAM_PREFIX`
#    VAR_OUTPUT_PREFIX1="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"
#    VAR_OUTPUT_PREFIX2="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"

        VAR_q=$MEA_BAM2WIG_PARAM_MIN_QUALITY     # min read quality [0]
        VAR_F=$MEA_BAM2WIG_PARAM_FILTERING_FLAG  # filtering flag [0]
        VAR_OUTPUT_TOTAL="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_total_F"$VAR_F"_q"$VAR_q"
        printProgress "[TrackHub Generate] started"

        meaCreateDir "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"
        touch "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt "$PARAM_OUTPUT_DIR"/Track_Hub/genomes.txt

        RPM_SCALING_FACTOR_tmp=$($MEA_BIN_SAMTOOLS view -c "$PARAM_BAM_PREFIX"_total_F"$VAR_F"_q"$VAR_q".bam)
        RPM_SCALING_FACTOR=$(echo "scale=25;1000000/$RPM_SCALING_FACTOR_tmp" | bc)
        echo "scaling factor: $RPM_SCALING_FACTOR"

        #ref
        awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
            $4 *= RPM_SCALE
            print $0;}' "$VAR_OUTPUT_PREFIX1".bedGraph > "$VAR_OUTPUT_PREFIX1"_tmp1.bedGraph
        grep -v "type" "$VAR_OUTPUT_PREFIX1"_tmp1.bedGraph > "$VAR_OUTPUT_PREFIX1"_tmp2.bedGraph
        sort -k1,1 -k2,2n "$VAR_OUTPUT_PREFIX1"_tmp2.bedGraph > "$VAR_OUTPUT_PREFIX1"_tmp3.bedGraph
        $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX1"_tmp3.bedGraph "$PARAM_CHROM_SIZES" "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_RPM.bw
        rm "$VAR_OUTPUT_PREFIX1"_tmp1.bedGraph "$VAR_OUTPUT_PREFIX1"_tmp2.bedGraph "$VAR_OUTPUT_PREFIX1"_tmp3.bedGraph
        #alt
        awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
            $4 *= RPM_SCALE
            print $0;}' "$VAR_OUTPUT_PREFIX2".bedGraph > "$VAR_OUTPUT_PREFIX2"_tmp1.bedGraph
        grep -v "type" "$VAR_OUTPUT_PREFIX2"_tmp1.bedGraph > "$VAR_OUTPUT_PREFIX2"_tmp2.bedGraph
        sort -k1,1 -k2,2n "$VAR_OUTPUT_PREFIX2"_tmp2.bedGraph > "$VAR_OUTPUT_PREFIX2"_tmp3.bedGraph
        $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX2"_tmp3.bedGraph "$PARAM_CHROM_SIZES" "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_RPM.bw
        rm "$VAR_OUTPUT_PREFIX2"_tmp1.bedGraph "$VAR_OUTPUT_PREFIX2"_tmp2.bedGraph "$VAR_OUTPUT_PREFIX2"_tmp3.bedGraph
        #total
        mv "$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_total.bedGraph "$VAR_OUTPUT_TOTAL".bedGraph
        awk -F "\t" -v RPM_SCALE="$RPM_SCALING_FACTOR" 'BEGIN {OFS="\t"; print $0} {
            $4 *= RPM_SCALE
        print $0;}' "$VAR_OUTPUT_TOTAL".bedGraph > "$VAR_OUTPUT_TOTAL"_RPM.bedGraph
        $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_TOTAL"_RPM.bedGraph "$PARAM_CHROM_SIZES" "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_total_F"$VAR_F"_q"$VAR_q"_RPM.bw

        ### code adapted from Aaron Bogutz, Louis Lefebvre lab (UBC)
        printf "genome "$MEA_BUILD"\ntrackDb "$MEA_BUILD"/trackDb.txt" > "$PARAM_OUTPUT_DIR"/Track_Hub/genomes.txt
        printf "hub <name>\nshortLabel <short name>\nlongLabel <Hub to display data at UCSC>\ngenomesFile genomes.txt\nemail <email>" > "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt
        printf "track %s\ncontainer multiWig\nshortLabel %s\nlongLabel %s\ntype bigWig\nvisibility full\nmaxHeightPixels 150:60:32\nconfigurable on\nautoScale on\naggregate transparentOverlay\nshowSubtrackColorOnUi on\npriority 1.0\n\n" $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
        printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" $VAR_OUTPUT_BASENAME"_total" $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME"_total" $VAR_OUTPUT_BASENAME"_total" $VAR_OUTPUT_BASENAME"_total.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
        printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1 $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1".bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
        printf "\ttrack %s\n\tparent %s\n\tshortLabel %s\n\tlongLabel %s\n\ttype bigWig\n\tbigDataUrl %s\n\tcolor 215,215,215\n\taltColor 225,225,225\n\tautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2 $VAR_OUTPUT_BASENAME $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2 $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2".bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
        $MEA_BIN_HUBCHECK -noTracks "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt

        printProgress "[TrackHub generate] finished successfully"
    #    if [ $MEA_DEBUG == 0 ]; then
    #        rm "$PARAM_OUTPUT_DIR"/*preProject*
    #        rm "$PARAM_OUTPUT_DIR"/*total.wig.gz
    #    fi
#}
    fi





if [ $MEA_USE_BISMARK = 1 ]; then
#function generateAllelicMethylTracks {
    VAR_OUTPUT_PREFIX_TOTAL="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_total
    VAR_BAM_PREFIX_1_2="$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
    VAR_OUTPUT_PREFIX_1_2="$PARAM_OUTPUT_DIR"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_"$PARAM_STRAIN2"
    
    mv "$VAR_OUTPUT_PREFIX1".bedGraph "$VAR_OUTPUT_PREFIX1"_coverage.bedGraph
    mv "$VAR_OUTPUT_PREFIX2".bedGraph "$VAR_OUTPUT_PREFIX2"_coverage.bedGraph
    mv "$VAR_OUTPUT_PREFIX1".bw "$VAR_OUTPUT_PREFIX1"_coverage.bw
    mv "$VAR_OUTPUT_PREFIX2".bw "$VAR_OUTPUT_PREFIX2"_coverage.bw
    mv "$VAR_OUTPUT_PREFIX1"_preProject.wig.gz "$VAR_OUTPUT_PREFIX1"_preProject_coverage.wig.gz
    mv "$VAR_OUTPUT_PREFIX2"_preProject.wig.gz "$VAR_OUTPUT_PREFIX2"_preProject_coverage.wig.gz
    mv "$VAR_OUTPUT_PREFIX1"_preProject.bedGraph "$VAR_OUTPUT_PREFIX1"_preProject_coverage.bedGraph
    mv "$VAR_OUTPUT_PREFIX2"_preProject.bedGraph "$VAR_OUTPUT_PREFIX2"_preProject_coverage.bedGraph
    mv "$VAR_OUTPUT_PREFIX1"_preProject.bw "$VAR_OUTPUT_PREFIX1"_preProject_coverage.bw
    mv "$VAR_OUTPUT_PREFIX2"_preProject.bw "$VAR_OUTPUT_PREFIX2"_preProject_coverage.bw
    mv "$VAR_OUTPUT_PREFIX_TOTAL".bedGraph "$VAR_OUTPUT_PREFIX_TOTAL"_coverage.bedGraph
    mv "$VAR_OUTPUT_PREFIX_TOTAL".bw "$VAR_OUTPUT_PREFIX_TOTAL"_coverage.bw
    mv "$VAR_OUTPUT_PREFIX_TOTAL".wig.gz "$VAR_OUTPUT_PREFIX_TOTAL"_coverage.wig.gz
    if [ $MEA_DEBUG != 0 ]; then
        mv "$VAR_OUTPUT_PREFIX1"_preProject.wig "$VAR_OUTPUT_PREFIX1"_preProject_coverage.wig
        mv "$VAR_OUTPUT_PREFIX2"_preProject.wig "$VAR_OUTPUT_PREFIX2"_preProject_coverage.wig
        mv "$VAR_OUTPUT_PREFIX_TOTAL".wig "$VAR_OUTPUT_PREFIX_TOTAL"_coverage.wig
    fi
    
    mergeTwoStrandMethylation "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN1"_preProject.CpG_report.txt "$VAR_OUTPUT_PREFIX1"_preProject.CpG_site_report.txt
    mergeTwoStrandMethylation "$PARAM_BAM_PREFIX"_"$PARAM_STRAIN2"_preProject.CpG_report.txt "$VAR_OUTPUT_PREFIX2"_preProject.CpG_site_report.txt
    
    convertMethylationToWig "$VAR_OUTPUT_PREFIX1"_preProject.CpG_site_report.txt "$VAR_OUTPUT_PREFIX1"_preProject_methyl.wig
    convertMethylationToWig "$VAR_OUTPUT_PREFIX2"_preProject.CpG_site_report.txt "$VAR_OUTPUT_PREFIX2"_preProject_methyl.wig
    
    projectToReferenceGenome "$VAR_OUTPUT_PREFIX1"_preProject_methyl.wig "$PARAM_REFMAP_FILE1" "$VAR_OUTPUT_PREFIX1"_methyl.bedGraph
    $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX1"_methyl.bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX1"_methyl.bw
    
    projectToReferenceGenome "$VAR_OUTPUT_PREFIX2"_preProject_methyl.wig "$PARAM_REFMAP_FILE2" "$VAR_OUTPUT_PREFIX2"_methyl.bedGraph
    $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX2"_methyl.bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX2"_methyl.bw
    
    cp "$PARAM_BAM_PREFIX"_total.CpG_report.txt "$VAR_OUTPUT_PREFIX_TOTAL".CpG_report.txt
    mergeTwoStrandMethylation "$VAR_OUTPUT_PREFIX_TOTAL".CpG_report.txt "$VAR_OUTPUT_PREFIX_TOTAL".CpG_site_report.txt
    convertMethylationToBedgraph "$VAR_OUTPUT_PREFIX_TOTAL".CpG_site_report.txt "$VAR_OUTPUT_PREFIX_TOTAL"_methyl.bedGraph
    $MEA_BIN_BEDGRAPH_TO_BW "$VAR_OUTPUT_PREFIX_TOTAL"_methyl.bedGraph "$PARAM_CHROM_SIZES" "$VAR_OUTPUT_PREFIX_TOTAL"_methyl.bw

    methyl2LOC "$VAR_OUTPUT_PREFIX1"_methyl.bedGraph
    methyl2LOC "$VAR_OUTPUT_PREFIX2"_methyl.bedGraph
    methyl2LOC "$VAR_OUTPUT_PREFIX_TOTAL"_methyl.bedGraph

    meaCreateDir "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"
    touch "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt "$PARAM_OUTPUT_DIR"/Track_Hub/genomes.txt

    cp "$VAR_OUTPUT_PREFIX1"_methyl_LOC.bw "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl_LOC.bw
    cp "$VAR_OUTPUT_PREFIX2"_methyl_LOC.bw "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl_LOC.bw
    cp "$VAR_OUTPUT_PREFIX_TOTAL"_methyl_LOC.bw "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_total_methyl_LOC.bw

    cp "$VAR_OUTPUT_PREFIX1"_methyl.bw "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl.bw
    cp "$VAR_OUTPUT_PREFIX2"_methyl.bw "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl.bw
    cp "$VAR_OUTPUT_PREFIX_TOTAL"_methyl.bw "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/"$VAR_OUTPUT_BASENAME"_total_methyl.bw



 ### code adapted from Aaron Bogutz, Louis Lefebvre lab (UBC)
    printf "genome "$MEA_BUILD"\ntrackDb "$MEA_BUILD"/trackDb.txt" > "$PARAM_OUTPUT_DIR"/Track_Hub/genomes.txt
    printf "hub <name>\nshortLabel <short name>\nlongLabel <Hub to display data at UCSC>\ngenomesFile genomes.txt\nemail <email>" > "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt

    printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 215,215,215\naltColor 225,225,225\nautoScale on\n\n" $VAR_OUTPUT_BASENAME"_total_methyl" $VAR_OUTPUT_BASENAME"_total_methyl" $VAR_OUTPUT_BASENAME"_total_methyl" $VAR_OUTPUT_BASENAME"_total_methyl.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
    printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 215,215,215\naltColor 225,225,225\nautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
    printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 215,215,215\naltColor 225,225,225\nautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
    printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 215,215,215\naltColor 225,225,225\nautoScale on\n\n" $VAR_OUTPUT_BASENAME"_total_methyl_LOC" $VAR_OUTPUT_BASENAME"_total_methyl_LOC" $VAR_OUTPUT_BASENAME"_total_methyl_LOC" $VAR_OUTPUT_BASENAME"_total_methyl_LOC.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
    printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 215,215,215\naltColor 225,225,225\nautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl_LOC" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl_LOC" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl_LOC" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN1"_methyl_LOC.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt
    printf "track %s\nshortLabel %s\nlongLabel %s\ntype bigWig\nbigDataUrl %s\ncolor 215,215,215\naltColor 225,225,225\nautoScale on\n\n" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl_LOC" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl_LOC" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl_LOC" $VAR_OUTPUT_BASENAME"_"$PARAM_STRAIN2"_methyl_LOC.bw" | tee -a "$PARAM_OUTPUT_DIR"/Track_Hub/"$MEA_BUILD"/trackDb.txt


    $MEA_BIN_HUBCHECK -noTracks "$PARAM_OUTPUT_DIR"/Track_Hub/hub.txt
    printProgress "[TrackHub generate] finished successfully"



#}
fi
