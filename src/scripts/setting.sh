#!/bin/bash

pushd `dirname $0` > /dev/null
MEA_DIR_TOOLS=`pwd -P` # get the full path to itself
popd > /dev/null


# MANUAL INSTALLATION
#source $MEA_DIR_TOOLS/mea.setting
#source $MEA_DIR_TOOLS/mea.config
# DOCKER INSTALLATION


if [ -f /mea-data/mea.config ]
then
    source /mea-data/mea.config
else
    source $MEA_DIR_TOOLS/mea.config
fi


if [ -f /mea-data/mea.setting ]
then
    source /mea-data/mea.setting
else
    source $MEA_DIR_TOOLS/mea.setting
fi



##############################################################################
#############   Module 0: interactive setting
##############################################################################

meaCreateDir $MEA_DIR_REFERENCES
function download_genome {
	local BUILD=$1
	cd mea-data/reference_genomes
	echo "Downloading "$BUILD".chrom.sizes from UCSC..."
	wget http://hgdownload-test.cse.ucsc.edu/goldenPath/"$BUILD"/bigZips/"$BUILD".chrom.sizes
    echo "Downloading "$BUILD".fa from UCSC"
	wget http://hgdownload-test.cse.ucsc.edu/goldenPath/"$BUILD"/bigZips/"$BUILD".2bit
	$MEA_DIR_TOOLS/twoBitToFa "$BUILD".2bit "$BUILD".fa && rm "$BUILD".2bit
	sed -i 's/>chr/>/g' "$BUILD".fa
	sed -i 's/chr//g' "$BUILD".chrom.sizes
    echo "Download complete. Setting new path variables"
	INPUT_BUILD="$BUILD"
	INPUT_REFERENCE_GENOME=""$MEA_DIR_REFERENCES"/"$BUILD".fa"
	INPUT_REFERENCE_CHROM_SIZES=""$MEA_DIR_REFERENCES"/"$BUILD".chrom.sizes"
}

#-----------------------------------------------------------------------------
# Modifies the setting file
#-----------------------------------------------------------------------------
function modifySetting {
    local PARAM_ANALYSIS_TYPE=$1
    local PARAM_ALIGNER_TYPE=$2
    local PARAM_REFERENCE_TYPE=$3
    
    cp $MEA_DIR_TOOLS/mea.setting $MEA_DIR_TOOLS/.mea.setting.old
    
    case $PARAM_ANALYSIS_TYPE in
        1 ) # in the case of ChIP-seq
            cat $MEA_DIR_TOOLS/mea.setting | sed -e "
                s/MEA_ANALYSIS_TYPE=.*/MEA_ANALYSIS_TYPE=$PARAM_ANALYSIS_TYPE/
                s/MEA_ALIGNER_TYPE_CHIP=.*/MEA_ALIGNER_TYPE_CHIP=$PARAM_ALIGNER_TYPE/
                s/MEA_REFERENCE_TYPE=.*/MEA_REFERENCE_TYPE=$PARAM_REFERENCE_TYPE/
            " > $MEA_DIR_TOOLS/.mea.setting.tmp
            ;;
        
        2 ) # in the case of RNA-seq
            cat $MEA_DIR_TOOLS/mea.setting | sed -e "
                s/MEA_ANALYSIS_TYPE=.*/MEA_ANALYSIS_TYPE=$PARAM_ANALYSIS_TYPE/
                s/MEA_ALIGNER_TYPE_RNA=.*/MEA_ALIGNER_TYPE_RNA=$PARAM_ALIGNER_TYPE/
                s/MEA_REFERENCE_TYPE=.*/MEA_REFERENCE_TYPE=$PARAM_REFERENCE_TYPE/
            " > $MEA_DIR_TOOLS/.mea.setting.tmp
            ;;
        
        3) # in the case of BS-seq
            cat $MEA_DIR_TOOLS/mea.setting | sed -e "
                s/MEA_ANALYSIS_TYPE=.*/MEA_ANALYSIS_TYPE=$PARAM_ANALYSIS_TYPE/
                s/MEA_ALIGNER_TYPE_BS=.*/MEA_ALIGNER_TYPE_BS=$PARAM_ALIGNER_TYPE/
                s/MEA_REFERENCE_TYPE=.*/MEA_REFERENCE_TYPE=$PARAM_REFERENCE_TYPE/
            " > $MEA_DIR_TOOLS/.mea.setting.tmp
            ;;
        
        * )
            echo "Invalid analysis type"
            exit 1
            ;;
    esac
    
    if [ $MEA_DEBUG = 0 ]
    then
		cp $MEA_DIR_TOOLS/.mea.setting.tmp /mea-data/mea.setting
        mv $MEA_DIR_TOOLS/.mea.setting.tmp $MEA_DIR_TOOLS/mea.setting
    else
        cp $MEA_DIR_TOOLS/.mea.setting.tmp $MEA_DIR_TOOLS/mea.setting
		cp $MEA_DIR_TOOLS/.mea.setting.tmp /mea-data/mea.setting
    fi
    
}

#-----------------------------------------------------------------------------
# Modifies the config file
#-----------------------------------------------------------------------------
function modifyConfig {
    local PARAM_ANALYSIS_TYPE=$1
    local PARAM_ALIGNER_TYPE=$2
    local PARAM_REFERENCE_TYPE=$3
    
    cp $MEA_DIR_TOOLS/mea.config $MEA_DIR_TOOLS/.mea.config.old
    
    case $PARAM_ANALYSIS_TYPE in
        1 ) # in the case of ChIP-seq
            case $PARAM_ALIGNER_TYPE in
                1 ) # in the case of BWA
                    cat $MEA_DIR_TOOLS/mea.config | sed -e "
                        s/MEA_USE_BWA=[01]/MEA_USE_BWA=1/
                        s/MEA_USE_BOWTIE1=[01]/MEA_USE_BOWTIE1=0/
                        s/MEA_USE_BOWTIE2=[01]/MEA_USE_BOWTIE2=0/
                        s/MEA_USE_BISMARK=[01]/MEA_USE_BISMARK=0/
                        s/MEA_USE_STAR=[01]/MEA_USE_STAR=0/
                        s/MEA_USE_TOPHAT2=[01]/MEA_USE_TOPHAT2=0/
                    " > $MEA_DIR_TOOLS/.mea.config.tmp
                    ;;
                
                2 ) # in the case of Bowtie1
                    cat $MEA_DIR_TOOLS/mea.config | sed -e "
                        s/MEA_USE_BWA=[01]/MEA_USE_BWA=0/
                        s/MEA_USE_BOWTIE1=[01]/MEA_USE_BOWTIE1=1/
                        s/MEA_USE_BOWTIE2=[01]/MEA_USE_BOWTIE2=0/
                        s/MEA_USE_BISMARK=[01]/MEA_USE_BISMARK=0/
                        s/MEA_USE_STAR=[01]/MEA_USE_STAR=0/
                        s/MEA_USE_TOPHAT2=[01]/MEA_USE_TOPHAT2=0/
                    " > $MEA_DIR_TOOLS/.mea.config.tmp
                    ;;
                
                3 ) # in the case of Bowtie2
                    cat $MEA_DIR_TOOLS/mea.config | sed -e "
                        s/MEA_USE_BWA=[01]/MEA_USE_BWA=0/
                        s/MEA_USE_BOWTIE1=[01]/MEA_USE_BOWTIE1=0/
                        s/MEA_USE_BOWTIE2=[01]/MEA_USE_BOWTIE2=1/
                        s/MEA_USE_BISMARK=[01]/MEA_USE_BISMARK=0/
                        s/MEA_USE_STAR=[01]/MEA_USE_STAR=0/
                        s/MEA_USE_TOPHAT2=[01]/MEA_USE_TOPHAT2=0/
                    " > $MEA_DIR_TOOLS/.mea.config.tmp
                    ;;
                
                * )
                    echo "Invalid aligner type"
                    exit 1
                    ;;
            esac
            ;;
        
        2 ) # in the case of RNA-seq
            case $PARAM_ALIGNER_TYPE in
                
                1 ) # in the case of STAR
                    cat $MEA_DIR_TOOLS/mea.config | sed -e "
                        s/MEA_USE_BWA=[01]/MEA_USE_BWA=0/
                        s/MEA_USE_BOWTIE1=[01]/MEA_USE_BOWTIE1=0/
                        s/MEA_USE_BOWTIE2=[01]/MEA_USE_BOWTIE2=0/
                        s/MEA_USE_BISMARK=[01]/MEA_USE_BISMARK=0/
                        s/MEA_USE_STAR=[01]/MEA_USE_STAR=1/
                        s/MEA_USE_TOPHAT2=[01]/MEA_USE_TOPHAT2=0/
                    " > $MEA_DIR_TOOLS/.mea.config.tmp
                    ;;
                
                2 ) # in the case of TOPHAT2
                    cat $MEA_DIR_TOOLS/mea.config | sed -e "
                        s/MEA_USE_BWA=[01]/MEA_USE_BWA=0/
                        s/MEA_USE_BOWTIE1=[01]/MEA_USE_BOWTIE1=0/
                        s/MEA_USE_BOWTIE2=[01]/MEA_USE_BOWTIE2=0/
                        s/MEA_USE_BISMARK=[01]/MEA_USE_BISMARK=0/
                        s/MEA_USE_STAR=[01]/MEA_USE_STAR=0/
                        s/MEA_USE_TOPHAT2=[01]/MEA_USE_TOPHAT2=1/
                    " > $MEA_DIR_TOOLS/.mea.config.tmp
                    ;;
                
                * )
                    echo "Invalid aligner type"
                    exit 1
                    ;;
            esac
            ;;
        
        3 ) # in the case of BS-seq
            case $PARAM_ALIGNER_TYPE in
                1 ) # in the case of Bismark
                    cat $MEA_DIR_TOOLS/mea.config | sed -e "
                        s/MEA_USE_BWA=[01]/MEA_USE_BWA=0/
                        s/MEA_USE_BOWTIE1=[01]/MEA_USE_BOWTIE1=0/
                        s/MEA_USE_BOWTIE2=[01]/MEA_USE_BOWTIE2=0/
                        s/MEA_USE_BISMARK=[01]/MEA_USE_BISMARK=1/
                        s/MEA_USE_STAR=[01]/MEA_USE_STAR=0/
                        s/MEA_USE_TOPHAT2=[01]/MEA_USE_TOPHAT2=0/
                    " > $MEA_DIR_TOOLS/.mea.config.tmp
                    ;;
                
                * )
                    echo "Invalid aligner type"
                    exit 1
                    ;;
            esac
            ;;
        
        * )
            echo "Invalid analysis type"
            exit 1
            ;;
    esac


#	Comment out 7 following lines JRA Dec2017    / 	nevermind, might need those!!
    if [ $MEA_DEBUG = 0 ]
    then
        mv $MEA_DIR_TOOLS/.mea.config.tmp $MEA_DIR_TOOLS/mea.config
    else
        cp $MEA_DIR_TOOLS/.mea.config.tmp $MEA_DIR_TOOLS/mea.config
        mv $MEA_DIR_TOOLS/.mea.config.tmp $MEA_DIR_TOOLS/.mea.config.tmp1
    fi
    
#    case $PARAM_METHOD_TYPE in
#        1 ) # in the case of concatenated method
#            cat $MEA_DIR_TOOLS/mea.config | sed -e "
#                s/MEA_USE_CONCATENATED_GENOME=[01]/MEA_USE_CONCATENATED_GENOME=1/
#            " > $MEA_DIR_TOOLS/.mea.config.tmp
#            ;;
#        
#        2 ) # in the case of separate method
#            cat $MEA_DIR_TOOLS/mea.config | sed -e "
#                s/MEA_USE_CONCATENATED_GENOME=[01]/MEA_USE_CONCATENATED_GENOME=0/
#            " > $MEA_DIR_TOOLS/.mea.config.tmp
#            ;;
#        
#        * )
#            echo "Invalid method type"
#            exit 1
#            ;;
#    esac

    case $PARAM_REFERENCE_TYPE in
        1 ) # in the case of do not change settings
        # DO NOTHING
            cat $MEA_DIR_TOOLS/mea.config | sed -e "
                s/MEA_REFERENCE_TYPE=[123]/MEA_REFERENCE_TYPE=1/
                s:MEA_REFERENCE_GENOME=.*:MEA_REFERENCE_GENOME="$MEA_REFERENCE_GENOME":
                s:MEA_REFERENCE_CHROM_SIZES=.*:MEA_REFERENCE_CHROM_SIZES="$MEA_REFERENCE_CHROM_SIZES":
                s:MEA_BUILD=.*:MEA_BUILD="$MEA_BUILD":
            " > $MEA_DIR_TOOLS/.mea.config.tmp
            ;;

        2 ) # in the case of download genome automatically
            # GET BUILD
            # DO REST AUTOMATICALLY
            cat $MEA_DIR_TOOLS/mea.config | sed -e "
                s/MEA_REFERENCE_TYPE=[123]/MEA_REFERENCE_TYPE=2/
                s:MEA_REFERENCE_GENOME=.*:MEA_REFERENCE_GENOME="$INPUT_REFERENCE_GENOME":
                s:MEA_REFERENCE_CHROM_SIZES=.*:MEA_REFERENCE_CHROM_SIZES="$INPUT_REFERENCE_CHROM_SIZES":
                s:MEA_BUILD=.*:MEA_BUILD="$INPUT_BUILD":
            " > $MEA_DIR_TOOLS/.mea.config.tmp
	    ;;

        3 ) # in the case of set paths manually
            # GET BUILD
            # GET .fa PATH
            # GET .sizes PATH
            cat $MEA_DIR_TOOLS/mea.config | sed -e "
                s/MEA_REFERENCE_TYPE=[123]/MEA_REFERENCE_TYPE=3/
                s:MEA_REFERENCE_GENOME=.*:MEA_REFERENCE_GENOME="$INPUT_REFERENCE_GENOME":
                s:MEA_REFERENCE_CHROM_SIZES=.*:MEA_REFERENCE_CHROM_SIZES="$INPUT_REFERENCE_CHROM_SIZES":
                s:MEA_BUILD=.*:MEA_BUILD="$INPUT_BUILD":
            " > $MEA_DIR_TOOLS/.mea.config.tmp
            ;;

        * )
            echo "Invalid method type"
            exit 1
            ;;
    esac
    
    if [ $MEA_DEBUG = 0 ]
    then
# MANUAL INSTALL
#        mv $MEA_DIR_TOOLS/.mea.config.tmp $MEA_DIR_TOOLS/mea.config
# DOCKER
        mv $MEA_DIR_TOOLS/.mea.config.tmp /mea-data/mea.config
    else
        cp $MEA_DIR_TOOLS/.mea.config.tmp /mea-data/mea.config
        mv $MEA_DIR_TOOLS/.mea.config.tmp $MEA_DIR_TOOLS/.mea.config.tmp2
    fi
}

#------------------------------------------------------------------------------------
# setting
#------------------------------------------------------------------------------------

echo
echo -n "Which type of analysis do you want to do?"
echo "$MEA_ANALYSIS_CHOICE"
echo -n "Press 1, 2 or 3 [current: $MEA_ANALYSIS_TYPE] : "

read INPUT_ANALYSIS_TYPE
if [ -z "$INPUT_ANALYSIS_TYPE" ]
then
    INPUT_ANALYSIS_TYPE=$MEA_ANALYSIS_TYPE
elif [[ ! "$INPUT_ANALYSIS_TYPE" =~ [1-3] ]]
then
    echo "Invalid input. The current status was held."
    INPUT_ANALYSIS_TYPE=$MEA_ANALYSIS_TYPE
fi

case $INPUT_ANALYSIS_TYPE in
    1 ) # in the case of ChIP-seq
        echo 
        echo -n "Which type of alignment tool do you want to use?"
        echo "$MEA_ALIGNER_CHOICE_CHIP"
        echo -n "Press 1, 2 or 3 [current: $MEA_ALIGNER_TYPE_CHIP] : "
        
        read INPUT_ALIGNER_TYPE
        if [ -z "$INPUT_ALIGNER_TYPE" ]
        then
            INPUT_ALIGNER_TYPE=$MEA_ALIGNER_TYPE_CHIP
        elif [[ ! "$INPUT_ALIGNER_TYPE" =~ [1-3] ]]
        then
            echo "Invalid input. The current status was held."
            INPUT_ALIGNER_TYPE=$MEA_ALIGNER_TYPE_CHIP
        fi
        ;;
    
    2 ) # in the case of RNA-seq
        echo 
        echo -n "Which type of alignment tool do you want to use?"
        echo "$MEA_ALIGNER_CHOICE_RNA"
        echo -n "Press 1, 2 or 3 [current: $MEA_ALIGNER_TYPE_RNA] : "
        
        read INPUT_ALIGNER_TYPE
        if [ -z "$INPUT_ALIGNER_TYPE" ]
        then
            INPUT_ALIGNER_TYPE=$MEA_ALIGNER_TYPE_RNA
        elif [[ ! "$INPUT_ALIGNER_TYPE" =~ [1-5] ]]
        then
            echo "Invalid input. The current status was held."
            INPUT_ALIGNER_TYPE=$MEA_ALIGNER_TYPE_RNA
        fi
        ;;
    
    3 ) # in the case of BS-seq
        echo 
        echo -n "Which type of alignment tool do you want to use?"
        echo "$MEA_ALIGNER_CHOICE_BS"
        echo -n "Press 1 [current: $MEA_ALIGNER_TYPE_BS] : "
        
        read INPUT_ALIGNER_TYPE
        if [ -z "$INPUT_ALIGNER_TYPE" ]
        then
            INPUT_ALIGNER_TYPE=$MEA_ALIGNER_TYPE_BS
        elif [[ ! "$INPUT_ALIGNER_TYPE" =~ [1] ]]
        then
            echo "Invalid input. The current status was held."
            INPUT_ALIGNER_TYPE=$MEA_ALIGNER_TYPE_BS
        fi
        ;;
    
    * )
        echo "Invalid analysis type"
        exit 1
        ;;
esac

#echo
#echo -n "Which type of alignment method do you want to use?"
#echo "$MEA_METHOD_CHOICE"
#echo -n "Press 1 or 2 [current: $MEA_METHOD_TYPE] : "

#read INPUT_METHOD_TYPE
#if [ -z "$INPUT_METHOD_TYPE" ]
#then
#    INPUT_METHOD_TYPE=$MEA_METHOD_TYPE
#elif [[ ! "$INPUT_METHOD_TYPE" =~ [1-3] ]]
#then
#    echo "Invalid input. The current status was held."
#    INPUT_METHOD_TYPE=$MEA_METHOD_TYPE
#fi


echo
echo -n "How would you like to specify the reference genome location?"
echo "$MEA_REFERENCE_CHOICE"
echo -n "Press 1, 2 or 3 [current: $MEA_REFERENCE_TYPE] : "

read INPUT_REFERENCE_TYPE
if [ -z "$INPUT_REFERENCE_TYPE" ]
then
    INPUT_REFERENCE_TYPE=$MEA_REFERENCE_TYPE
elif [[ ! "$INPUT_REFERENCE_TYPE" =~ [1-3] ]]
then
    echo "Invalid input. The current status was held."
    INPUT_REFERENCE_TYPE=$MEA_REFERENCE_TYPE
fi

case $INPUT_REFERENCE_TYPE in
    1 ) # in the case of do not change settings
        echo
    	echo "Current settings"
		echo "----------------"
		echo "build : "$MEA_BUILD""
		echo ".fa path : "$MEA_REFERENCE_GENOME""
		echo ".chrom.sizes path : "$MEA_REFERENCE_CHROM_SIZES""
	;;

    2 ) # in the case of automatically download
        echo 
        echo "Which genome would you like to download?"
        echo "Recent builds: mm10 ; hg19 ; ce11 ; dm6 ;rn6 ; sacCer3"
        echo "Older builds: mm9 ; h18 ; ce10 ; dm5 ; rn5 ; sacCer2 ; etc."
        echo -n "Please specify build code [current: $MEA_BUILD] : "

        read INPUT_BUILD
        if [ -z "$INPUT_BUILD" ]
        then
            INPUT_BUILD=$MEA_BUILD
        elif [[ ! "$INPUT_BUILD" =~ [1-9] ]]
        then
            echo "Invalid input. The current status was held."
        fi

	# check whether build already exists in mea-data/reference_genomes/ folder
	if [ -f mea-data/reference_genomes/"$INPUT_BUILD".fa ]; then
	# if present, modify reference genome paths to match user input build
   	    echo "File mea-data/reference_genomes/"$INPUT_BUILD".fa exists. Do not download new one."
	    echo "Setting new path variables to match genome build."
            INPUT_REFERENCE_GENOME=""$MEA_DIR_REFERENCES"/"$INPUT_BUILD".fa"
            INPUT_REFERENCE_CHROM_SIZES=""$MEA_DIR_REFERENCES"/"$INPUT_BUILD".chrom.sizes"
	else
	# if not present, download genome build and put it in mea-data/reference_genomes/
   	    echo "File "$BUILD".fa does not exist. Downloading it now..."
	    download_genome $INPUT_BUILD
	fi

	;;

    3 ) # in the case manually locate genome directory
        echo 
        echo -n "Please specify the complete path (.fa included) to the reference genome file : "

        read INPUT_REFERENCE_GENOME
        if [ -z "$INPUT_REFERENCE_GENOME" ]
        then
            INPUT_REFERENCE_GENOME=$MEA_REFERENCE_GENOME
        else
            echo "New reference genome file located : " $INPUT_REFERENCE_GENOME
	    ln $INPUT_REFERENCE_GENOME $MEA_DIR_REFERENCES/
	    INPUT_REFERENCE_GENOME="$MEA_DIR_REFERENCES"/$(basename $INPUT_REFERENCE_GENOME)
        fi
        echo 
        echo -n "Please specify the complete path (.sizes included) to the reference genome chromosome sizes file : "

        read INPUT_REFERENCE_CHROM_SIZES
        if [ -z "$INPUT_REFERENCE_CHROM_SIZES" ]
        then
            INPUT_REFERENCE_CHROM_SIZES=$MEA_REFERENCE_CHROM_SIZES
        else
            echo "New reference chromosome sizes file located : " $INPUT_REFERENCE_CHROM_SIZES
	    ln $INPUT_REFERENCE_CHROM_SIZES $MEA_DIR_REFERENCES/
	    INPUT_REFERENCE_CHROM_SIZES="$MEA_DIR_REFERENCES"/$(basename $INPUT_REFERENCE_CHROM_SIZES)
        fi
	echo 
	echo -n "Please specify build code [current: $MEA_BUILD] : "
        read INPUT_BUILD
        if [ -z "$INPUT_BUILD" ]
        then
            INPUT_BUILD=$MEA_BUILD
        elif [[ ! "$INPUT_BUILD" =~ [1-9] ]]
        then
            echo "Invalid input. The current status was held."
	    INPUT_BUILD=$MEA_BUILD
        fi

        ;;

    * )
        echo "Invalid reference genome type"
        exit 1
        ;;
esac

modifySetting "$INPUT_ANALYSIS_TYPE" "$INPUT_ALIGNER_TYPE" "$INPUT_REFERENCE_TYPE"
modifyConfig "$INPUT_ANALYSIS_TYPE" "$INPUT_ALIGNER_TYPE" "$INPUT_REFERENCE_TYPE"
