pushd `dirname $0` > /dev/null
AL_DIR=`pwd -P` # get the full path to itself
popd > /dev/null
AL_DIR_TOOLS="$AL_DIR"/bin

source "$AL_DIR_TOOLS"/alea.config

mkdir -p "$AL_DIR"/test-data/human
pushd "$AL_DIR"/test-data/human

wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/GRCh37-lite.fa.gz
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/GRCh37-lite.fa.fai
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/H3K36me3_fib01_A08393_1.fastq.gz
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/hg19.chrom.sizes
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/skin01_a20921.phased.vcf.gz
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/skin01_a20921.phased.vcf.gz.tbi
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/skin01_a20921.vcf.gz
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/human/skin01_haplotypes.tar.gz
wget -c http://mathgen.stats.ox.ac.uk/impute/ALL_1000G_phase1integrated_v3_impute.tgz
tar -xzf skin01_haplotypes.tar.gz
gunzip GRCh37-lite.fa.gz
gunzip skin01_a20921.vcf.gz


ls -go

popd

