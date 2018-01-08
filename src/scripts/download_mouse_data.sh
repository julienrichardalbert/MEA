pushd `dirname $0` > /dev/null
AL_DIR=`pwd -P` # get the full path to itself
popd > /dev/null
AL_DIR_TOOLS="$AL_DIR"/bin

source "$AL_DIR_TOOLS"/alea.config


mkdir -p "$AL_DIR"/test-data/mouse
pushd "$AL_DIR"/test-data/mouse

wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/mgp.v2.indels.chrname.vcf.gz
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/mgp.v2.indels.chrname.vcf.gz.tbi
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/mgp.v2.snps.chrname.vcf.gz
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/mgp.v2.snps.chrname.vcf.gz.tbi

wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/mm9_build37_mouse.fasta
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/mm9_build37_mouse.fasta.fai
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/mm9.fullchrom.sizes
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/C57BL6J.fasta.refmap
wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/H3K36me3.fastq.gz

wget -c ftp://ftp.bcgsc.ca/supplementary/ALEA/files/test-data/mouse/README.txt

#wget -c ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX160%2FSRX160404/SRR522267/SRR522267.sra
#"$AL_BIN_FASTQ_DUMP" -A H3K36me3 ./SRR522267.sra

ls -go

popd

