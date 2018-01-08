PARAM_FASTA1=$1
PARAM_FASTA2=$2
PARAM_FASTA_CONCAT=$3

#first make the chromosome names unique
cat $PARAM_FASTA1 | sed 's/>/>hap1_chr/g' > $PARAM_FASTA_CONCAT.1
cat $PARAM_FASTA2 | sed 's/>/>hap2_chr/g' > $PARAM_FASTA_CONCAT.2

# concatenate the two fasta files. (add a newline inbetween)
echo >> $PARAM_FASTA_CONCAT.1 
cat $PARAM_FASTA_CONCAT.2 >> $PARAM_FASTA_CONCAT.1 
mv $PARAM_FASTA_CONCAT.1 $PARAM_FASTA_CONCAT
rm -f $PARAM_FASTA_CONCAT.2

# create the index for the concatenated fasta files
samtools faidx $PARAM_FASTA_CONCAT


#set -x

### create fasta with correct chromosome names
#cat delete/C57BL6J2.fasta | sed 's/>C57BL6J_chr/>chr/g' > C57BL6J.fasta
#cat delete/CASTEiJ2.fasta | sed 's/>CASTEiJ_chr/>chr/g' > CASTEiJ.fasta
#cat /projects/edcc/Allele-specific/mouse/insilico_genomes/mm9/CASTEiJ/CASTEiJ.fasta.refmap | sed 's/>/>chr/g' | sed 's/chrMT/chrM/g' > CASTEiJ.fasta.refmap

#samtools faidx C57BL6J.fasta
#samtools faidx CASTEiJ.fasta

#samtools faidx C57BL6J_CASTEiJ.fasta
#bwa index C57BL6J_CASTEiJ.fasta > C57BL6J_CASTEiJ_bwaindex.log

### align reads to genome
#bwa aln -t 4 -l 25 C57BL6J_CASTEiJ.fasta ~/Allele-specific/mouse/celltypes/trophoblast/H3K36me3/CASTEiJ_C57BL6NJ/H3K36me3.nodup.fastq > H3K36me3_C57BL6J_CASTEiJ.sai
#bwa samse C57BL6J_CASTEiJ.fasta H3K36me3_C57BL6J_CASTEiJ.sai ~/Allele-specific/mouse/celltypes/trophoblast/H3K36me3/CASTEiJ_C57BL6NJ/H3K36me3.nodup.fastq > H3K36me3_C57BL6J_CASTEiJ.sam

### convert to bam
#samtools view -bt C57BL6J_CASTEiJ.fasta H3K36me3_C57BL6J_CASTEiJ.sam > H3K36me3_C57BL6J_CASTEiJ.bam

### sort by name
#samtools sort -n H3K36me3_C57BL6J_CASTEiJ.bam H3K36me3_C57BL6J_CASTEiJ.sortByName

### output headers
#samtools view -H H3K36me3_C57BL6J_CASTEiJ.sortByName.bam | awk '($0 ~ "C57BL6J") {print $0}' | sed 's/C57BL6J_chr/chr/g' > H3K36me3_C57BL6J.sortByName.sam
#samtools view -H H3K36me3_C57BL6J_CASTEiJ.sortByName.bam | awk '($0 ~ "CASTEiJ") {print $0}' | sed 's/CASTEiJ_chr/chr/g' > H3K36me3_CASTEiJ.sortByName.sam

### filter allele-specific
#samtools view H3K36me3_C57BL6J_CASTEiJ.sortByName.bam | awk '(($3 ~ "C57BL6J")&&($5 > 20)) {print $0}' >> H3K36me3_C57BL6J.sortByName.sam
#samtools view H3K36me3_C57BL6J_CASTEiJ.sortByName.bam | awk '(($3 ~ "CASTEiJ")&&($5 > 20)) {print $0}' >> H3K36me3_CASTEiJ.sortByName.sam
#samtools view H3K36me3_C57BL6J_CASTEiJ.sortByName.bam | awk '(($3 ~ "C57BL6J")&&($5 > 20)) {print $0}' | sed 's/C57BL6J_chr/chr/g' >> H3K36me3_C57BL6J.sortByName.sam
#samtools view H3K36me3_C57BL6J_CASTEiJ.sortByName.bam | awk '(($3 ~ "CASTEiJ")&&($5 > 20)) {print $0}' | sed 's/CASTEiJ_chr/chr/g' >> H3K36me3_CASTEiJ.sortByName.sam


### convert to bam
#samtools view -bt C57BL6J.fasta H3K36me3_C57BL6J.sortByName.sam > H3K36me3_C57BL6J.sortByName.bam
#samtools view -bt CASTEiJ.fasta H3K36me3_CASTEiJ.sortByName.sam > H3K36me3_CASTEiJ.sortByName.bam
#mv H3K36me3_C57BL6J.sortByName.bam H3K36me3_C57BL6J.sortByName.Qgt20.bam
#mv H3K36me3_CASTEiJ.sortByName.bam H3K36me3_CASTEiJ.sortByName.Qgt20.bam

#rm -f *.sam
