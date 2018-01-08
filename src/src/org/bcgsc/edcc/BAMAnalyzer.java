package org.bcgsc.edcc;

import htsjdk.samtools.SAMFileWriterImpl;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.broad.igv.feature.genome.FastaIndexedSequence;
import org.broad.igv.sam.Alignment;
import org.broad.igv.sam.PicardAlignment;
import org.broad.igv.sam.reader.SAMReader;

public class BAMAnalyzer {
    
//    String markName = "H3K36me3";
//    String refName1 = "mm9";
//    String refName2 = "CASTEiJ";
//    String xpDir = "q?"; // experiment directory e.g. q1-q2_gt_10 meaning q1 - q2 > 10
//    String inputBam1  = "/projects/edcc/Allele-specific/mouse/site_specific/" + markName + "/" + markName + "_" + refName1 + ".sort-n.bam";
//    String inputBam2  = "/projects/edcc/Allele-specific/mouse/site_specific/" + markName + "/" + markName + "_" + refName2 + ".sort-n.bam";
//    String outputBam1  = "/projects/edcc/Allele-specific/mouse/site_specific/" + markName + "/" + xpDir + "/" + markName + "_" + refName1 + ".filtered.bam";
//    String outputBam2  = "/projects/edcc/Allele-specific/mouse/site_specific/" + markName + "/" + xpDir + "/" + markName + "_" + refName2 + ".filtered.bam";
//    String inputFasta1  = "/projects/edcc/Allele-specific/mouse/processed/mm9_build37_mouse/mm9_build37_mouse.fasta";
//    String inputFasta2  = "/projects/edcc/Allele-specific/mouse/processed/mm9_" + refName2 + "/mm9_" + refName2 + ".fasta";
//  String outputStatsDir = "/projects/edcc/Allele-specific/mouse/site_specific/" + markName + "/" + markName + "_" + refName1 + "_" + refName2;
    
    /**
     * Filter the reads of two input bam files based on mapping quality.
     * Requires the input bam files to be sorted by read names and have identical reads.
     * @param inputBam1 filename of the first input bam file
     * @param inputBam2 filename of the second input bam file
     * @param inputFasta1 filename of the first reference genome file (.fasta)
     * @param inputFasta2 filename of the second reference genome file (.fasta)
     * @param outputBam1 output filename for the first filtered bam file
     * @param outputBam2 output filename for the second filtered bam file
     * @param outputStatsDir directory to output alignment statistics
     * @throws Exception
     */
    static public void filterReads(String inputBam1, String inputBam2, 
                                  String inputFasta1, String inputFasta2,
                                  String outputBam1, String outputBam2, 
                                  String outputStatsDir) throws Exception
    {
        //TODO: verify files exist, verify output writable
        //String snpPath   = "/projects/edcc/Allele-specific/mouse/VCF_SNPs/NCBIm37/mgp.v2.snps.annot.reformat.vcf.gz";
        //String indelPath = "/projects/edcc/Allele-specific/mouse/VCF_SNPs/NCBIm37/mgp.v2.indels.annot.reformat.vcf.gz";
        
        boolean bOutputStats = (outputStatsDir != null); // whether to create output stat files
        boolean bOutputBam = (outputBam1 != null && outputBam2 != null); // whether to create filtered output bam files
        
        boolean bPrintReadsInRange = false; // whether to show aligned reads in the range
        long rangeStart = 98503099;
        long rangeEnd  =  98503138;
        String rangeChr = "2";
        
        //chr2:
        
        // start analyzing
        System.out.println("\n[BAMAnalyzer.filterReads] started");

        // prepare BAM readers TODO: check for IO errors
        SAMReader samReader1 = new SAMReader(inputBam1, false);
        Iterator<PicardAlignment> alignIter1 = samReader1.iterator();

        SAMReader samReader2 = new SAMReader(inputBam2, false);
        Iterator<PicardAlignment> alignIter2 = samReader2.iterator();

        // prepare BAM writers TODO: check for IO errors
        SAMFileWriterImpl writer1 = null;
        SAMFileWriterImpl writer2 = null;
        if (bOutputBam) {
            File outFile1 = new File(outputBam1);
            outFile1.delete();
            writer1 = new BAMFileWriter(new BufferedOutputStream(new FileOutputStream(outFile1)), null);
            writer1.setHeader(samReader1.getHeader());
            System.out.println("\n[BAMAnalyzer.filterReads] writing filtered bam file at " + outputBam1);

            File outFile2 = new File(outputBam2);
            outFile2.delete();
            writer2 = new BAMFileWriter(new BufferedOutputStream(new FileOutputStream(outFile2)), null);
            writer2.setHeader(samReader2.getHeader());
            System.out.println("\n[BAMAnalyzer.filterReads] writing filtered bam file at " + outputBam2);
        }

        int count = 0; // number of reads processed so far
        int[][] mapQ = new int[256][256]; // #reads per parent mapping quality
        String[][] mapSeq = new String[256][256]; // sample mapped reads
        
        // prepare fasta readers (used for debug output stats)
        FastaIndexedSequence fastaSeq1 = inputFasta1 != null ? new FastaIndexedSequence(inputFasta1) : null;
        FastaIndexedSequence fastaSeq2 = inputFasta2 != null ? new FastaIndexedSequence(inputFasta2) : null;
        
        boolean bSkipQ0 = false; // skip reads with mapping quality 0.
        int minQualityDiff = 10;
        
        while (alignIter1.hasNext()) 
        {
            assert(alignIter2.hasNext());
            SamAlignment align1 = (SamAlignment)alignIter1.next();
            SamAlignment align2 = (SamAlignment)alignIter2.next();
            assert(align1.getReadName().equals(align2.getReadName()));
            assert(align1.getReadSequence().equals(align2.getReadSequence()));
            int q1 = align1.getMappingQuality();
            int q2 = align2.getMappingQuality();
            
            boolean bPassed = false;
            if (!bSkipQ0 || (q1 != 0 && q2 != 0))  
            {  
                //TODO: should we accept read1 when read2.chr=* ?
                if (q1 >= q2 + minQualityDiff) { // alignments which have a higher mapping quality to ref1
                    writer1.addAlignment(align1.getRecord());
                    bPassed = true;
                } else if (q2 >= q1 + minQualityDiff) { // alignments which have a higher mapping quality to ref2
                    writer2.addAlignment(align2.getRecord());
                    bPassed = true;
                }
                
                /*
                if (q1 == 37 && q2 == 37) 
                { // alignments where q1=q2=37
                    if (writer1 != null) {
                        writer1.addAlignment(align1.getRecord());
                        bPassed = true;
                    }
                    
                    if (writer2 != null) {
                        writer2.addAlignment(align2.getRecord());
                        bPassed = true;
                    }
                }
                */
            }
            
            if (bPassed && bPrintReadsInRange)
            {
                boolean b1 = isIntersectingInterval(align1, rangeChr, rangeStart, rangeEnd);
                boolean b2 = isIntersectingInterval(align2, rangeChr, rangeStart, rangeEnd);
                if (b1 || b2) {
                    System.out.printf("\n[%s] Q1:%2d %s:%d-%d  ", b1?"X":" ", q1, align1.getChr(), align1.getAlignmentStart(), align1.getAlignmentEnd());
                    System.out.printf("\n[%s] Q2:%2d %s:%d-%d\n", b2?"X":" ", q2, align2.getChr(), align2.getAlignmentStart(), align2.getAlignmentEnd());
                }
            }
            
            if (bOutputStats)
            {  // record the actual sequence it for only the first few alignments per quality.
                mapQ[q1][q2]++;
                if (mapQ[q1][q2] <= 3)
                {
                    if (mapSeq[q1][q2] == null) {
                        mapSeq[q1][q2] = "";
                    }
                    
                    mapSeq[q1][q2] +=     "read: " + align1.getReadSequence() + "\n";
                    if (fastaSeq1 != null) {
                        byte[] seqBytes = fastaSeq1.getSequence(align1.getChr(), align1.getAlignmentStart(), align1.getAlignmentEnd());
                        mapSeq[q1][q2] += "ref1: " + (seqBytes != null ? new String(seqBytes) : "")
                                + "  " + align1.getChr() + ":" + align1.getAlignmentStart() + "-" + align1.getAlignmentEnd() + "\n";
                    }
                    
                    if (fastaSeq2 != null) {
                        byte[] seqBytes2 = fastaSeq2.getSequence(align2.getChr(), align2.getAlignmentStart(), align2.getAlignmentEnd());
                        mapSeq[q1][q2] += "ref2: " + (seqBytes2 != null ? new String(seqBytes2) : "")
                                + "  " + align2.getChr() + ":" + align2.getAlignmentStart() + "-" + align2.getAlignmentEnd() + "\n";
                    }
                    
                    mapSeq[q1][q2] += "-----------------------------------------\n";
                }
            }
            
            if (count++ == 1000000) {
                System.out.print(".");
                count = 0;
            }
        }
        
        if (writer1 != null) {
            writer1.close();
        }
        
        if (writer2 != null) {
            writer2.close();
        }
        
        if (bOutputStats)
        {
            // output the comparison of mapping quality 
            String mapqPath = outputStatsDir + "/quality_matrix.txt";
            System.out.println("\n[BAMAnalyzer.filterReads] Started writing the mapping quality matrix at " + mapqPath);
            PrintWriter mapqWriter = new PrintWriter(new BufferedWriter(new FileWriter(mapqPath)));
            for (int i = 0; i < 256; ++i) {
                for (int j = 0; j < 256; ++j) {
                    mapqWriter.print(mapQ[i][j] + " ");
                }
                mapqWriter.println();
            }
            mapqWriter.close();
            System.out.println("\n[BAMAnalyzer.filterReads] Done writing the mapping quality matrix at " + mapqPath);

            // output sample reads per mapping quality
            String mapseqPath = outputStatsDir + "/sample_quality_reads.txt";
            System.out.println("\n[BAMAnalyzer.filterReads] Started writing the sample reads per quality at " + mapseqPath);
            PrintWriter mapSeqWriter = new PrintWriter(new BufferedWriter(new FileWriter(mapseqPath)));
            for (int i = 0; i < 256; ++i) {
                for (int j = 0; j < 256; ++j) {
                    if (mapSeq[i][j] != null) {
                        mapSeqWriter.println("("+ i + "," + j +"): ");
                        mapSeqWriter.println(mapSeq[i][j]);
                    }
                }
            }
            mapSeqWriter.close();
            System.out.println("\n[BAMAnalyzer.filterReads] Done writing the sample reads per quality at " + mapseqPath);
        }

        System.out.println("\n[BAMAnalyzer.filterReads] Done");
    }
    
    /**
     * checks if an aligned read intersects the specified interval.
     * Mainly used for debugging purpose.
     * 
     * @param align         aligned read
     * @param chr           chromosome name
     * @param rangeStart    interval start
     * @param rangeEnd      interval end
     * @return              true if the read intersects the interval
     */
    private static boolean isIntersectingInterval(SamAlignment align, String chr, long start, long end)
    {
        return  align.getChr().equals(chr)
             && (   (align.getAlignmentStart() >= start && align.getAlignmentStart() <= end)
                 || (align.getAlignmentEnd() >= start && align.getAlignmentEnd() <= end));
    }
    
    static public void viewBam(String inpath, String outpath) throws Exception 
    {
        //TODO skip igv and use net.sf.samtools.SAMRecord directly.
        
        //BAMFileReader bamreader = new BAMFileReader(new File(inpath));
        //bamreader.getSequenceNames();
        //bamreader.getHeader();
        //CloseableIterator<Alignment> bamiter = bamreader.query(chr, start, end, true);
        //SAMFileHeader samHeader = samReader.getHeader();

        //String bamPath2  = "/projects/edcc/Allele-specific/mouse/site_specific/H3K36me3/H3K36me3_CASTEiJ.sort-n.bam";
        //String statsPath = "/projects/edcc/Allele-specific/mouse/site_specific/H3K36me3/H3K36me3_mm9_CASTEiJ";
        
        //String bamPath1  = "/projects/edcc/Allele-specific/mouse/site_specific/H3K27me3_CB2/H3K27me3_CB2_mm9.sort-n.bam";
        //String bamPath2  = "/projects/edcc/Allele-specific/mouse/site_specific/H3K27me3_CB2/H3K27me3_CB2_CASTEiJ.sort-n.bam";
        //String statsPath = "/projects/edcc/Allele-specific/mouse/site_specific/H3K27me3_CB2/H3K27me3_CB2_mm9_CASTEiJ";
        
        SAMReader samReader = new SAMReader(inpath, false);
        Iterator<Alignment> alignIter = samReader.iterator();
        List<SamAlignment> alignList = new ArrayList<SamAlignment>();
        int counter = 0;
        
        while (alignIter.hasNext()) 
        {
            SamAlignment align = (SamAlignment)alignIter.next();
            alignList.add(align);
            
            int q = align.getMappingQuality();
            float score = align.getScore();
            
            boolean unmapped = align.getRecord().getReadUnmappedFlag();
            
            if (!unmapped) {
                System.out.println(align.getCigarString() + "\t" + q + "\t" + score);
                counter++ ;
            } else {
                System.out.println(align.getCigarString() + "\t" + q + "\t" + score);
            }
            
            if (align.getMate() != null) 
            {
                counter++; //WHY?!!
            }
            
            if (counter == 10)
            {
                counter = 0;
            }
        }

        if (outpath == null)
            return;
        
        boolean isBam = outpath.endsWith(".bam");

        File outFile = new File(outpath);
        outFile.delete();
        OutputStream outputStream = new BufferedOutputStream(new FileOutputStream(outFile));
        SAMWriter writer = new SAMWriter(samReader.getHeader());
        writer.writeToStream(outputStream, alignList, isBam);
        
        //checkRecordsMatch(shi.alignments, outFile, inpath);
    }


}
