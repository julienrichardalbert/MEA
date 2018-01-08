package org.bcgsc.edcc;

import htsjdk.tribble.Feature;

import java.util.List;

import org.broad.igv.feature.genome.FastaIndexedSequence;
import org.broad.igv.feature.genome.FastaUtils;
import org.broad.igv.feature.genome.Genome;
import org.broad.igv.track.Track;
import org.broad.igv.track.TrackLoader;
import org.broad.igv.util.ResourceLocator;
import org.broad.igv.variant.VariantTrack;
import org.broad.igv.variant.vcf.VCFGenotype;
import org.broad.igv.variant.vcf.VCFVariant;

/**
 * Example
 * 
 * @author Hamid Younesy
 * @date 2013-03-20
 */
public class VCFAnalyzer {

    static boolean s_bIsVerbose = false;
    
    /**
     * 
     * @param refPath
     * @param vcfPath
     * @param strainName
     * @param outPath
     * @param chr
     * @throws Exception
     */
    public static void createInsilicoGenome(String refPath, String vcfPath, String strainName, String chr, String outPath) throws Exception
    {
        //TODO: create a fasta index if .fai doesn't exist
        //FastaUtils.createIndexFile(refPath, refPath + ".fai");
        
        FastaIndexedSequence sequence = new FastaIndexedSequence(refPath);
        List<String> chrNames = sequence.getChromosomeNames();
        
        Genome genome = null;
        VariantTrack variantTrack = (VariantTrack) (new TrackLoader()).load(new ResourceLocator(vcfPath), genome).get(0);

        if (outPath == null) {
            outPath = refPath + strainName + ".fasta";
        }
        
        ReferenceMapper refMap = new ReferenceMapper();
        
        FastaWriter fastaWriter = new FastaWriter(outPath);
        for (String chrName : chrNames)
        {
            if (chr != null) {
                if (!chr.equals(chrName)) {
                    continue; // skip
                }
            }
            
            long iNumRefBases = 0;   // total number of reference bases processed
            long iNumIndelBases = 0; // number of bases inserted or deleted

            System.out.println("\n[ALEA] Processing chr=" + chrName + " length = " + sequence.getChromosomeLength(chrName));

            refMap.addChr(chrName);
            
            int chrLength = sequence.getChromosomeLength(chrName);
            
            iNumRefBases += chrLength;
            
            //TODO: couldn't get this from FastaIndexedSequence. need to parse it directly from the input .fai file. 
            int iNumRowBases = 80;//chrName.equals("MT") ? 70 : 80;
            fastaWriter.startChromosome(chrName, iNumRowBases);
            
            int chunkSize = 1 << 20; // size of the interval chunks to read from the reference
            int refStart = 0;   // start position of the current interval chunk in the reference game
            int refOffset = 0;  // offset from start of chunk of bytes that are already processed
            
            int lastVCFStart = -1;
            
            int lastProgress = 0;
            while (refStart < chrLength)
            {
                // make sure refEnd doesn't exceed chromosome length
                int refEnd = Math.min(refStart + chunkSize, chrLength);

                if (refStart - lastProgress  > chrLength / 20) {
                    System.out.printf("[%d%%]", ((long)refStart) * 100 / chrLength);
                    lastProgress = refStart;
                }

                // read reference sequence to a buffer
                byte[] seqBuffer = sequence.getSequence(chrName, refStart, refEnd);

                // get the vcf features within the current reference interval
                List<Feature> featuresList = variantTrack.getFeatures(chrName, refStart, refEnd - 1);
                
                for (Feature f : featuresList)
                {
                    VCFVariant vf = ((VCFVariant) f);
                    VCFGenotype genotype = (VCFGenotype)(vf.getGenotype(strainName));
                    String genotypeString = genotype.getGenotypeString();
                    
                    if (genotypeString.equals(".")) {
                        continue; // skip. no variant for this sample
                    }
                    
                    if (vf.getStart() == lastVCFStart) {
                        continue; // skip. this variant was already processed (in last chunk)
                    }
                    lastVCFStart = vf.getStart();
                    
                    assert(vf.getStart() >= refStart + refOffset);
                    
                    if (vf.getStart() > refStart + refOffset) {
                        // write the sequence up to this variant
                        int len = vf.getStart() - (refStart + refOffset);
                        refMap.mapInterval(refStart + refOffset, fastaWriter.getNumWrittenBases(), len);
                        fastaWriter.writeBases(seqBuffer, refOffset, len);
                        refOffset = (vf.getStart() - refStart); //variantContext.getType() == VariantContext.Type.INDEL
                    }
                    
                    // write the variant
                    int variantLength = genotypeString.indexOf('/');
                    assert(variantLength >= 0);
                    if (variantLength > 0) {
                        fastaWriter.writeBases(genotypeString.getBytes(), 0, variantLength);
                    }
                    
                    // also skip the mutated reference region
                    refOffset += (vf.getEnd() - vf.getStart()); 

                    String referenceString = vf.getVariantContext().getReference().getBaseString();
                    iNumIndelBases += variantLength - referenceString.length();
                    
                    if (s_bIsVerbose) 
                    {
                        System.out.println("[" + chrName + ":" + refStart + "-" + refEnd + "]");
                        System.out.print(vf.toString() + "\t");
                        System.out.print("[" + referenceString + "]\t");
                        System.out.print("[" + genotypeString + "]\t"); //TTTG/TTTG   or  .  or  /
                        System.out.println();
                    }
                }
                
                // done with the variants in this chunk.
                if (refStart + refOffset < refEnd)
                {
                    //write what's remaining of the sequence
                    int len = refEnd - (refStart + refOffset);
                    refMap.mapInterval(refStart + refOffset, fastaWriter.getNumWrittenBases(), len);
                    fastaWriter.writeBases(seqBuffer, refOffset, len);
                    refStart = refEnd;
                } else {
                    // the last variant went over the current chunk.
                    refStart += refOffset;
                }
                
                refOffset = 0;
            }
            System.gc();
            
            System.out.println("\n[ALEA] Read: " + iNumRefBases 
                                + "  Wrote: " + fastaWriter.getNumWrittenBases()
                                + "  Difference: " + iNumIndelBases);
        } //for (String chrName : chrNames)
        fastaWriter.close();
        
        FastaUtils.createIndexFile(outPath, outPath + ".fai");
        refMap.writeToFile(outPath + ".refmap");
    }

    /**
     * Outputs the sequence for the query interval to console.  Coordinates are 0 based.
     * Outputs the list of chromosomes when null passed to chr.
     * @param refPath   path to reference fasta file
     * @param chr       chromosome name
     * @param start     interval start
     * @param end       interval end
     * @throws Exception
     */
    public static void viewReference(String refPath, String chr, int start, int end) throws Exception 
    {

        //FastaUtils.createIndexFile(fastaPath, fastaPath + ".fai");
        FastaIndexedSequence sequence = new FastaIndexedSequence(refPath);
        
        if (chr == null)
        {
            List<String> chrNames = sequence.getChromosomeNames();
            for (String chrName : chrNames) {
                System.out.println(chrName + " length = " + sequence.getChromosomeLength(chrName));
            }
        } else {
            byte[] seqBytes = sequence.getSequence(chr, start, end);
            if (seqBytes != null) {
                System.out.println(new String(seqBytes));
            }
        }
    }
    
    public static String getSequence(String refPath, String chr, int start, int end) throws Exception 
    {
        FastaIndexedSequence sequence = new FastaIndexedSequence(refPath);
        byte[] seqBytes = sequence.getSequence(chr, start, end);
        if (seqBytes != null) {
            return new String(seqBytes);
        }
        return null;
    }

    /**
     * 
     * @param vcfPath
     * @param sampleName
     * @param chr
     * @param start
     * @param end
     * @throws Exception
     */
    public static void viewVCF(String vcfPath, String sampleName, String chr, int start, int end) throws Exception 
    {
        // TestUtils.createIndex(vcfPath);
        Genome genome = null;
        List<Track> trackList = (new TrackLoader()).load(new ResourceLocator(vcfPath), genome);
        System.out.println(trackList.size());
        Track newTrack = trackList.get(0);
        VariantTrack variantTrack = (VariantTrack) newTrack;

        if (chr == null) {
            chr = "1";
        }
        if (start == -1) {
            start = 3000000;
        }
        
        if (end == -1) {
            end = start + 10000;
        }
        
        String sampleNames[] = {"C57BL6NJ", "BALBcJ", "CBAJ", "AKR", "LPJ", "NODShiLtJ", "NZO", "WSBEiJ",
                                "C3HHeJ", "CASTEiJ", "PWKPhJ", "DBA2J", "FVB_NJ", "129S1", "Spretus", "AJ"};
        if (sampleName == null) {
            sampleName = sampleNames[2];
        }
        
        List<Feature> featuresList = variantTrack.getFeatures(chr, start, end);
        
        System.out.println(variantTrack.getSample()); //small.snps.vcf
        System.out.println("[" + chr + ":" + start + "-" + end + "]  #features = " + featuresList.size());
        
        for (Feature f : featuresList) {
            // System.out.println(((VCFVariant)f).getSampleNames().toString());
            VCFVariant vf = ((VCFVariant) f);
            
            VCFGenotype genotype = null;
            if (sampleName != null) {
                genotype = (VCFGenotype)(vf.getGenotype(sampleName));
                if (genotype != null) {
                    try {
                    if (genotype.getGenotypeString().equals(".")) {
                        // no variant
                        continue;
                    }
                    } catch (java.lang.NullPointerException e) {
                        genotype = null;
                    }
                }
            }

            //System.out.println(vf.getVariantContext().getSource()); //Unknown
            System.out.print(vf.toString() + "\t"); //VCFVariant[1:3000748-3000748]
            System.out.print(vf.getType() + "\t"); //SNP or INDEL 
            System.out.print(vf.getVariantContext().getReference().getBaseString() + "\t"); //CTTACA  same as getDisplayString() ?
            
            if (genotype != null) {
                System.out.print(genotype.getGenotypeString() + "\t"); //TTTG/TTTG   or  .  or  /
            }
            System.out.print(vf.getVariantContext().getAlleles() + "\t"); //[GGG*, CAGGG, -, G]
            

            //System.out.print(vf.getAttributes().toString()); // the INFO field in vcf: {AF1=0, AC=2, VDB=0.0001, DP=551, DP4=[209, 93, 0, 1], AC1=0, MQ=46, AN=32}
            //System.out.println(vf.getVariantContext().getGenotypes().get(0).getGenotypeString()); //CAGGG/CAGGG
            //System.out.println(vf.getVariantContext().getAlleleStringWithRefPadding(alleleList.get(0)));
            
            //List<Allele> alleleList = genotype.getAlleles();
            //System.out.print(((alleleList.size() > 1) ? new String(alleleList.get(0).getBases()) : alleleList.toString()) + "\t");
            //System.out.print(genotype.getAttributes()); //{SP=0, FI=1}
            
            // read the reference sequence
            //if (sequence != null) {
            //    byte[] seqBytes = sequence.getSequence(chr, vf.getStart(), vf.getEnd());
            //    System.out.print("ref:[" + (seqBytes != null ? new String(seqBytes) : "") + "]\t");
            //}
            
            System.out.println();
        }
        
        System.out.println("Done");
    }

}
