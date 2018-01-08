package org.bcgsc.edcc;

import java.io.*;

/**
 * Example
 * 
 * @author Hamid Younesy
 * @date 2013-03-20
 */
public class Alea {

    static boolean s_bIsVerbose = false;
    static boolean s_bIsViewing = false;
    
    String m_Chr = null;
    int    m_Start = -1;
    int    m_End = -1;

    /**
     * 
     * @param args
     * @throws IOException
     */
    public static void main(String[] args) throws IOException
    {/*
        if (true)
        {
            int mb = 1024*1024;
            int size = 30 * mb;
            Runtime runtime = Runtime.getRuntime();
            System.out.println("Used Memory:" + (runtime.totalMemory() - runtime.freeMemory()) / mb);
            System.out.println("Free Memory:" + runtime.freeMemory() / mb);
            long startTime = System.currentTimeMillis();
            byte[] b1 = new byte[size];
            byte[] b2 = new byte[size];
            byte[] b3 = new byte[size];
            long endTime   = System.currentTimeMillis();
            System.out.println("Allocation time: " + (endTime - startTime));
            System.out.println("Used Memory:" + (runtime.totalMemory() - runtime.freeMemory()) / mb);
            System.out.println("Free Memory:" + runtime.freeMemory() / mb);
            startTime = System.currentTimeMillis();
            for (int i = 0; i < size; ++i)
            {
                b3[i] = (byte)(128.0 * (b1[i] - b2[i]) / (b1[i] + b2[i] + 1));
            }
            endTime   = System.currentTimeMillis();
            System.out.println("Iteration time: " + (endTime - startTime));
            return;
        }
        
        */
        String sInputFasta  = null; // input .fasta file
        String sInputVCF    = null; // input .vcf or .vcf.gz file
        String sStrainName  = null; // sample name. e.g. "CBAJ"
        String sOutputFasta = null; // output .fasta file
        
        String sInputBam1   = null;
        String sInputBam2   = null;
//        String sInputFasta1 = null;
//        String sInputFasta2 = null;
//        String sOutputBam1  = null;
//        String sOutputBam2  = null;
//        String sOutputStatsDir = null;
        
        String sInputWig    = null;
        String sInputRefMap = null;
        String sOutputBedGraph = null;
        
        boolean bInsilico     = false;
        boolean bFilter       = false;
        boolean bProject      = false;
        
        // used with -view
        String sChrom       = null;
        String sStartPos    = null;
        String sEndPos      = null;

        // First get file
        for (String a : args) {
            if (a.equals("-verbose")) {
                s_bIsVerbose = true;
            } else if (a.equals("view")) {
                s_bIsViewing = true;
            } else if (a.equals("insilico")) {
                bInsilico = true;
            } else if (a.equals("filter")) {
                bFilter = true;
            } else if (a.equals("project")) {
                bProject = true;
            } else {
                String tokens[] = a.split("=");
                if (tokens.length == 2) {
                    String key = tokens[0].toLowerCase();
//                    if (key.equals("genome")) {
//                        genomeFile = tokens[1];
//                    } else if (key.equals("ref")) {
//                        ref = tokens[1];
//                    } else if (key.equals("vcf")) {
//                        vcf = tokens[1];
//                    } else if (key.equals("sample")) {
//                        sample = tokens[1];
//                    } else if (key.equals("out")) {
//                        out = tokens[1];
                    if (key.equals("--chr")) {
                        sChrom = tokens[1];
                    } else if (key.equals("--start")) {
                        sStartPos = tokens[1];
                    } else if (key.equals("--end")) {
                        sEndPos = tokens[1];
                    } else if (key.equals("--input-fasta")) {
                        sInputFasta = tokens[1];
                    } else if (key.equals("--input-vcf")) {
                        sInputVCF = tokens[1];
                    } else if (key.equals("--strain")) {
                        sStrainName = tokens[1];
                    } else if (key.equals("--output-fasta")) {
                        sOutputFasta = tokens[1];
                    } else if (key.equals("--input-bam1")) {
                        sInputBam1 = tokens[1];
                    } else if (key.equals("--input-bam2")) {
                        sInputBam2 = tokens[1];
//                    } else if (key.equals("--input-fasta1")) {
//                        sInputFasta1 = tokens[1];
//                    } else if (key.equals("--input-fasta2")) {
//                        sInputFasta2 = tokens[1];
//                    } else if (key.equals("--output-bam1")) {
//                        sOutputBam1 = tokens[1];
//                    } else if (key.equals("--output-bam2")) {
//                        sOutputBam2 = tokens[1];
//                    } else if (key.equals("--output-stats-dir")) {
//                        sOutputStatsDir = tokens[1];
                    } else if (key.equals("--input-wig")) {
                        sInputWig = tokens[1];
                    } else if (key.equals("--input-refmap")) {
                        sInputRefMap = tokens[1];
                    } else if (key.equals("--output-bedgraph")) {
                        sOutputBedGraph = tokens[1];
                    }
                }
            }
        }
        
        try {
            if (bInsilico) {
                if (checkForParameters(sInputFasta != null && sInputVCF != null && sStrainName != null)) {
                    VCFAnalyzer.createInsilicoGenome(sInputFasta, sInputVCF, sStrainName, sChrom, sOutputFasta);
                }
            } if (bFilter) {
                if (checkForParameters(sInputBam1 != null && sInputBam2 != null)) {
                    System.out.println("[ALEA] Error: Unsupported function: filter");
                    //BAMAnalyzer.filterReads(sInputBam1, sInputBam2, sInputFasta1, sInputFasta2, sOutputBam1, sOutputBam2, sOutputStatsDir);
                }
            } else if (bProject) {
                if (checkForParameters(sInputWig != null && sInputRefMap != null && sOutputBedGraph != null)) {
                    WigAnalyzer.projectWig(sInputWig, sInputRefMap, sOutputBedGraph);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        //if (genomeFile != null) {
        //    genome = IgvTools.loadGenome(genomeFile, true);
        //}
        
        if (s_bIsViewing)
        {
            int iStart = -1, iEnd = -1;

            try {
                iStart = sStartPos != null ? Integer.parseInt(sStartPos) : -1;
                iEnd = sEndPos != null ? Integer.parseInt(sEndPos) : -1;
            } catch (Exception e) {
                e.printStackTrace();
            }
            
            if (sInputFasta != null) {
                try {
                    VCFAnalyzer.viewReference(sInputFasta, sChrom, iStart, iEnd);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
    
            if (sInputVCF != null) {
                try {
                    VCFAnalyzer.viewVCF(sInputVCF, sStrainName, sChrom, iStart, iEnd);
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }
    
    static boolean checkForParameters(boolean condition)
    {
        if (!condition) {
            //TODO create a more informative error
            System.out.println("[ALEA] Error: Insufficient parameters.");
        }
        return condition;
    }
}
