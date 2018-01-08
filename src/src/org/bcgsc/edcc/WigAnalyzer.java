package org.bcgsc.edcc;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import org.broad.igv.track.TrackType;
import org.broad.igv.tools.parsers.DataConsumer;
import org.broad.igv.tools.parsers.ToolsWiggleParser;

/**
 * Utility functions to analyze wig files
 * @author hyounesy
 *
 */
public class WigAnalyzer 
{
    
    /**
     * Implements DataConsumer interface. Projects data back to reference before writing as a bedgraph 
     * inspired from org.broad.igv.tools.converters.WigToBed
     * @author hyounesy
     */
    static class WigMapper implements DataConsumer
    {
        PrintWriter bedGraphWriter;
        ReferenceMapper refMap;

        String currChr = null;
        int currStart = -1;
        int currEnd = -1;
        float currData = -1;
        
        WigMapper(String refMapPath, String outBedPath)
        {
            refMap = new ReferenceMapper();
            refMap.readFromFile(refMapPath);
            try {
                bedGraphWriter = new PrintWriter(new BufferedWriter(new FileWriter(outBedPath)));
                System.out.println("[ALEA:WigMapper] started writing projected track at " + outBedPath);
                bedGraphWriter.println("track type=bedGraph");
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        
        @Override
        public void setType(String type) {
            // TODO Auto-generated method stub
            //System.out.println("setType "+type);
            
        }

        @Override
        public void addData(String chr, int start, int end, float[] data, String name)
        {
            if (!refMap.setCurrentChr(chr)) {
                System.out.println("[ALEA:WigMapper] Error: unable to find chromosome " + chr);
            }
            
            for (int i = start; i < end; i++)
            {
                long pos1 = refMap.getPositionRef1(start);
                if (pos1 != -1) {
                    addSingleData(chr, (int)pos1, data[0]);
                    //bedWriter.println(chr + "\t" + pos2 + "\t" + (pos2+1) + "\t" + data[0]);
                }
            }
            // TODO Auto-generated method stub
        }
        
        private void addSingleData(String chr, int pos, float data)
        {
            if (currChr != null && currChr.equals(chr) && pos == currEnd && data == currData)
            {
                currEnd++;
                return;
            }
            
            if (currChr != null && currStart >= 0)
            {
                bedGraphWriter.println(currChr + "\t" + currStart + "\t" + currEnd + "\t" + currData);
            }
            
            currChr = chr;
            currStart = pos;
            currEnd = currStart + 1;
            currData = data;
            
        }

        @Override
        public void parsingComplete() {
            addSingleData("", -1, 0);
            bedGraphWriter.close();
            System.out.println("[ALEA:WigMapper] done writing projected track.");
        }

        @Override
        public void setTrackParameters(TrackType trackType, String trackLine,
                String[] trackNames) {
            // TODO Auto-generated method stub
            //System.out.println("setTrackParameters");
        }

        @Override
        public void setTrackParameters(TrackType trackType, String trackLine,
                String[] trackNames, boolean b) {
            // TODO Auto-generated method stub
            //System.out.println("setTrackParameters");
        }

        @Override
        public void setSortTolerance(int tolerance) {
            // TODO Auto-generated method stub
            //System.out.println("setSortTolerance");
        }

        @Override
        public void setAttribute(String key, String value) {
            // TODO Auto-generated method stub
            //System.out.println("setAttribute");
        }
    }
    
//  String markName = "H3K36me3";
//  String refName = "CASTEiJ";
//  String xpDir = "q37";//"q1-q2_gt_10"; // experiment directory e.g. q1-q2_gt_10 meaning q1 - q2 > 10
//  String wigPath     = "/projects/edcc/Allele-specific/mouse/site_specific/" + markName + "/" + xpDir + "/" + markName + "_" + refName + ".filtered.sorted.nodup.wig";
//  String outBedPath  = "/projects/edcc/Allele-specific/mouse/site_specific/" + markName + "/" + xpDir + "/" + markName + "_" + refName + ".filtered.sorted.nodup.projected.bedGraph";
//  String refMapPath = "/projects/edcc/Allele-specific/mouse/processed/mm9_" + refName + "/mm9_" + refName + ".fasta.refmap";

    /**
     * Projects a wig file to reference genome using the .refmap index
     * @param inputWig
     * @param inputRefMap
     * @param outputBedGraph
     */
    static public void projectWig(String inputWig, String inputRefMap, String outputBedGraph)
    {
        //package org.broad.igv.tools.converters.WigToBed;                
        //ReferenceMapper refMap = new ReferenceMapper();
        //refMap.readFromFile(refMapPath);
        //TrackLoader trackLoader = new TrackLoader();
        //ResourceLocator locator = new ResourceLocator(wigPath);
        
        WigMapper wigMapper = new WigMapper(inputRefMap, outputBedGraph);
        ToolsWiggleParser wigParser = new ToolsWiggleParser(inputWig, wigMapper, null);
        try {
            wigParser.parse();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        //WiggleDataset wigDataSet = wigParser.parse();
        //TrackProperties props = wigDataSet.getTrackProperties();
        
    }
}
