package org.bcgsc.edcc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Map;
import java.util.TreeMap;

/**
 * Keeps the mapping of intervals between two reference genomes. Used for faster projection back to reference genome.
 * @author hyounesy
 *
 */
public class ReferenceMapper 
{

    /**
     * Stores the mapping between two reference genomes for one chromosome
     */
	static class ChrMap
	{
	    String          m_ChrName; // chromosome name
		ArrayList<Long> m_PosRef1; // starting position of each interval in reference genome 1
		ArrayList<Long> m_PosRef2; // starting position of each interval in reference genome 2
		ArrayList<Long> m_Lengths; // length of each interval
	};
	
	// http://docs.oracle.com/javase/6/docs/api/java/util/TreeMap.html
	// The map is sorted according to the natural ordering of its keys ...
    TreeMap<String, ChrMap> m_ReferenceMap = new TreeMap<String, ChrMap>();
    
    String	m_CurrChrName = null;
    ChrMap 	m_CurrChrMap  = null;

    /**
     * Adds a new chromosome and sets the current chomosome to it.
     * @param chrName  chromosome name
     */
	void addChr(String chrName)
	{
	    if (!setCurrentChr(chrName))
	    {
            m_CurrChrName = chrName;
	        m_CurrChrMap = new ChrMap();
	        m_CurrChrMap.m_ChrName = chrName;
	        m_CurrChrMap.m_PosRef1 = new ArrayList<Long>();
	        m_CurrChrMap.m_PosRef2 = new ArrayList<Long>();
	        m_CurrChrMap.m_Lengths = new ArrayList<Long>();
	        m_ReferenceMap.put(chrName, m_CurrChrMap);
	    }
	}
    
	/**
	 * Sets current chromosome.
     * @param chrName  chromosome name
     * @return true if the chromosome name was found
	 */
    boolean setCurrentChr(String chrName)
	{
        if (m_CurrChrName == null || !m_CurrChrName.equals(chrName))
        {
    		m_CurrChrName = chrName;
    		m_CurrChrMap = m_ReferenceMap.get(chrName);
        }
		return m_CurrChrMap != null;
	}
	
    /**
     * Adds a new mapping from an interval in reference1 to reference2.
     * Requirement: these new positions must be AFTER the last one.
     * @param pos1  starting position in reference1
     * @param pos2  starting position in reference2
     * @param len   length of the interval
     */
    void mapInterval(long pos1, long pos2, long len)
    {
        assert(m_CurrChrMap != null);
        int lastIndex = m_CurrChrMap.m_PosRef1.size() - 1;
        if (lastIndex >= 0)
        {
            long prev1 = m_CurrChrMap.m_PosRef1.get(lastIndex).longValue();
            long prev2 = m_CurrChrMap.m_PosRef2.get(lastIndex).longValue();
            long prevlen = m_CurrChrMap.m_Lengths.get(lastIndex).longValue();
            assert(prev1 + prevlen <= pos1);
            assert(prev2 + prevlen <= pos2);
            
            //if ((prev1 + prevlen == pos1) && (prev2 + prevlen == pos2)) {
            // this interval is right after the previous one. so just increment it
            if ((pos1 - prev1) == (pos2 - prev2)) {
                //m_CurrChrMap.m_Lengths.set(lastIndex, prevlen + len);
                m_CurrChrMap.m_Lengths.set(lastIndex, (pos1 - prev1) + len);
                return;
            } 
        }
        
        m_CurrChrMap.m_PosRef1.add(pos1);
        m_CurrChrMap.m_PosRef2.add(pos2);
        m_CurrChrMap.m_Lengths.add(len);
    }

    /**
     * Maps a position from reference2 to reference1.
     * @param pos2 position in reference2
     * @return  mapped position in reference1
     */
    long getPositionRef1(long pos2)
    {
        int index = Collections.binarySearch(m_CurrChrMap.m_PosRef2, pos2);
        // the index of the search key, if it is contained in the list; 
        // otherwise, (-(insertion point) - 1): 
        //   The insertion point is defined as the point at which the key would be inserted into the list:
        //     - the index of the first element greater than the key, 
        //     - or list.size() if all elements in the list are less than the specified key.
        //       Note that this guarantees that the return value will be >= 0 if and only if the key is found.

        if (index < 0) { // no exact match
            index = -index - 2;
        }
        if (index >= 0 && index < m_CurrChrMap.m_PosRef2.size())
        {
            long start2 = m_CurrChrMap.m_PosRef2.get(index).longValue();
            long start1 = m_CurrChrMap.m_PosRef1.get(index).longValue();
            long len = m_CurrChrMap.m_Lengths.get(index).longValue();
            if (start2 <= pos2 && start2 + len > pos2) 
            {
                return start1 + (pos2 - start2);
            }
            return -1;
        }
        
        return -1;
    }

    /**
     * Maps a position from reference1 to reference2.
     * @param pos1 position in reference1
     * @return  mapped position in reference2
     */
    long getPositionRef2(long pos1)
    {
        //TODO
        //Collections.binarySearch(m_CurrChrMap.m_PosRef1, pos1)
        return -1;
    }
    
    /**
     * Writes the mapping to a .refmap file.
     * @param filename   output filename
     * @return  true if successful
     */
    boolean writeToFile(String filename)
    {
        try {
            System.out.println("\n[ALEA:ReferenceMapper] Creating refmap file at " + filename);
            PrintWriter mapWriter = new PrintWriter(new BufferedWriter(new FileWriter(filename)));
            for (Map.Entry<String, ChrMap> entry : m_ReferenceMap.entrySet()) 
            {
                mapWriter.println(">" + entry.getKey());
                ChrMap chrMap = entry.getValue();
                int size = chrMap.m_Lengths.size();
                
                for (int i = 0; i < size; ++i) {
                    mapWriter.print(chrMap.m_PosRef1.get(i).longValue() + "\t");
                    mapWriter.print(chrMap.m_PosRef2.get(i).longValue() + "\t");
                    mapWriter.println(chrMap.m_Lengths.get(i).longValue());
                }
            }
            mapWriter.close();
            System.out.println("\n[ALEA:ReferenceMapper] Done creating refmap file.");
            return true;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return false;
    }
    
    /**
     * Reads the mapping from a .refmap file.
     * @param filename   input file path
     * @return true if successful
     */
    boolean readFromFile(String filename)
    {
        try {
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            System.out.println("\n[ALEA:ReferenceMapper] Reading refmap file at " + filename);
            String strLine = null;
            while ((strLine = reader.readLine()) != null)
            {
                if (strLine.startsWith(">"))
                {// add the new chromosome
                    addChr(strLine.substring(1, strLine.length()));
                } else {
                    // read tab delimited pos1, pos2, len
                    String tokens[] = strLine.split("\t");
                    assert(tokens != null && tokens.length == 3); //TODO error check
                    mapInterval(Long.parseLong(tokens[0]), Long.parseLong(tokens[1]), Long.parseLong(tokens[2]));
                }
            }
            reader.close();
            System.out.println("\n[ALEA:ReferenceMapper] Done reading refmap file.");
            return true;
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        return false;
    }
}
