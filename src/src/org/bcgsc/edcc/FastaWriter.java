package org.bcgsc.edcc;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;


/**
 * Simple class to output a fasta file using a buffered writer.
 * 
 * @author Hamid Younesy
 * @date 2013-03-20
 */
class FastaWriter
{
    BufferedOutputStream m_OutStream = null;
    int m_iNumBasesPerLine = 0;  // number of bases to write at each line (e.g. 60 or 80)
    int m_iCurrLineBases   = 0;  // the number of bases written in the current line
    long m_iNumWrittenBases = 0; // total number of bases written so far
    long m_iLineNumber = 1;     // current line number
    byte[] m_NL = {'\n'};       // new line character
          
    /**
     * Constructor
     * @param outpath output filename
     */
    FastaWriter(String outpath)
    {
        try {
            System.out.println("\n[FastaWriter] Creating " + outpath);
            FileOutputStream fos = new FileOutputStream(new File(outpath));
            m_OutStream = new BufferedOutputStream(fos);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * closes the output stream writer
     */
    void close()
    {
        try {
            System.out.println("\n[FastaWriter] Writing to file...");
            m_OutStream.close();
            System.out.println("\n[FastaWriter] Done");
        } catch (Exception e1) {
            // TODO Auto-generated catch block
            e1.printStackTrace();
        }
        m_OutStream = null;
    }

    /**
     * @return number of bases written for current chromosome
     */
    public long getNumWrittenBases() {
        return m_iNumWrittenBases;
    }
    
    /**
     * @return current output file line number.
     */
    public long getLineNumber() {
        return m_iLineNumber;
    }

    /**
     * Starts writing the bases for a new chromosome
     * @param chrName   chromosome name
     * @param numBasesPerLine  number of bases to write per each line
     */
    void startChromosome(String chrName, int numBasesPerLine)
    {
        try {
            if (m_iCurrLineBases != 0) {
                m_OutStream.write(m_NL);
        		m_iLineNumber++;
                m_iCurrLineBases = 0;
                m_iNumWrittenBases = 0;
            }
            
            m_OutStream.write((">" + chrName + "\n").getBytes());
            m_iLineNumber++;
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        m_iNumBasesPerLine = numBasesPerLine;
    }
    
    /**
     * Writes a sequence of bases
     * @param buffer    input byte array
     * @param offset    offset from the start of the buffer
     * @param len       number of characters to write starting from offset
     */
    void writeBases(byte[] buffer, int offset, int len)
    {
        try {
            m_iNumWrittenBases += len;

            while (m_iNumBasesPerLine - m_iCurrLineBases <= len)
            {
                m_OutStream.write(buffer, offset, m_iNumBasesPerLine - m_iCurrLineBases);
                m_OutStream.write(m_NL);
                m_iLineNumber++;
                
                len -= (m_iNumBasesPerLine - m_iCurrLineBases);
                offset +=  (m_iNumBasesPerLine - m_iCurrLineBases);
                m_iCurrLineBases = 0;
            }

            if (len > 0) {
                m_OutStream.write(buffer, offset, len);
                m_iCurrLineBases += len;
            }
        } catch (Exception e) {
           e.printStackTrace();
        }
    }
}
