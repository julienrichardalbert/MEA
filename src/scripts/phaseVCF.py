#!/usr/bin/python

import sys
import os

def phaseVCF(inputHaps, inputVCF, outputVCF):
    print 'creating dictionary of the phased SNPs from ' + inputHaps
    hapsFile = open(inputHaps,'rU')
    chromDict={} # {chrom, {position,snp info}}
    
    for line in hapsFile:
        line=line.strip() #remove \n
        #e.g. line: 'X rs2377582 77582 A G 0 1'
        hap = line.split(' ')
        
        hapsDict=chromDict.get(hap[0])
        if not hapsDict: # first time adding the chromosome
            chromDict[hap[0]] = {}
            hapsDict = chromDict[hap[0]]
            print 'Adding chromosome ' + hap[0]
        
        hapsDict[hap[2]]=hap
        
    print 'total: ' + str(len(chromDict))
    
    hapsFile.close();
    

    print 'Reading vcf file from ' + inputVCF
    iFormat = 8 #index for the Format field
    iSample = 9 #index for the sample field
    iNumSamples = 1
    # process the vcf line by line
    # vcf 4.0 format: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
    vcfFile = open(inputVCF,'rU')
    vcfOut = open(outputVCF,'w')
    
    #todo: make sure that the input VCF has exactly 9 fields. otherwise this will not work.
    
    dictProgress={}
    for line in vcfFile:
        line=line.strip() #remove \n
        if line.startswith('##'):
            vcfOut.write(line + '\n');
            #meta information lines. write as is
            pass
            #print line
        elif line.startswith('#'):
            # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	JOC224
            # 0         1       2       3       4       5       6       7       8       9
            header=line.split('\t')
            assert(iFormat==header.index('FORMAT'))
            header[iSample]='hap1'
            header.append('hap2')
            vcfOut.write('\t'.join(header) + '\n')
            pass
        else:
            #e.g.: 1	103250	.	G	A	172.77	.	AC=1;AF=0.500	GT:AD:DP:GQ:PL	0/1:9,12:20:46:201,0,46
            variant=line.split('\t')
            assert(len(variant) == iSample + iNumSamples)
            
            if not dictProgress.get(variant[0]):
                print 'Processing chromosome '+variant[0]
                dictProgress[variant[0]] = True
            
            hapsDict=chromDict.get(variant[0])
            found_phased=False
            if hapsDict:
                hap=hapsDict.get(variant[1])
                if hap: # found phased haplotype
                    assert(hap[3]==variant[3] and hap[4] == variant[4])
                    #print line,'\n', hap
                    found_phased=True
                    sample = variant[iSample].split(':')
                    sample[0]=hap[5]+'/'+hap[5]
                    variant[iSample]=':'.join(sample) # update the field for hap1
                    sample[0]=hap[6]+'/'+hap[6]
                    variant.append(':'.join(sample))  # add a second field for hap2
                    
                    # print indels
                    if len(hap[3]) > 1 or len(hap[4]) > 1 :
                        print '\t'.join(variant)
                    
                    #write two fields, one for each allele
                    vcfOut.write('\t'.join(variant) +'\n')
                    #print '\n'
                    
                else:
                    pass #unphased. 
            else:
                pass #entire chromosome / contig did not exist in the phased hap file. Still, won't hurt adding the homozygous snps
            
            if not found_phased:
                #unphased. now check if it is homozygous
                sample = variant[iSample].split(':')
                hap = sample[0].split('/')
                if hap[0]==hap[1]:
                    # yes, a homozygous variant
                    variant.append(variant[iSample]) # add a second identical field for hap2
                    vcfOut.write('\t'.join(variant) +'\n')

    vcfOut.close();
    vcfFile.close();




def main():
    if len(sys.argv) > 3:
        phaseVCF(sys.argv[1], sys.argv[2], sys.argv[3])
        sys.exit(0)
    else:
        print 'Usage: phaseVCF inputHaps inputVCF outputVCF '
        sys.exit(1)

# Standard boilerplate to call the main() function to begin the program.
if __name__ == '__main__':
    main()