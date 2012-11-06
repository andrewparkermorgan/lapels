'''
A module for storing meta data of a (reference) genome

Created on Nov 5, 2012

@author: Shunping Huang
'''

import os
import xml.etree.ElementTree as ET
from modtools import alias

__all__ = ['defaultXML', 'MetaData']

mm9_xml = '''
<genome>
<name>mm9</name>
<alias>NCBI37</alias>
<alias>NCBI Build 37</alias>
<alias>MGSCv37</alias>
<source>
  <file><url>ftp://ftp-mouse.sanger.ac.uk/ref/NCBIM37_um.fa</url></file>
</source>
<source>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr1.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr2.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr3.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr4.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr5.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr6.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr7.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr8.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr9.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr10.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr11.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr12.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr13.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr14.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr15.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr16.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr17.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr18.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chr19.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chrX.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chrY.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/chrM.fa.gz</url></file>
</source>
<source>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr1.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr2.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr3.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr4.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr5.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr6.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr7.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr8.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr9.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr10.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr11.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr12.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr13.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr14.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr15.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr16.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr17.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr18.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chr19.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chrX.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chrY.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.37.1/Assembled_chromosomes/mm_ref_chrM.fa.gz</url></file>
</source>
<chromosome><name>1</name> <alias>chr1</alias> <length>197195432</length></chromosome>
<chromosome><name>2</name> <alias>chr2</alias> <length>181748087</length></chromosome>
<chromosome><name>3</name> <alias>chr3</alias> <length>159599783</length></chromosome>
<chromosome><name>4</name> <alias>chr4</alias> <length>155630120</length></chromosome>
<chromosome><name>5</name> <alias>chr5</alias> <length>152537259</length></chromosome>
<chromosome><name>6</name> <alias>chr6</alias> <length>149517037</length></chromosome>
<chromosome><name>7</name> <alias>chr7</alias> <length>152524553</length></chromosome>
<chromosome><name>8</name> <alias>chr8</alias> <length>131738871</length></chromosome>
<chromosome><name>9</name> <alias>chr9</alias> <length>124076172</length></chromosome>
<chromosome><name>10</name><alias>chr10</alias><length>129993255</length></chromosome>
<chromosome><name>11</name><alias>chr11</alias><length>121843856</length></chromosome>
<chromosome><name>12</name><alias>chr12</alias><length>121257530</length></chromosome>
<chromosome><name>13</name><alias>chr13</alias><length>120284312</length></chromosome>
<chromosome><name>14</name><alias>chr14</alias><length>125194864</length></chromosome>
<chromosome><name>15</name><alias>chr15</alias><length>103494974</length></chromosome>
<chromosome><name>16</name><alias>chr16</alias> <length>98319150</length></chromosome>
<chromosome><name>17</name><alias>chr17</alias> <length>95272651</length></chromosome>
<chromosome><name>18</name><alias>chr18</alias> <length>90772031</length></chromosome>
<chromosome><name>19</name><alias>chr19</alias> <length>61342430</length></chromosome>
<chromosome><name>X</name><alias>chrX</alias>  <length>166650296</length></chromosome>
<chromosome><name>Y</name><alias>chrY</alias>   <length>15902555</length></chromosome>
<chromosome><name>M</name><alias>chrM</alias><alias>MT</alias><alias>chrMT</alias><length>16299</length></chromosome>
</genome>
'''


mm10_xml = '''
<genome>
<name>mm10</name>
<alias>GRCm38</alias>
<alias>Genome Reference Consortium Mouse Build 38</alias>
<source>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr1.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr2.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr3.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr4.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr5.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr6.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr7.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr8.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr9.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr10.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr11.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr12.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr13.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr14.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr15.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr16.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr17.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr18.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chr19.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chrX.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chrY.fa.gz</url></file>
  <file><url>http://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes/chrM.fa.gz</url></file>
</source>
<source>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr1.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr2.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr3.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr4.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr5.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr6.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr7.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr8.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr9.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr10.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr11.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr12.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr13.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr14.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr15.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr16.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr17.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr18.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chr19.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chrX.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chrY.fa.gz</url></file>
  <file><url>ftp://ftp.ncbi.nih.gov/genomes/M_musculus/ARCHIVE/BUILD.38.1/Assembled_chromosomes/seq/mm_ref_GRCm38_chrM.fa.gz</url></file>      
</source>
<chromosome><name>1</name> <alias>chr1</alias> <length>195471971</length></chromosome>
<chromosome><name>2</name> <alias>chr2</alias> <length>182113224</length></chromosome>
<chromosome><name>3</name> <alias>chr3</alias> <length>160039680</length></chromosome>
<chromosome><name>4</name> <alias>chr4</alias> <length>156508116</length></chromosome>
<chromosome><name>5</name> <alias>chr5</alias> <length>151834684</length></chromosome>
<chromosome><name>6</name> <alias>chr6</alias> <length>149736546</length></chromosome>
<chromosome><name>7</name> <alias>chr7</alias> <length>145441459</length></chromosome>
<chromosome><name>8</name> <alias>chr8</alias> <length>129401213</length></chromosome>
<chromosome><name>9</name> <alias>chr9</alias> <length>124595110</length></chromosome>
<chromosome><name>10</name><alias>chr10</alias><length>130694993</length></chromosome>
<chromosome><name>11</name><alias>chr11</alias><length>122082543</length></chromosome>
<chromosome><name>12</name><alias>chr12</alias><length>120129022</length></chromosome>
<chromosome><name>13</name><alias>chr13</alias><length>120421639</length></chromosome>
<chromosome><name>14</name><alias>chr14</alias><length>124902244</length></chromosome>
<chromosome><name>15</name><alias>chr15</alias><length>104043685</length></chromosome>
<chromosome><name>16</name><alias>chr16</alias> <length>98207768</length></chromosome>
<chromosome><name>17</name><alias>chr17</alias> <length>94987271</length></chromosome>
<chromosome><name>18</name><alias>chr18</alias> <length>90702639</length></chromosome>
<chromosome><name>19</name><alias>chr19</alias> <length>61431566</length></chromosome>
<chromosome><name>X</name><alias>chrX</alias>  <length>171031299</length></chromosome>
<chromosome><name>Y</name><alias>chrY</alias>   <length>91744698</length></chromosome>
<chromosome><name>M</name><alias>chrM</alias><alias>MT</alias><alias>chrMT</alias><length>16299</length></chromosome>
</genome>
'''

defaultXML = {'mm9': mm9_xml, 'mm10': mm10_xml}

#open('mm9.xml','wb').write(mm9_xml)
#mm9tree = ET.ElementTree(file='mm9.xml')
#mm9genome = mm9tree.getroot()
#chroms = mm9genome.findall('chromosome') 
#mm9len = dict()
#for chrom in chroms:
#    mm9len[chrom.find('name').text] = int(chrom.find('length').text) 
##mm9tree.write('mm9.out.xml')
#
#open('mm10.xml','wb').write(mm10_xml)
#mm10tree = ET.ElementTree(file='mm10.xml')
#mm10genome = mm9tree.getroot()
#chroms = mm10genome.findall('chromosome')
##mm10tree.write('mm10.out.xml')
# 
#mm10len = dict()
#for chrom in chroms:
#    mm10len[chrom.find('name').text] = int(chrom.find('length').text) 
#
#for chrom in sorted(mm9len.keys()):
#    print chrom, mm9len[chrom], mm10len[chrom], mm10len[chrom]-mm9len[chrom], (mm10len[chrom]-mm9len[chrom])*100.0/mm9len[chrom] 


class MetaData:
        
    def __init__(self, name=None, fileName = None):
        if name is not None:
            if name in defaultXML.keys():
                self.load(name)
            elif os.path.isfile(name+'.xml'):
                self.loadFromFile(name+'.xml')
            else:
                raise ValueError("Cannot find meta data for '%s'" % name)
        elif fileName is not None:
            self.loadFromFile(fileName)
    
    
    def load(self, name):
        '''Load a default genome'''
        open(name+'.xml', 'wb').write(defaultXML[name])                
        self.loadFromFile(name+'.xml')        
    
    
    def loadFromFile(self, fileName):
        '''Load meta data from an external XML file'''
        tree = ET.ElementTree(file=fileName)
        root = tree.getroot()
        self.chromNames = [chrom.find('name').text for chrom in root.findall('chromosome')]
        self.chromLengths = dict([(chrom.find('name').text, 
                                   int(chrom.find('length').text)) 
                                  for chrom in root.findall('chromosome')])
#        print(self.chromNames)
#        print(self.chromLengths)
        chromClasses = []
        for chrom in root.findall('chromosome'):
            chromClasses.append([])
            chromClasses[-1].append(chrom.find('name').text)
            for tag in chrom.findall('alias'): 
                chromClasses[-1].append(tag.text)            
#        print(chromClasses)
        self.chromAliases = alias.Alias(chromClasses)
                
    
    def getChromNames(self):
        return self.chromNames
    
    
    def getChromAliases(self):
        return self.chromAliases
    
    
    def getChromLengths(self):
        return self.chromLengths

    
    def getChromLength(self, chrom):
        basicName = self.chromAliases.getBasicName(chrom)        
        return self.chromLengths.get(basicName, None)
    
    
    def verify(self, chrom, fastaFileName):
        # TODO: examine a fasta file to see if it matches the metadata.
        # So far only the length of chromosome should be checked.
        #
        # Please note that the fastaFile may contain one or more chromosomes,
        # and you need to check all chromosomes in it.
        #
        # Use pysam to open/index/read fasta file, such as:
        #    pysam.Fastafile
        #    pysam.faidx
        #    fasta.fetch
        pass

    

if __name__ == '__main__':
    genome = MetaData('mm9')    
    print(genome.getChromLength('chr1'))
    