'''
The module of mod file. 

Created on Sep 20, 2012

@author: Shunping Huang
'''

import gc
import pysam
import gzip
import logging
from modtools import posmap
from modtools import metadata


VERSION = '0.1.0'

__all__ = ['Mod', 'VERSION']

class Mod:
    '''The class for parsing a piece of a mod file from the same chromosome.'''
    
    def __init__(self, fileName):
        self.logger = logging.getLogger('mod')
        fp = gzip.open(fileName, 'rb')
        self.header = dict()
        for line in fp:
            line = line.strip()
            if len(line) == 0:
                continue
            elif line.startswith('#'):
                tup = line[1:].split('=')
                assert len(tup) == 2
                self.header[tup[0]] = tup[1]
            else:
                break
        fp.close()
#        print(header)
#        # Compress with bgzip
#        if not fileName.endswith('.gz'): 
#            if not os.path.isfile(fileName+'.gz'):
#                pysam.tabix_compress(fileName, fileName+'.gz')
#            fileName += '.gz'
#                
#        # Build tabix index
#        if not os.path.isfile(fileName+'.tbi'):
#            pysam.tabix_index(fileName, force=True, seq_col=1, start_col=2, 
#                              end_col=2, meta_char='#', zerobased=True)
                        
        self.tabix = pysam.Tabixfile(fileName)
        self.fileName = fileName
        self.chroms = self.tabix.contigs
        self.chrom = -1
        try:
            self.meta = metadata.MetaData(self.header['reference'])
        except KeyError:
            pass        
    
        
    def load(self, chrom):
        '''Load data from an iterator and do conversion of integer if needed.'''                    
        if self.chrom != chrom: # chrom not loaded                                                                                            
            self.chrom = chrom
            # Reset posmap, seq, and data
            self.posmap = None
            self.seq = None                            
            self.data = []        
            
            if chrom not in self.chroms:                    
                self.logger.warning("chromosome '%s' not found in MOD", chrom)
                return
        
            gc.disable()
            append = self.data.append            
            for line in self.tabix.fetch(reference=chrom):
                try:
                    cols = line.split('\t')
                    cols[2] = int(cols[2]) # Convert positions to integers.
                    cols[-1] = cols[-1].rstrip()
                    append(cols)                        
                except:
                    print(line)
                    print(cols)
                    raise Exception("ERROR!! at line %d" % len(self.data))
            gc.enable()
            
        assert len(self.data) > 0 
        self.logger.info("%d line(s) found in MOD" % len(self.data))


    def buildPosMap(self, chromLen):
        '''Build the position mapping instance.'''
        assert self.data is not None        
        gc.disable()
        # Read rows at the same position
        data = self.data
        chrom = self.chrom        
        nRows = len(data)
        maps = []

        self.logger.info("[%s]: building position map ...", chrom)
        # Current position in reference/new genome coordinate
        refPos = 0
        newPos = 0
        
        if len(data) > 0:        
            varPos = data[0][2]
    
            # Rows in data[startIdx:endIdx] have the same position
            startIdx = 0
            endIdx = 0
            for i in range(nRows+1):
                if i < nRows and data[i][2] == varPos:
                    endIdx+=1
                    continue
    
                if (refPos > varPos):
                    raise ValueError("Position not in order at line %d" %(i+1))
                
                # Fill 'M's in the gap.
                if refPos < varPos:       
                    maps.append((chrom, refPos, chrom, newPos, varPos-refPos, '+'))
                    newPos += varPos - refPos
                    refPos = varPos
    
                subSegs=[(1, 'm')]
                for j in range(startIdx,endIdx):
                    tup = data[j]
                    if tup[0] == 's':
                        subSegs[0] = (1, 's', tup[3])
                    elif tup[0] == 'i':
                        subSegs.append((len(tup[3]), 'i', tup[3]))
                    elif tup[0] == 'd':
                        subSegs[0] = (1, 'd')
                    else:
                        raise ValueError("Unknown operation %s" % tup[0])
    
                for seg in subSegs:
                    segLen = seg[0]
                    segType = seg[1]
                    if segType == 'm':
                        maps.append((chrom, refPos, chrom, newPos, segLen, '+'))
                        refPos += segLen
                        newPos += segLen
                    elif segType == 's':
                        maps.append((chrom, refPos, chrom, newPos, segLen, '+'))
                        refPos += segLen
                        newPos += segLen
                    elif segType == 'i':
                        # Insertion
                        # Set the ref position to the preceding ref position
                        maps.append((chrom, -refPos+1, chrom, newPos, segLen, '+'))
                        newPos += segLen
                    elif segType == 'd':
                        # Deletion
                        # Set the new position to the preceding new position
                        maps.append((chrom, refPos, chrom, -newPos+1, segLen, '+'))
                        refPos += segLen
                    else:
                        raise ValueError("Unknown operation %s" % segType)
    
                startIdx = endIdx
                endIdx += 1
    
                if i<len(data):
                    varPos = data[i][2]

#        assert refPos <= refLens[chrom]
        if refPos > chromLen:
            raise ValueError("Variant position %d out of reference boundary"
                             % refPos)

        if refPos < chromLen:
            maps.append((chrom, refPos, chrom, newPos, chromLen-refPos, '+'))
 
        assert len(maps) > 0

        # Compress consecutive matches or deletions.
        compressed = []
        buf = maps[0]
        for r in maps[1:]:            
            if (buf[0] == r[0] and buf[2] == r[2] and buf[5] == r[5] and 
                ((buf[1] >= 0 and r[1] >= 0 and buf[3] >= 0 and r[3] >= 0 and 
                  (buf[1]-buf[3]) == (r[1]-r[3])) or 
                 (buf[3] < 0 and buf[3] == r[3]))):
                buf=(buf[0],buf[1],buf[2],buf[3],buf[4]+r[4],buf[5])
            else:
                compressed.append(buf)
                buf = r
        compressed.append(buf)

        self.posmap = posmap.PosMap()
        self.posmap.load(compressed, isConverted=True)

        gc.enable()      


    def getPosMap(self, chrom, chromLen = None):
        '''Return a PosMap instance of the current mod instance.'''        
        if self.chrom != chrom:
            self.load(chrom)
        
        if chromLen is None:
            chromLen = self.meta.getChromLength(chrom)
            
        if self.posmap is None:
            self.buildPosMap(chromLen)
        assert self.posmap is not None
        return self.posmap


    def buildSeq(self, fasta, chrom, fastaChroms):
        '''Build the sequence based on the mod data and reference sequences.'''
        assert chrom == self.chrom
        
        data = self.data        
        assert data is not None
        
        self.logger.info("[%s]: building sequence ...", chrom)
        
        meta = self.meta 
        basicName = meta.chromAliases.getBasicName(chrom)
        fastaChrom = meta.chromAliases.getMatchedAlias(basicName, fastaChroms)        
        if fastaChrom is None:
            raise ValueError("Chromosome '%s' not found in FASTA. " % chrom +
                             "Possible names: %s. " % 
                             ','.join(sorted(fastaChroms)))
        
        # If no content in MOD for this chromosome
        if len(data) == 0:
            self.seq = fasta.fetch(reference=fastaChrom, start=0) 
            return 
        
        gc.disable()        
        
        chromLen = meta.getChromLength(chrom)       
        nRows = len(data)
        seqs = []

        # Current position in reference/new genome coordinate
        refPos = 0
        newPos = 0
        varPos = data[0][2]

        # Rows in data[startIdx:endIdx] have the same position
        startIdx = 0
        endIdx = 0
        for i in range(nRows+1):
            if i < nRows and data[i][2] == varPos:
                endIdx+=1
                continue

            #assert refPos <= varPos
            if (refPos > varPos):
                raise ValueError("Position not in order at line %d" %(i+1))

            # Fill 'M's in the gap.
            if refPos < varPos:                
                seqs.append(fasta.fetch(reference=fastaChrom, start=refPos, 
                                        end=varPos))                
                newPos += varPos - refPos
                refPos = varPos                

            subSegs=[(1, 'm')]
            for j in range(startIdx,endIdx):
                tup = data[j]
                if tup[0] == 's':
                    subSegs[0] = (1, 's', tup[3])
                elif tup[0] == 'i':
                    subSegs.append((len(tup[3]), 'i', tup[3]))
                elif tup[0] == 'd':
                    subSegs[0] = (1, 'd')
                else:
                    raise ValueError("Unknown operation %s" % tup[0])

            for seg in subSegs:
                segLen = seg[0]
                segType = seg[1]
                if segType == 'm':
                    seqs.append(fasta.fetch(reference=fastaChrom, start=refPos,
                                            end=refPos+segLen))
                    refPos += segLen
                    newPos += segLen
                elif segType == 's':
                    seqs.append(seg[2][-1])
                    refPos += segLen
                    newPos += segLen
                elif segType == 'i':
                    seqs.append(seg[2])
                    newPos += segLen
                elif segType == 'd':
                    refPos += segLen
                else:
                    raise ValueError("Unknown operation %s" % segType)

            startIdx = endIdx
            endIdx += 1

            if i<len(data):
                varPos = data[i][2]

        #assert refPos <= refLens[chrom]
        if refPos > chromLen:
            raise ValueError("Variant position out of reference boundary")

        if refPos < chromLen:
            seqs.append(fasta.fetch(reference=fastaChrom, start=refPos))

        self.seq = ''.join(seqs)
        gc.enable()


    def getSeq(self, fasta, chrom, fastaChroms):
        '''Return the sequence of the in silico chromosome'''
        if self.chrom != chrom:
            self.load(chrom)
            
        if self.seq is None:
            self.buildSeq(fasta, chrom, fastaChroms)
        assert len(self.seq) > 0
        return self.seq
    
        
    
#    def buildRawCigar(self, fasta, chrom, chromMap, chromLens):
#        '''Build the sequence based on the mod data and reference sequences.'''
#        
#        M = 'M'
#        I = 'I'
#        S = 'S'
#        D = 'D'
#        
#        assert chrom == self.chrom
#        
#        data = self.data        
#        assert data is not None
#                
#        fastaChrom = getOutChrom(chromMap, chrom)
#        if fastaChrom not in chromLens.keys():
#            raise ValueError("Chromosome '%s' not found in FASTA. " % 
#                             fastaChrom +
#                             "Possible names: %s. " % 
#                             ','.join(sorted(chromLens.keys())) +
#                             "Chromosome name mapping may be used.\n")
#                
#        # If no content in MOD for this chromosome
#        if len(data) == 0:
#            self.seq = fasta.fetch(reference=fastaChrom, start=0) 
#            return 
#        
#        gc.disable()
#        
#        chromLen = chromLens[fastaChrom]        
#        nRows = len(data)
#        seqs = []
#
#        # Current position in reference/new genome coordinate
#        refPos = 0
#        newPos = 0
#        varPos = data[0][2]
#
#        # Rows in data[startIdx:endIdx] have the same position
#        startIdx = 0
#        endIdx = 0
#        for i in range(nRows+1):
#            if i < nRows and data[i][2] == varPos:
#                endIdx+=1
#                continue
#
#            #assert refPos <= varPos
#            if (refPos > varPos):
#                raise ValueError("Position not in order at line %d" %(i+1))
#
#            # Fill 'M's in the gap.
#            if refPos < varPos:                
#                #seqs.append(fasta.fetch(reference=fastaChrom, start=refPos, 
#                #                        end=varPos))
#                for k in range(refPos,varPos):
#                    seqs.append(M)
#                                    
#                newPos += varPos - refPos
#                refPos = varPos                
#
#            subSegs=[(1, 'm')]
#            for j in range(startIdx,endIdx):
#                tup = data[j]
#                if tup[0] == 's':
#                    subSegs[0] = (1, 's', tup[3])
#                elif tup[0] == 'i':
#                    subSegs.append((len(tup[3]), 'i', tup[3]))                        
#                elif tup[0] == 'd':
#                    subSegs[0] = (1, 'd')
#                else:
#                    raise ValueError("Unknown operation %s" % tup[0])
#
#            for seg in subSegs:
#                segLen = seg[0]
#                segType = seg[1]
#                if segType == 'm':                    
#                    #seqs.append(fasta.fetch(reference=fastaChrom, start=refPos,
#                    #                        end=refPos+segLen))
#                    for k in range(segLen):
#                        seqs.append(M)
#                    refPos += segLen
#                    newPos += segLen
#                elif segType == 's':
#                    #seqs.append(seg[2][-1])                    
#                    seqs.append(S)
#                    refPos += segLen
#                    newPos += segLen
#                elif segType == 'i':
#                    #seqs.append(seg[2])
#                    for k in range(len(seg[2])):
#                        seqs.append(I)
#                    newPos += segLen
#                elif segType == 'd':
#                    refPos += segLen
#                else:
#                    raise ValueError("Unknown operation %s" % segType)
#
#            startIdx = endIdx
#            endIdx += 1
#
#            if i<len(data):
#                varPos = data[i][2]
#
#        #assert refPos <= refLens[chrom]
#        if refPos > chromLen:
#            raise ValueError("Variant position out of reference boundary")
#
#        if refPos < chromLen:
#            #seqs.append(fasta.fetch(reference=fastaChrom, start=refPos))
#            for k in range(refPos, chromLen):
#                seqs.append(M)
#
#        gc.enable()
#        return ''.join(seqs)

