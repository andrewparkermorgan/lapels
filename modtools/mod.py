'''
The module of mod file. 

Created on Sep 20, 2012

@author: Shunping Huang
'''

import gc
from modtools import posmap


VERSION = '0.0.2'

class Mod:
    '''The class for parsing a piece of a mod file from the same chromosome.'''

    def __init__(self, chrom, maxLen, dataIter=None):
        '''Initialize with the data from one chromosome.
        chrom: the template chromosome id
        maxLen: the template chromosome length
        dataIter: an iterator of tuples extracted from a mod file. Each tuple
            should have 4 fields.
        '''

        self.seq = None
        self.posmap = None
        self.chrom = chrom
        self.maxLen = maxLen
        self.data = None
        if dataIter is not None:
            self.load(dataIter)


    def load(self, dataIter, isConverted=False):
        '''Load data from an iterator and do conversion of integer if needed.'''
        assert dataIter is not None
        gc.disable()
        self.data = []
        append = self.data.append
        if isConverted:
            for row in dataIter:
                append(row)
        else:
            for row in dataIter:
                row[2] = int(row[2]) # Convert positions to integers.
                append(row)
        gc.enable()


    def buildPosMap(self):
        '''Build the position mapping instance.'''
        assert self.data is not None
        gc.disable()
        # Read rows at the same position
        data = self.data
        chrom = self.chrom
        nRows = len(data)
        maps = []

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

        #assert refPos <= refLens[chrom]
        if refPos > self.maxLen:
            raise ValueError("Variant position %d out of reference boundary"
                             % refPos)

        if refPos < self.maxLen:
            maps.append((chrom, refPos, chrom, newPos, self.maxLen-refPos, '+'))
 
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


    def getPosMap(self):
        '''Return a PosMap instance of the current mod instance.'''
        if self.posmap is None:
            self.buildPosMap()
        assert self.posmap is not None
        return self.posmap


    def buildSeq(self, refSeqs):
        '''Build the sequence based on the mod data and reference sequences.'''
        assert self.data is not None
        gc.disable()
        # Read rows at the same position
        data = self.data
        chrom = self.chrom
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
                seqs.append(refSeqs[chrom][refPos:varPos])
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
                    seqs.append(refSeqs[chrom][refPos:refPos+segLen])
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
        if refPos > self.maxLen:
            raise ValueError("Variant position out of reference boundary")

        if refPos < self.maxLen:
            seqs.append(refSeqs[chrom][refPos:self.maxLen])

        self.seq = ''.join(seqs)
        gc.enable()


    def getSeq(self, refSeqs):
        '''Return the sequence of the in silico chromosome'''
        if self.seq is None:
            self.buildSeq(refSeqs)
        assert len(self.seq) > 0
        return self.seq

