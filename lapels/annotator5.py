'''
The module of annotate reads that aligned to in silico (pseudo) genome, and
map them back to reference coordinate.


Created on Sep 1, 2012

@author: Shunping Huang
'''

import bisect
import re
import logging

import cigarutils as cu
import cigarbuilder
import regionutils

from lapels.utils import log

VERSION = '0.0.5'
TESTING = False # For unit test: return the annotated read data once set to 1.
VERBOSITY = 0   # Print out detail log during the process.



def logRegion(reg):    
    log("  op: %d.\n" % reg[0])
    log("  region start: %d, end: %d, mpos: %d.\n" 
        % (reg[2], reg[3], reg[4]))
    log("  region cigar: '%s'.\n" % cu.toString(reg[1]))
    if len(reg) > 5:
        log("  annotation: s: %d, i: %d, d: %d.\n" % (reg[5], reg[6], reg[7]))  
                            
    
def getReadOffset(rseq, pos):
    '''
    Given the position in the alignment coordinate, return its offset in read
    '''
    readPos = rseq.pos    
    readLen = rseq.qlen
    cigar = rseq.cigar

    if (pos < readPos):        
        raise ValueError("position underflows")

    offset = 0
    curPos = readPos
    for op, length in cigar:        
        if op == 0 or op == 7 or op == 8:     # Match
            if (pos < curPos + length):       # (a) pos in this match         
                ret = offset + pos - curPos
                if ret < readLen:
                    return ret
                else:
                    msg = "cigar '%s' and length '%s' conflict in read '%s'"
                    cigarStr = cu.toString(cigar)
                    raise ValueError(msg % (cigarStr, readLen, rseq.qname))
            else:
                offset += length              # (b) pos not in this match
                curPos += length
        elif op == 1:                         # Insertion
            offset += length
        elif op == 2 or op == 3:              # Deletion or Splicing junction
            curPos += length
        else:
            raise NotImplementedError("unknown op '%s' in read '%s'" 
                                      % (op,rseq.qname))
        
        # In the case of offset == readLen, there may be an overflow        
        if offset > readLen:
            msg = "cigar '%s' and length '%s' conflict in read '%s'"
            cigarStr = cu.toString(cigar)
            raise ValueError(msg % (cigarStr, readLen, rseq.qname))
          
        if curPos > pos:
            raise ValueError('position in deletion or splicing junction')
        
    raise ValueError('position overflows')



def getTargetRegions(rseq):
    '''Given a read, return tuples of region data in the target coordinate.'''
    ret = []
    newpos = rseq.pos
    cigar = rseq.cigar
    for op, length in cigar:
        if op == 0 or op == 7 or op == 8: ## Match(M_1)
            assert length > 0
            ret.append((op, None, newpos, newpos+length-1))
            newpos += length
        elif op == 1:                     ## Insertion(I1)
            ret.append((op, None, newpos, newpos-1))
        elif op == 2 or op == 3:          ## Deletion(D_1)/Splice junction(N_1)
            ret.append((op, None, newpos, newpos+length-1))
            newpos += length
    if newpos != rseq.aend:
        raise ValueError("cigar %s conflicts with read region %d-%d." 
                         % (str(rseq.cigar), rseq.pos, rseq.aend))
    return ret


    
class Annotator:
    '''The class for annotating reads'''
    
    def __init__(self, chrom, chromLen, mod, inBam, nReads=None, 
                 tagPrefixes=None, outBam = None):
        self.logger = logging.getLogger('annotator')            
        self.mod = mod
        self.chrom = chrom          #chrom in mod
        self.chromLen = chromLen
        self.nReads = nReads
        self.inBam = inBam
        self.tagPrefixes = tagPrefixes
        self.outBam = outBam
        
    
    def setTag(self, tags, key, value):
        '''Set read tag'''        
        if key is not None:
            tags[key] = value
    
    
    def parseTargetRegion(self, region, rseq):
        '''
        Parse a read region in the target coordinate and get data in ref.
        
        input region: a tuple of (op, None, start, end)
        output region: a tuple of (op, newCigar, newStart, newEnd, newPos,
                                   nSNPs, nInsertions, nDeletions)
        '''
                        
        data = self.data
        posmap = self.posmap
        modKeys = self.modKeys  # All variant positions in reference coordinate 
        
        nSNPs = 0
        nInsertions = 0
        nDeletions = 0
                                       
        # The region's exact start and end in the target coordinate.
        op = region[0]
        tstart = region[2]
        tend = region[3]      
        if tstart > tend:
            assert op == 1  # It must be an insertion region(I_1)
            return (op,)
                    
        # Find the region's fuzzy boundaries in the reference coordinate.
        # The position will be imprecise if it falls in an Insertion(I_0).        
        rstart = posmap.bmap((self.chrom, tstart))
        rend = posmap.bmap((self.chrom, tend))
        
        # Detect the case of translocation between two chromosomes.
        if rstart[0] != self.chrom or rend[0] != self.chrom:
            raise NotImplementedError("cannot parse this region.")
        rstart = abs(rstart[1])
        rend = abs(rend[1])
        # Detect the case of translocation/duplication/inversion.
        if rstart > rend:
            raise NotImplementedError("cannot parse this region.")
                
        if VERBOSITY > 1:            
            log("T: %d-%d; R: %d-%d\n" %(tstart, tend, rstart, rend))
                        
        # Initialize the new attributes of the region
        nstart = rend + 1   # The ref position of the first M or D.          
        nend = -1           # The ref position of the last M or D.
        npos = rend + 1     # The ref position of the first M (and M only).                     
        ncigar = []         # The new cigar in the reference coordinate                
        
        # The next processing position in ref and target coordinate
        rpos = rstart               
        tpos = posmap.fmap((self.chrom, rstart))[1]  # tpos <= tstart
        if tpos < 0:        # it falls in a deletion(D_0)
            tpos = -tpos + 1
        
        # The searching boundaries in mod data.
        # Variants in the region will be in data[lo:hi]
        lo = bisect.bisect_left(modKeys, rstart)
        hi = bisect.bisect_right(modKeys, rend)
                            
        if lo < hi: # There are some variants in the region            
            # The boundaries of the processing block in mod data.                    
            # Rows in data[startIdx:endIdx] will have the same ref position.
            startIdx = lo
            endIdx = lo            
            # The next processing variant position in reference coordinate.
            vpos = data[lo][2]
            for i in range(lo, hi+1):  # The plus 1 trick.
                if tpos > tend:        # Already reach the region's end.
                    break
                
                # Find the rows that share the same ref position.
                if i < hi and data[i][2] == vpos:  # The plus 1 trick
                    endIdx += 1
                    continue
                
                #assert rpos <= vpos
                if (rpos > vpos):
                    raise ValueError("position not sorted in MOD at line %d." % 
                                     (i+1))
                            
                if rpos < vpos:     
                    # Fill M's to the head or the gap between sub-regions.
                    if tpos >= tstart and tpos <= tend:
                        ncigar.append((0, vpos-rpos))                                                                         
                        nstart = min(nstart, rpos)
                        nend = max(nend, vpos-1)
                        npos = min(npos, rpos)                                              
                    tpos += vpos - rpos
                    rpos = vpos
                
                if VERBOSITY > 1:
                    for j in range(startIdx, endIdx):
                        log(','.join(map(str,data[j])))
                        log('\n')
                
                # Handle the sub-regions for each variant position.
                subRegs=[('m',1)]
                for j in range(startIdx,endIdx):
                    tup = data[j]                    
                    if tup[0] == 's' and subRegs[0] != 'd': # 'd' overrides 's'
                        subRegs[0] = ('s', tup[3])
                    elif tup[0] == 'i':
                        subRegs.append(('i', tup[3]))
                    elif tup[0] == 'd':
                        subRegs[0] = ('d')
                    else:
                        raise NotImplementedError("unknown op '%s'" % tup[0])
    
                for reg in subRegs:
                    segType = reg[0]                                  
                    if segType == 'm':      # Match
                        if tpos >= tstart and tpos <= tend:
                            ncigar.append((0, 1))                                
                            nstart = min(nstart, rpos)
                            nend = max(nend, rpos)
                            npos = min(npos, rpos)
                        rpos += 1
                        tpos += 1
                    elif segType == 's':    # Substitution
                        if tpos >= tstart and tpos <= tend:                  
                            ncigar.append((0, 1))                                
                            nstart = min(nstart, rpos)
                            nend = max(nend, rpos)
                            npos = min(npos, rpos)
                            
                            if op != 2 and op != 3: # D_1 or N_1                                                                
                                rbase = rseq.seq[getReadOffset(rseq, tpos)]                            
                                if VERBOSITY > 1:                                                        
                                    log("SNP found at %d: %s\n" 
                                             % (tpos, reg[1]))                                
                                    log("Read base: %s\n" % rbase)
                                if rbase == reg[1][-1]:
                                    nSNPs += 1
                        rpos += 1
                        tpos += 1
                    elif segType == 'i':    # Insertion         
                        # TODO: Support multiple insertions in a line.
                        segLen = len(reg[1])                            
                        tmax = min(tpos + segLen, tend + 1)                        
                        if tpos > tstart:
                            ncigar.append((1, tmax - tpos))
                        else:
                            if tmax > tstart:
                                ncigar.append((1, tmax - tstart))                                                                                                                                    
                        tpos = tmax                                
                    elif segType == 'd':    # Deletion
                        if tpos > tstart and tpos <= tend:                                                  
                            ncigar.append((2, 1))
                            # NO assignment of 'npos' here: npos is for M only.
                            nstart = min(nstart, rpos)
                            nend = max(nend, rpos)
                        rpos += 1
                    else:
                        raise NotImplementedError("unknown op '%s'" % segType)                                                                
                    if tpos > tend:
                        break
                                
                startIdx = endIdx
                endIdx += 1                
                if i < hi:  # The plus 1 trick
                    vpos = data[i][2]
                    
        #assert rpos <= refLens[chrom]
        if rpos > rend + 1:
            raise ValueError("variant position out of boundary")

        # Fill M's to the region end if bases are not enough
        if rpos < rend + 1:
            ncigar.append((0, rend-rpos+1))
            npos = min(npos, rpos)
            nstart = min(nstart, rpos)
            nend = max(nend, rend)
        
        # Count the number of bases for insertions and deletions in a region.
        ncigar = cu.simplify(ncigar)
        for cig in ncigar:
            if cig[0] == 1:     # I_0
                nInsertions += cig[1]
            elif cig[0] == 2:   # D_0
                nDeletions += cig[1]
        
        if nstart > rend:            
            # If no previous assignment of nstart/nend, then this region 
            # contains merely I_0 regions.  
            # Assign them according to the anchor position (the last preceding 
            # position in reference), so that the gap of regions can be filled 
            # properly later.            
            
            # start is one-base after the anchor position, so that when it 
            # subtracts the previous end will yield 0.
            nstart = rpos
            
            # end is the position of the anchor point, so that when it is 
            # subtracted by the next start will yield 0. 
            nend = rpos - 1
            
        # Assign -1 if the region contains no M_0.
        # TODO: assign unmapped properties for this reads
        if npos > rend:
            npos = -1
        
        return (op, ncigar, nstart, nend, npos, nSNPs, nInsertions, nDeletions)
    
        
    def execute(self):
        '''The driver method for the module'''
        self.logger.info("%d read(s) found in BAM", self.nReads)
                
        self.data = self.mod.data
        self.modKeys = [tup[2] for tup in self.data]            
        self.posmap = self.mod.getPosMap(self.chrom, self.chromLen)        
            
        count = 0
        count2 = 0
        if TESTING:
            results = []
        if self.nReads == 0:
            return 0
        if not TESTING:            
            self.logger.info("progress: %3d%%" % (count*100/self.nReads))
                     
        for rseq in self.inBam:      # Annotate each reads
            if VERBOSITY > 1:       
                log("read name: %s\n" % rseq.qname)                
                log("t. alignment pos: %d\n" % rseq.pos)
                log("t. alignment cigar: '%s'\n\n" % cu.toString(rseq.cigar))                
            rseq.cigar = cu.simplify(rseq.cigar)   # Simplify the cigar first.
            regions = []
            tregs = getTargetRegions(rseq)
            for idx, treg in enumerate(tregs):   
                if treg[0] == 0 or treg[0] == 7 or treg[0] == 8: # Match
                    if VERBOSITY > 1:
                        log("process match\n")
                        log("t. region cigar(%d): %s\n" 
                            % (idx, cu.toString([rseq.cigar[idx]])))
                    rreg = self.parseTargetRegion(treg, rseq)  
                    if VERBOSITY > 1:                     
                        logRegion(rreg)                                                                  
                elif treg[0] == 2 or treg[0] == 3: # Deletion/Splice junction
                    rreg = (treg[0], )                    
                else:                    
                    rreg = (1, [rseq.cigar[idx]], 0, -1, -1) # Insertion
                regions.append(rreg)
            
            if VERBOSITY > 1:                                                                
                log("\nafter parsing match regions:\n")
                log('\n'.join(map(str,regions)))
                log("\n\n")
                
            nRegions = len(regions)
            assert nRegions == len(rseq.cigar)
                        
            for idx in range(nRegions):
                op = regions[idx][0]
                if op == 2 or op == 3:   # Handle deletions and splicing                    
                    if VERBOSITY > 1:                                    
                        log("process deletion/splicing junction\n")
                        log("t. region cigar(%d): %s\n" 
                            % (idx, cu.toString([rseq.cigar[idx]])))                    
                    if ((idx > 0) and idx < (nRegions - 1) and 
                        (regions[idx-1][0] == 0) and (regions[idx+1][0] == 0)):
                        delta = regions[idx+1][2] - regions[idx-1][3] - 1
                        assert delta >= 0
                        if delta == 0:                                                   
                            rreg = (op, [], 
                                    regions[idx-1][3]+1, 
                                    regions[idx+1][2]-1, 
                                    -1)
                        else:                                                                                                       
                            rreg = (op, [(op, delta)], 
                                    regions[idx-1][3]+1, 
                                    regions[idx+1][2]-1, 
                                    -1)
                        if VERBOSITY > 1:
                            logRegion(rreg)           
                                
                    else:                        
                        rreg = self.parseTargetRegion(tregs[idx], rseq)                        
                        if VERBOSITY > 1:                                        
                            logRegion(rreg)                             
                        rreg = regionutils.modifyRegion(rreg)                                            
                        if VERBOSITY > 1:                    
                            log("after modification\n")  
                            logRegion(rreg)                        
                    regions[idx] = rreg                             
            
            if VERBOSITY > 1:                                                                
                log("\nafter parsing deletions/splicing junction regions:\n")
                log('\n'.join(map(str,regions)))  
                log('\n\n')              
            
            cb = cigarbuilder.CigarBuilder()
#            cigar=cb.build(regions)
            for reg in regions:
                cb.append(reg)
            cigar = cb.cigar
            
            # Fix the MIDM pattern
            if re.match('.*\d*I,\d*D', cu.toString(cigar)) is not None:
                if VERBOSITY > 1:
                    log("fix MIDM pattern: %s\n" % cu.toString(cigar))                
#                print("%d,%d,%d" %(nSNPs, nInsertions, nDeletions))
#                print(regions)           
                for i in range(nRegions):
                    if tregs[i][0] == 1:                        
                        length = rseq.cigar[i][1]
                        assert length > 0
                                                                                              
                        offset = getReadOffset(rseq, tregs[i][3]) + 1
                        ins = (rseq.seq[offset:offset+length])
                        if VERBOSITY > 1:
                            log('seq in insertion: %s\n' % ins)
                        
                        if i > 0:
                            assert regions[i-1][3] >= 0                            
                            loKey = regions[i-1][3]
                        elif i < nRegions - 1 and regions[i+1][2] >= length:
                            loKey = regions[i+1][2] - length
                        else:
                            loKey = -1
                                                    
                        if i < nRegions - 1:
                            assert regions[i+1][2] >= 0
                            hiKey = regions[i+1][2] 
                        elif i > 0 and regions[i-1][3] >= 0:
                            hiKey = regions[i-1][3] + length
                        else:
                            hiKey = -1
                              
#                        loKey = regions[i-1][3]                      
#                        hiKey = regions[i+1][2]
                      
                        lo = bisect.bisect_left(self.modKeys, loKey)
                        hi = bisect.bisect_right(self.modKeys, hiKey)
                        if VERBOSITY > 1:
                            log('variants from %d-%d\n' % (loKey,hiKey))                          
                            for j in range(lo,hi):
                                log('%s\n' % str(self.data[j]))
                                                                                            
                        isMatched = False
                        matchStart = -1
                        pivot = 0                                                                                                                   
                        for j in range(lo,hi):
                            if self.data[j][0] != 'd':
                                continue                                           
                            if matchStart == -1:
                                matchStart = int(self.data[j][2])                            
                            if ins[pivot] == self.data[j][3]:
                                pivot += 1
                            else:
                                break
                            if j == hi - 1 or pivot >= length:
                                isMatched = True
                                break
                        
                        if isMatched:
                            if VERBOSITY > 1:
                                log('insertion matches gap from left\n')
                                log("before:\n")
                                logRegion(regions[i])
                            if pivot < length:                                                                                                              
                                regions[i] = (0, [(0, pivot),(1,length-pivot)], 
                                              matchStart, 
                                              matchStart + pivot - 1,
                                              matchStart)
                            else:
                                regions[i] = (0, [(0, pivot)], matchStart, 
                                              matchStart + pivot - 1,
                                              matchStart)
                            if VERBOSITY > 1:
                                log("after:\n")
                                logRegion(regions[i])
                                log("\n")
                        else:                                                
                            pivot = length - 1
                            for j in range(hi-1,lo-1,-1):
                                if self.data[j][0] != 'd':
                                    continue                                                                                                                         
                                if ins[pivot] == self.data[j][3]:
                                    matchStart = int(self.data[j][2])
                                    pivot -= 1
                                else:
                                    break
                                if j == lo or pivot <= 0:
                                    isMatched = True
                                    break
                        
                            if isMatched:
                                if VERBOSITY > 1:
                                    log('insertion matches gap from right\n')
                                    log("before:\n")
                                    logRegion(regions[i])
                                if pivot >= 0:                                                                                            
                                    regions[i] = (0, [(1, pivot+1), (0, length-1-pivot)], 
                                                  matchStart,
                                                  matchStart + length - pivot - 2,
                                                  matchStart)
                                else:
                                    regions[i] = (0, [(0, length-1-pivot)], 
                                                  matchStart,
                                                  matchStart + length - pivot - 2,
                                                  matchStart)
                                if VERBOSITY > 1:
                                    log("after:\n")
                                    logRegion(regions[i])
                                    log("\n")
                            else:
                                if VERBOSITY > 1:
                                    log('insertion not matches\n')                  
                            
                cb = cigarbuilder.CigarBuilder()
#                cigar=cb.build(regions)
                for reg in regions:
                    cb.append(reg)
                cigar = cb.cigar
#                print(cu.toString(cigar))
            
            
            if VERBOSITY > 1:                                                                
                log("\nafter fixing special pattern:\n")
                log('\n'.join(map(str,regions)))  
                log('\n\n')              
                
            nSNPs = 0
            nInsertions = 0
            nDeletions = 0
            for reg in regions:
                if len(reg) > 5:
                    nSNPs += reg[5]
                    nInsertions += reg[6]
                    nDeletions += reg[7]
                                                        
            ## Set tags
            tags = dict(rseq.tags)
            if self.tagPrefixes is not None:
                self.setTag(tags, self.tagPrefixes[0]+'0', nSNPs)            
                self.setTag(tags, self.tagPrefixes[1]+'0', nInsertions)
                self.setTag(tags, self.tagPrefixes[2]+'0', nDeletions)            
                self.setTag(tags, 'OC', cu.toString(rseq.cigar).translate(None,','))
                self.setTag(tags, 'OM', tags['NM'])            
                del tags['NM']  ## Delete the old 'NM' tag.
                  
            if nSNPs != 0 or nInsertions != 0 or nDeletions != 0:
                count2 += 1
                                          
            rseq.tags = [(key, tags[key]) for key in sorted(tags.keys())] 
                        
            ## Set pos to be the first M            
            pos = -1
            for reg in regions:
                if reg[4] >= 0:
                    pos = reg[4]
                    break
            rseq.pos = pos
            
            ## Set cigar
            rseq.cigar = cu.simplify(cigar)

            if self.outBam is not None:                
                self.outBam.write(rseq)
        
            if TESTING:
                results.append((rseq.cigar, rseq.pos, nSNPs, nInsertions, 
                                nDeletions))
        
            if VERBOSITY > 1:                   
                log("output read pos: %d\n" % rseq.pos)
                log("output read cigar: '%s'\n" 
                    % cu.toString(rseq.cigar))
                log("output annotation: s: %d, i: %d, d: %d.\n========================\n\n" % 
                    (nSNPs, nInsertions, nDeletions))
            count += 1
            
            if not TESTING and count % 100000 == 0:
                self.logger.info("progress: %3d%%" % (count*100/self.nReads))        
                
        if not TESTING:
            self.logger.info("progress: %3d%%" % (count*100/self.nReads))
                                
        if TESTING:              
            return results
        
        self.logger.info("%d read(s) written to file", count)
        self.logger.info("%d read(s) have variants", count2)
        return count
