'''
The module of annotate reads that aligned to in silico (pseudo) genome, and
map them back to reference coordinate.


Created on Sep 1, 2012

@author: Shunping Huang
'''

import bisect

import cigarutils as cu
import cigarbuilder
import regionutils

from lapels.utils import log

VERSION = '0.0.5'
TESTING = False # For unit test: return the annotated read data once set to 1.
VERBOSITY = 0   # Print out detail log during the process.


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
            self.log("T: %d-%d; R: %d-%d\n" %(tstart, tend, rstart, rend))
                        
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
                        self.log(','.join(map(str,data[j])))
                        self.log('\n')
                
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
                                    self.log("SNP found at %d: %s\n" 
                                             % (tpos, reg[1]))                                
                                    self.log("Read base: %s\n" % rbase)
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
        
        self.data = self.mod.data
        self.modKeys = [tup[2] for tup in self.data]            
        self.posmap = self.mod.getPosMap(self.chrom, self.chromLen)        
            
        count = 0
        count2 = 0
        if TESTING:
            results = []
        if self.nReads == 0:
            return 0
        if VERBOSITY > 0 and not TESTING:            
            log("progress: %3d%%\r" % (count*100/self.nReads), 1, True)
                     
        for rseq in self.inBam:      # Annotate each reads
            if VERBOSITY > 1:                
                log("%s\n" % rseq.qname)                
                log("tgt read pos: %d.\n" % rseq.pos)
                log("tgt read cigar: '%s'.\n\n" % cu.toString(rseq.cigar))
                
            rseq.cigar = cu.simplify(rseq.cigar)   # Simplify the cigar first.
            regions = []
            rend = -1
            for idx, treg in enumerate(getTargetRegions(rseq)):
                if VERBOSITY > 1:
                    log("tgt reg cigar(%d): %s.\n" 
                        % (idx, cu.toString([rseq.cigar[idx]])))                 
                rreg = self.parseTargetRegion(treg, rseq)
                
                if rreg[0] != 1:                # Not insertions (I_1)
                    if VERBOSITY > 1:                                        
                        log("rop: %d.\n" % rreg[0])
                        log("ref reg start: %d, end: %d, pos: %d.\n" 
                            % (rreg[2], rreg[3], rreg[4]))
                        log("ref reg cigar: '%s'.\n" 
                            % cu.toString(rreg[1]))
                        log("s: %d, i: %d, d: %d.\n" 
                            % (rreg[5], rreg[6], rreg[7]))                    
                                                  
                    rreg = regionutils.modifyRegion(rreg)
                    
                    if VERBOSITY > 1:                      
                        log("fixed start: %d, end: %d, pos: %d.\n" 
                            % (rreg[2], rreg[3], rreg[4]))
                        log("fixed cigar: '%s'.\n" 
                            % cu.toString(rreg[1]))
                        log("s: %d, i: %d, d: %d.\n\n" 
                            % (rreg[5], rreg[6], rreg[7]))                                            
                    if rend < rreg[3]:
                        rend = rreg[3]
                else:
                    rreg = (1, [rseq.cigar[idx]], rend+1, rend, -1)
                regions.append(rreg)
            
            #print regions
            assert len(regions) == len(rseq.cigar)
            
#            maxIdx = len(regions)
#            ## I1/M1/(N1,D1)
#            leftIdx = -1
#            leftPos = -1             
#            rightPos = -1
#            rightIdx = -1
#            lastInsertionEnd = -2
#            for idx, reg in enumerate(regions):
#                ## Get an M, the left bound will be one base after
#                if reg[0] != 1:
#                    if reg[0] == 0:                   #M_1
#                        leftIdx = idx
#                        leftPos = regions[leftIdx][3] + 1
#                    ## Get a D or N, left bound will include it 
#                    elif reg[0] == 2 and leftPos < 0: #D_1
#                        leftIdx = idx
#                        leftPos = regions[leftIdx][2]
#                    elif reg[0] == 3 and leftPos < 0: #N_1
#                        leftIdx = idx
#                        leftPos = regions[leftIdx][2]                
#                else:                    
##                    assert reg[0] == 1
#                    ## Get an I
#                    ## Find the right position bound of this insertion
#                    if rightIdx < idx:
#                        rightPos = -1
#                        rightIdx = idx + 1                                                
#                        while rightIdx < maxIdx:
#                            ## Got an M, the right bound will be one base before
#                            if regions[rightIdx][0] == 0:        #M_1
#                                #print(regions[rightIdx])    
#                                rightPos = regions[rightIdx][2] - 1                              
#                                break
#                            ## Got a D or N, the right bound will include it
#                            elif regions[rightIdx][0] == 2:      #D_1
#                                rightPos = regions[rightIdx][3]
#                            elif regions[rightIdx][0] == 3:      #N_1
#                                rightPos = regions[rightIdx][3]
#                            rightIdx += 1 
#                    assert rightIdx > idx
#                    if VERBOSITY > 1: 
#                        log("-check tgt reg cigar(%d) from %d @ %d to %d @ %d\n" 
#                            %(idx, leftPos, leftIdx, rightPos, rightIdx)) 
#                    
#                    ## Shrink the search region if a previous I is in the region 
#                    leftPos = max(leftPos, lastInsertionEnd+1)
#                                        
#                    ## Either left bound or right bound is missing.
#                    ## FIX ME!!! may try to search the upstream or downstream
#                    ## for D_0 if only one bound is missing.
#                    cigar = [rseq.cigar[idx]]
#                    if leftPos < 0:
#                        if VERBOSITY > 1:                                      
#                            log("?? left bound not found.\n")
#                        ## Put it before any read position
#                        regions[idx]=(reg[0], cigar, -1, -2, -1)
#                    elif rightPos < 0:
#                        if VERBOSITY > 1:
#                            log("?? right bound not found.\n")
#                        ## Put it after any read position
#                        regions[idx]=(reg[0], cigar, rend+1, rend, -1)
#                    else:                                                                    
#                        if leftPos <= rightPos:
#                            if VERBOSITY > 1:
#                                log("!! some gap between left and right.\n")
#                                
#                            ## FIXME: find the correct position in D_0                                 
##                            lo = bisect.bisect_left(self.modKeys, leftPos)
##                            hi = bisect.bisect_right(self.modKeys, rightPos)                            
##                            if VERBOSITY > 1:
##                                for j in range(lo,hi):
##                                    log(str(self.data[j]))
##                                    log("\n")
#                                                                                                               
#                            ## do the match here
#                            ## data[lo:hi] to see if D_0 is there.
#                            ## inserted bases in the read (how to get this?)
#                            ## return position that match
#                                                        
##                            regions[idx]=(reg[0], leftPos, 
##                                          leftPos + cigar[0][1] - 1, 
##                                          leftPos, [(0,cigar[0][1])])
#                            regions[idx]=(reg[0], cigar, leftPos, leftPos-1, -1)
#                        else:
#                            ## Insert into an M_0
#                            assert leftPos - rightPos == 1
#                            if VERBOSITY > 1:
#                                log("## no gap between left and right.\n")
#                            regions[idx] =(reg[0], cigar, leftPos, rightPos, -1)
#                    lastInsertionEnd = regions[idx][3]
#                    if VERBOSITY > 1:
#                        log("tgt reg start: %d, end: %d, pos: %d.\n" 
#                            % (regions[idx][2], regions[idx][3], regions[idx][4]))
#                        log("tgt reg cigar: '%s'.\n" 
#                            % cu.toString(regions[idx][1]))
#
            cb = cigarbuilder.CigarBuilder()
#            cigar=cb.build(regions)
            for reg in regions:
                cb.append(reg)
            cigar = cb.cigar
            
            
            if VERBOSITY > 1:                                                                
                log("\nregions:\n")
                log('\n'.join(map(str,regions)))
                log("\n\nread cigar: '%s' ; %s\n\n" 
                    % (cu.toString(cigar), str(cigar)))
                                                                                                                                    
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
                log("s: %d, i: %d, d: %d.\n" % (nSNPs, nInsertions, nDeletions))
                log("new read pos: %d.\n" % rseq.pos)
                log("new read cigar: '%s'.\n========================\n\n" 
                    % cu.toString(rseq.cigar))
            count += 1
            
            if VERBOSITY > 0 and not TESTING:                
                if count % 10000 == 0:
                    log("progress: %3d%%\r" % (count*100/self.nReads), 1, True)        
                
        if VERBOSITY > 0 and not TESTING:
            log("progress: %3d%%\n" % (count*100/self.nReads), 1, True)
                                
        if TESTING:              
            return results
        log("%d read(s) have variants\n" % count2, 1, True)
        return count
