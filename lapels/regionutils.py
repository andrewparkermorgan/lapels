'''
The module of manipulating regions. 

Created on Oct 8, 2012

@author: Shunping Huang
'''

import cigarutils


def checkRegion(region):
    '''sanity check for a region'''
    assert len(region) >= 4 
    cigar = region[1]
    nRefBases = 0
    for op, length in cigar:
        if op == 0 or op == 2 or op == 3:
            nRefBases += length

    if region[2]+ nRefBases - 1 != region[3]:
        raise ValueError("Error: cigar '%s' conflicts with read region %d-%d."
                         % (cigarutils.toString(cigar), region[2], region[3]))

        
def makeReadRegion(rop, sCigar, start, end):
    cigar = cigarutils.make(sCigar)
    region = (rop, cigar, start, end)
    checkRegion(region)
    return region


def isRightTo(r1, r2):
    '''return true if r1 is right to r2'''
    checkRegion(r1)
    checkRegion(r2)    
    if r2[2] <= r2[3]:          ## r2 is not an insertion
        return r1[2] > r2[3]
    else:                       ## r2 is an insertion
        if r1[2] <= r1[3]:      ## r1 is not an insertion
            return r1[2] > r2[3]
        else:                   ## r1 is an insertion
            return r1[3] > r2[3]
    

def isLeftTo(r1, r2):
    '''return true if r1 is left to r2'''
    checkRegion(r1)
    checkRegion(r2)
    if r2[2] <= r2[3]:
        return r1[3] < r2[2]
    else:
        if r1[2] <= r1[3]:
            return r1[3] < r2[2]
        else:
            return r1[2] < r2[2]


def modifyRegion(region):
    ''' Modify a region according to the region op (rop). 
    A region is a tuple of (rop, cigar, start, end, pos, #SNPs, #Ins, #Del)
    '''    
    rop = region[0]
    cigar = region[1]    
    if rop == 0:                     ## M_1
        ret = region    
    elif rop == 2:                   ## D_1
        for i, (op, length) in enumerate(cigar):
            if op == 0:         ## M_0
                cigar[i] = (2, length)
            elif op == 1:       ## I_0
                cigar[i] = (-1, length)
            elif op == 2:       ## D_0
                cigar[i] = (2, length)
        
        cigar = cigarutils.simplify(cigar)
        ret = (rop, cigar, region[2], region[3], -1, 0, 0, 0)                
        
    elif rop == 3:                   ## N_1
        for i, (op, length) in enumerate(cigar):
            if op == 0:         ## M_0
                cigar[i] = (3, length)
            elif op == 1:       ## I_0
                cigar[i] = (-1, length)
            elif op == 2:       ## D_0
                cigar[i] = (3, length)  #D_0 in N_1 should be N
                #cigar[i] = (2, length)                        
        cigar = cigarutils.simplify(cigar)            
        ret = (rop, cigar, region[2], region[3], -1, 0, 0, 0)
    
    else:    
        raise ValueError("Error: unknown op type %d" % rop)
                    
    return ret


