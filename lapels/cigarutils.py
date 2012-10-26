'''
The module of cigar utilities. 

Created on Oct 7, 2012

@author: Shunping Huang
'''

fMap = {'M':0, 'I':1, 'D':2, 'N': 3}
bMap = dict([(tup[1],tup[0]) for tup in fMap.items()]) 


def toString(cigar):
    '''Convert a cigar to its string representation.'''    
    segs = []        
    for op, length in cigar:        
        if length <= 0:
            continue
        #assert length > 0
        try :        
            segs.append("%d%s" % (length, bMap[op]))
        except KeyError:
            raise NotImplementedError("unknown op %d in cigar." % op)                    
    return ','.join(segs)


def make(sCigar):
    '''Create a cigar from a string'''
    cigar = []
    if len(sCigar) > 0:        
        for seg in sCigar.split(','):
            try:        
                op=fMap[seg[-1]]
                length = int(seg[:-1])
                assert length > 0
                cigar.append((op,length))
            except KeyError:
                raise NotImplementedError("unknown type %s in cigar." % seg[-1])
    return cigar


def sub(cigar, cpos=0, start=None, end=None):
    '''Get a sub cigar by start and end''' 
    assert cpos >= 0
#    print(start,end)
    if (start is None and end == cpos - 1) or (start == cpos and end == None)\
            or (len(cigar) == 0):
        return []
    
    idx1 = 0
    offset1 = 0
    idx2 = len(cigar) - 1
    offset2 = cigar[-1][1] - 1
    pos = cpos
    iBuffer = []
       
    if start is not None:
        if start < pos:
            raise ValueError("start position underflows.")        
        for op, length in cigar:
            if op == 0 or op == 7 or op == 8 or op == 2 or op == 3:
                if start < pos + length:
                    offset1 = start - pos                    
                    break
                else:
                    pos += length
                    idx1 += 1
                    iBuffer = []
            else:
                if op == 1:         ##store the insertion before
                    iBuffer.append((op, length))
                idx1 += 1
        
#        print start, end, pos
        ## start position not overflows                        
        if end == pos - 1 and start == pos:
#            print start, end, pos
            return []
         
        if idx1 >= len(cigar):
            raise ValueError("start position overflows.")
        
            
    if end is not None:
        if end < pos:
            raise ValueError("end position underflows.")            
            
        idx2 = idx1
        for op, length in cigar[idx1:]:
            if op == 0 or op == 7 or op == 8 or op == 2 or op == 3:
                if end < pos + length:
                    offset2 = end - pos                    
                    break
                else:
                    pos += length
                    idx2 += 1
            else:
                idx2 += 1
        if idx2 >= len(cigar):
            raise ValueError("end position overflows.")
    
    ret = []
    ## If the cigar is ...IM..., and the region start at the first base of M,
    ## then I is included.
    if len(iBuffer) > 0 and offset1 == 0: 
        ret.extend(iBuffer)
    if idx1 == idx2:
        assert offset1 <= offset2
        ret.append((cigar[idx1][0],offset2-offset1+1))
    else:
        ret.append((cigar[idx1][0],cigar[idx1][1]-offset1))
        for cig in cigar[idx1+1:idx2]:
            ret.append(cig)
        ret.append((cigar[idx2][0],offset2+1))
    return ret
    #return (idx1,offset1,idx2,offset2)



def replace(cigar, cpos, regions):
    '''Replace a sorted list of regions in a cigar.'''
    ret = []
    last = None
    for reg in regions:        
        cig = sub(cigar, cpos, last, reg[2]-1)
        ret.extend(cig)
#        print "what to append",reg[3]
        #ret.append(reg[3])
        ret.extend(reg[1])
        last = reg[3] + 1
    try:
        cig = sub(cigar, cpos, last, None)
    except ValueError:
        cig = []
    ret.extend(cig)    
    return ret
            

def simplify(cigar):
    '''To simplify a cigar by combining adjacent regions of the same type.'''    
    ## Find the first non-negative op region.
    i = 0
    while i < len(cigar):
        if cigar[i][0] >= 0:
            break
        i += 1
    else:
        return []
    
    ret = []
    buf = cigar[i]
    for cig in cigar[i+1:]:        
        if cig[0] < 0: ## ignore negative op region.
            continue
        if buf[0] == cig[0]:
            buf = (buf[0], buf[1]+cig[1])
        else:
            ret.append(buf)
            buf = cig
    ret.append(buf)
    return ret
        