
DEL='d' #Deletion
INS='i' #Insertion
SUB='s' #Substitution/SNP
INV='v' #Inversion
DUP='u' #Duplication
TRN='t' #Translocation

#Inversion not supported yet.
VARIANT_TYPES=set([DEL,SUB,INS,INV,DUP,TRN])

class Variant:
    def __init__(self, chrom, offset, length, vtype, extra=None):        
        self.start = (str(chrom), int(offset)) #offset is 0-based
        self.length = int(length)
        self.type = vtype
        self.extra = extra
        
    def __str__(self):        
        return "%s:%d-%d:%s" % (self.start[0], self.start[1], 
                                self.start[1]+self.length-1, self.type)
    
    def ragic(self):
        return "%d%s" %(self.length,self.type)

    
class SNP(Variant):
    def __init__(self, chrom, offset, allele):
        Variant.__init__(self, chrom, offset, 1, SUB, allele)

    def __str__(self):
        return "%s:%d:%s:%s" % (self.start[0], self.start[1], self.type, 
                                self.extra)


class Deletion(Variant):    
    def __init__(self, chrom, offset, seq):
        Variant.__init__(self, chrom, offset, len(seq), DEL, seq)
        
    def __str__(self):
        return "%s:%d-%d:%s" % (self.start[0], self.start[1], 
                                self.start[1]+self.length-1, self.type)

    
class Insertion(Variant):
    def __init__(self, chrom, offset, seq):
        Variant.__init__(self, chrom, offset, len(seq), INS, seq)

    def __str__(self):
        return "%s:%d:%s:%s" % (self.start[0], self.start[1], self.type, 
                                self.extra)


class Inversion(Variant):
    def __init__(self, chrom, offset, length):
        Variant.__init__(self, chrom, offset, length, INV, None)


class Duplication(Variant):
    def __init__(self, chrom, offset, length, times):
        Variant.__init__(self, chrom, offset, length, DUP, None)
        self.times = times


class Translocation:
    def __init__(self, chrom1, offset1, length1, chrom2, offset2, length2):
        self.r1 = Variant(chrom1, offset1, length1)
        self.r2 = Variant(chrom2, offset2, length2)



def variantFactory(chrom, offset, length, vtype, extra=None):              
#    if type == MAT:
#        assert extra == None and length > 0
#        return Match(chrom, offset, length)
    if vtype == SUB:
        assert length == 1 and len(extra)>=1
        return SNP(chrom, offset, extra)
    if vtype == DEL:
        assert length == len(extra) and length > 0
        return Deletion(chrom, offset, extra)
    if vtype == INS:
        assert length == len(extra) and length > 0
        return Insertion(chrom, offset, extra)
    raise ValueError("Unknown variant type '%s'" % vtype)



#def parseVariant(chrom, offset, ref, alt):
#    if ref == alt:
#        if len(ref) == 1:
#            return SNP(chrom, offset, "%s/%s" % (ref,alt))
#        return None
#        
#    refLen = len(ref)
#    altLen = len(alt)
#    if refLen == 1 and altLen == 1:
#        return SNP(chrom, offset, "%s/%s" % (ref,alt))
#    
#    left = 0
#    while left < refLen and left < altLen:
#        if ref[left] != alt[left]:
#            break
#        left += 1
#    else:
#        if left == refLen:
#            return Insertion(chrom, offset+left-1, alt[left:])
#        else: #left == altLen
#            #self.assertEqual(str pos,left    
#            return Deletion(chrom, offset+left, ref[left:])
#            #return Deletion(chrom, offset+left, len(ref[left:]))
#    
#    right = 0
#    while left <= refLen-1-right and left <= altLen-1-right:
#        if ref[refLen-1-right] != alt[altLen-1-right]:
#            break
#        right += 1
#    else:
#        if left > refLen-1-right:
#            return Insertion(chrom, offset+left-1, alt[left:altLen-right])
#        else: #left > altLen-1-right
#            return Deletion(chrom, offset+left, ref[left:refLen-right])
#            #return Deletion(chrom, offset+left, len(ref[left:refLen-right]))
#    
#    if left+right == refLen-1 and left+right == altLen-1:
#        return SNP(chrom, offset+left, "%s/%s" % (ref[left],alt[left]))
#    else:
#        raise ValueError("Cannot parse variant. (ref:%s,alt:%s)" % (ref,alt))



def parseVariant(chrom, offset, ref, alt):
    if ref == alt:
        if len(ref) == 1:
            return SNP(chrom, offset, "%s/%s" % (ref,alt))
        return None

    refLen = len(ref)
    altLen = len(alt)
    
    if refLen == 1 and altLen == 1:
        return SNP(chrom, offset, "%s/%s" % (ref,alt))
    
    assert ref[0] ==alt[0]  # The first base of ref and alt should be the same
            
    left = 1 # Skip the first position
    right = 0
    
    # Match from the right
    while left <= refLen-1-right and left <= altLen-1-right:
        if ref[refLen-1-right] != alt[altLen-1-right]:
            break
        right += 1
    else:
        if left > refLen-1-right:
            return Insertion(chrom, offset+left-1, alt[left:altLen-right])
        else: #left > altLen-1-right
            return Deletion(chrom, offset+left, ref[left:refLen-right])
            #return Deletion(chrom, offset+left, len(ref[left:refLen-right]))
    
    # Match from the left
    while left <= refLen-1-right and left <= altLen-1-right:
        if ref[left] != alt[left]:
            break
        left += 1
    else:
        if left > refLen-1-right:
            return Insertion(chrom, offset+left-1, alt[left:altLen-right])
        else: #left > altLen-1-right        
            return Deletion(chrom, offset+left, ref[left:refLen-right])
            #return Deletion(chrom, offset+left, len(ref[left:]))    
    
    if left+right == refLen-1 and left+right == altLen-1:
        return SNP(chrom, offset+left, "%s/%s" % (ref[left],alt[left]))
    else:
        raise ValueError("Cannot parse variant. (ref:%s,alt:%s)" % (ref,alt))
    
    

if __name__ == '__main__':
    print(parseVariant(8,1,'TACACACACACACACACAAACACACACACACACACAC','TACAC'))
#    print(parseVariant2(8,1,'TACACACACACACACACAAACACACACACACACACAC','TACAC'))
#    print(parseVariant(1,100,'ACA','ACCTAA'))
#    print(parseVariant(1,100,'ACTA','ACTTATA'))
#    print(parseVariant(1,100,'ACCA','ACCTATA'))
#    print(parseVariant(1,100,'ATAT','AT'))    
#    print(parseVariant(1,100,'ATACT','ATA'))
#    print(parseVariant(1,100,'ATCT','ATTAAT'))
