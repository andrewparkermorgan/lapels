import collections
import pysam
import os
import gzip

CHR = 0
POS = 1
REF = 3
ALT = 4
FMT = 8

GT = 'GT'
REF_ALIAS ='.'
ALT_FS = ','        
FS = '\t'    
FMT_FS = ':'

VERBOSITY = 0
    
__all__ = ['parseFormat', 'getGenotype', 'parseGenotype', 
           'VCFIterator', 'VCFReader']

def parseFormat(fmt, data):
    formatFields = fmt.split(FMT_FS)
    dataFields = data.split(FMT_FS)
    
    ret = dict()
    for i in range(min(len(formatFields), len(dataFields))):
        ret[formatFields[i]] = dataFields[i]
    return ret


def getGenotype(genotypeStr):
    if '/' in genotypeStr:
        return genotypeStr.split('/')
    if '|' in genotypeStr:
        return genotypeStr.split('|')
    return [genotypeStr]


def parseGenotype(ref, alt, genotype):
    altFields = alt.split(ALT_FS)
    if genotype == '0':
        return ref
    else:
        return altFields[int(genotype) - 1]


class VCFIterator(collections.Iterator):
    
    def __init__(self, parent, fetched):
        self.parent = parent
        self.fetched = fetched
    
    
    def __iter__(self):
        return self
    
    
    def next(self):            
        try:
            while True:
                line = self.fetched.next().rstrip()             
                                  
                fields = line.split(FS)
                if len(fields) != self.parent.nColumns:
                    raise ValueError("Number of columns not consistent. (%s)" %
                                     line)                
                                
                genotypes = []
                for i in self.parent.sampleIndexes:
                    fmtFields = parseFormat(fields[FMT], fields[i])                
                    genotype = getGenotype(fmtFields[GT])
                    if len(set(genotype)) > 1:
                        if VERBOSITY > 1:
                            print("Hets found in %s:%s of sample %s (%s). " % 
                                  (fields[CHR], fields[POS], 
                                   self.parent.samples[i], fmtFields[GT]) + 
                                  "Use the first allele.")
                    genotype = genotype[0]
                    if genotype == REF_ALIAS:
                        genotype = '0'                
                    genotypes.append(genotype) 
                
                # For only one sample, append a dumb ref genotype to compare
                if len(self.parent.sampleIndexes) == 1:
                    genotypes.append('0')
                                                                
                isVariant = False
                for i in range(len(genotypes)-1):
                    if genotypes[i] != genotypes[i+1]:
                        isVariant = True
                        break
                
                if not isVariant:
                    continue
                
                # Remove the dumb ref genotype
                if len(self.parent.sampleIndexes) == 1:
                    del genotypes[-1]
                            
                ret = []
                ret.append(fields[CHR])
                ret.append(int(fields[POS])-1)  # convert to 0-based
                ret.append(fields[REF])
                ret.extend([parseGenotype(fields[REF], fields[ALT], genotype) for genotype in genotypes])
                return ret
            
        except StopIteration:            
            raise StopIteration()


class VCFReader():            
    def __init__(self, fileName, samples):        
        self.samples = samples
        self.sampleIndexes = []
        self.nColumns = 0      
        
        # Compress with bgzip
        if not fileName.endswith('.gz'):
            if not os.path.isfile(fileName+'.gz'):
                pysam.tabix_compress(fileName, fileName+'.gz')
            fileName += '.gz'
        
        # Build tabix index
        if not os.path.isfile(fileName+'.tbi'): 
            pysam.tabix_index(fileName, preset='vcf')                                 
                 
        nLines = 0        
        fp = gzip.open(fileName, 'r')
        line = fp.readline()        
        while line:
            nLines += 1          
            if line.startswith('##'):
                line = fp.readline()                                    
            elif line.startswith('#'):  # Header line
                break
            else:
                line = None        # Content line, no header line found
        else:
            raise ValueError("Header not found.")
        
        # Get the column index of selected samples        
        headers = line[1:].rstrip().split(FS)
        self.nColumns = len(headers)        
        if self.nColumns <= 9:
            raise ValueError("Not enough columns in header.")                
        
        for name in self.samples:
            if name in headers[9:]:
                self.sampleIndexes.append(headers.index(name))
            else:
                raise ValueError("Sample %s not found in header." % name)
        
        self.tabix = pysam.Tabixfile(fileName)
        self.chroms = self.tabix.contigs
        self.fileName = fileName

                        
    def fetch(self, region):
        return VCFIterator(self, self.tabix.fetch(region=region))
                        



