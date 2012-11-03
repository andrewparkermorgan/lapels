'''
The module of position mapping. 

Created on Sep 20, 2012

@author: Shunping Huang
'''

import bisect
import gc

__all__ = ['PosMap']

class PosMap:
    def __init__(self, dataIter=None):
        '''Initialize a position map.
        dataIter: an iterator of tuples extracted from a mod file. Each tuple
            should have 4 fields.
        '''
            
        self.fvals = None  #forward mapping
        self.bvals = None  #backward mapping        
        self.fkeys = []
        self.bkeys = []        
        self.data = None
        if dataIter is not None:
            self.loadData(dataIter)  
        
    
    def load(self, dataIter, isConverted=False):
        assert dataIter is not None
        gc.disable()                
        self.data = []
        append = self.data.append
        if isConverted:
            for row in dataIter:
                append(((row[0], row[1]), (row[2], row[3]), row[4], row[5]))
        else:
            for row in dataIter:
                append(((row[0], int(row[1])), (row[2], int(row[3])), int(row[4]), row[5]))
        self.build()
        gc.enable()


    def build(self):                
        data = self.data
        self.fvals = sorted(data, key=lambda tup: tup[0])
        self.fkeys = [tup[0] for tup in self.fvals]
        self.bvals = sorted(data, key=lambda tup: tup[1])
        self.bkeys = [tup[1] for tup in self.bvals]        
    
    
    ##FIX ME!!!! Overlapping regions make the index from bisect not correct. 
    def fmap(self, pos):
        '''Mapping a position from reference to in silico genome.'''
        assert self.fkeys is not None
        assert pos[1] >= 0            
        i = bisect.bisect_right(self.fkeys, pos) - 1
        if i < 0 or self.fvals[i][0][1] < 0:
            raise ValueError("Error: Reference position %d underflows." % pos[1])                   
        #print(self.fvals[i])
        refpos = self.fvals[i][0]
        newpos = self.fvals[i][1]        
        length = self.fvals[i][2]
        direction = self.fvals[i][3][0]      
        if pos[0] == refpos[0]:
            if pos[1] >= refpos[1] and pos[1] < refpos[1] + length:
                ## May need to consider direction!!!
                ## Deletion, return the inverse of newest preceding position.                                    
                if newpos[1] < 0:       
                    return newpos
                                        
                if direction == '+':    ## Regular region
                    d = newpos[1] - refpos[1]         
                    return (newpos[0], pos[1] + d)                                    
                else:                   ## Inverted region
                    s = newpos[1] + refpos[1] + length - 1                
                    return (newpos[0], s - pos[1])                    
            else:
                raise ValueError("Error: Reference position %d overflows." % pos[1])
        else:
            raise ValueError("Error: Reference chromosome %s not found." % pos[0])
    
    
    ##FIX ME!!!! Overlapping regions make the index from bisect not correct.
    def bmap(self, pos):
        '''Mapping a position from in silico genome to reference'''
        assert self.fkeys is not None
        assert pos[1] >= 0        
        i=bisect.bisect_right(self.bkeys, pos)-1
        if i < 0 or self.bvals[i][1][1] < 0:
            raise ValueError("Error: In silico position %d underflows." % pos[1])                            
        #print(self.bvals[i])
        refpos = self.bvals[i][0]
        newpos = self.bvals[i][1]
        length = self.bvals[i][2]
        direction = self.bvals[i][3][0]
        if pos[0] == newpos[0]:
            if pos[1] >= newpos[1] and pos[1] < newpos[1] + length:
                ## May need to consider direction!!!
                ## Insertion, return the inverse of newest preceding position.
                if refpos[1] < 0:                 
                    return refpos
                
                if direction == '+':    #regular direction
                    d = refpos[1] - newpos[1]                
                    return (refpos[0], pos[1] + d)                                
                else:                   #reversed direction
                    s = refpos[1] + newpos[1] + length - 1                
                    return (refpos[0], s - pos[1])
            else:
                raise ValueError("Error: In silico position %d overflows." % pos[1])                
        else:
            raise ValueError("Error: In silico chromosome not found.")


    def toCSV(self):
        out = []
        append = out.append
        join = '\t'.join
        for row in self.data:
            append(join((row[0][0], str(row[0][1]), row[1][0], str(row[1][1]), str(row[2]), row[3])))
        return '\n'.join(out)


