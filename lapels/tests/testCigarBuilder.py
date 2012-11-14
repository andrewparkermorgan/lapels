'''
Created on Oct 8, 2012

@author: Shunping Huang
'''

from lapels import cigarutils
from lapels import cigarbuilder
from lapels.regionutils import makeReadRegion

import unittest

class TestCigarBuilder(unittest.TestCase):
    
    def test1(self):
        regions = []
        add = regions.append
        add(makeReadRegion(1, '2I', -1, -2))## an inserted I_1 before any M_1    
        add(makeReadRegion(0, '4M', 0, 3))
        add(makeReadRegion(1, '2I', -1, -1))  ## an inserted I_1 after an M_1
        add(makeReadRegion(1, '3I', -1, -1))  ## an inserted I_1 after an inserted I_1
                                            ## a gap
        add(makeReadRegion(1, '1M', 5, 5))  ## a match I_1 after is an inserted I_1
        add(makeReadRegion(1, '4I', 6, 5))  ## an inserted I_1 after a match I_1   
                                            ## a gap here 
        add(makeReadRegion(0, '2M', 7, 8))    
        add(makeReadRegion(0, '2I', 9, 8))  ## an M_1 that is an insertion  
        add(makeReadRegion(1, '3I', 9, 8))    
        add(makeReadRegion(0, '4I', 9, 8))  ## an M_1 that is an insertion
        add(makeReadRegion(1, '5I', 9, 8))
                                            ## a gap here
        add(makeReadRegion(0, '5M', 10, 14))    
        cb = cigarbuilder.CigarBuilder()
        cigar = cb.build(regions)
        print(cigar)
        self.assertEqual(cigar, [(1, 2), (0, 4), (1, 2), (1, 3), (2, 1), 
                                 (0, 1), (1, 4), (2, 1), (0, 2), (1, 2), 
                                 (1, 3), (1, 4), (1, 5), (2, 1), (0, 5)])
        self.assertEqual(cigarutils.toString(cigar), '2I,4M,2I,3I,1D,1M,4I,1D,2M,2I,3I,4I,5I,1D,5M')        

        
    def test2(self):
        regions = []
        add = regions.append
        add(makeReadRegion(1, '1I', -1, -2))
        add(makeReadRegion(2, '10I', 5, 4)) ## D_1 of I_0
        add(makeReadRegion(1, '2I', 5, 4))
        add(makeReadRegion(1, '3I', 5, 4))
        add(makeReadRegion(2, '20I', 5, 4,)) ## D_1 of I_0
        add(makeReadRegion(1, '4I', 5, 4 ))
        add(makeReadRegion(1, '4I', 7, 6 ))        
        cb = cigarbuilder.CigarBuilder()
        cigar = cb.build(regions)
        self.assertEqual(cigar, [(1, 1), (1, 2), (1, 3), (1, 4), (1, 20), (2, 2), (1, 4)])
        
    
    def test3(self):
        regions = []
        add = regions.append
        add(makeReadRegion(0, '25M', 10, 34))
        add(makeReadRegion(1, '1I', 35, 34))
        add(makeReadRegion(3, '1N', 35, 35))
        add(makeReadRegion(0, '74M', 36, 109))
        cb = cigarbuilder.CigarBuilder()
        cigar = cb.build(regions)        
        self.assertEqual(cigarutils.simplify(cigar),[(0, 25), (1, 1), (3, 1), (0, 74)])
    
    
    def test4(self):
        regions = []
        add = regions.append
        add(makeReadRegion(0, '25M', 10, 34))
        add(makeReadRegion(1, '1M,1I', 35, 35))        
        add(makeReadRegion(3, '1N', 35, 35))
        add(makeReadRegion(0, '74M', 36, 109))
        cb = cigarbuilder.CigarBuilder()
        cigar = cb.build(regions)
        self.assertEqual(cigarutils.simplify(cigar),[(0,26),(1,1),(0,74)])        
    
    
    def test5(self):
        regions = []
        add = regions.append
        add(makeReadRegion(0, '25M', 10, 34))            
        add(makeReadRegion(3, '1N', 35, 35))
        add(makeReadRegion(1, '1M,1I', 35, 35))
        add(makeReadRegion(0, '74M', 36, 109))
        cb = cigarbuilder.CigarBuilder()
        cigar = cb.build(regions)
        self.assertEqual(cigarutils.simplify(cigar),[(0,26),(1,1),(0,74)])   
        
        
    
if __name__ == '__main__':
    unittest.main()
