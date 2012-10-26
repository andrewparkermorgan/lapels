'''
Created on Oct 8, 2012

@author: Shunping Huang
'''


from lapels.cigarutils import make, toString, sub, replace

import unittest

class TestSubCigar(unittest.TestCase):
    def setUp(self):
        self.cigar1 = [(0,15),(1,5),(2,10),(0,5)]
        self.cigar2 = [(1,10),(0,15),(1,5),(2,10),(0,5),(1,10)]
    
    
    def testToString(self):
        self.assertEqual(toString(self.cigar1),'15M,5I,10D,5M')
        self.assertEqual(toString(self.cigar2),'10I,15M,5I,10D,5M,10I')
        
        
    def testMakeCigar(self):
        self.assertEqual(make(toString(self.cigar1)), self.cigar1)
        self.assertEqual(make(toString(self.cigar2)), self.cigar2)
        
                
    def testSubCigar(self):        
        self.assertEqual(sub(self.cigar1,20),self.cigar1)
        self.assertEqual(sub(self.cigar1,20,20,49),self.cigar1)
        self.assertEqual(sub(self.cigar1,20,22,26),[(0,5)])
        self.assertEqual(sub(self.cigar1,20,38,38),[(2,1)])
        self.assertEqual(sub(self.cigar1,20,21,45),[(0,14),(1,5),(2,10),(0,1)])
        self.assertEqual(sub(self.cigar1,20,21),[(0,14),(1,5),(2,10),(0,5)])
        self.assertEqual(sub(self.cigar1,20,None,35),[(0,15),(1,5),(2,1)])        
        self.assertRaises(ValueError,sub,self.cigar1,20,10,20)
        self.assertRaises(ValueError,sub,self.cigar1,20,20,100)
        
        self.assertEqual(sub(self.cigar2,20),self.cigar2)
        self.assertEqual(sub(self.cigar2,20,20,49),[(1,10),(0,15),(1,5),(2,10),(0,5)])
        self.assertEqual(sub(self.cigar2,20,22,26),[(0,5)])
        self.assertEqual(sub(self.cigar2,20,38,38),[(2,1)])
        self.assertEqual(sub(self.cigar2,20,21,45),[(0,14),(1,5),(2,10),(0,1)])        
        self.assertEqual(sub(self.cigar2,20,21),[(0,14),(1,5),(2,10),(0,5),(1,10)])                
        self.assertEqual(sub(self.cigar2,20,None,35),[(1,10),(0,15),(1,5),(2,1)])                
        self.assertRaises(ValueError,sub,self.cigar2,20,10,20)
        self.assertRaises(ValueError,sub,self.cigar2,20,20,100)

        self.assertEqual(sub([(1,10)], 5, 5, 4), [])
        self.assertRaises(ValueError,sub,[(1,10)], 5, 4, 5)
        self.assertRaises(ValueError,sub,[(1,10)], 5, 5, 5)
        self.assertRaises(ValueError,sub,[(1,10)], 5, 6, 5)        
        self.assertRaises(ValueError,sub,[(1,10)], 5, 4, 4)
        self.assertRaises(ValueError,sub,[(1,10)], 5, 4, 6)                
        self.assertRaises(ValueError,sub,[(1,10)], 5, 6, 4)
        
        
    def testReplaceCigar(self):
        ##insert before an insertion        
        self.assertEqual(replace(self.cigar1,0,[(1,[(1,10)],15,14),(1,[(1,12)],15,14)]), 
                         [(0, 15), (1, 10), (1, 12), (1, 5), (2, 10), (0,5)])
        ##insert after an deletion
        self.assertEqual(replace(self.cigar1,0,[(1,[(1,10)],25,24)]), 
                         [(0, 15), (1, 5), (2, 10), (1, 10), (0, 5)])          
        ##insert in a match
        self.assertEqual(replace(self.cigar1,0,[(1,[(1,10)],26,25)]),
                         [(0, 15), (1, 5), (2, 10), (0, 1), (1, 10), (0, 4)])  
        ##insert in a match
        self.assertEqual(replace(self.cigar1,0,[(1,[(0,2)],16,17),(1,[(0,3)],20,22)]),
                         [(0, 15), (1, 5), (2, 1), (0,2), (2,2), (0,3), (2,2), (0, 5)])  
#        


if __name__ == '__main__':
    unittest.main()