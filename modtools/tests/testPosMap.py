'''
Created on Oct 11, 2012

@author: Shunping Huang
'''


import unittest
from modtools import posmap

class TestMap(unittest.TestCase):
    '''Test Case 1
         10M  5D     10I    10D    10M    5I     10D   10M
    ref: 0-9, 10-14, -14  , 15-24, 25-34, -34  , 35-44, 45-54
    new: 0-9, -9   , 10-19, -19  , 20-29, 30-34, -34  , 35-44 
    '''       
    def setUp(self):
        data = [('1', 0, '1', 0, 10, '+'), 
                ('1', 10, '1', -9, 5, '+'), 
                ('1', -14, '1', 10, 10, '+'), 
                ('1', 15, '1', -19, 10, '+'), 
                ('1', 25, '1', 20, 10, '+'), 
                ('1', -34, '1', 30, 5, '+'), 
                ('1', 35, '1', -34, 10, '+'), 
                ('1', 45, '1', 35, 10, '+')]
        
        self.posmap = posmap.PosMap()
        self.posmap.load(data, isConverted = True)
    
      
    def test_load(self):
        self.assertEqual(self.posmap.data[1][0], ('1',10))
        self.assertEqual(self.posmap.data[1][1], ('1',-9))
        self.assertEqual(self.posmap.data[1][2], 5)
        self.assertEqual(self.posmap.data[1][3], '+')
        
        self.assertEqual(self.posmap.data[7][0], ('1',45))
        self.assertEqual(self.posmap.data[7][1], ('1',35))
        self.assertEqual(self.posmap.data[7][2], 10)
        self.assertEqual(self.posmap.data[7][3], '+')
    
    
    def test_buildMap(self):
        self.posmap.build()
        self.assertEqual(self.posmap.fkeys, [('1', -34), ('1', -14), ('1', 0), 
                                             ('1', 10), ('1', 15), ('1', 25), 
                                             ('1', 35), ('1', 45)])
    
        self.assertEqual(self.posmap.bkeys, [('1', -34), ('1', -19), ('1', -9), 
                                             ('1', 0), ('1', 10), ('1', 20), 
                                             ('1', 30), ('1', 35)])    
    
    
    def test_fmap(self):
        self.posmap.build()
        self.assertEqual(self.posmap.fmap(('1',0)), ('1',0))
        self.assertEqual(self.posmap.fmap(('1',5)), ('1',5))
        self.assertEqual(self.posmap.fmap(('1',9)), ('1',9))
        self.assertEqual(self.posmap.fmap(('1',10)), ('1',-9))
        self.assertEqual(self.posmap.fmap(('1',12)), ('1',-9))
        self.assertEqual(self.posmap.fmap(('1',14)), ('1',-9))
        self.assertEqual(self.posmap.fmap(('1',15)), ('1',-19))
        self.assertEqual(self.posmap.fmap(('1',20)), ('1',-19))
        self.assertEqual(self.posmap.fmap(('1',24)), ('1',-19))
        self.assertEqual(self.posmap.fmap(('1',25)), ('1',20))
        self.assertEqual(self.posmap.fmap(('1',30)), ('1',25))
        self.assertEqual(self.posmap.fmap(('1',34)), ('1',29))
        self.assertEqual(self.posmap.fmap(('1',35)), ('1',-34))
        self.assertEqual(self.posmap.fmap(('1',40)), ('1',-34))
        self.assertEqual(self.posmap.fmap(('1',44)), ('1',-34))
        self.assertEqual(self.posmap.fmap(('1',45)), ('1',35))
        self.assertEqual(self.posmap.fmap(('1',50)), ('1',40))
        self.assertEqual(self.posmap.fmap(('1',54)), ('1',44))
        
        self.assertRaises(ValueError, self.posmap.fmap, ('0',0))
        self.assertRaises(AssertionError, self.posmap.fmap, ('1',-1))
        self.assertRaises(ValueError, self.posmap.fmap, ('1',55))
    
    
    def test_bmap(self):
        self.posmap.build()
        self.assertEqual(self.posmap.bmap(('1',0)), ('1',0))
        self.assertEqual(self.posmap.bmap(('1',2)), ('1',2))
        self.assertEqual(self.posmap.bmap(('1',5)), ('1',5))
        self.assertEqual(self.posmap.bmap(('1',9)), ('1',9))
        self.assertEqual(self.posmap.bmap(('1',10)), ('1',-14))
        self.assertEqual(self.posmap.bmap(('1',15)), ('1',-14))
        self.assertEqual(self.posmap.bmap(('1',19)), ('1',-14))
        self.assertEqual(self.posmap.bmap(('1',20)), ('1',25))
        self.assertEqual(self.posmap.bmap(('1',25)), ('1',30))
        self.assertEqual(self.posmap.bmap(('1',29)), ('1',34))        
        self.assertEqual(self.posmap.bmap(('1',30)), ('1',-34))
        self.assertEqual(self.posmap.bmap(('1',32)), ('1',-34))
        self.assertEqual(self.posmap.bmap(('1',34)), ('1',-34))
        self.assertEqual(self.posmap.bmap(('1',35)), ('1',45))
        self.assertEqual(self.posmap.bmap(('1',40)), ('1',50))
        self.assertEqual(self.posmap.bmap(('1',44)), ('1',54))
                
        self.assertRaises(ValueError, self.posmap.bmap, ('0',0))
        self.assertRaises(AssertionError, self.posmap.bmap, ('1',-1))
        self.assertRaises(ValueError, self.posmap.bmap, ('1',45))
              
        
if __name__ == '__main__':    
    unittest.main()
