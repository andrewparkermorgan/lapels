'''
Created on Oct 11, 2012

@author: Shunping Huang
'''

import unittest
import StringIO
import csv
from modtools import mod

class TestMod1(unittest.TestCase):
    '''Test Case 1
         10M  5D     10I    10D    10M    5I     10D   10M
    ref: 0-9, 10-14, -14  , 15-24, 25-34, -34  , 35-44, 45-54
    new: 0-9, -9   , 10-19, -19  , 20-29, 30-34, -34  , 35-44
    '''

    def setUp(self):
        modFile = StringIO.StringIO('''d\t1\t10\t1
d\t1\t11\t2
d\t1\t12\t3
d\t1\t13\t4
d\t1\t14\t5
i\t1\t14\tabcdefghij
d\t1\t15\t1
d\t1\t16\t2
d\t1\t17\t3
d\t1\t18\t4
d\t1\t19\t5
d\t1\t20\t6
d\t1\t21\t7
d\t1\t22\t8
d\t1\t23\t9
d\t1\t24\t0
i\t1\t34\tabcde
d\t1\t35\t1
d\t1\t36\t2
d\t1\t37\t3
d\t1\t38\t4
d\t1\t39\t5
d\t1\t40\t6
d\t1\t41\t7
d\t1\t42\t8
d\t1\t43\t9
d\t1\t44\t0
''')

        fp =  csv.reader(modFile, delimiter='\t')
        self.mod = mod.Mod('1', 55, fp)


    def test_buildMap(self):
        self.mod.buildPosMap()
        self.assertEqual(self.mod.posmap.data, [(('1', 0), ('1', 0), 10, '+'),
                                                (('1', 10), ('1', -9), 5, '+'),
                                                (('1', -14), ('1', 10), 10, '+'),
                                                (('1', 15), ('1', -19), 10, '+'),
                                                (('1', 25), ('1', 20), 10, '+'),
                                                (('1', -34), ('1', 30), 5, '+'),
                                                (('1', 35), ('1', -34), 10, '+'),
                                                (('1', 45), ('1', 35), 10, '+')])


    def test_getPosMap(self):

        posmap = self.mod.getPosMap()
        self.assertEqual(posmap.data, [(('1', 0), ('1', 0), 10, '+'),
                                       (('1', 10), ('1', -9), 5, '+'),
                                       (('1', -14), ('1', 10), 10, '+'),
                                       (('1', 15), ('1', -19), 10, '+'),
                                       (('1', 25), ('1', 20), 10, '+'),
                                       (('1', -34), ('1', 30), 5, '+'),
                                       (('1', 35), ('1', -34), 10, '+'),
                                       (('1', 45), ('1', 35), 10, '+')])


    def test_buildSeqs(self):
        refSeqs=dict()
        refSeqs['1']=''.join(['ABCDEFGHIJK']*5)
        self.assertEqual(len(refSeqs['1']), 55)
        self.mod.buildSeq(refSeqs)
        self.assertEqual(len(self.mod.seq), 45)
        self.assertEqual(self.mod.seq,'ABCDEFGHIJabcdefghijDEFGHIJKABabcdeBCDEFGHIJK')


    def test_getSeq(self):
        refSeqs=dict()
        refSeqs['1']=''.join(['ABCDEFGHIJK']*5)
        seq = self.mod.getSeq(refSeqs)
        self.assertEqual(len(seq), 45)
        self.assertEqual(seq,'ABCDEFGHIJabcdefghijDEFGHIJKABabcdeBCDEFGHIJK')



if __name__ == '__main__':
    unittest.main()
