'''
Test cases for class AnnotateCommand

The construction of the test case are driven by the fact that the target cigar
only has three types of regions: M(M_0), D(D_0), and I(I_0). For a region mapped 
to a target genome, its start and end position will always be in M_0 or I_0, 
because only M_0 and I_0 have target positions. Thus, D_0 will be sandwiched by 
M_0 and I_0. 

There nine types of sandwiches: 
M, MDM, MDI, (M)I, (M)IDM, (M)IDI, (D)I, (D)IDM, (D)IDI

A read can starts at the beginning or the middle of a M/I region, while
it can ends at the middle or end of a M/I region.

The tests will generally verify: the new cigar, the new start and end in 
reference coordinate, the number of snps, insertions, deletions of the target 
being observed. 


Created on Oct 3, 2012

@author: Shunping Huang
'''

import unittest
import StringIO
import tempfile
import os
import pysam

from lapels import annotator as annot
from lapels import cigarutils
from modtools import mod



polish = lambda x: cigarutils.toString(cigarutils.simplify(x))



class Read:
    '''Class for simulating reads from a bam file'''
    def __init__(self, start, end, cigar=None, qlen=None):
        self.qname = 'unknown'
        self.pos = start
        self.aend = end     #one base after the actual end
        self.tags = dict()        
        if cigar is None:
            self.cigar = [(0, self.aend - self.pos)]
        else:
            self.cigar = cigar
        if qlen is None:
            self.qlen = 0
            for op,length in self.cigar:
                if op == 0 or op == 7 or op == 8 or op == 1:
                    self.qlen += length
        else:
            self.qlen = qlen



class TestGetReadOffset(unittest.TestCase):
    '''Test class for getReadOffset() method '''
    def setUp(self):
        pass

    
    def test1(self):
        r = Read(10, 50, [(0,10),(1,5),(0,10),(2,10),(0,10)])
        self.assertRaisesRegexp(ValueError, 'underflows', annot.getReadOffset, r, 1)
        self.assertEquals(annot.getReadOffset(r, 10), 0)
        self.assertEquals(annot.getReadOffset(r, 19), 9)
        self.assertEquals(annot.getReadOffset(r, 20), 15)
        self.assertEquals(annot.getReadOffset(r, 29), 24)
        self.assertRaisesRegexp(ValueError, 'deletion', annot.getReadOffset, r, 30)
        self.assertRaisesRegexp(ValueError, 'deletion', annot.getReadOffset, r, 39)
        self.assertEquals(annot.getReadOffset(r, 40), 25)
        self.assertEquals(annot.getReadOffset(r, 49), 34)
        self.assertRaisesRegexp(ValueError, 'overflows', annot.getReadOffset, r, 50)


    def test2(self):
        # qlen is set wrongly on purpose
        r = Read(10, 50, [(0,10),(1,5),(0,10),(2,10),(0,10)], 30)  
        self.assertRaisesRegexp(ValueError, 'underflows', annot.getReadOffset, r, 1)
        self.assertEquals(annot.getReadOffset(r, 10), 0)
        self.assertEquals(annot.getReadOffset(r, 19), 9)
        self.assertEquals(annot.getReadOffset(r, 20), 15)
        self.assertEquals(annot.getReadOffset(r, 29), 24)
        self.assertRaisesRegexp(ValueError, 'deletion', annot.getReadOffset, r, 30)
        self.assertRaisesRegexp(ValueError, 'deletion', annot.getReadOffset, r, 39)
        self.assertEquals(annot.getReadOffset(r, 40), 25)
        self.assertEquals(annot.getReadOffset(r, 44), 29)
        self.assertRaisesRegexp(ValueError, 'conflict', annot.getReadOffset, r, 45)
        self.assertRaisesRegexp(ValueError, 'conflict', annot.getReadOffset, r, 49)
        self.assertRaisesRegexp(ValueError, 'conflict', annot.getReadOffset, r, 50)                

        
    
class TestAnnotator(unittest.TestCase):    
    ''' Test class for Annotator '''
            
    def setUp(self):
        annot.TESTING = 1
        annot.VERBOSITY = 1
       

    def batchTestHelper(self, modFile, pool, refLens):        
        tmpName = tempfile.mkstemp('.tsv')[1]
        tmpfp = open(tmpName, 'wb')
        for line in modFile:
            tmpfp.write(line)
        tmpfp.close()
        pysam.tabix_index(tmpName, force=True, seq_col=1, start_col=2, end_col=2, 
                      meta_char='#', zerobased=True)
        tmpName += '.gz'
        modFile.close()
        
        self.chromoID = '1'
        self.modobj = mod.Mod(tmpName)
        self.modobj.load(self.chromoID)
        
        for tup in pool:       
            bamIter=[Read(tup[0], tup[1]+1, tup[2]) for tup in pool]        
                                   
        a = annot.Annotator(self.chromoID, refLens[self.chromoID],
                                self.modobj, bamIter)
        results = a.execute()
        
        for i,res in enumerate(results):            
            self.assertEqual(polish(res[0]),pool[i][3])
            self.assertEqual(res[1], pool[i][4])
            self.assertEqual(res[2], pool[i][5])
            self.assertEqual(res[3], pool[i][6])
            self.assertEqual(res[4], pool[i][7])
        
        os.remove(tmpName)
        os.remove(tmpName+'.tbi')
        
              
                
    def test1(self):
        '''Test case for (D)I, (D)ID*M, and (M)ID*M
            
             10M |  5D   |  10I  |  10D  |  10M  |  5I   |  10D  |  10M
    Ref   :  0-9 | 10-14 | -14   | 15-24 | 25-34 | -34   | 35-44 | 45-54
    Tgt   :  0-9 | -9    | 10-19 | -19   | 20-29 | 30-34 | -34   | 35-44
    read0 :   == 
    read1 :                 ==
    read2 :               =======
    read3 :                  ================
    read4 :               =======================
    read5 :                                  ========================
    read6 :                               ===============================
    read7 : =============================================================
    read8 :                 ==***************==***********===========
    read9 :                 ==...............==...........============
    read10:   ==**********==***************====
    '''

        annot.LOG = 1
        refLens = {'1':55}        
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
        pool = [ 
                (2, 6, None, '5M', 2, 0, 0, 0),
                 (12, 16, None, '5I', -1, 0, 5, 0),
                 (10, 19, None, '10I', -1, 0, 10 ,0),
                 (13, 22, None, '7I,10D,3M', 25, 0, 7, 10),
                 (10, 29, None, '10I,10D,10M', 25, 0, 10, 10),
                 (23, 37, None, '7M,5I,10D,3M', 28, 0, 5, 10),
                 (20, 44, None, '10M,5I,10D,10M', 25, 0, 5, 10),
                 (0, 44, None, '10M,5D,10I,10D,10M,5I,10D,10M', 0, 0, 15, 25),
                 (13, 37, [(0,4), (2,5), (0, 6), (2, 7), (0,3)], '4I,12D,6M,12D,3M', 27, 0 , 4 ,0), 
                 (13, 37, [(0,4), (3,5), (0, 6), (3, 7), (0,3)], '4I,12N,6M,12N,3M', 27, 0 , 4 ,0), 
                 #(13, 37, [(0,4), (3,5), (0, 6), (3, 7), (0,3)], '4I,12N,6M,2N,10D,3M', 27, 0 , 4 ,0),
                 (2, 28, [(0,4),(2,4),(0,5),(2,5),(0,9)], '4M,9D,5I,10D,9M', 2, 0, 5, 0)
                ]
    
        self.batchTestHelper(modFile, pool, refLens)            


    def test2(self):
        '''Test case for M, MDM, MDI
              10M |  5D   |  10M  |  10D  | 10I
    Ref   :  0-9 | 10-14 | 15-24 | 25-34 | -34
    Tgt   :  0-9 | -9    | 10-19 | -19   | 20-29
    Read0 :   ==
    Read1 :   ================
    Read2 :  =====================
    Read3 :                  ================
    Read4 :               =======================             
    '''
                    
        refLens = {'1':35}
        modFile = StringIO.StringIO('''d\t1\t10\t1
d\t1\t11\t2
d\t1\t12\t3
d\t1\t13\t4
d\t1\t14\t5        
d\t1\t25\t1
d\t1\t26\t2
d\t1\t27\t3
d\t1\t28\t4
d\t1\t29\t5
d\t1\t30\t6
d\t1\t31\t7
d\t1\t32\t8
d\t1\t33\t9
d\t1\t34\t0
i\t1\t34\tabcdefghij
''')
                
        pool = [ (2, 6, None, '5M', 2, 0, 0, 0),
                 (3, 12, None, '7M,5D,3M', 3, 0, 0, 5),
                 (0, 19, None, '10M,5D,10M', 0, 0, 0, 5),
                 (13, 22, None, '7M,10D,3I', 18, 0, 3, 10),
                 (10, 29, None, '10M,10D,10I', 15, 0, 10, 10),
                ]
        
        self.batchTestHelper(modFile, pool, refLens)


    def test3(self):
        '''Test case for (M)I, (M)IDI, (D)IDI
             10M |  10I  |  10D  |  10I  |  5D   | 10I
    Ref   :  0-9 | -9    | 10-19 | -19   | 20-24 | -24
    Tgt   :  0-9 | 10-19 | -19   | 20-29 | -29   | 30-39
    Read1 :         ===
    Read2 :       =======
    Read3 :           ===============
    Read4 :       =======================
    Read5 :                           ================
    Read6 :                       =======================
    Read7 :    =======================================
    '''
        refLens = {'1':40}
        modFile = StringIO.StringIO('''i\t1\t9\tabcdefghij
d\t1\t10\t1
d\t1\t11\t2
d\t1\t12\t3
d\t1\t13\t4
d\t1\t14\t5
d\t1\t15\t6
d\t1\t16\t7
d\t1\t17\t8
d\t1\t18\t9
d\t1\t19\t0
i\t1\t19\tabcdefghij
d\t1\t20\t1
d\t1\t21\t2
d\t1\t22\t3
d\t1\t23\t4
d\t1\t24\t5
i\t1\t24\tabcdefghij
''')        
        pool = [(12, 16, None, '5I', -1, 0, 5, 0),
                (10, 19, None, '10I', -1, 0, 10, 0),
                (15, 24, None, '5I,10D,5I', -1, 0, 10, 10),
                (10, 29, None, '10I,10D,10I', -1, 0, 20, 10),
                (25, 34, None, '5I,5D,5I', -1, 0, 10, 5),
                (20, 39, None, '10I,5D,10I', -1, 0, 20, 5),
                (5, 34, None, '5M,10I,10D,10I,5D,5I', 5, 0, 25, 15)
                ]
        self.batchTestHelper(modFile, pool, refLens)



class TestAnnotator2(unittest.TestCase):    
    '''
    Test case for insertions/deletion/splicing junction in read
        
            10M |  10I  | 10M   | 5D    | 10M   |  5I  | 5D    | 5I   | 10M
    Ref   : 0-9 |  -9   | 10-19 | 20-24 | 25-34 | -34  | 35-39 | -39  | 40-49
    Tgt   : 0-9 | 10-19 | 20-29 | -29   | 30-39 | 40-44| -44   | 45-49| 50-59
    Read1 : =^=
    Read2 :    =^=
    Read3 :        =^=
    Read4 :            =^=
    Read5 :                    =^=
    Read6 :                                    =^=
    Read7 :                                       =^=
    Read8 :                                           =^=    
    Read9:                                                       =^=
    Read10:                                                           =^=    
    '''
            
    def setUp(self):
        annot.TESTING = 1
        annot.VERBOSITY = 1
        annot.LOG = 1
        
        self.refLens = {'1':50}
        self.modFile = StringIO.StringIO('''i\t1\t9\tabcdefghij
d\t1\t20\t1
d\t1\t21\t2
d\t1\t22\t3
d\t1\t23\t4
d\t1\t24\t5
i\t1\t34\tklmno
d\t1\t35\t6
d\t1\t36\t7
d\t1\t37\t8
d\t1\t38\t9
d\t1\t39\t0
i\t1\t39\tpqrst
''')
    
    
    def batchTestHelper(self, modFile, pool, refLens):                
        tmpName = tempfile.mkstemp('.tsv')[1]
        tmpfp = open(tmpName, 'wb')
        for line in modFile:
            tmpfp.write(line)
        tmpfp.close()
        pysam.tabix_index(tmpName, force=True, seq_col=1, start_col=2, end_col=2, 
                      meta_char='#', zerobased=True)
        tmpName += '.gz'
        modFile.close()
        
        self.chromoID = '1'
        self.modobj = mod.Mod(tmpName)
        self.modobj.load(self.chromoID)
        
        for tup in pool:       
            bamIter=[Read(tup[0], tup[1]+1, tup[2]) for tup in pool]        
                                   
        a = annot.Annotator(self.chromoID, refLens[self.chromoID],
                                self.modobj, bamIter)
        results = a.execute()
        
        for i,res in enumerate(results):            
            self.assertEqual(polish(res[0]),pool[i][3])
            self.assertEqual(res[1], pool[i][4])
            self.assertEqual(res[2], pool[i][5])
            self.assertEqual(res[3], pool[i][6])
            self.assertEqual(res[4], pool[i][7])
        
        os.remove(tmpName)
        os.remove(tmpName+'.tbi')
        
                            
    def test4(self):
        cigar = [(0,2),(1,1),(0,2)]  #MIM
        
        pool = [(2,5,cigar,'2M,1I,2M', 2, 0, 0, 0),
                (8,11,cigar,'2M,3I', 8, 0, 2, 0),
                (12,15,cigar,'5I', -1, 0, 4, 0),
                (18,21,cigar,'3I,2M', 10, 0, 2, 0),
                (28,31,cigar,'2M,1I,5D,2M', 18, 0, 0, 0), #########
                (38,41,cigar,'2M,3I', 33, 0, 2, 0),
                (41,44,cigar,'5I', -1, 0, 4, 0),
                (43,46,cigar,'3I,5D,2I', -1, 0, 4, 0),  ########
                (45,48,cigar,'5I', -1, 0, 4, 0),
                (48,51,cigar,'3I,2M', 40, 0, 2, 0),     
                ]
                
        self.batchTestHelper(self.modFile, pool, self.refLens)


    def test5(self):
        cigar = [(0,1),(2,1),(1,1),(2,1),(0,1)] #MDIDM
        pool = [
                (2,5,cigar,'1M,1D,1I,1D,1M', 2, 0, 0, 0),
                (8,11,cigar,'1M,1D,2I', 8, 0, 1, 0),
                (12,15,cigar,'3I', -1, 0, 2, 0),
                (18,21,cigar,'2I,1D,1M', 11, 0, 1, 0),
                (28,31,cigar,'1M,1D,1I,6D,1M', 18, 0, 0, 0), #########
                (38,41,cigar,'1M,1D,2I', 33, 0, 1, 0),
                (41,44,cigar,'3I', -1, 0, 2, 0),
                (43,46,cigar,'2I,5D,1I', -1, 0, 2, 0),  ########
                (45,48,cigar,'3I', -1, 0, 2, 0),
                (48,51,cigar,'2I,1D,1M', 41, 0, 1, 0),     
                ]
        self.batchTestHelper(self.modFile, pool, self.refLens)
    


#    def test5alt(self):
#        cigar = [(0,1),(2,1),(1,1),(2,1),(0,1)] #MDIDM
#        pool = [
#                (2,5,cigar,'1M,2D,1I,1M', 2, 0, 0, 0),
#                (8,11,cigar,'1M,1D,2I', 8, 0, 1, 0),
#                (12,15,cigar,'3I', -1, 0, 2, 0),
#                (18,21,cigar,'1I,1D,1I,1M', 11, 0, 1, 0),
#                (28,31,cigar,'1M,7D,1I,1M', 18, 0, 0, 0), #########
#                (38,41,cigar,'1M,1D,2I', 33, 0, 1, 0),
#                (41,44,cigar,'3I', -1, 0, 2, 0),
#                (43,46,cigar,'1I,5D,2I', -1, 0, 2, 0),  ########
#                (45,48,cigar,'3I', -1, 0, 2, 0),
#                (48,51,cigar,'1I,1D,1I,1M', 41, 0, 1, 0),     
#                ]
#        self.batchTestHelper(self.modFile, pool, self.refLens)

        
    def test6(self):
        cigar = [(2,2),(1,1),(2,1),(0,1)]
        pool = [(2,5,cigar,'2D,1I,1D,1M', 5, 0, 0, 0),
                (8,11,cigar,'2D,2I', -1, 0, 1, 0),
                (12,15,cigar,'2I', -1, 0, 1, 0),
                (18,21,cigar,'1I,1D,1M', 11, 0, 0, 0),
                (28,31,cigar,'2D,1I,6D,1M', 26, 0, 0, 0), #########
                (38,41,cigar,'2D,2I', -1, 0, 1, 0),
                (41,44,cigar,'2I', -1, 0, 1, 0),
                (43,46,cigar,'1I,5D,1I', -1, 0, 1, 0),  ########
                (45,48,cigar,'2I', -1, 0, 1, 0),
                (48,51,cigar,'1I,1D,1M', 41, 0, 0, 0),     
                ]
        self.batchTestHelper(self.modFile, pool, self.refLens)


        
#        cigar = [(0,1),(2,1),(1,1),(2,2)]
#        pool = [(2,5,cigar,'2M,1I,2M', 2, 0, 0, 0),
#                (8,11,cigar,'2M,3I', 8, 0, 2, 0),
#                (12,15,cigar,'5I', -1, 0, 4, 0),
#                (18,21,cigar,'3I,2M', 10, 0, 2, 0),
#                (28,31,cigar,'2M,1I,5D,2M', 18, 0, 0, 0), #########
#                (38,41,cigar,'2M,3I', 33, 0, 2, 0),
#                (41,44,cigar,'5I', -1, 0, 4, 0),
#                (43,46,cigar,'3I,5D,2I', -1, 0, 4, 0),  ########
#                (45,48,cigar,'5I', -1, 0, 4, 0),
#                (48,51,cigar,'3I,2M', 40, 0, 2, 0),     
#                ]
#        
#        #cigar = [(1,1),(0,4)]
#        pool = [(2,5,cigar,'2M,1I,2M', 2, 0, 0, 0),
#                (8,11,cigar,'2M,3I', 8, 0, 2, 0),
#                (12,15,cigar,'5I', -1, 0, 4, 0),
#                (18,21,cigar,'3I,2M', 10, 0, 2, 0),
#                (28,31,cigar,'2M,1I,5D,2M', 18, 0, 0, 0), #########
#                (38,41,cigar,'2M,3I', 33, 0, 2, 0),
#                (41,44,cigar,'5I', -1, 0, 4, 0),
#                (43,46,cigar,'3I,5D,2I', -1, 0, 4, 0),  ########
#                (45,48,cigar,'5I', -1, 0, 4, 0),
#                (48,51,cigar,'3I,2M', 40, 0, 2, 0),     
#                ]
#        
#        #cigar = [(0,4),(1,1)]
#        pool = [(2,5,cigar,'2M,1I,2M', 2, 0, 0, 0),
#                (8,11,cigar,'2M,3I', 8, 0, 2, 0),
#                (12,15,cigar,'5I', -1, 0, 4, 0),
#                (18,21,cigar,'3I,2M', 10, 0, 2, 0),
#                (28,31,cigar,'2M,1I,5D,2M', 18, 0, 0, 0), #########
#                (38,41,cigar,'2M,3I', 33, 0, 2, 0),
#                (41,44,cigar,'5I', -1, 0, 4, 0),
#                (43,46,cigar,'3I,5D,2I', -1, 0, 4, 0),  ########
#                (45,48,cigar,'5I', -1, 0, 4, 0),
#                (48,51,cigar,'3I,2M', 40, 0, 2, 0),     
#                ]
#        
#        self.batchTestHelper(modFile, pool, refLens)

        
if __name__ == '__main__':
    unittest.main()
