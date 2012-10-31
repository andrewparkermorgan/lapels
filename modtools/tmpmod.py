'''
Created on Oct 30, 2012

@author: Shunping Huang
'''

import tempfile
import gzip
import pysam


def getTabixMod(filename):
    '''Unzip a mod file, use bgzip to rezip it, and and build tabix index.'''
    modfp = gzip.open(filename, 'rb')
    fd, tmpName = tempfile.mkstemp('.tsv')
    tmpfp = open(tmpName, 'wb')
    tmpfp.writelines(modfp)
    tmpfp.close()
    modfp.close()    
    pysam.tabix_index(tmpName, force=True, seq_col=1, start_col=2, end_col=2, 
                      meta_char='#', zerobased=True)
    return tmpName+'.gz'

#print(getTabixMod("../data/B.mod"))