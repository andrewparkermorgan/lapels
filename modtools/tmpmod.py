'''
Convert a mod format to a tabix accessible format.

Created on Oct 30, 2012

@author: Shunping Huang
'''

import logging
import gzip
import tempfile
import pysam


__all__ = ['getTabixMod']


def getTabixMod(filename):
    '''Unzip a mod file, use bgzip to rezip it, and and build tabix index.'''
    logger = logging.getLogger('tmpmod') 
    logger.info('extracting MOD file ...')   
    modfp = gzip.open(filename, 'rb')
    tmpName = tempfile.mkstemp('.tsv')[1]    
    tmpfp = open(tmpName, 'wb')
    tmpfp.writelines(modfp)
    tmpfp.close()
    modfp.close()    
    pysam.tabix_index(tmpName, force=True, seq_col=1, start_col=2, end_col=2, 
                      meta_char='#', zerobased=True)
    tmpName += '.gz'
    logger.info('temporary file %s created', tmpName)
    return tmpName

#print(getTabixMod("../data/B.mod"))