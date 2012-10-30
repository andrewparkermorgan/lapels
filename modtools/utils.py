'''
Created on Oct 29, 2012

@author: Shunping Huang
'''

#chromMap = dict([(str(i), 'chr'+str(i)) for i in range(1,20)] + 
#                [('X', 'chrX'), ('Y', 'chrY'), 
#                 ('M', 'chrM'), ('MT',' chrM')])

import sys
from time import localtime, strftime

def log(s, verbosity=2, showtime=False):     
    if verbosity == 1:
        if showtime:
            msg = "[%s] %s" % (strftime("%Y/%m/%d %H:%M:%S", localtime()), s)
        else:
            msg = s
        sys.stdout.write(msg)
        sys.stdout.flush()
    elif verbosity > 1:
        if showtime:
            msg = "[%s] %s" % (strftime("%Y/%m/%d %H:%M:%S", localtime()), s)
        else:
            msg = s
        sys.stderr.write(msg)
        
        
def buildChromMap(fp):
    chromMap = {}
    if fp is not None:
        for line in fp:
            tmp = line.rstrip().split('\t')
            if tmp[0] in chromMap.keys():
                raise ValueError("Duplicated keys in chromosome mappings.")
            else:
                try:
                    chromMap[tmp[0]] = tmp[1]
                except IndexError:
                    raise ValueError('Mapping format not correct.')
    return chromMap


def getOutChrom(chromMap, inChrom, prefix=None):
    assert chromMap is not None
    if inChrom in chromMap.keys():
        return chromMap[inChrom]
    elif prefix is not None:
        return prefix + inChrom
    else:
        return inChrom


                