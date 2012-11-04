'''
Miscellaneous utilities

Created on Oct 29, 2012

@author: Shunping Huang
'''


import os
import sys
import argparse as ap
from time import localtime, strftime

__all__ = ['log', 'validChromList', 'readableFile', 'writableFile'] 


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
        

def validChromList(s):
    chromList = []
    chromSet = set()                    
    for chrom in s.split(','):
        if chrom in chromSet:
            raise ap.ArgumentTypeError("Duplicated chromosome '%s' found."
                                       % chrom)
        else:
            chromSet.add(chrom)
            chromList.append(chrom)            
    return chromList


def readableFile(fileName):
    if os.path.isfile(fileName) and os.access(fileName, os.R_OK):
        return fileName
    else:
        raise ap.ArgumentTypeError("Cannot read file '%s'." % fileName)


def writableFile(fileName):
    if os.access(os.path.dirname(fileName), os.W_OK):
        return fileName
    else:        
        raise ap.ArgumentTypeError("Cannot write file '%s'." % fileName)

            
#def buildChromMap(fileName):
#    assert fileName is not None
#    fp = open(fileName, 'r')
#    chromMap = {}
#    if fp is not None:
#        for line in fp:
#            tmp = line.rstrip().split('\t')
#            if tmp[0] in chromMap.keys():
#                raise ValueError("Duplicated keys in chromosome mappings.")
#            else:
#                try:
#                    chromMap[tmp[0]] = tmp[1]
#                except IndexError:
#                    raise ValueError('Mapping format not correct.')
#    fp.close()
#    return chromMap
#
#
#def getOutChrom(chromMap, inChrom, prefix=None):
#    assert chromMap is not None
#    if inChrom in chromMap.keys():
#        return chromMap[inChrom]
#    elif prefix is not None:
#        return prefix + inChrom
#    else:
#        return inChrom
                