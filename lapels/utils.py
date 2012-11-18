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
            
                