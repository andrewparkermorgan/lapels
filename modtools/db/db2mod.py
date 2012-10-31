#! /bin/env python 
'''
Created on Oct 28, 2012

@author: Shunping Huang
'''

import csv
import gc
import argparse as ap
import gzip
import os 
from time import localtime,strftime

import dbutils as dbu
from modtools.variants import parseVariant, INS, DEL
from modtools.utils import buildChromMap,getOutChrom,log
from modtools import version

DESC = 'A DB to MOD converter.'
__version__ = '0.0.2'
VERBOSITY = 1

#dbu.db = "/playpen/data/newgenes.db"
#outfile = "./test.mod"
#ref = 'mm9'
#sample = 'A_J'
#chroms = []


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
 

if __name__ == '__main__':
    p = ap.ArgumentParser(description=DESC, 
                          formatter_class = ap.RawTextHelpFormatter)
    # Optional arguments
    group = p.add_mutually_exclusive_group()    
    group.add_argument("-q", dest='quiet', action='store_true',
                        help='quiet mode')
    group.add_argument('-v', dest='verbosity', action="store_const", const=2,
                        default=1, help='verbose mode')                    
    p.add_argument('-c', metavar='chromList', dest='chroms', 
                   type=validChromList, default = [],
                   help='a comma-separated list of chromosomes in VCF' +
                        ' (default: all)')
    p.add_argument('--map', metavar='chromMap', dest='cmap', 
                   type=readableFile, default = None,
                   help='the file of chromosome name mapping from VCF to MOD' +
                        ' (default: none)')
    p.add_argument('-o', metavar='mod', dest='mod', 
                   type=writableFile, default=None, 
                   help='the output mod file'\
                        +' (default: <sample>.mod)')
    # Required arguments
    p.add_argument('ref', help='reference name')
    p.add_argument('sample', help='requested sample name in VCF')
    p.add_argument('db', type=readableFile, help='database location')                    
        
    args = p.parse_args()
    
    if args.quiet:
        VERBOSITY = 0
    else:            
        VERBOSITY = args.verbosity
        
    if args.mod is None:                                
        modfp = gzip.open(args.sample + '.mod', 'wb')
    else:
        modfp = gzip.open(args.mod, 'wb')        
    
    if args.cmap is None:        
        chromMap = {}
    else:
        chromMap = buildChromMap(args.cmap)
                        
    chroms = args.chroms        
    sample = args.sample
    ref = args.ref    
    db = args.db    
    
    if VERBOSITY > 0:
        log("from %s to %s\n" %(ref, sample), 1 ,True)
        log("input DB file: %s\n" % db, 1, True)                
        log("output MOD file: %s\n" % modfp.name, 1, True)
                                
    if len(chroms) == 0:                    
        chroms = [str(i) for i in range(1,20)] + ['X','Y','M']
    
    modfp.write("#version=%s\n" % version.__mod_version__)
    modfp.write("#date=%s\n" % strftime("%Y%m%d",localtime()))
    modfp.write("#reference=%s\n" % ref)
    modfp.write("#sample=%s\n" % sample)    
    out = csv.writer(modfp, delimiter='\t',lineterminator='\n')
    
    for chrom in chroms:  # for each chromosome
        gc.disable()
        nSub = 0            
        nIns = 0
        nDel = 0
        pool = []
        modChrom = getOutChrom(chromMap, chrom)
        
        if VERBOSITY > 0:
            log("processing chromosome '%s' in db\n" % 
                chrom, 1, True)  
        
        # Read SNPs        
        snps = dbu.readSNPsFromDB(db, dbu.chromMap[chrom], 
                                  dbu.strainMap[sample])        
        nSub = len(snps)
        for snp in snps: 
            pool.append(('s', modChrom, snp[0], "%s/%s" % (snp[1],snp[2])))
        
        if VERBOSITY > 0:    
            log("%d SNP(s) read from DB\n" % nSub, 1, True)
        del snps

        ## Read indels    
        indels = dbu.readIndelsFromDB(db, dbu.chromMap[chrom], 
                                      dbu.strainMap[sample])        
        nIndels = len(indels)
        for indel in indels:
            v=parseVariant(modChrom, indel[0], indel[1], indel[2])            
            if v.type == INS:            
                pool.append(('i', modChrom, v.start[1], v.extra))
                nIns += v.length                                           
            elif v.type == DEL:
                # Change non-atomic deletions to atomic
                for j in range(v.length): 
                    pool.append(('d', modChrom, v.start[1]+j, v.extra[j]))
                nDel += v.length                                        
            else:
                raise ValueError("Unknown variant type: %s" % v.type)
        
        if VERBOSITY > 0:
            log("%d indel(s) read from DB\n" % nIndels, 1, True)
        del indels                    
        
        pool=sorted(set(pool), key = lambda tup: tup[2])            
        out.writerows(pool)
        if VERBOSITY > 0:            
            log("%d line(s) written to MOD\n" % len(pool), 1, True)
            log("SNPs: %d base(s)\n" % nSub, 1, True)
            log("Insertions: %d base(s)\n" % nIns, 1, True)
            log("Deletions: %d base(s)\n" % nDel, 1, True)
            
        del pool
        gc.enable()                        
                        
    modfp.close()
    
    if VERBOSITY > 0:
        log("All Done!\n", 1, True)