#! /bin/env python 
'''
Created on Oct 28, 2012

@author: Shunping Huang
'''

import csv
import gc
import argparse as ap
import gzip
import logging
from time import localtime,strftime

import dbutils as dbu
from modtools.variants import parseVariant, INS, DEL
from modtools.utils import readableFile, writableFile, validChromList
from modtools import version
from modtools import metadata

DESC = 'A DB to MOD converter.'
__version__ = '0.1.0'
VERBOSITY = 1
logger = None


def initLogger():
    global logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    formatter = logging.Formatter("[%(asctime)s] %(name)-10s: %(levelname)s: %(message)s",
                                  "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    logger.addHandler(ch)

 
if __name__ == '__main__':
    initLogger()
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
                   help='a comma-separated list of chromosomes in output' +
                        ' (default: all)')  
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
        logger.setLevel(logging.CRITICAL)
    elif args.verbosity == 2:                    
        logger.setLevel(logging.DEBUG)
            
    if args.mod is None:                                
        modfp = gzip.open(args.sample + '.mod', 'wb')
    else:
        modfp = gzip.open(args.mod, 'wb')        
                        
    chroms = args.chroms        
    sample = args.sample
    ref = args.ref    
    db = args.db
    
    meta = metadata.MetaData(ref)
    chromAliases = meta.chromAliases
    
    logger.info("from %s to %s", ref, sample)
    logger.info("input DB file: %s", db)                
    logger.info("output MOD file: %s", modfp.name)
                                
    if len(chroms) == 0:                    
        chroms = [str(i) for i in range(1,20)] + ['X','Y','M']
    
    modfp.write("#version=%s\n" % version.__mod_version__)
    modfp.write("#date=%s\n" % strftime("%Y%m%d",localtime()))
    modfp.write("#reference=%s\n" % ref)
    modfp.write("#sample=%s\n" % sample)    
    out = csv.writer(modfp, delimiter='\t',lineterminator='\n')
    
    for modChrom in chroms:  # for each chromosome
        gc.disable()
        nSub = 0            
        nIns = 0
        nDel = 0
        pool = []
        
        chrom = chromAliases.getBasicName(modChrom)
        
        
        logger.info("processing chromosome '%s' in db", chrom)  
        
        # Read SNPs        
        snps = dbu.readSNPsFromDB(db, dbu.chromMap[chrom], 
                                  dbu.strainMap[sample])        
        nSub = len(snps)
        for snp in snps: 
            pool.append(('s', modChrom, snp[0], "%s/%s" % (snp[1],snp[2])))
        
            
        logger.info("%d SNP(s) read from DB", nSub)
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
        
        
        logger.info("%d indel(s) read from DB", nIndels)
        del indels                    
        
        pool=sorted(set(pool), key = lambda tup: tup[2])                                    
        out.writerows(pool)
        
        logger.info("%d line(s) written to MOD", len(pool))
        if len(pool) > 0:
            logger.info("SNPs: %d base(s)", nSub)
            logger.info("Insertions: %d base(s)", nIns)
            logger.info("Deletions: %d base(s)", nDel)
                    
        del pool
        gc.enable()                        
                        
    modfp.close()
    
    logger.info("All Done!")