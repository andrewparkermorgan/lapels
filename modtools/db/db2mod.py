#! /bin/env python 
'''
Created on Oct 28, 2012

@author: Shunping Huang
'''

import pysam
import csv
import gc
import argparse

import dbutils as dbu
from modtools.variants import parseVariant, INS, DEL
from modtools.utils import buildChromMap,getOutChrom,log

DESC = 'A DB to MOD converter.'
__version__ = '0.0.2'
VERBOSITY = 1

#dbu.db = "/playpen/data/newgenes.db"
#outfile = "./test.mod"
#ref = 'mm9'
#sample = 'A_J'
#chroms = []


def csvChromList(s):
    #delchars = ''.join(c for c in map(chr, range(256)) if not c.isdigit())    
    chroms = set()
    for chrom in s.split(','):
        #chrom = tmp[:-1].translate(None, delchars)+tmp[-1]
        chroms.add(chrom)            
    return chroms


if __name__ == '__main__':    
    ## Parse arguments
    p = argparse.ArgumentParser(description=DESC, 
                                formatter_class = argparse.RawTextHelpFormatter)
    group = p.add_mutually_exclusive_group()    
    group.add_argument("-q", dest='quiet', action='store_true',
                       help='quiet mode')
    group.add_argument('-v', dest='verbosity', action="store_const", const=2,
                       default=1, help='verbose mode')        
    p.add_argument('-c', metavar='chromList', dest='chroms', type=csvChromList,
                   default = set(),
                   help='a comma-separated list of chromosomes in DB' +
                        ' (default: all)')
    p.add_argument('-o', metavar='mod', dest='modfp', type=argparse.FileType('w'), 
                   default=None, help='the output mod file'\
                        +' (default: sample.mod)')    
    p.add_argument('--map', metavar='chromMap', dest='mapfp', type=argparse.FileType('r'),
                   default = None,
                   help='the file of chromosome name mapping from DB to MOD' +
                        ' (default: none)')
    p.add_argument('ref', help='reference name')
    p.add_argument('sample', help='requested sample name in DB')
    p.add_argument('db', type=argparse.FileType('r'), help='database location')                    
        
    args = p.parse_args()
    
    
    if args.modfp is None:                
        outfile = args.sample + '.mod'            
        modfp = open(outfile, 'wb')
    else:
        outfile = args.modfp.name
        modfp = args.modfp
        
    chromMap = buildChromMap(args.mapfp)
             
    chroms = list(args.chroms)
    sample = args.sample
    ref = args.ref
    
    if args.quiet:
        VERBOSITY = 0
    else:            
        VERBOSITY = args.verbosity
    
    args.db.close()
    db = args.db.name    
    
    if VERBOSITY > 0:
        log("from %s to %s\n" %(ref, sample), 1 ,True)
        log("input DB file: %s\n" % db, 1, True)                
        log("output MOD file: %s\n" % outfile, 1, True)
                        
    out = csv.writer(modfp, delimiter='\t',lineterminator='\n')
    
    if len(chroms) == 0:                    
        chroms = [str(i) for i in range(1,20)] + ['X','Y','M']
    
    modfp.write("#ver=0.1\n")
    modfp.write("#ref=%s\n" % ref)
    modfp.write("#sample=%s\n" % sample)
    for chrom in chroms:  # for each chromosome
        gc.disable()
        pool = []
        modChrom = getOutChrom(chromMap, chrom)
        
        if VERBOSITY > 0:
            log("processing chromosome %s in db\n" % 
                chrom, 1, True)  
        
        # Read SNPs        
        snps = dbu.readSNPsFromDB(db, dbu.chromMap[chrom], 
                                  dbu.strainMap[sample])        
        nSNPs = len(snps)
        for snp in snps: 
            pool.append(('s', modChrom, snp[0], "%s/%s" % (snp[1],snp[2])))
        
        if VERBOSITY > 0:    
            log("%d SNP(s) read from DB\n" % nSNPs, 1, True)
        del snps

        ## Read indels    
        indels = dbu.readIndelsFromDB(db, dbu.chromMap[chrom], 
                                      dbu.strainMap[sample])        
        nIndels = len(indels)
        for indel in indels:                  
            v=parseVariant(modChrom, indel[0], indel[1], indel[2])            
            if v.type == INS:            
                pool.append(('i', modChrom, v.start[1], v.extra))                                           
            elif v.type == DEL:
                # Change non-atomic deletions to atomic
                for j in range(v.length): 
                    pool.append(('d', modChrom, v.start[1]+j, v.extra[j]))                                        
            else:
                raise ValueError("Unknown variant type: %s" % v.type)
        
        if VERBOSITY > 0:
            log("%d indel(s) read from DB\n" % nIndels, 1, True)
        del indels                    
        
        pool=sorted(set(pool), key = lambda tup: tup[2])            
        out.writerows(pool)
        if VERBOSITY > 0:
            log("%d unique atomic variants written to MOD\n" % len(pool), 1, True)
        del pool
        gc.enable()                        
                    
    
    modfp.close()
    
    if VERBOSITY > 0:
        log("building tabix index for MOD\n", 1, True)
    pysam.tabix_index(outfile, force=True, seq_col=1, start_col=2, end_col=2, 
                      meta_char='#', zerobased=True)
    
    if VERBOSITY > 0:
        log("All Done!\n", 1, True)