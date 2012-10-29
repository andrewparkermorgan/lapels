#! /bin/env python 
'''
Created on Oct 28, 2012

@author: Shunping Huang
'''

import pysam
import csv
import gc
import sys
import argparse
from time import localtime, strftime

import dbutils as dbu
from modtools.variants import parseVariant, INS, DEL

DESC = 'A DB to MOD converter.'
__version__ = '0.0.1'
VERBOSITY = 1

chromMap = dict([(str(i), 'chr'+str(i)) for i in range(20)] + 
                [('X', 'chrX'), ('Y', 'chrY'), 
                 ('M', 'chrM'), ('MT',' chrM')])
                

#dbu.db = "/playpen/data/newgenes.db"
#outfile = "./test.mod"
#ref = 'mm9'
#sample = 'A_J'
#chroms = []


def csvChromList(s):
    delchars = ''.join(c for c in map(chr, range(256)) if not c.isdigit())    
    chroms = set()
    for tmp in s.split(','):
        chrom = tmp[:-1].translate(None, delchars)+tmp[-1]
        chroms.add(chrom)            
    return chroms


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
                   help='a comma-separated list of chromosomes' +
                        ' (default: all)')
    p.add_argument('-o', metavar='mod', dest='mod', type=argparse.FileType('w'), 
                   default=None, help='the output mod file'\
                        +' (default: sample.mod)')
    
    p.add_argument('ref', help='reference name')
    p.add_argument('sample', help='requested sample name in DB')
    p.add_argument('db', type=argparse.FileType('r'), help='database location')                    
        
    args = p.parse_args()
    
    
    if args.mod is None:                
        outfile = args.sample + '.mod'            
        mod = open(outfile, 'wb')
    else:
        outfile = args.mod.name
        mod = args.mod
        
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
                        
    out = csv.writer(mod, delimiter='\t',lineterminator='\n')
    
    if len(chroms) == 0:                    
        chroms = [str(i) for i in range(1,20)] + ['X','Y','M']
    
    mod.write("#ver=0.1\n")
    mod.write("#ref=%s\n" % ref)
    mod.write("#sample=%s\n" % sample)
    for chrom in chroms:  # for each chromosome
        gc.disable()
        pool = []
        
        if VERBOSITY > 0:
            log("processing chromosome %s in db\n" % 
                chrom, 1, True)  
        
        # Read SNPs        
        snps = dbu.readSNPsFromDB(db, dbu.chromMap[chrom], 
                                  dbu.strainMap[sample])        
        nSNPs = len(snps)
        for snp in snps: 
            try:
                newChrom = chromMap[chrom]
            except KeyError:
                newChrom = chrom               
            pool.append(('s', newChrom, snp[0], "%s/%s" % (snp[1],snp[2])))
        
        if VERBOSITY > 0:    
            log("%d SNP(s) read from DB\n" % nSNPs, 1, True)
        del snps

        ## Read indels    
        indels = dbu.readIndelsFromDB(db, dbu.chromMap[chrom], 
                                      dbu.strainMap[sample])        
        nIndels = len(indels)
        for indel in indels:
            try:
                newChrom = chromMap[chrom]
            except KeyError:
                newChrom = chrom            
            v=parseVariant(newChrom, indel[0], indel[1], indel[2])            
            if v.type == INS:            
                pool.append(('i', newChrom, v.start[1], v.extra))                                           
            elif v.type == DEL:
                # Change non-atomic deletions to atomic
                for i in range(v.length): 
                    pool.append(('d', newChrom, v.start[1]+i, v.extra[i]))                                        
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
                    
    
    mod.close()
    
    if VERBOSITY > 0:
        log("building tabix index for MOD\n", 1, True)
    pysam.tabix_index(outfile, force=True, seq_col=1, start_col=2, end_col=2, 
                      meta_char='#', zerobased=True)
    
    if VERBOSITY > 0:
        log("All Done!\n", 1, True)