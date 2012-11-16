'''
Created on Nov 13, 2012

@author: Shunping Huang
'''

__all__ = ['fixmate']

import gc
import pysam
import logging
logger = logging.getLogger() 

MAX_SEGMENTS_PER_HIT = 2


def process(reads, rname):
    '''Process a set of reads with the same read name'''
    nReads = len(reads)    
    assert nReads > 0
         
    tags = [dict(r.tags) for r in reads]        
    hits = dict()
    for i in range(nReads):       
        hi = int(tags[i].get('HI',0))
        if hi in hits.keys():
            hits[hi].append(i)
        else:
            hits[hi] = [i]
    
#    if logger.isEnabledFor(logging.DEBUG):
#        logger.debug(reads[0].qname)
#        for k,v in hits.items():
#            logger.debug('%s -> %s' %(k, str(v)))
    
    try:
        NH = len(hits)
#        if logger.isEnabledFor(logging.DEBUG):
#            logger.debug("before:")
#            for i in range(len(reads)):
#                logger.debug(str(reads[i]))
#                logger.debug(reads[i].tlen)
                                        
        # Update NH, HI, mate chrom and position, and insertion size                
        for HI,k in enumerate(sorted(hits.keys())):
            hit = hits[k]
            n = len(hit) # Number of segments in a hit
            assert n > 0
            
            if n > MAX_SEGMENTS_PER_HIT:
                logger.warning("%d reads in a hit" % n)
                raise ValueError("%d reads in a hit" % n)
            
            for j in hit:
                tag = tags[j]            
                if NH != tag.get('NH', None):
#                    logger.debug("replacing NH=%d with %d" % (tag['NH'], NH))
                    tag['NH'] = NH
                if HI != tag.get('HI', 0):
#                    logger.debug("replacing NH=%d with %d" % (tag['NH'], NH))
                    if NH > 1:
                        tag['HI'] = HI
                    else:
                        del tag['HI']
                                                            
            # Single-end reads
            if n == 1:
                reads[hit[0]].mrnm = -1
                reads[hit[0]].mpos = -1
                reads[hit[0]].tlen = 0
            # Multi-end reads
            else:
                for j in range(n):
                    cur = reads[hit[j]]
                    nxt = reads[hit[(j+1)%n]]                            
                    cur.mrnm = nxt.tid
                    cur.mpos = nxt.pos
                            
                    if cur.tid == nxt.tid:                    
                        if cur.pos < nxt.pos:
                            cur.tlen = nxt.aend - cur.pos 
                        else:
                            cur.tlen = -(cur.aend - nxt.pos)
                    else:
                        cur.tlen = 0                        
        
        # Update CC and CP
        prev = dict()
        for HI,k in enumerate(sorted(hits.keys())):
            if HI > 0 and HI < NH:
                for j in hits[k]:
                    idx = prev[reads[j].is_read1]                    
                    try:
                        if reads[idx].tid == reads[j].tid:
                            tags[idx]['CC'] = '='
                        else:
                            tags[idx]['CC'] = rname(reads[j].tid)
                    except KeyError:
                        raise ValueError('No previous read1 or read2 is found')
                    tags[idx]['CP'] = reads[j].pos + 1 # 1-based CP
                prev=dict()
                    
            for j in hits[k]:      
                is_read1 = reads[j].is_read1
                if is_read1 in prev:                                    
                    raise ValueError("More than one read1 or read2 is found.")
                else:
                    prev[is_read1] = j
                    
            if HI == NH - 1:
                for j in hits[k]:
                    if 'CC' in tags[j].keys():
                        del tags[j]['CC']
                    if 'CP' in tags[j].keys():
                        del tags[j]['CP']
        
        # Update tags
        for i in range(nReads):
            reads[i].tags = [(k,v) for k,v in tags[i].items()]
            
    except ValueError, err:
        logger.warning('%s : %s', str(err), reads[0].qname)        
        return 0
        
#    if logger.isEnabledFor(logging.DEBUG):
#        logger.debug("after:")
#        for i in range(len(reads)):
#            logger.debug(str(reads[i]))
#            logger.debug(reads[i].tlen)
            
    return nReads


def fixmate(infile, outfile):
    inbam = pysam.Samfile(infile, 'rb')
    outbam = pysam.Samfile(outfile, 'wb', header=inbam.header, 
                           referencenames=inbam.references)
    qname = None    
    nTotal = 0
    nFixed = 0
    count = 0;        
    reads = []    
    gc.disable()    
    for rseq in inbam.fetch(until_eof=True):
        nTotal += 1        
        if qname is None or qname == rseq.qname:
            qname = rseq.qname
            reads.append(rseq)            
        else:            
            count = process(reads, inbam.getrname) 
            if count > 0:
                for r in reads:
                    outbam.write(r)
                nFixed += count                
            qname = rseq.qname
            del reads
            reads = [rseq]
        if nTotal % 200000 == 0:        
            logger.info('%d read(s) fixed' % nTotal)
            gc.enable()
            gc.disable()
    
    count = process(reads, inbam.getrname)
    if count > 0:
        for r in reads:
            outbam.write(r)
        nFixed += count
    logger.info('%d read(s) processed' % nTotal)
    logger.info('%d read(s) fixed and written to file' % nFixed)
        
    inbam.close()
    outbam.close()
    return (nTotal, nFixed)