'''
The module of building cigar from regions

Created on Oct 8, 2012

@author: Shunping Huang
'''

import cigarutils
import regionutils

class CigarBuilder():
    def __init__(self):
        self.pend = -1
        self.cigar = []
            
        
    def append(self, region):
        if region[0] != 1:
            if self.pend >= 0:            
                delta = region[2] - 1 - self.pend
                if delta < 0:
                    print region
                    print self.pend                        
                assert delta >= 0
                if delta > 0:
                    self.cigar.append((2, delta))  ##Insert deletions to gaps
        if region[3] >= 0:
            self.pend = region[3]
        self.cigar.extend(region[1])
        
    
            
    def build(self, regions):                
        maxIdx = len(regions)
        idx1 = 0  ## Index of the current non-I_1 region
        idx2 = 0  ## Index of the current I_1 region
        while idx1 < maxIdx:
            if regions[idx1][0] != 1:
                break
            idx1 += 1
        while idx2 < maxIdx:
            if regions[idx2][0] == 1:
                break
            idx2 += 1
                  
        while idx1 < maxIdx and idx2 < maxIdx:
            if regions[idx1][0] == 0:   ## M_1                
                ## append I_1 before this M_1
                while idx1 > idx2:
                    assert not regionutils.isRightTo(regions[idx2],regions[idx1])                                        
                    self.append(regions[idx2])
                    ## next idx2
                    idx2 += 1                    
                    while idx2 < maxIdx:
                        if regions[idx2][0] == 1:
                            break
                        idx2 += 1

                assert idx1 < idx2                
                ## append M_1
                self.append(regions[idx1])
                idx1 += 1
                ## next idx2                
                while idx1 < maxIdx:
                    if regions[idx1][0] != 1:
                        break
                    idx1 += 1               
            else:                       ## D_1/N_1
                ## append the I_1 left on the left
                while idx2 < maxIdx:                    
                    if regionutils.isLeftTo(regions[idx2], regions[idx1]):
                        self.append(regions[idx2])
                        ## next idx2
                        idx2 += 1
                        while idx2 < maxIdx:
                            if regions[idx2][0] == 1:
                                break
                            idx2 += 1
                    else:
                        break                                
                ## get the overlapping I_1           
                overlap = []
                while idx2 < maxIdx:                    
                    if regionutils.isRightTo(regions[idx2], regions[idx1]):
                        break
                    else:
                        overlap.append(regions[idx2])
                        idx2 += 1
                        while idx2 < maxIdx:
                            if regions[idx2][0] == 1:
                                break
                            idx2 += 1
#                print(regions[idx1])
#                print(overlap)
                ## replace the current one with overlapping I_1
                replaced = cigarutils.replace(regions[idx1][1], 
                                              regions[idx1][2], 
                                              overlap)
#                print "replaced",replaced                
                self.append((regions[idx1][0], replaced, 
                             regions[idx1][2], regions[idx1][3]))                                
                ## next idx1
                idx1 += 1
                while idx1 < maxIdx:
                    if regions[idx1][0] != 1:
                        break
                    idx1 += 1
                        
        for i in range(idx1, maxIdx):
            if regions[i][0] != 1:
                self.append(regions[i])    
        for i in range(idx2, maxIdx):
            if regions[i][0] == 1:
                self.append(regions[i])
                
        return self.cigar