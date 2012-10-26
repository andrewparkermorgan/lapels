'''
The module of reading blocks

Created on Sep 1, 2012

@author: Shunping Huang
'''


import collections

class BlockReader(collections.Iterator):
    def __init__(self, fp, isBlockDelimiter=lambda x: x.startswith('>'), 
                 isSkipped=lambda x: x.startswith('##')):
        #self.fp = open(fileName,'rb')
        self.fp = fp.__iter__()
        self.isBlockDelimiter = isBlockDelimiter
        self.isSkipped = isSkipped        
        self.buffer = None
        self.isBlockOpen = False
        self.isBlockEnd = False
        self.isFileEnd = False


    def __iter__(self):
        return self


    def next(self):        
        if self.isBlockEnd:                 ## The current block is over.            
            raise StopIteration()
        else:                        
            if self.buffer is None:         ## Read from file
                try:
                    line = self.fp.next()            
                    if self.isSkipped is not None:
                        while self.isSkipped(line):
                            line = self.fp.next()
                except StopIteration:
                    line = None
            else:    
                line = self.buffer          ## Read from buffer
                self.buffer = None                                            

            if line is None:                ## The end of file
                self.isFileEnd = True
                self.isBlockEnd = True                
                raise StopIteration()
            
            #print(line)
            if self.isBlockDelimiter(line): ## Block delimiter line
                if self.isBlockOpen:        ## If an existing block is open.
                    self.isBlockEnd = True                    
                    self.buffer = line                    
                    raise StopIteration()
                else:                       ## The block is not open
                    self.isBlockOpen = True
                    return line
            else:                           ## Common data line
                self.isBlockOpen = True
                return line



    def peek(self):
        if self.buffer is None:
            try:
                line = self.fp.next()            
                if self.isSkipped is not None:
                    while self.isSkipped(line):
                        line = self.fp.next()                
            except StopIteration:
                line = None
            self.buffer = line
        return self.buffer
    
    

    def nextBlock(self):
        if not self.isBlockEnd:
            while not self.isBlockEnd:
                try:
                    self.next()
                except StopIteration:
                    pass

        if self.isFileEnd:
            return False
        else:
            self.isBlockOpen = False
            self.isBlockEnd = False
            return True


    def readline(self):
        try:
            return self.next()
        except StopIteration:
            return ""


