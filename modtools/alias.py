'''
A module for handling with aliases of a name.
  For example, the sequence name in Mitochondria can be M, MT, chrM, and etc. 

Created on Nov 4, 2012

@author: Shunping Huang
'''


class Alias:
    def __init__(self, aliasClasses=None):
        if aliasClasses is not None:        
            self.load(aliasClasses)
        
    
    def load(self, aliasClasses):
        '''Load alias from the equivalence class of alias'''
        
        # The first one in a class is considered the basic name of this class
        self.basicNames = [tup[0] for tup in aliasClasses]
        
        # A dict with basic names as keys and aliases as values
        self.aliases = dict([(tup[0],list(tup)) for tup in aliasClasses])
        
        # A dict with aliases as keys and basic names as values        
        self.basics = {}
        for k in self.basicNames:            
            for vi in self.aliases[k]:
                if vi in self.basics.keys():
                    raise KeyError("Duplicated definition of alias '%s' of '%s" 
                                   % (vi, k))
                self.basics[vi] = k
    
    
    def loadFromFile(self, fileName):
        '''Load from a file of which lines are comma-separated'''
        
        fp = open(fileName, 'rb')
        aliasClasses = []        
        for line in fp:
            stripped = line.strip()
            if stripped == "":  # Skip any blank line
                continue
            aliasClasses.append(stripped.split(','))        
        self.load(aliasClasses)
    
    
    def getBasics(self):        
        return self.basics
    
    
    def getAliases(self):
        return self.aliases
    
    
    def getBasicName(self, aliasName):        
        return self.basics.get(aliasName, aliasName)

    
    def getAliasNames(self, basicName):
        return self.aliases.get(basicName, [basicName])
    
    
    def getMatchedAlias(self, basicName, matchList):
        '''Fine the first alias occurring in the list for a basic name'''
        
        matchSet = set(matchList)
        aliases = self.getAliasNames(basicName)   
        for alias in aliases:
            if alias in matchSet:
                return alias
        return None
    
    
#chromClasses = [('1', 'chr1'),
#                ('2', 'chr2'),
#                ('3', 'chr3'),
#                ('4', 'chr4'),
#                ('5', 'chr5'),
#                ('6', 'chr6'),
#                ('7', 'chr7'),
#                ('8', 'chr8'),
#                ('9', 'chr9'),
#                ('10', 'chr10'),
#                ('11', 'chr11'),
#                ('12', 'chr12'),
#                ('13', 'chr13'),
#                ('14', 'chr14'),
#                ('15', 'chr15'),
#                ('16', 'chr16'),
#                ('17', 'chr17'),
#                ('18', 'chr18'),
#                ('19', 'chr19'),
#                ('X', 'chrX'),
#                ('Y', 'chrY'),
#                ('M', 'MT', 'chrM', 'chrMT'),
#                ]
#chromAliases = Alias(chromClasses)
#print(chromAliases.getBasics())
#print(chromAliases.getAliases())
    
