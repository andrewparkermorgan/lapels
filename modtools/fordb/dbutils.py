import sqlite3

db="/playpen/data/newgenes.db"

chromMap = dict([(str(i),i) for i in range(1,20)]+[('X',20),('Y',21),('M',22)])
strainMap = {'A':'aAJ',
           'A/J': 'aAJ',
           'AJ': 'aA_J',
           'B':'aC57BL',
           'C57BL/6J':'aC57BL',
           'B6':'aC57BL',
           'C':'a129S1',
           '129S1/SvImJ':'a129S1',
           '129':'a129S1',
           '129S1':'a129S1',
           'D':'aNOD',
           'NOD/ShiLtJ':'aNOD',
           'NOD/LtJ':'aNOD',
           'NOD':'aNOD',
           'E':'aNZO',
           'NZO/HlLtJ':'aNZO',
           'NZO':'aNZO',
           'F':'aCAST',
           'CAST/EiJ':'aCAST',
           'CAST':'aCAST',
           'G':'aPWK',
           'PWK/PhJ':'aPWK',
           'PWK':'aPWK',
           'H':'aWSB',
           'WSB/EiJ':'aWSB',
           'WSB':'aWSB'}


def readIndelsFromDB(db, chromoID, sample):
    assert sample.startswith('a')
    db = sqlite3.connect(db)
    db.row_factory = sqlite3.Row
    cursor = db.cursor()
    ## Exclude ref. allele(0).
    query = "select position, alleles, %s from SangerIndels where chromoID = %d and %s <> 0"
    cursor.execute(query % (sample, chromoID, sample))
    answers = cursor.fetchall()
    ## From 1-based position in DB to 0-based
    response = [(p-1, a.split(';')[0], a.split(';')[i]) for p, a, i in answers]
    db.commit()
    db.close()
    return response


def readSNPsFromDB(db, chromoID, sample):
    assert sample.startswith('a')
    db = sqlite3.connect(db)
    db.row_factory = sqlite3.Row
    cursor = db.cursor()
    ## Exclude ref. allele(0) and hets(-1).
    query = "select position, alleles, %s from SangerSNP where chromoID = %d and %s <> 0 and %s <> -1"
    cursor.execute(query % (sample, chromoID, sample, sample))
    answers = cursor.fetchall()
    ## From 1-based position in DB to 0-based
    response = [(p-1, a[0], a[i]) for p, a, i in answers]
    db.commit()
    db.close()
    return response


if __name__ == '__main__':
    indels=readIndelsFromDB(db, 1, strainMap['CAST'])
    snps=readSNPsFromDB(db, 1, strainMap['CAST'])
    
    print("indels:")
    for i in range(10):
        print(indels[i])
    print("snps:")
    for i in range(10):
        print(snps[i])
    
