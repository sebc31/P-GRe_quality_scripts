

#### FUNCTIONS

def newSet(file):
    setToCreate = set()
    with open(file) as fileToRead:
        for line in fileToRead:
            setToCreate.add(line.replace('\n', ''))
    return setToCreate

def initializeDict(dict):
    dict['Completness'] = {}
    dict['Completness']['Copy'] = 0
    dict['Completness']['Fragment'] = 0
    dict['Completness']['Fragments'] = 0
    dict['Completness']['Fragment or degraded copy'] = 0
    dict['Type'] = {}
    dict['Type']['Chimeric pseudogene'] = 0
    dict['Type']['Duplicated pseudogene'] = 0
    dict['Type']['(Iso)retropseudogene'] = 0
    dict['Type']['Retropseudogene'] = 0
    dict['Type']['Unknown'] = 0

def addToDic(dic, completness, type):
    dic['Completness'][completness] += 1
    dic['Type'][type] += 1

def dictRatioCompute(dic):
    totalType = sum(dic['Type'].values())
    totalCompletness = sum(dic['Completness'].values())
    print('----- Completness:')
    for completnessKind, completnessNb in dic['Completness'].items():
        print(completnessKind, round(completnessNb / totalCompletness * 100, 2), completnessNb)
    print('----- Type:')
    for typeKind, typeNb in dic['Type'].items():
        print(typeKind, round(typeNb / totalType * 100, 2), typeNb)
    print('')

#### STEP1: import the different sets of pseudogenes
PGcPG = newSet('PG_overlapping_kPG.id')  # set of P-GRe pseudogenes overlapping TAIR10 pseudogenes by at least 60%
PGcTE = newSet('PG_overlapping_kTE.id')  # set of P-GRe pseudogenes overlapping TAIR10 transposable elements (TEs) by at
                                         # least 60%
PGdcPGdcdaTE = newSet('ukPG.id')  # set of P-GRe pseudogenes that don't cover any TAIR10 pseudogenes or any TAIR10
                               # transposable elements, and that don't align with any transposable element, aka
                               # "unknown" pseudogenes

#### STEP2: retrieve the number of each kind of completness and type for each set
PGcPGComputation = {}
initializeDict(PGcPGComputation)
PGcTEComputation = {}
initializeDict(PGcTEComputation)
PGdcPGdcdaTEComputation = {}
initializeDict(PGdcPGdcdaTEComputation)
totalComputation = {}
initializeDict(totalComputation)
with open('pseudogenes.info') as info:
    next(info)  # Skip header
    for line in info:
        pg, completness, type, parent = line.split('\t')
        if pg in PGcPG:
            addToDic(PGcPGComputation, completness, type)
        if pg in PGcTE:
            addToDic(PGcTEComputation, completness, type)
        if pg in PGdcPGdcdaTE:
            addToDic(PGdcPGdcdaTEComputation, completness, type)
        addToDic(totalComputation, completness, type)

#### STEP3: print the % for each set
print('----------------- P-GRe pseudogenes overlapping known pseudogene -----------------')
dictRatioCompute(PGcPGComputation)
print('--------------------- P-GRe pseudogenes overlapping known TE ---------------------')
dictRatioCompute(PGcTEComputation)
print('----- P-GRe predictions not overlapping known PG or TE, not aligning with TE -----')
dictRatioCompute(PGdcPGdcdaTEComputation)
print('------------------------------ Total predictions ---------------------------------')
dictRatioCompute(totalComputation)

print('/!\ Please note that P-GRe classify as "Chimeric pseudogene" pseudogenes that have more than one parent gene.\n'
      'Having more than one parent can indeed be due to true chimeric case, but also to sequence divergence. /!\\')
