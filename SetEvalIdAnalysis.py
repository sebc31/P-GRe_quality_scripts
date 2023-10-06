'''
This script was used to create the different set of predicted pseudogenes based on different criterion:
- Do they cover (by at least 60%) a known pseudogene ?
- Do they cover (by at least 60%) a knwon transposable element (TE) gene ?
- Do their sequences align with a known TE sequence ?
The last part also computs the mean E-value for some of this sets based on the local alignments made by P-GRe

Files used and origins:
pseudogenes.gff is the GFF file returned by P-GRe
TAIR10_GFF3_genes_transposons.gff is the TAIR10 GFF file downloaded from Ensembl Plant with supplementary informations
 concerning TEs
TAIR10.pseudogenes.gff is a TAIR10 GFF file. It was modified with the grep command on a Linux system to keep only
 pseudogenes-related informations.
idFromBlastTEG.id is the id of predicted pseudogenes that aligned with a known (annotated) TE gene sequences from
 TAIR10. It was obtained from the blastCDSvsTEG.blast file, which is the ouput of BLAST alignment of all predicted
 pseudogenes against all TE gene sequences
proteome_vs_masked_genome_tblastn.blast is a P-GRe working file. It is the output of the very first local alignment made
 by P-GRe (see README.md on P-GRe GitHub page)

Output of this script is directed to the strandard output. A copy of the result is saved on the SetsEvalIdAnalysis.txt
 file
'''

########################################################################################################################
####################################################### FUNCTIONS ######################################################
#"######################################################################################################################

def splitLine(line):
    line = line.split('\t')
    chr = line[0]
    structure = line[2]
    left = min(int(line[3]), int(line[4]))
    right = max(int(line[3]), int(line[4]))
    id = line[8].split(';')[0].replace('ID=','')
    return chr, structure, left, right, id

def getCovered60(covDic, line, structure_list, multiple):
    global positionsDic
    chr, structure, left, right, id = splitLine(line)
    if structure in structure_list:
        length = right - left + 1
        for i in range(0, len(positionsDic[chr]['id'])):
            if (right > positionsDic[chr]['left'][i] and right <= positionsDic[chr]['right'][i]) \
                    or (left >= positionsDic[chr]['left'][i] and left < positionsDic[chr]['right'][i]) \
                    or (left >= positionsDic[chr]['left'][i] and right <= positionsDic[chr]['right'][i]) \
                    or (left <= positionsDic[chr]['left'][i] and right >= positionsDic[chr]['right'][i]):
                if id not in covDic:
                    covDic[id] = {}
                    covDic[id]['cov'] = 0
                    covDic[id]['pg'] = []
                upper_left = max(left, positionsDic[chr]['left'][i])
                lower_right = min(right, positionsDic[chr]['right'][i])
                overlap = lower_right - upper_left + 1
                coverage = overlap / length
                if multiple:
                    covDic[id]['cov'] += coverage
                else:
                    covDic[id]['cov'] = max(coverage, covDic[id]['cov'])
                covDic[id]['pg'].append(positionsDic[chr]['id'][i])

def cov60Stats(covDic):
    res = set()
    nb_res = 0
    for ids in covDic:
        if covDic[ids]['cov'] >= 0.60:
            nb_res += 1
            for pgs in covDic[ids]['pg']:
                res.add(pgs)
    return res, nb_res

def saveRes(fileName, nb_res, PG, name):
    resFile = open(fileName, 'w')
    print(str(nb_res) + ' known ' + name + ' are overlapped by ' + str(len(PG)) + \
          ' unique predicted pseudogenes.')
    for ids in PG:
        resFile.write(ids + '\n')

def appendEvalDic(setToChange):
    global evalDic, eval, aln_length, identity
    evalDic[setToChange]['mean_Eval'].append(eval)
    evalDic[setToChange]['mean_aln_len'].append(aln_length)
    evalDic[setToChange]['mean_id'].append(identity)

########################################################################################################################
################################################### P-GRe predictions ##################################################
########################################################################################################################

# Get results from P-GRe output
positionsDic = {}
predictionsDic = {}
all_PG = set()
with open('pseudogenes.gff') as PG:
    for line in PG:
        chr, structure, left, right, id = splitLine(line)
        if structure == 'pseudogene':
            if chr not in positionsDic:
                positionsDic[chr] = {}
                positionsDic[chr]['left'] = []
                positionsDic[chr]['right'] = []
                positionsDic[chr]['id'] = []
            positionsDic[chr]['left'].append(left)
            positionsDic[chr]['right'].append(right)
            positionsDic[chr]['id'].append(id)
            all_PG.add(id)
            parent = line.split('\t')[8].split(';')[-1].replace('Parent_gene=','').replace('\n','').split(',')
            for parents in parent:
                if parents not in predictionsDic:
                    predictionsDic[parents] = {}
                    predictionsDic[parents]['left'] = []
                    predictionsDic[parents]['right'] = []
                    predictionsDic[parents]['pgre_id'] = []
                predictionsDic[parents]['left'].append(left)
                predictionsDic[parents]['right'].append(right)
                predictionsDic[parents]['pgre_id'].append(id)

########################################################################################################################
######################################### Covered known TE and PG ######################################################
########################################################################################################################

##### TE #####
covDicTE = {}
with open('TAIR10_GFF3_genes_transposons.gff') as TE:
    for line in TE:
        getCovered60(covDicTE, line, ['transposable_element', 'transposable_element_gene'], False)
PG_that_covers_TE, nb_res = cov60Stats(covDicTE)
saveRes('PG_overlapping_kTE.id', nb_res, PG_that_covers_TE, 'TE')

##### PG #####
covDicPG = {}
with open('TAIR10.pseudogenes.gff') as PG:
    for line in PG:
        getCovered60(covDicPG, line, ['pseudogene'], True)
PG_that_covers_PG, nb_res = cov60Stats(covDicPG)
saveRes('PG_overlapping_kPG.id', nb_res, PG_that_covers_PG, 'PG')

########################################################################################################################
################################################ SET CONSTRUCTIONS #####################################################
########################################################################################################################

PG_that_aligns_TE = set()
with open('idFromBlastTEG.id') as TEBl:
    for line in TEBl:
        id = line.replace('\n','')
        if id in all_PG:
            PG_that_aligns_TE.add(id)
print(str(len(PG_that_aligns_TE)) + ' unique predicted pseudogenes aligns with known TE sequences.', end = ' ')

PG_that_aligns_and_cover_TE = PG_that_aligns_TE.intersection(PG_that_covers_TE)
print(str(len(PG_that_aligns_and_cover_TE)) + ' of them also cover known TE.')

PG_that_aligns_and_dont_cover_TE = PG_that_aligns_TE.difference(PG_that_covers_TE)
print(str(len(PG_that_aligns_and_dont_cover_TE)) + ' of them don t cover known TE.')

PG_that_dont_align_or_cover_TE = all_PG.difference(PG_that_aligns_TE).difference(PG_that_covers_TE)
print(str(len(PG_that_dont_align_or_cover_TE)) + ' unique predicted pseudogenes don t align or cover any known TE.')

unknown_PG = PG_that_dont_align_or_cover_TE.difference(PG_that_covers_PG)
print(str(len(unknown_PG)) + ' unique predicted pseudogenes might be new.')

########################################################################################################################
########################################### TO SORT BLAST RESULTS ! ####################################################
############################## I know, this should be easier and integrated in P-GRe...#################################

evalDic = {}
# c = covers, a = aligns, k = known, d = don't
set_list = ['PG_c_kPG', 'PG_c_kTE', 'PG_dca_kTE', 'PG_dc_kPG_dcda_kTE']
for sets in set_list:
    evalDic[sets] = {}
    evalDic[sets]['max_Eval'] = 1
    evalDic[sets]['mean_Eval'] = []
    evalDic[sets]['mean_aln_len'] = []
    evalDic[sets]['mean_id'] = []

with open('proteome_vs_masked_genome_tblastn.blast') as bl:
    for line in bl:
        toks = line.split('\t')
        id = toks[0]
        identity = float(toks[2])
        aln_length = int(toks[3])
        left = min(int(toks[8]), int(toks[9]))
        right = max(int(toks[8]), int(toks[9]))
        eval = float(toks[10])
        if id in predictionsDic:
            for i in range(0, len(predictionsDic[id]['left'])):
                if (right > predictionsDic[id]['left'][i] and right <= predictionsDic[id]['right'][i]) \
                    or (left >= predictionsDic[id]['left'][i] and left < predictionsDic[id]['right'][i]) \
                    or (left >= predictionsDic[id]['left'][i] and right <= predictionsDic[id]['right'][i]) \
                    or (left <= predictionsDic[id]['left'][i] and right >= predictionsDic[id]['right'][i]):
                    if predictionsDic[id]['pgre_id'][i] in PG_that_covers_PG:
                        appendEvalDic('PG_c_kPG')
                    if predictionsDic[id]['pgre_id'][i] in PG_that_covers_TE:
                        appendEvalDic('PG_c_kTE')
                    if predictionsDic[id]['pgre_id'][i] in PG_that_aligns_and_dont_cover_TE:
                        appendEvalDic('PG_dca_kTE')
                    if predictionsDic[id]['pgre_id'][i] in unknown_PG:
                        appendEvalDic('PG_dc_kPG_dcda_kTE')
for sets in evalDic:
    print('\n-------------------------------------------\n' + sets + '\n-------------------------------------------')
    evalDic[sets]['max_Eval'] = max(evalDic[sets]['mean_Eval'])
    print('Maximum E-value: ' + str(evalDic[sets]['max_Eval']))
    for attributes in evalDic[sets]:
        if 'mean_' in attributes:
            evalDic[sets][attributes] = sum(evalDic[sets][attributes]) / len(evalDic[sets][attributes])
            print(attributes + ': ' + str(evalDic[sets][attributes]))

##### Export new PG (non-TE associated and not covering known PG) for more analysis
ukPG = open('ukPG.id', 'w')
for pg in unknown_PG:
    ukPG.write(pg + '\n')
