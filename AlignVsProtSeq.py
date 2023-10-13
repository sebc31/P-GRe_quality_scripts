'''
This script was used to compare the alignments between pseudogenes' protein sequences predicted by P-GRe and their TAIR10 
equivalent. Only pseudogenes common to P-GRe and TAIR10 are compared (i.e. TAIR10 pseudogenes that were overlapped by at 
least 60% of their length by predictions).

Files used and origins:
common.pseudoprot.PGRE.faa and common.pseudoprot.TAIR10.faa are the FASTA sequences of the protein encoded by the 
 pseudogenes which were both found by P-GRe and annotated by TAIR10. For the TAIR10 sequences, they were generated by 
 assembling the genomic sequences annotated as "pseudogenic_exon" for each pseudogene and translated in silico. Note 
 that despite the name "pseudogenic_exon", these type of sequences are most likely rather pseudo-CDS.
Arabidopsis_thaliana.TAIR10.pep.all.fa is the proteome sequences file, downloaded from Ensembl Plant

Output of this script is directed to the strandard output. A copy of the result is saved on the AlignVsProtSeq.txt file
'''

import math
import re
import Bio.SeqIO
import statistics
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import pairwise2
from Bio.Align import substitution_matrices
from pprint import pprint

########################################################################################################################
####################################################### FUNCTIONS ######################################################
########################################################################################################################

def localScore(dic,prot):
    global prtDic, blosum62
    alignement = pairwise2.align.globalds( \
        prtDic[prot]['seq'], \
        dic[prot]['parent_seq'], \
        blosum62, -0.5, -0.1, one_alignment_only=True, score_only=True, penalize_end_gaps=False)
    return alignement

def getLen(file):
    global prtDic
    for sequence in list(Bio.SeqIO.parse(file, 'fasta')):
        id = sequence.id.split("|")[0]
        prtDic[id] = {}
        prtDic[id]['seq'] = sequence.seq
        prtDic[id]['len'] = len(sequence.seq)

def stat(blastres):
    res = {}
    with open(blastres) as file:
        for line in file:
            line = line.split("\t")
            if line[0] not in res:
                if "|" in line[0]: line[0]=line[0].split("|")[0]
                res[line[0]] = {}
                res[line[0]]['nb_id'] = 0
                res[line[0]]['percent_id'] = 0.0
                res[line[0]]['coverage_aln'] = 0.0
                res[line[0]]['coverage_hit'] = 0.0
                res[line[0]]['total_hits_len'] = set()
                res[line[0]]['total_hits_len_w_overlap'] = 0
                res[line[0]]['left'] = 0
                res[line[0]]['global_aln_score'] = 0.0
                res[line[0]]['right'] = 0
                res[line[0]]['parent_seq'] = ""
                res[line[0]]['parent_len'] = 0
                res[line[0]]['total_aln_len'] = 0
                res[line[0]]['dist_from_start'] = math.inf
                res[line[0]]['dist_from_end'] = math.inf
                res[line[0]]['bit_sum'] = 0
            res[line[0]]['parent_len'] = prtDic[line[1]]['len']
            res[line[0]]['parent_seq'] = prtDic[line[1]]['seq']
            res[line[0]]['nb_id'] += int(float(line[2]) / 100 * int(line[3]))
            res[line[0]]['left'] = min(int(line[8]), res[line[0]]['left'])
            res[line[0]]['right'] = max(int(line[9]), res[line[0]]['right'])
            for i in range(int(line[6]),int(line[7])+1):  # avoid hits overlap
                res[line[0]]['total_hits_len'].add(i)
            res[line[0]]['total_hits_len_w_overlap'] += int(line[3])
            res[line[0]]['dist_from_start'] = min(int(line[8]), res[line[0]]['dist_from_start'])
            res[line[0]]['dist_from_end'] = min(abs(res[line[0]]['parent_len'] - int(line[9])), res[line[0]]['dist_from_end'])
            res[line[0]]['bit_sum'] += float(line[11].replace('\n', ''))
    for result in res:
        res[result]['total_hits_len']=len(res[result]['total_hits_len'])
        res[result]['total_aln_len'] = res[result]['right'] - res[result]['left']
        res[result]['coverage_aln'] = res[result]['total_aln_len'] / res[result]['parent_len']
        res[result]['coverage_hit'] = res[result]['total_hits_len'] / res[result]['parent_len']
        res[result]['percent_id'] = res[result]['nb_id'] / res[result]['total_hits_len_w_overlap']
    return res

def printStat(dic):
    percent_id=[]
    coverage_aln=[]
    coverage_hit = []
    dist_from_start=[]
    dist_from_end=[]
    bit_sum=[]
    global_score=[]
    for id in dic:
        percent_id.append(dic[id]['percent_id'])
        coverage_aln.append(dic[id]['coverage_aln'])
        coverage_hit.append(dic[id]['coverage_hit'])
        dist_from_start.append(dic[id]['dist_from_start'])
        dist_from_end.append(dic[id]['dist_from_end'])
        bit_sum.append(dic[id]['bit_sum'])
        global_score.append(dic[id]['global_aln_score'])
    print('hits % of identity: '+ str(round(float(statistics.mean(percent_id)),2)))
    print('coverage hit: ' + str(round(float(statistics.mean(coverage_hit)), 2)))
    print('coverage aln: ' + str(round(float(statistics.mean(coverage_aln)), 2)))
    print('distance from start: ' + str(round(float(statistics.mean(dist_from_start)), 2)))
    print('distance from end: ' + str(round(float(statistics.mean(dist_from_end)), 2)))
    print('bit sum: ' + str(round(float(statistics.mean(bit_sum)), 2)))
    print('global score: ' + str(round(float(statistics.mean(global_score)), 2)))


########################################################################################################################
################################################### LOCAL ALIGNMENTS ###################################################
########################################################################################################################

# Retrieve the knwon pseudogene overlapped by the P-GRe predictions (previously added in the prediction headers via
# Linux command lines)
pseudoMerg={}
for sequence in list(Bio.SeqIO.parse("common.pseudoprot.PGRE.faa", 'fasta')):
    header=sequence.id.split("|")
    if header[1] not in pseudoMerg:
        pseudoMerg[header[1]]=[]
    pseudoMerg[header[1]].append(header[0].replace("_PGRe",""))

# Retrieve proteins length to compute coverage later
prtDic={}
getLen("../P-GRe/tmp/Arabidopsis_thaliana.TAIR10.pep.all.fa")
getLen("common.pseudoprot.PGRe.faa")
getLen("common.pseudoprot.TAIR10.faa")

# Local alignments: P-GRe predictions vs A. thaliana proteome / TAIR10 pseudogenes vs A. thaliana proteome. Path to BLAST
# binary as been modified for confidentiality
cmd_blastp = NcbiblastpCommandline(cmd="blastp", query="common.pseudoprot.PGRe.faa", out="PGRe.tsv", outfmt=6, db="Arabidopsis_thaliana.TAIR10.pep.all.fa", max_target_seqs=1)
cmd_blastp()
cmd_blastp = NcbiblastpCommandline(cmd="blastp", query="common.pseudoprot.TAIR10.faa", out="TAIR10.tsv", outfmt=6, db="Arabidopsis_thaliana.TAIR10.pep.all.fa", max_target_seqs=1)
cmd_blastp()

# Retrieve data from the BLAST results and do some small calculations
pgrestat=stat("PGRe.tsv")
tairstat=stat("TAIR10.tsv")

########################################################################################################################
############################################### SEMI-GLOBAL ALIGNMENTS #################################################
########################################################################################################################

# Semi-global alignment is performed between each pseudogene sequences and the sequence of the protein that aligned best
# during the local alignment step ("parent"). Parent id and sequence where previously stocked and associated with each 
# pseudogenes by the stat() function.

blosum62 = substitution_matrices.load('BLOSUM62')
i = 0; sum = 0

for prot in pgrestat:
    pgrestat[prot]['global_aln_score'] = localScore(pgrestat,prot)

for prot in tairstat:
    tairstat[prot]['global_aln_score'] = localScore(tairstat,prot)
    sum += tairstat[prot]['global_aln_score']

# In the case where a P-GRe prediction overlaps several pseudogenes annotated by TAIR10, each pseudogene annotated by
# TAIR10 is assigned as a score the sum of the scores of these pseudogenes. The reverse is not done. Note that this 
# artificially boosts the score of pseudogenes annotated by TAIR10.
for prot in pseudoMerg:
    bit_sum = 0
    global_sum = 0
    if len(pseudoMerg[prot]) > 1:
        for pg in pseudoMerg[prot]:
            pg_id=pg+"_TAIR10"
            if pg_id in tairstat:
                bit_sum += tairstat[pg_id]['bit_sum']
                global_sum += tairstat[pg_id]['global_aln_score']
        for pg in pseudoMerg[prot]:
            pg_id = pg + "_TAIR10"
            if pg_id in tairstat:
                tairstat[pg_id]['bit_sum'] = bit_sum
                tairstat[pg_id]['global_aln_score'] = global_sum

# Each pseudogenes annotated by P-GRe is compared to its TAIR10 counterpart, and the number times the semi-global
# alignment score is greater for P-GRe or for TAIR10 is counted
pgre_win = 0
tair_win = 0
total = 0
for prot in tairstat:
    pgre_id = prot.replace('_TAIR10','_PGRe')
    if pgre_id in pgrestat:
        total += 1
        if pgrestat[pgre_id]['global_aln_score']>tairstat[prot]['global_aln_score']:
            pgre_win += 1
        elif pgrestat[pgre_id]['global_aln_score']<tairstat[prot]['global_aln_score']:
            tair_win += 1

########################################################################################################################
###################################################### STAT OUTPUT #####################################################
########################################################################################################################

print("-----P-GRe-----")
printStat(pgrestat)
print("global score win: "+str(pgre_win)+" ("+str(pgre_win/total*100)+"%)")
print("-----TAIR10-----")
printStat(tairstat)
print("global score win: "+str(tair_win)+" ("+str(tair_win/total*100)+"%)")
