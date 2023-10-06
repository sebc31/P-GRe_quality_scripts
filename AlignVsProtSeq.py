# -*- coding: utf-8 -*-

import math
import re

import Bio.SeqIO
import statistics
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio import pairwise2
from Bio.Align import substitution_matrices
from pprint import pprint

#Fun
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

#See where P-GRe has merged pseudogenes, so we can sum bit-score for TAIR10 "splitted" pseudogenes
pseudoMerg={}
for sequence in list(Bio.SeqIO.parse("common.pseudoprot.PGRE.faa", 'fasta')):
    header=sequence.id.split("|")
    if header[1] not in pseudoMerg:
        pseudoMerg[header[1]]=[]
    pseudoMerg[header[1]].append(header[0].replace("_PGRe",""))

#Retrieve prot len
prtDic={}
getLen("../P-GRe/tmp/Arabidopsis_thaliana.TAIR10.pep.all.fa")
getLen("common.pseudoprot.PGRe.faa")
getLen("common.pseudoprot.TAIR10.faa")

#BLAST
cmd_blastp = NcbiblastpCommandline(cmd="C:\\Users\\sebastien.cabanac\\AppData\\Local\\NBCI\\blast-BLAST_VERSION+\\bin\\blastp.exe", query="common.pseudoprot.PGRe.faa", out="PGRe.tsv", outfmt=6, db="Arabidopsis_thaliana.TAIR10.pep.all.fa", max_target_seqs=1)
cmd_blastp()
cmd_blastp = NcbiblastpCommandline(cmd="C:\\Users\\sebastien.cabanac\\AppData\\Local\\NBCI\\blast-BLAST_VERSION+\\bin\\blastp.exe", query="common.pseudoprot.TAIR10.faa", out="TAIR10.tsv", outfmt=6, db="Arabidopsis_thaliana.TAIR10.pep.all.fa", max_target_seqs=1)
cmd_blastp()

#Compute stat
pgrestat=stat("PGRe.tsv")
tairstat=stat("TAIR10.tsv")

#Alignement semi-global
blosum62=substitution_matrices.load('BLOSUM62')
i=0; sum=0
for prot in pgrestat:
    pgrestat[prot]['global_aln_score'] = localScore(pgrestat,prot)

for prot in tairstat:
    tairstat[prot]['global_aln_score'] = localScore(tairstat,prot)
    sum += tairstat[prot]['global_aln_score']

#Sum bit-score and global score for splitted pseudogene in ref
print(tairstat['AT1G34797.1_TAIR10'])
print(tairstat['AT1G34803.1_TAIR10'])
print(tairstat['AT1G34808.1_TAIR10'])

print(pgrestat['AT1G34797.1_PGRe'])
print(pgrestat['AT1G34803.1_PGRe'])
print(pgrestat['AT1G34808.1_PGRe'])

for prot in pseudoMerg:
    bit_sum=0
    global_sum=0
    if len(pseudoMerg[prot])>1:
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

#Compare two by two
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

print("-----P-GRe-----")
printStat(pgrestat)
print("global score win: "+str(pgre_win)+" ("+str(pgre_win/total*100)+"%)")
print("-----TAIR10-----")
printStat(tairstat)
print("global score win: "+str(tair_win)+" ("+str(tair_win/total*100)+"%)")
