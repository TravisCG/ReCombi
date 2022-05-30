#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import gzip
import sys

def getGT(info, field):
    tokens = info.split(":")
    gtcol = -1
    for i in range(len(tokens)):
        if tokens[i] == "GT":
            gtcol = i
            break
    if gtcol == -1:
        return 3
    gt = field.split(":")[gtcol]
    if len(gt) < 3:
        return 3
    if gt[0] == gt[2] and gt[0] == '0':
        return 0
    if (gt[0] == '0' and gt[2] == '1') or (gt[0] == '1' and gt[2] == '0'):
        return 1
    if gt[0] == '1' and gt[2] == '1':
        return 2
    return 3

def recombipos(prevkm, actualkm, sub):
    (h,w) = sub.shape
    for i in range(w):
        for syb in range(h):
            print(sub[syb,i], end=" ")

# Normalizing class. K-means some times swap the cluster IDs
# this method makes the clusters similar to the previous one
def normclasses(prevkm, actualkm):
    diff = 0
    negative = list()
    for i in range(len(prevkm)):
        if prevkm[i] != actualkm[i]:
            diff += 1
        if actualkm[i] == 1:   #FIXME can you make it branchless?
            negative.append(0)
        else:
            negative.append(1)

    if diff > len(prevkm) / 2:
        return negative
    return actualkm

def processMatrix(matrix, pmatrix, poslist, chrom, window, sybnames):
    m       = np.array(matrix)
    prevkm  = list()
    counter = 0
    sub     = np.empty([len(sybnames), window])
    pca     = PCA(n_components = 2)

    for i in range(len(matrix)):
        if pmatrix[i][0] == 1 and (pmatrix[i][1] == 0 or pmatrix[i][1] == 2):
            for j in range(len(sybnames)):
                sub[j,counter] = m[i,j]
            counter += 1
            if counter == window:
                # we have enough mutations, lets use it
                pr = pca.fit_transform(sub)
                km = KMeans(n_clusters = 2, random_state = 0).fit(pr)
                classes = normclasses(prevkm, km.labels_)
                recombipos(prevkm, classes, sub)
                # prepare for the next cycle
                counter = 0
                sub = np.empty([len(sybnames), window])
                prevkm = classes

##### Main ######
motherid = sys.argv[2] #T1419022
fatherid = sys.argv[3] #T2319031
window   = 200

if len(sys.argv) == 5:
    window   = int(sys.argv[4])

mother = 0
father = 0

matrix    = list()
pmatrix   = list()
prevchrom = ''
poslist   = list()
sybnames  = list()

for vcfline in gzip.open(sys.argv[1]):
    line = vcfline.decode('utf-8').rstrip().split("\t")

    if line[0].startswith('##'):
        continue

    if line[0][0] == '#':
        for i in range(9,len(line)):
            if line[i] == motherid:
                mother = i
            elif line[i] == fatherid:
                father = i
            else:
                sybnames.append(line[i])
        continue

    chrom = line[0]
    if prevchrom != chrom and prevchrom != '':
        processMatrix(matrix, pmatrix, poslist, prevchrom, window, sybnames)
        matrix = list()
        poslist = list()
        pmatrix = list()

    pos   = int(line[1])
    info  = line[8]
    if len(line[3]) > 1 or len(line[4]) > 1:
        continue

    if getGT(info, line[mother]) != 3 and getGT(info, line[father]) != 3:
        row = [0] * (len(line) - 9 - 2)
        prow = [0,0]
        ri = 0
        for i in range(9, len(line)):
            gt = getGT(line[8], line[i])
            if i == mother:
                prow[0] = gt
                continue
            if i == father:
                prow[1] = gt
                continue
            row[ri] = gt
            ri += 1
        matrix.append(row)
        pmatrix.append(prow)
        poslist.append(pos)

    prevchrom = chrom

processMatrix(matrix, pmatrix, poslist, prevchrom, window, sybnames)
