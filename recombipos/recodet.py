#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import gzip
import sys

motherid = sys.argv[2]
fatherid = sys.argv[3]
window   = 200

if len(sys.argv) == 5:
    window   = int(sys.argv[4])

mother = 0
father = 0

def getGT(info, field):
    tokens = info.split(":")
    gtcol = -1
    for i in range(len(tokens)):
        if tokens[i] == "GT":
            gtcol = i
            break
    if gtcol == -1:
        return 0
    gt = field.split(":")[gtcol]
    if len(gt) < 3:
        return 0
    if gt[0] == gt[2] and gt[0] == '0':
        return 0
    if gt[0] == '0' and gt[2] != '0':
        return 1
    return 2

def diffcounter(kmlabels, prevkm):
    dc = 0
    for j in range(len(prevkm)):
        if kmlabels[j] != prevkm[j]:
            dc += 1
    return dc

def recombipos(classes, prevkm, selected, sub, window):
    prevgroupmates = list()
    for i in range(len(prevkm)):
        if prevkm[i] == prevkm[selected]:
            prevgroupmates.append(i)
    thisgroupmates = list()
    for i in range(len(classes)):
        if classes[i] == classes[selected]:
            thisgroupmates.append(i)
    for i in range(window):
        same = 0
        same2 = 0
        for p in prevgroupmates:
            if sub[p,i] == sub[selected,i]:
                same += 1
        for t in thisgroupmates:
            if sub[t,i] == sub[selected,i]:
                same2 += 1
        ratio1 = float(same) / float(len(prevgroupmates))
        ratio2 = float(same2) / float(len(thisgroupmates))
        if ratio2 > ratio1:
            return i
    return window

def processMatrix(matrix, poslist, chrom, window, sybnames):
    m = np.array(matrix)
    prevkm = list()
    for i in range(len(matrix)-window):
        sub = m[i:i+window,].T
        pca = PCA(n_components=2)
        pr  = pca.fit_transform(sub)
        km  = KMeans(n_clusters = 2, random_state = 0).fit(pr)
        # KMeans random, so the group ids can be the opposite
        dc = diffcounter(km.labels_, prevkm)

        if dc > 0 and dc < len(prevkm):
            # we have a recombination, but in which siblings?
            if dc > len(prevkm) / 2:
                classes = [int(((x-0.5)*-2)/2.0+0.5) for x in km.labels_]
            else:
                classes = km.labels_
            for j in range(len(classes)):
                if classes[j] != prevkm[j]:
                    #print("Recombination in",j,"at",chrom,poslist[i],sybnames[j])
                    # where is the exact position of the recombination?
                    index = recombipos(classes, prevkm, j, sub, window)
                    print(chrom, poslist[i+index], sybnames[j])
        prevkm = km.labels_

matrix = list()
prevchrom = ''
poslist = list()
sybnames = list()

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
        processMatrix(matrix, poslist, prevchrom, window, sybnames)
        matrix = list()
        poslist = list()

    pos   = int(line[1])
    info  = line[8]
    if len(line[3]) > 1 or len(line[4]) > 1:
        continue

    if getGT(info, line[mother]) == 1 and getGT(info, line[father]) == 0:
        row = [0] * (len(line) - 9 - 2)
        ri = 0
        for i in range(9, len(line)):
            gt = getGT(line[8], line[i])
            if i == mother or i == father:
                continue
            row[ri] = gt
            ri += 1
        matrix.append(row)
        poslist.append(pos)

    prevchrom = chrom

processMatrix(matrix, poslist, prevchrom, window, sybnames)
