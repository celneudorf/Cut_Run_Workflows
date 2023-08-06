
import mappy
import sys
import numpy
import os
ActualGenomeKmers=sys.argv[4]

kmerSet=set()
with open(sys.argv[2]) as f:
    for line in f:
        a=line.strip().split('\t')
        kmer=a[0]
        kmerSet.add(kmer)


stepsize = 10000
posDict={}
for name,seq,q in mappy.fastx_read(sys.argv[1]):
    if name not in posDict:
        posDict[name]={}
    print(name)
    for index in  range(0,len(seq),stepsize):
        print(name,index)
        GenomekmerSet=set()
        posDict[name][index]=[0,0]
        subSeq=seq[index:index+stepsize]
        for subIndex in range(0,len(subSeq),1):
            kmer=subSeq[subIndex:subIndex+25]
            GenomekmerSet.add(kmer)
            if kmer in kmerSet:
                posDict[name][index][0]+=1
