import sys
import numpy
import scipy.stats

enrichments=[]
counter=0
total=0
for line in open(sys.argv[1]):
    total+=1
    a=line.strip().split('\t')
    enrichment=float(a[1])
    genomeCount = int(a[2])
    enrichments.append(enrichment)
    if enrichment>3:
        counter+=1
        if genomeCount == 1:
            print(a[0])

#print(numpy.median(enrichments),scipy.stats.median_abs_deviation(enrichments))
