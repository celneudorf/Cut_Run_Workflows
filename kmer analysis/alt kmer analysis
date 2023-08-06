import mappy
import sys
import numpy
enrichedKmerDict={} #reading in kmers from txt
with open(sys.argv[2]) as f:
    for line in f:
        a=line.strip().split('\t')
        kmer=a[0]
        kmerEnrichment=float(a[1])
        enrichedKmerDict[kmer]=kmerEnrichment
fastaKmers={} #reading in Jellyfish (tab) dump of genome fasta
with open(sys.argv[3]) as f:
    for line in f:
        a=line.strip().split('\t')
        kmer=a[0]
        count=int(a[1])
        fastaKmers[kmer] = count
posDict={}
for name,seq,q in mappy.fastx_read(sys.argv[1]): #reads in x fasta
    if name not in posDict:
        posDict[name]={} #adding X chromosome into dictionary
    print(name)
    for index in  range(0,len(seq),50000): #taking whole sequence from fasta, breaking area into 50kb chunks
        print(name,index)
        posDict[name][index]=0 #name is chromosome index is
        subSeq=seq[index:index+50000] #chunking into 50kb
        for subIndex in range(0,len(subSeq),1):
            kmer=subSeq[subIndex:subIndex+25] #pulling every single 25 kmer in 50kb chunk, #kmer multiplicity, kmers counted in more stops , #coumt by the number of kmers and by count of enrichment, can we discount$
            if kmer in enrichedKmerDict.keys(): #checking to see if 25 kmer set
                print('yes')
                fastaCountForKmer = fastaKmers[kmer]
                kmerEnrichment=enrichedKmerDict[kmer]
                posDict[name][index]+= (kmerEnrichment / fastaCountForKmer)
counter=0
out=open(sys.argv[4],'w')
out.write('track type=wiggle_0\n')
for name in posDict:
    out.write('fixedStep chrom=%s start=1 step=50000 span=50000\n' %(name))
    for index in posDict[name]:
        counter+=1
#        out.write(name+'\t'+str(index)+'\t'+str(index+50000)+'\tbin'+str(counter)+'\t'+str(posDict[name][index])+'\n')
        out.write(str(posDict[name][index])+'\n')
out.close()

