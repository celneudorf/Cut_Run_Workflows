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
for name,seq,q in mappy.fastx_read(sys.argv[1]): #opening fasta file
    if name not in posDict:
        posDict[name]={}                         #creating dictionary of positions 
    print(name)
    for index in range(0,len(seq),stepsize):     #iterating through step size 
        print(name,index)
        GenomekmerSet=set()
        posDict[name][index]=[0,0]               
        subSeq=seq[index:index+stepsize] 
        for subIndex in range(0,len(subSeq),1):  
            kmer=subSeq[subIndex:subIndex+25]     
            GenomekmerSet.add(kmer)              #adding 25 bp kmer to kmer set
            if kmer in kmerSet:                  #counting how many times 25 bp kmer in window
                posDict[name][index][0]+=1       
        out=open('kmers.fasta','w')
        for kmer in GenomekmerSet:
            out.write(f'>1\n{kmer}\n')           #writing out fasta file with kmer and count
        out.close()
        
        os.system(f"jellyfish query {ActualGenomeKmers} -s kmers.fasta > query_output.tab") #querying unique kmers 
        print('finished jellyfish')
        for line in open("query_output.tab"):
            GenomeLine=line.strip().split(' ')
            count=int(GenomeLine[1])
            if count==1:
                posDict[name][index][1]+=1       #counting how many times unique kmer in genome
        if posDict[name][index][1]>0:           
            print(posDict[name][index][0],posDict[name][index][1],posDict[name][index][0]/posDict[name][index][1])
            
counter=0
out=open(sys.argv[3],'w')                        #writing out wigfile 
out.write('track type=wiggle_0\n') 
for name in posDict:
    out.write(f'fixedStep chrom={name} start=1 step={stepsize} span={stepsize}\n')
    for index in posDict[name]:                  
        counter+=1
        if posDict[name][index][1]>0:           #if kmer in reference genome
            out.write(str(posDict[name][index][0]/posDict[name][index][1])+'\n') 
        else:
            out.write('0\n')
out.close()
