import sys
import argparse
import numpy
import mappy
import os

parser = argparse.ArgumentParser()

parser.add_argument('--actualgenomeKmers', '-a', type=str, action='store', help='actual genome kmer counts from reference genome as .jf file')
parser.add_argument('--genomeFastq', '-g', type=str, action='store', help='genome Fastq')
parser.add_argument('--cutRunFastq', '-c', type=str, action='store', help='cutrun Fastq')
parser.add_argument('--GenomeKmers', '-G', type=str, action='store', help='genome kmers (as jf file)')
parser.add_argument('--CutRunKmers', '-C', type=str, action='store', help='cutrun kmers (as tab file)')
parser.add_argument('--KmerEnrichmentFile', '-K', type=str, action='store', help='Kmer Enrichment File')

args = parser.parse_args()

ActualGenomeKmers = args.actualgenomeKmers
GenomeFastq = args.genomeFastq
CutRunFastq = args.cutRunFastq
GenomeKmers = args.GenomeKmers
CutRunKmers = args.CutRunKmers
outfile = args.KmerEnrichmentFile

###############################################################################
                    ##########counting lines ######
###############################################################################


####Counting lines for genome fastq (control)
print('reading genome fastq lines')
if GenomeFastq.isdigit(): 
    lineCounter1=int(GenomeFastq) #will usually be false, not sure when this is actually something that is true?
else:
    lineCounter1=0
    for name,seq,q in mappy.fastx_read(GenomeFastq): #name, sequence, quality score using mappy
        lineCounter1+=1                              #add 1 to line counter for every sequence 
print(lineCounter1, 'reads in genome Fastq')

####Counting lines for cut and run experimental fastq 
print('reading cutrun fastq lines')
if CutRunFastq.isdigit():
    lineCounter2=int(CutRunFastq)
else:
    lineCounter2=0
    for name,seq,q in mappy.fastx_read(CutRunFastq):
        lineCounter2+=1
print(lineCounter2, 'reads in cutRun Fastq')

###############################################################################
##                              Calculating Enrichment                      ##
###############################################################################

outmain=open(outfile,'w') #writing out a .tab file

i = 0
cutrunKmerList=[]
cutrunCount={} #{kmer: count}
counter = 0 
cutRunKmersLines = 0 #keeps count of reads 
kmerList = [] #kmer cut and run count

#counting lines from cutrun kmer count.tab
print('counting lines in cutrun data')

for line in open(CutRunKmers): #for line in cutrun kmers fasta file
    cutRunKmersLines += 1      #counts lines

print('reading cutrun kmers')
for line in open(CutRunKmers): #for line in cutrun kmers fasta file
    cutrunLine=line.strip().split('\t')
    kmer=cutrunLine[0] 
    cutrunCount[kmer]=int(cutrunLine[1])
    counter+=1 #counting kmers
    kmerList.append((kmer,cutrunCount)) 
    if counter%1000000==0 or counter==cutRunKmersLines: #includes end of file with extra counts
        print(counter)
        out = open('kmer.fasta', 'w')
        for kmer1,cutrunCount in kmerList:
            out.write('>'+str(cutrunCount[kmer1])+'\n'+kmer1+'\n')
        out.close()

#subsetting kmer count of kmers only found in the reference genome
        os.system(f"jellyfish query {ActualGenomeKmers} -s kmer.fasta > query_output.tab") 
        print('finished jellyfish')

#Create Genome Dictionary from query_output.tab
        GenomeDict={} #kmers only found in the actual reference genome with their counts from the cut and run exp data
        for line in open("query_output.tab"):
            GenomeLine=line.strip().split(' ')
            GenomeCount=int(controlLine[1])
            Genomekmer=controlLine[0]
            GenomeDict[Genomekmer]=GenomeCount #adding actual genome kmer to genome count 

#Jellyfish query kmers found in illumina cut and run control
        print('starting jellyfish')
        os.system(f"jellyfish query {GenomeKmers} -s kmer.fasta > query_output.tab") #only looking at kmers found also in the control 
        print('finished jellyfish')

#parsing query output of kmers only found in control 
        for line in open("query_output.tab"):
            controlLine=line.strip().split(' ')
            controlCount=int(controlLine[1])
            controlkmer=controlLine[0]

#calculating enrichment by (cut and run kmer/reads)/(cut and run kmers found in illumina genomic/reads)
            ## only calculate if we have a non-zero controlCount
            if controlCount != 0:
                enrichment= (cutrunCount[controlkmer]/lineCounter2)/(controlCount/lineCounter1)
                outmain.write(controlkmer+'\t'+str(enrichment)+'\t'+str(GenomeDict[controlkmer])+'\n')
        kmerList=[]
        cutrunCount={}
        print('finished', counter)

