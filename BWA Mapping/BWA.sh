set -eux #if something fails stop running script
sample_name=$1
r1_fastq=$2
r2_fastq=$3
ref_file=$4
output_bw_file=$5 #year_month_date_sample_controlorcutrun

java -jar ~/Downloads/Trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 10 \
   $r1_fastq \
   $r2_fastq \
   R1.paired.fastq.gz R1.unpaired.fastq.gz \
   R2.paired.fastq.gz R2.unpaired.fastq.gz \
   CROP:150 ILLUMINACLIP:/home/k2so/Downloads/Trimmomatic/Trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:20 \
   LEADING:3 TRAILING:3 MINLEN:36

export PATH="/home/k2so/bwa:$PATH"

bwa mem \
    -t 40 \
    -k 50 \
    -c 1000000 \
    $ref_file \
    R1.paired.fastq.gz \
    R2.paired.fastq.gz \
    > ${sample_name}.bwa.sam


samtools view \
    -F 2308 \
    -b \
    ${sample_name}.bwa.sam \
    > ${sample_name}.bam

samtools sort \
    -@10 \
    ${sample_name}.bam \
    > ${sample_name}.sorted.bam

samtools index ${sample_name}.sorted.bam

export PATH="/home/k2so/.local/bin:$PATH"

bamCoverage \
    -b ${sample_name}.sorted.bam \
    -o ${output_bw_file}.bw
