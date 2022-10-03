# RNA-Seq-Tutorial-Step-by-Step
## FastQC
```
cd /scratch/jingliu/RNA-Seq
mkdir fastqc

# Run Fastqc on all the files using a wildcard (*)
###### fastqc.sbatch #######
#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

module load fastqc

fastqc -o fastqc/ -f fastq ../FRG_RNAseq/rawdata/*.fq.gz

# -o -> Create all output files in the specified output directory
# -f -> Bypasses the normal sequence file format detection and forces the program to use the specified format.Valid formats are bam,sam,bam_mapped,sam_mapped and fastq
# *.fq.gz -> Will run FASTQC for all files in the current folder that ends with .fq.gz
```

## SolexaQA++

```
mkdir trimmed_data

#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=ALL

cd rawdata

module load solexaqa

for file in {1..4}; do
{
SolexaQA++ dynamictrim OC${file}_L_1.fq.gz OC${file}_L_2.fq.gz \
-h 20 \
-d trimmed_data/
};
done

for file in {1..4}; do
{
SolexaQA++ dynamictrim OG${file}_B_1.fq.gz OG${file}_B_2.fq.gz \
-h 20 \
-d trimmed_data/
};
done

for file in {1..4}; do
{
SolexaQA++ dynamictrim OG${file}_L_1.fq.gz OG${file}_L_2.fq.gz \
-h 20 \
-d trimmed_data/
};
done

for file in {1..4}; do
{
SolexaQA++ dynamictrim OG${file}_B_1.fq.gz OG${file}_B_2.fq.gz \
-h 20 \
-d trimmed_data/
};
done

```

### length sorted data
```
mkdir length_sorted

#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

cd rawdata

module load solexaqa

for file in {1..4}; do
{
SolexaQA++ lengthsort trimmed_data/OC${file}_L_1.fq.trimmed.gz trimmed_data/OC${file}_L_2.fq.trimmed.gz \
-l 60 \
-d length_sorted/
};
done

for file in {1..4}; do
{
SolexaQA++ lengthsort trimmed_data/OC${file}_B_1.fq.trimmed.gz trimmed_data/OC${file}_B_2.fq.trimmed.gz \
-l 60 \
-d length_sorted/
};
done

for file in {1..4}; do
{
SolexaQA++ lengthsort trimmed_data/OG${file}_B_1.fq.trimmed.gz trimmed_data/OG${file}_B_2.fq.trimmed.gz \
-l 60 \
-d length_sorted/
};
done

for file in {1..4}; do
{
SolexaQA++ lengthsort trimmed_data/OG${file}_L_1.fq.trimmed.gz trimmed_data/OG${file}_L_2.fq.trimmed.gz \
-l 60 \
-d length_sorted/
};
done
```
# Build a genome index/database 
```
cd mus_reference/
mkdir hisat2_index

cd mus_reference/
mkdir hisat2_index

# hisat2.sbatch
#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

module load hisat2
hisat2-build GCF_000001635.27_GRCm39_genomic.fna hisat2_index/GRCm39_genomic
```
# Extract splicesites from gtf file
Here we use `hisat2_extract_splice_sites.py` to extract splicesites from Gff file https://github.com/DaehwanKimLab/hisat2/blob/master/hisat2_extract_splice_sites.py

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

module load hisat2
python hisat2_extract_splice_sites.py /mus_reference/genomic.gff > splicesites.txt
```
# Run Alighments
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

module load hisat2
cd rawdata/length_sorted

for i in {1..4};do
{
hisat2 -p 12 --mp 6,6 --score-min L,0,-0.2 --dta-cufflinks --rna-strandness RF \
-x ../../mus_reference/hisat2_index/GRCm39_genomic \
-1 OG${i}_B_1.fq.trimmed.paired.gz \
-2 OG${i}_B_2.fq.trimmed.paired.gz \
-S ../sam_files/OG${i}_B.sam
};
done

for i in {1..4};do
{
hisat2 -p 12 --mp 6,6 --score-min L,0,-0.2 --dta-cufflinks --rna-strandness RF \
-x ../../mus_reference/hisat2_index/GRCm39_genomic \
-1 OG${i}_L_1.fq.trimmed.paired.gz \
-2 OG${i}_L_2.fq.trimmed.paired.gz \
-S ../sam_files/OG${i}_L.sam
};
done

for i in {1..4};do
{
hisat2 -p 12 --mp 6,6 --score-min L,0,-0.2 --dta-cufflinks --rna-strandness RF \
-x ../../mus_reference/hisat2_index/GRCm39_genomic \
-1 OC${i}_B_1.fq.trimmed.paired.gz \
-2 OC${i}_B_2.fq.trimmed.paired.gz \
-S ../sam_files/OC${i}_B.sam
};
done

for i in {1..4};do
{
hisat2 -p 12 --mp 6,6 --score-min L,0,-0.2 --dta-cufflinks --rna-strandness RF \
-x ../../mus_reference/hisat2_index/GRCm39_genomic \
-1 OC${i}_L_1.fq.trimmed.paired.gz \
-2 OC${i}_L_2.fq.trimmed.paired.gz \
-S ../sam_files/OC${i}_L.sam
};
done

```

# Use Samtools to convert to coordinate sorted bam files

## Fisrt convert SAM file to BAM file
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

module load samtools
cd rawdata/sam_files/
for i in {1..4};do
{
samtools view -S -b OC{i}.sam > ../bam_files/OC{i}.bam
samtools sort -@ 12 -o /projects/dsn001/Rivera_HIC_Collaboration_human/198.2.192.40:2129/CGTPJLL220126-RNA/Bam_files_stranded/HC03028${i}.sorted.bam /projects/dsn001/Rivera_HIC_Collaboration_human/198.2.192.40:2129/CGTPJLL220126-RNA/Bam_files_stranded/HC03028${i}.bam
};
done
```



