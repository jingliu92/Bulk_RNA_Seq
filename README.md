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

1. Convert SAM file to BAM file
2. Sort Bam file
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

for i in {1..4};do samtools
{
 view -S -b OC${i}_L.sam > ../bam_files/OC${i}_L.bam
 view -S -b OC${i}_B.sam > ../bam_files/OC${i}_B.bam
 view -S -b OG${i}_L.sam > ../bam_files/OG${i}_L.bam
 view -S -b OG${i}_B.sam > ../bam_files/OG${i}_B.bam
};
done

for i in {1..4};do
{
sort -@ 12 -o ../bam_files/OC${i}_L.sorted.bam ../bam_files/OC${i}_L.bam
sort -@ 12 -o ../bam_files/OC${i}_B.sorted.bam ../bam_files/OC${i}_B.bam
sort -@ 12 -o ../bam_files/OG${i}_L.sorted.bam ../bam_files/OG${i}_L.bam
sort -@ 12 -o ../bam_files/OG${i}_B.sorted.bam ../bam_files/$OG{i}_B.bam
};
done

```
# Use stringtie for transcriptome assembly
## Assemble and merge transcripts. Create gtf file for each sample

```
#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

module load stringtie

cd rawdata/stringtie_gtf/

for i in {1..4};do stringtie
{
 -o OC${i}_L.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OC${i}_L.sorted.bam
 -o OC${i}_B.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OC${i}_B.sorted.bam
 -o OG${i}_L.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OG${i}_L.sorted.bam
 -o OG${i}_B.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OG${i}_B.sorted.bam
};
done

```

### Merge gtf files
- Create a `gtf_list.txt` which include the full path of the .gtf file that generated from above step
```
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC1_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC2_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC3_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC4_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG1_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG2_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG3_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG4_L.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC1_B.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC2_B.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC3_B.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OC4_B.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG1_B.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG2_B.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG3_B.gtf
/scratch/jingliu/FRG_RNAseq/rawdata/stringtie_gtf/OG4_B.gtf
```



stringtie --merge \
-p 12 \
-G ../../mus_reference/genomic.gff \
-o stringtie_merged.gtf gtf_list.txt 

stringtie -e -B -p 16 \
-G stringtie_merged.gtf \
-o ../final_gtf/OC${i}_L/OC${i}_L.gtf \
--rf ../bam_files/OC${i}_L.sorted.bam
};
done

