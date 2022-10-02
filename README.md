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

mkdir trimmed_data

#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=ALL

cd ../FRG_RNAseq/rawdata

module load solexaqa

SolexaQA++ dynamictrim forward.fq.gz reverse.fq.gz -h 20 -d ../trimmed_data
SolexaQA++ lengthsort /trimmed_data/forward.fq.trimmed.gz /trimmed_data/reverse.fq.trimmed.gz \
-l 60 -d ../rawdata/length_sorted

for file in *.fq.gz; do
{
SolexaQA++ dynamictrim ${file}__1.fq.gz ${file}__1.fq.gz

-x /home/usrname/DGE_Virtual_Oct2020/human_reference/hisat2_index/GRCh38.p12.genome \
-U ${file} \
--threads 4 \
-S /home/usrname/DGE_Virtual_Oct2020/hisat2_alignments_subset/${file}.sam; done

for i in {94..95};do
{
#SolexaQA++ dynamictrim HC03028${i}-HC03028${i}_combined_R1.fastq.gz HC03028${i}-HC03028${i}_combined_R2.fastq.gz -h 20 -d /projects/dsn001/Rivera_HIC_Collaboration_human/198.2.192.40:2129/CGTPJLL220126-RNA/trimmed_data/
