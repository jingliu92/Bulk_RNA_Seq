# RNA-Seq-Tutorial-Step-by-Step
## Install required packages
1. Create Env
```
conda create -n "RNA-seq"
conda activate RNA-seq
conda deactivate
```
2. Install fastq-dump for NCBI data downloading
```
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -vxzf sratoolkit.tar.gz
export PATH=$PWD/sratoolkit.3.2.1-ubuntu64/bin:$PATH
which fastq-dump
fastq-dump --outdir ./ SRR33981241 (1:50PM started)
```
3. Install FastQC for quality control
```
conda install bioconda::fastqc
```
4. 

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
#SBATCH --mail-type=All

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

for i in {1..4};do
{
samtools view -S -b OC${i}_L.sam > ../bam_files/OC${i}_L.bam
samtools view -S -b OC${i}_B.sam > ../bam_files/OC${i}_B.bam
samtools view -S -b OG${i}_L.sam > ../bam_files/OG${i}_L.bam
samtools view -S -b OG${i}_B.sam > ../bam_files/OG${i}_B.bam
};
done

for i in {1..4};do
{
samtools sort -@ 12 -o ../bam_files/OC${i}_L.sorted.bam ../bam_files/OC${i}_L.bam
samtools sort -@ 12 -o ../bam_files/OC${i}_B.sorted.bam ../bam_files/OC${i}_B.bam
samtools sort -@ 12 -o ../bam_files/OG${i}_L.sorted.bam ../bam_files/OG${i}_L.bam
samtools sort -@ 12 -o ../bam_files/OG${i}_B.sorted.bam ../bam_files/$OG{i}_B.bam
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

for i in {1..4};do
{
stringtie -o OC${i}_L.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OC${i}_L.sorted.bam
stringtie -o OC${i}_B.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OC${i}_B.sorted.bam
stringtie -o OG${i}_L.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OG${i}_L.sorted.bam
stringtie -o OG${i}_B.gtf -p 16 -G ../../mus_reference/genomic.gff --rf ../bam_files/OG${i}_B.sorted.bam
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
stringtie --merge -p 12 -G ../../mus_reference/genomic.gff -o stringtie_merged.gtf gtf_list.txt

for i in {1..4};do
{
stringtie -e -B -p 16 -G stringtie_merged.gtf -o ../final_gtf/OC${i}_L/OC${i}_L.gtf --rf ../bam_files/OC${i}_L.sorted.bam 
stringtie -e -B -p 16 -G stringtie_merged.gtf -o ../final_gtf/OC${i}_B/OC${i}_B.gtf --rf ../bam_files/OC${i}_B.sorted.bam
stringtie -e -B -p 16 -G stringtie_merged.gtf -o ../final_gtf/OG${i}_L/OG${i}_L.gtf --rf ../bam_files/OG${i}_L.sorted.bam
stringtie -e -B -p 16 -G stringtie_merged.gtf -o ../final_gtf/OG${i}_B/OG${i}_B.gtf --rf ../bam_files/OG${i}_B.sorted.bam
};
done

```
## Generate gene and transcript count tables
- Use `prepDE.py` to generate a count matrix for genes and transcripts

  Here I used a python script that downloaded online (https://ccb.jhu.edu/software/stringtie/dl/prepDE.py) to generate count tables for DE analysis
- Make a sample list that include sampleID and full path
  
```
OC1_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC1_L/OC1_L.gtf
OC2_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC2_L/OC2_L.gtf
OC3_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC3_L/OC3_L.gtf
OC4_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC4_L/OC4_L.gtf
OG1_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG1_L/OG1_L.gtf
OG2_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG2_L/OG2_L.gtf
OG3_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG3_L/OG3_L.gtf
OG4_L /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG4_L/OG4_L.gtf
OC1_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC1_B/OC1_B.gtf
OC2_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC2_B/OC2_B.gtf
OC3_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC3_B/OC3_B.gtf
OC4_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OC4_B/OC4_B.gtf
OG1_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG1_B/OG1_B.gtf
OG2_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG2_B/OG2_B.gtf
OG3_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG3_B/OG3_B.gtf
OG4_B /scratch/jingliu/FRG_RNAseq/rawdata/final_gtf/OG4_B/OG4_B.gtf
```
```
#!/bin/bash
#SBATCH -p batch
#SBATCH -t 120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mail-user=jing.liu12@okstate.edu
#SBATCH --mail-type=end

module load python

cd rawdata/final_gtf/

python PrepDE.py -i sample_list.txt
```
# Got the gene count table for DE analysis


# GO term enrichment
## Among our list of DE genes, which GO terms are enriched?

Transcripts that are longer or more highly expressed give more statistical power for detecting differential expression between samples

The same bias holds true for GO categories: categories with predominantly highly expressed or long genes are more likely to be found to be over-represented within the DEG.

`GOseq`:

1. determine DEG
2. quantify likelihood of DE as a function of gene length (–> weight)
3. statistical test of each GO category’s significance taking the DE probability into account
Manual work-around for mouse:

```
#biocLite("org.Mm.eg.db")
#biocLite("goseq")
#biocLite("biomaRt")
library(biomaRt)
library(org.Mm.eg.db)
library(goseq)
library(geneLenDataBase)

```
## Determine DEG
```{r}
# edger R analysis
# Import gene expression maxtrix and remove duplicate rows
```{r}
rm(list=ls()) #Clears the existing environment before we start downstream analysis, if there is one. 

# Import all sample gene count matrix remove duplicate rows.
gene_all<- read.csv('RNA-Seq/gene_count_matrix.csv',  sep = ',',check.names = FALSE,header = T)

# Find duplicate rows: 
gene_all$gene_id[duplicated(gene_all$gene_id)]

# remove duplicate IDS (here we use "SYMBOL", but it should be whatever was selected as keyType)
gene_all = gene_all[!gene_all(ids[c("SYMBOL")]),]

# Remove duplicate rows (delete lower gene count rows)
```

```{r}
# Liver_NCD gene expression
```{r}

rm(list=ls())
# Liver_NCD gene expression
liver_NCD <- read.csv('RNA-Seq/Raw_data/liver_NCD.csv',  sep = ',',row.names = 1,check.names = FALSE,header = T)

# Remove rowSums ==0
liver_NCD <- liver_NCD %>% 
  filter(if_any(where(is.numeric),~.x >0))

# Specify the grouping, make sure that the order of the samples in the expression matrix and the order of the grouping are the same, with control being the first
group_liver_NCD <- rep(c('NCD', 'NCD+CUR'), times= c(4,4))
```

# Creating a DGEList object and filter low count genes
```{r}
dgelist <- DGEList(counts = liver_NCD, remove.zeros = TRUE,group = group_liver_NCD)

#filter low count genes by using CPM (recommended)
keep1 <- rowSums(cpm(dgelist) > 1 ) >= 2
summary(keep1)

# We can also do this using the {edgeR} command filterByExpr().
keep <- filterByExpr(dgelist) 
summary(keep)

dgelist <- dgelist[keep1, , keep.lib.sizes = FALSE]
```

## Normalization
### We now need to normalise the library sizes between the samples in DGEList. This is done to minimise bias towards highly expressed genes. The {edgeR} function calcNormFactors() normalises the library sizes by finding a set of scaling factors for the library sizes that minimizes the log-fold changes between the samples for most genes. This function can be assigned directly to the DGEList object. By printing the normalisation factors for the samples (dgelist$samples$norm.factors) before and after assigning calcNormFactors() to the DGEList object, we can see the effect of this command.

```{r}
dgelist$samples$norm.factors
dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')
dgelist_norm <- calcNormFactors(dgelist)
dgelist_norm$samples$norm.factors

# Export normalized gene counts
write.table(dgelist_norm$counts, 'RNA-Seq/Processed_data/liver_NCD_norm.csv', sep = ',',col.names = NA, quote = FALSE)
```

## Differential analysis
### First approach: exact test between the two groups
```{r}
# Estimate dispersion
dge.exact <- estimateCommonDisp(dgelist_norm)
dge.exact
plotBCV(dge.exact)

# Perform the test
dgeExactTest <- exactTest(dge.exact)
dgeExactTest

# p values are corrected with the function topTags :
resExactTest <- topTags(dgeExactTest, n = nrow(dgeExactTest$table))
head(resExactTest$table)

#p-value and (BH) adjusted p-value distribution can be assessed with:
par(mfrow = c(1,2))
hist(resExactTest$table$PValue, xlab = "p-value", main = "raw p-values")
hist(resExactTest$table$FDR, xlab = "p-value", main = "adjusted p-values")

# And finally, genes with a FDR smaller than 5% and a log Fold Change larger than 1 or smaller than -1 are extracted:

selectedET <- resExactTest$table$FDR < 0.05 & abs(resExactTest$table$logFC) > 2
selectedET <- resExactTest$table[selectedET, ]
nrow(selectedET)
```
