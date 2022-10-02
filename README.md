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
