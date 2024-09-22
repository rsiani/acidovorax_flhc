#!/bin/bash
#
#SBATCH --job-name=kallistoQuant
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=roberto.siani
#SBATCH --partition=cpu_p
#SBATCH --qos cpu_short
#SBATCH --time=05:00:00
#SBATCH --nice=10000

source $HOME/.bashrc
conda activate kallisto


echo "kallisto quantification"
mkdir kall17_$(date +"%d%m")
mkdir kall19_$(date +"%d%m")
mkdir kall21_$(date +"%d%m")
mkdir kall23_$(date +"%d%m")
mkdir kall25_$(date +"%d%m")
mkdir kall27_$(date +"%d%m")
mkdir kall31_$(date +"%d%m")
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_k17.idx -o kall17_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_k17.idx -o kall17_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_k19.idx -o kall19_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_k19.idx -o kall19_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_k21.idx -o kall21_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_k21.idx -o kall21_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_k23.idx -o kall23_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_k23.idx -o kall23_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_k25.idx -o kall25_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_k25.idx -o kall25_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_k27.idx -o kall27_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_k27.idx -o kall27_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_k31.idx -o kall31_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_k31.idx -o kall31_$(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 8 --plaintext ${i%%_L001_R1_001.fastq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
echo "All done!!!"

