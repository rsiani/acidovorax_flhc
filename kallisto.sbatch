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
mkdir kall29_$(date +"%d%m")
mkdir kall31_$(date +"%d%m")
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_17.idx -o kall17_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_17.idx -o kall17_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_19.idx -o kall19_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_19.idx -o kall19_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_21.idx -o kall21_$
(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fas
tq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001
.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_21.idx -o kall21_$
(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang  ${i%%_L001_R1_001.fa
stq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_00
1.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_23.idx -o kall23_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang  ${i%%_L001_R1_001.fas
tq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001
.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_23.idx -o kall23_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang  ${i%%_L001_R1_001.fas
tq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001
.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_25.idx -o kall25_$
(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fas
tq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001
.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_25.idx -o kall25_$
(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fas
tq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001
.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_27.idx -o kall27_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_27.idx -o kall27_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_29.idx -o kall29_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_29.idx -o kall29_$(
date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fast
q.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001.
fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR140.list | while read i; do kallisto quant -i ~/0723_extras/LR140_31.idx -o kall31_$
(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fas
tq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001
.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
cat ../LR124.list | while read i; do kallisto quant -i ~/0723_extras/LR124_31.idx -o kall31_$
(date +"%d%m")/${i%%_L001_R1_001.fastq.gz} -t 16 -b 99 --single-overhang ${i%%_L001_R1_001.fas
tq.gz}_L001_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L002_R*_001.fastq.gz ${i%%_L001_R1_001
.fastq.gz}_L003_R*_001.fastq.gz ${i%%_L001_R1_001.fastq.gz}_L004_R*_001.fastq.gz; done
echo "All done!!!"

