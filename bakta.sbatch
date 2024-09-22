#!/bin/bash

#SBATCH -J emapper
#SBATCH -p normal_q
#SBATCH -q normal
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 01-00
#SBATCH --mail-user=roberto.siani
#SBATCH --mail-type=ALL

source $HOME/.bashrc
conda activate bakta

echo "Started bakta annotation"

bakta LR140.fa --db ~/db --output bakta_results40/  --prefix LR140 --genus Acidovorax --species delafieldii --strain LjRoot140 --gram "?" --locus 'LR140'  -t 8 
bakta LR124.fa --db ~/db --output bakta_results24  --prefix LR124 --genus Acidovorax --species delafieldii --strain LjRoot124 --gram "?"  --locus 'LR140' -t 8

