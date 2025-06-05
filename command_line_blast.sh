#!/bin/bash
#SBATCH --job-name="blastn"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G 
database=$1 
query=$2 
format=$3

makeblastdb -in $database -dbtype nucl -out tmp
blastn -db tmp -query $query -outfmt $format > blast_out.txt
