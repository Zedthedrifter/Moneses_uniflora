#!/bin/bash
#SBATCH --job-name="Moneses uniflora"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G #ask for 32G memory

#also run as interactive job on slurm
#a piloting project. later if I have other jobs that can run in parallel I'll put them to queue

function download_data {
accession=$1
indir=$2
#get FTP paths
r1=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp" | awk -F'\t' '{print $2}'|awk -F';' '{print $1}')
r1=${r1/'fastq_ftp'/''}#potential bug
r2=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${accession}&result=read_run&fields=fastq_ftp" | awk -F'\t' '{print $2}'|awk -F';' '{print $2}')
echo 'downloading read files from' $r1 $r2
#download with wget
wget -P $indir $r1 #could probably use -o to rename the files and reduce the number of variables? perhaps next time
wget -P $indir $r2
}

function quality_control {
r1=$1
r2=$2

#fastqc output checking
fastqc $r1 $r2 -o QC
#no adapter detected
#do we need to trim anything?

#trim just in case
function trim {
adapters=/home/zchen/apps/conda/share/trimmomatic-0.39-2/adapters/TruSeq3-PE.fa #need to specify the full path to the adapter files
trimmomatic PE $r1 $r2 ${r1/fastq/trimmed.fastq} ${r1/fastq/unpaired.fastq} ${r2/fastq/trimmed.fastq} ${r2/fastq/unpaired.fastq} \
ILLUMINACLIP:${adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
/
}
trim #run once is enough
#check QC again after trimming 
fastqc ${r1/fastq/trimmed.fastq} ${r2/fastq/trimmed.fastq} -o QC
}

#remove the plastmid reads
function remove_chlor {
#reference genome: NC_000932.1 Thale cress chloroplast reference genome
#update: use Vaccinium macrocarpon reference genome. a common reference for Ericaceae
ref=$1
r1=$2
r2=$3
strain=$4
#alignment with bowtie2
bowtie2-build -f $ref chlo_ref
bowtie2-align-s --wrapper basic-0 -x chlo_ref -p 20 --threads 16 --phred33 --very-sensitive --quiet --time --rg-id ${strain} --rg SM:${strain} --rg PL:'ILLUMINA' -S reads.bam -1 ${r1} -2 ${r2}
#split mapped vs unmapped reads; output the bam files into directory 'mapping'
samtools view -b -f 4 -o mapping/${name}_unmapped.reads.bam reads.bam #non plasmid reads
samtools view -b -F 4 -o mapping/${name}_mapped.reads.bam reads.bam #plasmid reads
#output as paired reads
samtools fastq -1 ${r1/trimmed.fastq/trimmed.nc.fastq} -2 ${r2/trimmed.fastq/trimmed.nc.fastq} -0 /dev/null -s /dev/null -n mapping/${name}_unmapped.reads.bam
samtools fastq -1 ${r1/trimmed.fastq/trimmed.pl.fastq} -2 ${r2/trimmed.fastq/trimmed.pl.fastq} -0 /dev/null -s /dev/null -n mapping/${name}_mapped.reads.bam
#can do assembly for chloroplast and nuclear genome separately
#when using Thale cress chloroplast reference genome
#before: 11558738 pairs
#after:  11497885 pairs
#not reducing too much apparently
#when using cranberry reference genome: 
#after:  11464275 pairs, not reduced too much...
}

#===================================================================
#de novo assembly of nuclear genome

#abyss
function abyss_assembly { #memory heavy? killed over and over
name=$1
r1=$indir/${name}_1.trimmed.nc.fastq.gz 
r2=$indir/${name}_2.trimmed.nc.fastq.gz

#abyss-pe k=63 name=$name in='${r1} ${r2}' l=50 n=10 #the input format does not support regular expression
#type the file names directly
abyss-pe k=63 name=$name in='/home/zchen/maternity_cover/moneses_uniflora_202505/inputs/ERR7756292_1.trimmed.nc.fastq.gz /home/zchen/maternity_cover/moneses_uniflora_202505/inputs/ERR7756292_2.trimmed.nc.fastq.gz' l=50 n=10 j=16
echo 'done'
}

#Use if coverage >5X. --only-assembler Mode 
function spades {

name=$1
r1=$indir/${name}_1.trimmed.nc.fastq.gz 
r2=$indir/${name}_2.trimmed.nc.fastq.gz
outdir=/home/zchen/maternity_cover/moneses_uniflora_202505/careful_assembly_Norway
contigdir=/home/zchen/maternity_cover/moneses_uniflora_202505/Norway_assembly_contigs

#spades.py --only-assembler -1 $r1 -2 $r2 -o $outdir -k 31,55 --careful -t 32 -m 128 #try larger kmer, quest 32 GB memory with --careful
#mv $outdir/contigs.fasta $contigdir/contigs.1.fasta
#
function contig_iter {

  trust_ctg=$1
  k1=$2
  k2=$3
  nb=$4
  
  seqtk seq -L 1000 $contigdir/$trust_ctg > $contigdir/${trust_ctg/fasta/filtered.fasta}
  spades.py --only-assembler -1 $r1 -2 $r2 -o $outdir -k $k1,$k2 --trusted-contigs $contigdir/${trust_ctg/fasta/filtered.fasta} --careful -t 32 -m 128
  assembly_evaluation $outdir contigs.fasta Norway_assembly_quast
  mv Norway_assembly_quast/report.html Norway_assembly_quast/report.${nb}.html 
  mv $outdir/contigs.fasta $contigdir/contigs.${nb}.fasta
  }  

#contig_iter contigs.1.fasta 23 65 2
#contig_iter contigs.2.fasta 69 75 3
contig_iter contigs.3.fasta 45 61 4
contig_iter contigs.4.fasta 25 35 5
contig_iter contigs.5.fasta 39 49 6
contig_iter contigs.6.fasta 29 59 7
}

#===========================================================================

#de novo assembly evaluation
function assembly_evaluation {

indir=$1
contigs=$2
outdir=$3

quast.py $indir/$contigs -o $outdir  
}


#===========================================================================

#SSR detection
function MISA_SSR {

indir=$1
contigs=$2
#remove reads <500 bp
seqtk seq -L 500 $indir/$contigs > $indir/${contigs/fasta/filtered.fasta}
#MISA is installed by source code directly
misa.pl $indir/${contigs/fasta/filtered.fasta} 
}
#===========================================================================

function main {
accession="ERX5335444" #the Norway sample
name="ERR5555402"
indir='/home/zchen/maternity_cover/moneses_uniflora_202505/inputs'
#download data (run once)
#download_data $accession $indir #comment out after downloading once, if other analysis needs to be run again

#quality control
r1=inputs/${name}_1.fastq.gz #has to write the accesstion number 
r2=inputs/${name}_2.fastq.gz
#quality_control $r1 $r2
#quality looks good

#remove/seperate plasmid reads
#remove_chlor $indir/Vaccinium_macrocarpon_chloroplast_ref.fasta ${r1/fastq/trimmed.fastq} ${r2/fastq/trimmed.fastq} norway


#de novo assembly 
#use nc reads only
#spades $name #include assessment into the assembly steps

#assembly assessment
#assembly_evaluation Norway_assembly_contigs contigs.2.fasta Norway_assembly_quast

#SSR detection
MISA_SSR /home/zchen/maternity_cover/moneses_uniflora_202505/Norway_assembly_contigs/ contigs.7.fasta 

}

main
