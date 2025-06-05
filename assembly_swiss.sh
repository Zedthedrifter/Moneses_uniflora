#!/bin/bash
#SBATCH --job-name="Moneses uniflora"
#SBATCH --export=ALL
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G #ask for 128G memory

#also run as interactive job on slurm
#a piloting project. later if I have other jobs that can run in parallel I'll put them to queue

#===================================================================
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

#===================================================================
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

#===================================================================
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
contig_dir=Swiss_assembly_contigs

#spades.py --only-assembler -1 $r1 -2 $r2 -o Swiss_assembly -k 23,31 --meta -t 16 -m 32 #try smaller kmers, also requested 32GB memory when submitting queue, less thread
#spades.py --only-assembler -1 $r1 -2 $r2 -o Swiss_assembly_4k -k 31,55,65,75 --meta -t 16 -m 32 #try larger kmer, quest 32 GB memory. k31 finished. 
#spades.py --only-assembler -1 $r1 -2 $r2 -o careful_assembly_Swiss -k 31,55 --careful -t 16 -m 32 #try larger kmer, quest 32 GB memory with --careful this is actually worse than without --careful
#try interative method: 
#at least the assembly improves each time and builds upon previous assembly
#filter out long contigs from previous assembly to facilitate downstream assembly
#seqtk seq -L 1000 $contig_dir/contigs.3.fasta > $contig_dir/contigs.3.filtered.fasta
#spades.py --only-assembler -1 $r1 -2 $r2 -o careful_assembly_Swiss -k 23,65 --trusted-contigs $contig_dir/contigs.3.filtered.fasta --careful -t 16 -m 32 #continue from contig.3 --> contigs.4.fasta
#continue from contigs.4.fasta
#seqtk seq -L 1000 $contig_dir/contigs.4.fasta > $contig_dir/contigs.4.filtered.fasta
#spades.py --only-assembler -1 $r1 -2 $r2 -o careful_assembly_Swiss -k 69,75 --trusted-contigs $contig_dir/contigs.4.filtered.fasta --careful -t 32 -m 128
#continue from contigs.6.fasta
#seqtk seq -L 1000 $contig_dir/contigs.6.fasta > $contig_dir/contigs.6.filtered.fasta
#spades.py --only-assembler -1 $r1 -2 $r2 -o careful_assembly_Swiss -k 45,61 --trusted-contigs $contig_dir/contigs.6.filtered.fasta --careful -t 32 -m 128
#continue from contigs.7.fasta
seqtk seq -L 1000 $contig_dir/contigs.7.fasta > $contig_dir/contigs.7.filtered.fasta
spades.py --only-assembler -1 $r1 -2 $r2 -o careful_assembly_Swiss -k 25,35 --trusted-contigs $contig_dir/contigs.7.filtered.fasta --careful -t 32 -m 128
}

#try assembly with megahit
function megahit_assembly {

name=$1
outdir=$2
r1=$indir/${name}_1.trimmed.nc.fastq.gz 
r2=$indir/${name}_2.trimmed.nc.fastq.gz

megahit \
  -1 $r1 \
  -2 $r2 \
  -o $outdir \
  --k-list 21,41,61,81 \  
  --min-count 2 \          
  --prune-level 1 \        
  --prune-depth 1 \
  --no-mercy \             
  --memory 0.5 \           
  --low-local-ratio	0.1 \  
  -t 16       
/

# Reduced k-mer range for low coverage
# Lower threshold for rare k-mers
# Moderate pruning to retain fragmented contigs
# Disable "mercy" k-mers (improves sensitivity)
# Use 50% of available RAM (adjust if needed)
#Retains contigs with local low-coverage regions
echo 'done'
}
#===================================================================
#a de novo assembler specialized in plant genome assembly
function assembly_platanus {

name=$1
indir=$2
r1=$indir/${name}_1.trimmed.fastq.gz 
r2=$indir/${name}_2.trimmed.fastq.gz
prefix=$3

#god it's great now I've configured everything
platanus_allee assemble -f $r1 $r2 -o $prefix -t 24 -m 128 2>assemble.log
mv ${prefix}_contig.fa Swiss_assembly_contigs

##phasing -- this module won't work on the currently available platanus version if your genome is longer than ~600Mbp. Email me for the bug-free version (lcampos@rbge.org.uk)
#platanus_allee phase -o ${prefix}-phased -c ${prefix}_contig.fa -IP1 $r1 $r2 -mapper minimap2 -t 32  2>phase.log

#platanus consensus
#platanus_allee consensus -o ${prefix}-consensus -c ${prefix}-phased_allPhaseBlock.fa -IP1 $r1 $r2 -t 72 2>consensus.log

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
function init {

mkdir BUSCO_results Swiss_assembly_contigs inputs
}
#===========================================================================

function main {
accession="ERX7324453"
name="ERR7756292"
indir='/home/zchen/maternity_cover/moneses_uniflora_202505/inputs'

#init
#download data (run once)
#download_data $accession $indir #comment out after downloading once, if other analysis needs to be run again

#quality control
r1=$indir/${name}_1.fastq.gz
r2=$indir/${name}_2.fastq.gz
#final outputs
contigdir=/home/zchen/maternity_cover/moneses_uniflora_202505/Swiss_assembly_contigs
finalctg=contigs.8.fasta
buscodir=BUSCO_results
#quality_control $r1 $r2
#quality looks good

#remove/seperate plasmid reads
#remove_chlor $indir/Thale_cress_chloroplast_ref.fasta ${r1/fastq/trimmed.fastq} ${r2/fastq/trimmed.fastq} swiss
#remove_chlor $indir/Vaccinium_macrocarpon_chloroplast_ref.fasta ${r1/fastq/trimmed.fastq} ${r2/fastq/trimmed.fastq} swiss
#-----------------------------------------------------------
#de novo assembly 
#use nc reads only: spades, megahit --> compare N50 L50
#abyss_assembly $name #pass
#spades $name 
#megahit_assembly $name megahit_swiss #need to make sure the output directory is created freshly each time, otherwise, error
#assembly_platanus $name $indir mu1
#-----------------------------------------------------------
#assembly assessment
#assembly_evaluation Swiss_assembly contigs.fasta Swiss_assembly_quast #N50=1437
#assembly_evaluation Swiss_assembly_contigs mu1_contig.fa Swiss_assembly_quast 
#assembly_evaluation megahit_swiss final.contigs.fa megahit_swiss_quast #N50=1058
#assembly_evaluation Swiss_assembly contigs.fasta Swiss_assembly_quast #N50=760, total length 3970014, largest: 84731

#completeness assessment using BUSCO
busco -i $contigdir/$finalctg -o $buscodir -m genome -c 32 -l eudicotyledons_odb12 -r

#SSR detection
#MISA_SSR $contigdir $finalctg #files will be saved to the same directory as the input file


}

main
