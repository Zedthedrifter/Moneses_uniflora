#!/usr/bin/env python3

from os import system
import subprocess as sbp
import argparse
from Bio import SeqIO
from Bio import Align
from Bio.SeqRecord import SeqRecord

###############################################
#standard input: ./extract_flanked.py xx.filtered.fasta xx.filtered.fasta.misa prefix desired_flank_size
###############################################

def read_SSR_to_dict(misa):
  with open(misa) as handle:
    keys=[x for x in next(handle).strip('\n').split('\t') if x !='']
    misa_dict={f'SSR_{i}': {i:j for i,j in zip(keys,[x for x in line.strip('\n').split('\t') if x !=''])}for i,line in enumerate(handle)}
    misa_dict={k:{i:int(j) if i in ['size', 'start', 'end'] else j for i,j in v.items() } for k,v in misa_dict.items()} #make sure int entries are actually numbers
  return(misa_dict)  

def extract_contigs(contigs,misa_dict):
    hits=set([v['ID'] for v in misa_dict.values()])
    print(f'Number of SSR containing sequences: {len(hits)}')
    seqs={k:v for k,v in SeqIO.to_dict(SeqIO.parse(contigs,'fasta')).items() if k in hits} #extract only the hit contigs
    print(f'Number of SSR containing sequences extracted: {len(seqs)}')
    return (seqs)

def extract_SSR(misa_dict,seqs,flank):
  count=0
  ssrs={}
  for k,v in misa_dict.items():
    if v['start']<flank:
      print(f'Discard {k} due to no enough space for upstream flank: {v["start"]} < {flank}')
    elif v['end']+flank > len(seqs[v['ID']].seq):
      print(f'Discard {k} due to no enough space for downstream flank: {v["end"]+flank}>{len(seqs[v["ID"]].seq)}')
    else:
      count+=1
      ssrs[k]=v
      ssrs[k]['seq']=seqs[v['ID']].seq[v['start']-flank:v["end"]+flank]
  print(f"{count} SSRs have over {flank} bp flanking regions and are selected")
  return (ssrs)

def write_SSR(ssrs,prefix,flank):
  records=[SeqRecord(v['seq'],id=k,description=v['SSR'],name=f"{v['ID']}_{v['start']}_{v['end']}") for k,v in ssrs.items()]
  SeqIO.write(records,f"{prefix}_SSR_flk_{flank}.fasta",'fasta')

def main(contigs,misa,prefix,flank=150):
  misa_dict=read_SSR_to_dict(misa)
  seqs=extract_contigs(contigs,misa_dict)
  ssrs=extract_SSR(misa_dict,seqs,int(flank))
  write_SSR(ssrs,prefix,flank)
  #print(misa_dict)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='run kDNA annotation')
  parser.add_argument('contigs', help='Name of the genome assembly, filtered.fasta', metavar='contigs')
  parser.add_argument('misa', help='Name of the misa output', metavar='misa')
  parser.add_argument('prefix', help='Prefix of the SSR output', metavar='prefix')
  parser.add_argument('flank', help='Required length of flanking region of the SSR', metavar='flank')
  options = parser.parse_args()
  
  main(options.contigs, options.misa,options.prefix, options.flank)
