#!/usr/bin/env python3

from Bio import SeqIO
import os, argparse, sys, re, csv
from os import listdir
from os.path import isfile, join

parser = argparse.ArgumentParser(
                    prog='calc_polyQ',
                    description='PolyQ distribution',
                    epilog='calc some polyQ stats in proteins')
parser.add_argument('-i','--indir',)
parser.add_argument('-o', '--outdir')
parser.add_argument('-v', '--verbose',
                    action='store_true')
args = parser.parse_args()

with open(os.path.join(args.outdir,f'summary.tsv'),'w') as sumout:
    sumwriter = csv.writer(sumout, delimiter='\t')
    sumwriter.writerow(['Species','TotalProteins','polyQ10_count','polyQ10_freq'])
    for file in listdir(args.indir):
        if file.endswith('.fa') or file.endswith('fasta'):
            (name) = file
            m = re.search('(\S+)\.pep\.fa',file)
            if m:
                name = m.group(1)
            else:
                name = file.split("_")[0]
            if args.verbose:
                print(name)            
            with open(os.path.join(args.outdir,f'{name}.polyQ.tsv'),'w') as spout:
                polyQwriter = csv.writer(spout, delimiter='\t')
                polyQwriter.writerow(['Protein','PolyQlength'])
                pepcount = 0
                polyQcount = 0
                for record in SeqIO.parse(os.path.join(args.indir,file), "fasta"):
                    seqstr = str(record.seq)
                    qLens = []
                    pepcount += 1
                    for m in re.finditer('(Q{10,})',seqstr):
                        qLens.append(len(m.group(1)))
                    qLens.sort(reverse=True) # get longest for individual gene
                    if len(qLens):
                        polyQwriter.writerow([record.id,qLens[0]])
                        polyQcount += 1
                if pepcount > 0:
                    sumwriter.writerow([name,pepcount,polyQcount,
                                    f'{(polyQcount/pepcount):.4f}'])