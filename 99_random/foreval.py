#!/usr/bin/env python
'''
Creates joint alignments of forward and reverse reads
'''

import subprocess
import argparse
import glob
import code

try:
  from Bio import SeqIO
except ImportError:
  print 'Please install biopython first'
  exit(1)
import math

# Function for changing the names iof input files
def replname(file_name, suffix):
  base = file_name.split('.')[0]
  return base+suffix
def delete(delete_list):
  for x in delete_list:
    subprocess.check_call("rm "+x, shell = True)
# Align forward and reverse reads

class Searcher:
  def hmmsearch(self, hmm, read_files):
    suff = ['for', 'rev']
    counter = 0
    for read_file in read_files:
      splt = read_file.split('/')
      out_dir = "./"+splt[len(splt)-1].split('.')[0]+"_"+suff[counter]+".txt"
      if read_file.endswith(('.fq.gz','fastq.gz')):
        subprocess.check_call(["/bin/bash", "-c", " hmmsearch -E "+ args.evalue +" --domtblout "+out_dir+" "+hmm+" <(getorf -sequence <(awk '{print \">\" substr($0,2);getline;print;getline;getline}' <(zcat "+ read_file +" | sed 's/:/$/g')) -outseq >(cat) -minsize 98 2>/dev/null) 2>&1 > /dev/null"])
      elif read_file.endswith(('fa', 'fasta', 'fna', 'faa')):
        subprocess.check_call(["/bin/bash", "-c", " hmmsearch -E "+ args.evalue +" --domtblout "+out_dir+" "+hmm+" <(getorf -sequence <(sed 's/:/$/g' " +read_file+ ") -outseq >(cat) -minsize 98 2>/dev/null) 2>&1 > /dev/null"])
      counter += 1
  def fxtract(self, sequence_titles, sequences_list):
    counter = 0
    suff = ['_for', '_rev']
    for sequences in sequences_list:
      output_file_name = replname(sequences.split('/')[len(sequences.split('/')) -1], suff[counter]+'.faa')
      subprocess.check_call("fxtract -H -f " +sequence_titles+ " " +sequences+ " > "+output_file_name, shell=True)
      counter += 1
  def HMM_hit_extract(self, table_list, read_list):
    read_counter = 1
    counter = 0
    orf_for_names = []
    orf_rev_names = []
    for_names = []
    rev_names = []
    read_names = []
    for_sequences = []
    rev_sequences = []
    parsed_fasta_file = []
    for line in open(table_list[0], 'r'):
      if line.startswith('#'):
        continue
      else:
        orf_for_names.append(line.split(' ')[0])
        for_names.append(line.split(' ')[0].replace('$', ':')[:-2])
    for line in open(table_list[1], 'r'):
      if line.startswith('#'):
        continue
      else:
        orf_rev_names.append(line.split(' ')[0])
        rev_names.append(line.split(' ')[0].replace('$', ':')[:-2])
    r = set(for_names)
    for name in rev_names:
      if name in r:
        read_names.append(name)
    output_file_name = replname(read_list[0].split('/')[len(read_list[0].split('/')) - 1], '_names.txt')
    out_name_file = open(output_file_name, 'w')
    for name in read_names:
      out_name_file.write(name+'\n')
    out_name_file.close()

    counter += 1
    read_counter += 1
    orf_n = [orf_for_names, orf_rev_names]
    return orf_n

  def orf_extractor(self, orf_read_files, extract_names):
    counter = 0
    for read_file in orf_read_files:
      hits = extract_names[counter]
      write_list = []
      out_dir = "./"+read_file.split('.')[0]+".orf"
      subprocess.check_call(["/bin/bash", "-c", "getorf -sequence <(sed 's/:/$/g' "+read_file+") -outseq "+out_dir+" -minsize 98 2>/dev/null"])
      parsed_fasta = SeqIO.parse(open(out_dir), "fasta")
      for fasta in parsed_fasta:
        name, sequence = fasta.id, fasta.seq.tostring()
        if name in set(hits):
          write_list.append(fasta)
      out_dir_new = out_dir.replace('.orf', '_hits.orf')
      open_outdir = open(out_dir_new, 'w')
      for fasta in write_list:
        SeqIO.write(fasta, open_outdir, "fasta")
      counter += 1

  def fasta_parser(self, to_splits, hmm):
    sequence_list = []
    for to_split in to_splits:
      out_dir = to_split.replace('.orf', '.sto')
      subprocess.check_call('hmmalign --trim -o '+out_dir+' ' +hmm+ ' ' +to_split, shell = True)
      subprocess.check_call('seqmagick convert '+out_dir+' '+out_dir.replace('.sto', '.fa'), shell = True)
      sequences = list(SeqIO.parse(open(out_dir.replace('.sto', '.fa')), 'fasta'))
      rm_pos = []
      for fasta in sequences:
        pos = 0
        for n in list(fasta.seq):
          if n.islower():
            if pos not in rm_pos:
              rm_pos.append(pos)
          pos += 1
        continue
      rm_pos = sorted(rm_pos, reverse = True)
      for fasta in sequences:
        s = list(fasta.seq)
        for p in rm_pos:
          del s[p]
        fasta.seq = ''.join(s)
      sequence_list.append(sequences)
    for f in list(sequence_list[0]):
      for r in list(sequence_list[1]):
        p_seq = []
        if f.id[:-2] == r.id[:-2]:
          for i, j in zip(list(f.seq), list(r.seq)):
            if i == j:
              p_seq.append(i)
            elif i.isalpha() and i != j:
              p_seq.append(i)
            elif j.isalpha() and j != i:
              p_seq.append(j)
            elif i.isalpha() and j.isalpha() and i != j:
              p_seq.append('.')
        else:
          continue
        startend = []
        ind = 0
        for n in list(f.seq[::-1]):
          if n.isalpha():
            startend.append(len(f.seq)-ind)
            break
          else:
            ind += 1
        ind = 0
        for n in list(r.seq):
          if n.isalpha():
            startend.append(ind)
            break
          else:
            ind += 1
        print '>'+f.id[:-2]
        out = [''.join(p_seq[:startend[0]]),''.join(p_seq[startend[0]:startend[1]]).replace('-','.'),''.join(p_seq[startend[1]:])]
        print ''.join(out)




  def pynast(self, read_files, gg_db):
    suff = ['for', 'rev']
    sequence_list = []
    counter = 0
    for read_file in read_files:
      splt = read_file.split('/')
      out_dir = "./"+splt[len(splt)-1].split('.')[0]+"_"+suff[counter]+".aln.fasta"
      subprocess.check_call('pynast -l 0 -p 0.5 -i ' +read_file+ ' -t ' +gg_db+ ' -a ' +out_dir,  shell=True)
      counter +=1
      sequences = list(SeqIO.parse(open(out_dir), 'fasta'))
      sequence_list.append(sequences)
    for f in list(sequence_list[0]):
      for r in list(sequence_list[1]):
        p_seq = []
        if f.id == r.id:
          for i,j in zip(list(f.seq),list(r.seq)):
            if i == j:
              p_seq.append(i)
            elif i.isalpha() and j == '-' and i != j:
              p_seq.append(i)
            elif j.isalpha() and i == '-' and j != i:
              p_seq.append(j)
            elif i.isalpha() and j.isalpha() and i != j:
              p_seq.append('.')
            else:
              print 'Error in aligning encountered'
              exit(1)
        else:
          continue
        startend = []
        ind = 0
        for n in list(f.seq[::-1]):
          if n.isalpha():
            startend.append(len(f.seq)-ind)
            break
          else:
            ind += 1
        ind = 0
        for n in list(r.seq):
          if n.isalpha():
            startend.append(ind)
            break
          else:
            ind += 1
        print '>'+f.id
        print startend
        if startend[0] < startend[1]:
          out = [''.join(p_seq[:startend[0]]),''.join(p_seq[startend[0]:startend[1]]).replace('-','.'),''.join(p_seq[startend[1]:])]
          print ''.join(out)
        elif startend[0] > startend[1]:
          startend = []
          ind = 0
          for n in list(r.seq[::-1]):
            if n.isalpha():
              startend.append(len(f.seq)-ind)
              break
            else:
              ind += 1
          ind = 0
          for n in list(f.seq):
            if n.isalpha():
              startend.append(ind)
              break
            else:
              ind += 1
          out = [''.join(p_seq[:startend[0]]),''.join(p_seq[startend[1]:startend[0]]).replace('-','.'),''.join(p_seq[startend[1]:])]
          print ''.join(out)

# Reads in the HMM model, forward and reverse reads
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Creates joint alignments of forward and reverse reads.')
    parser.add_argument('-f','--forward_reads', help='Forward reads to be aligned', required = True)
    parser.add_argument('-r','--reverse_reads', help='Reverse reads to be aligned', required = True)
    parser.add_argument('-m','--hmm', help='HMM model to align to (protein)')
    parser.add_argument('-g','--gg_alignment', help='GG database to align to (dna)')
    parser.add_argument('-t','--type', help='dna or prot', choices = ['dna','prot'], required = True)
    parser.add_argument('-e','--evalue', help='Cutoff for hmmsearch')
    args=parser.parse_args()

    # Create list of reads based or input
    reads_list = [args.forward_reads, args.reverse_reads]
    table_list = []
    counter = 0
    suff = ['for', 'rev']
    for read_file in reads_list:
      splt = read_file.split('/')
      table_list.append("./"+splt[len(splt)-1].split('.')[0]+"_"+suff[counter]+".txt")
      counter += 1
    base = args.forward_reads.split('/')[len(args.forward_reads.split('/'))-1].split('.')[0]

    if args.type == 'dna':
      Searcher().pynast(reads_list, args.gg_alignment)
      delete([base+'_for*', base+'_rev*', '*pynast*'])
      exit(1)
    Searcher().hmmsearch(args.hmm, reads_list)
    orf_names = Searcher().HMM_hit_extract(table_list, reads_list)
    Searcher().fxtract(base+'_names.txt', reads_list)
    orf_seq_list = [base+'_for.faa', base+'_rev.faa']
    Searcher().orf_extractor(orf_seq_list, orf_names)
    aln_list = [base+'_for_hits.orf', base+'_rev_hits.orf']
    Searcher().fasta_parser(aln_list, args.hmm)
    delete([base+'_for*', base+'_rev*', base+'_names.txt'])
    exit(1)

