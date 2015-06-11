#!/usr/bin/env python
import sys
import random
import os
try:
  from Bio import SeqIO
except ImportError:
  print """

fasta_windower failed to load Biopython. Please install.

  """
  exit(1)


in_file = sys.argv[1]

out = []
record = 0

for sequence in SeqIO.parse(open(in_file, "rU"), "fasta"):
  ### Defining the window
  # Define start of chop
  start = 0
  # Define end of chop
  finish = int(sys.argv[2])
  # Define increment length
  inc = finish/2

  ### Reset naming scheme
  # Contig name
  fragment = 0
  # Sequence name
  record += 1

  ### Times to iterate over sequence
  times = len(sequence.seq)/ int(sys.argv[2]) *2 + 1

  while times > 0:
    fragment += 1

    random_number = str(random.randint(1, 1000000))
    if len(sequence.seq[start:finish]) == int(sys.argv[2]):
      out.append('>r_'+random_number+'Fragment_from_'+sequence.id+'_frag_'+ str(fragment)+'_fasta\n')
      out.append(sequence.seq[start:finish] + '\n')
    else:
      out.append('>r_'+random_number+'_Fragment_from_'+sequence.id+'_frag_'+ str(fragment)+'_fasta\n')
      out.append(sequence.seq[start:finish] + '\n')
      break
    start += inc
    finish += inc
    times -=1


with open(os.path.basename(in_file).split('.')[0]+'.windowed.fa', 'w') as o:
  for x in out:
    o.write(str(x))
