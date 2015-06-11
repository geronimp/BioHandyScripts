#!/usr/bin/env python
# Returns the number of reads that have an AVERAGE quality that sits above a specified phred cutoff
# assess_quality.py [FILE_NAME] [CUTOFF]

import sys
try:
  from Bio import SeqIO
except ImportError:
  print 'Please install Biopython first'
  exit(1)

fastq_file = sys.argv[1]

reads_pass = 0

for record in SeqIO.parse(open(fastq_file, 'r'), "fastq"):
  q = record.letter_annotations["phred_quality"]
  score = sum(q) / float(len(q))
  if score > float(sys.argv[2]):
    reads_pass += 1

print reads_pass
