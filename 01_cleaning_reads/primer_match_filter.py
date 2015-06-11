#!/usr/bin/env python
# Returns the number of reads that have an AVERAGE quality that sits above a specified phred cutoff
# assess_quality.py [FILE_NAME] [OUTOUT_FILE]

import sys
try:
    from Bio import SeqIO
except ImportError:
    print 'Please install Biopython first'
    exit(1)

fastq_file = sys.argv[1]
forward={'AGAGTTTGATCATGGCTCAG', 'AGAGTTTGATCCTGGCTCAG'}
reverse={'GAATTACCGCGGCGGCTG', 'GAATTACCGCGGCTGCTG', 'GTATTACCGCGGCGGCTG', 'GTATTACCGCGGCTGCTG'}
with open(sys.argv[2], 'w') as out:
  for record in SeqIO.parse(open(fastq_file, 'r'), "fastq"):
    if str(record.seq)[:18] in reverse:
      SeqIO.write(record, out, "fastq")
    elif str(record.seq)[:20] in forward:
      SeqIO.write(record, out, "fastq")
    else:
      continue
