#!/usr/bin/env python
## returns reads without ambiguous bases (Ns)
## e.g. usage: clean_Ns.py <INPUT.fa> <OUTPUT.fa>

import sys
try:
  from Bio import SeqIO
except ImportError:
  print 'Install Biopython'
  exit(0)

class CleanNs:
  def __init__(self):
    pass

  def collect_clean_reads(self, reads, output):
    with open(output, 'w') as out:
      for record in list(SeqIO.parse(open(reads), 'fastq')):
        if 'N' not in record.seq:
          SeqIO.write(record, out, "fastq")

  def main(self, reads, output):
    self.collect_clean_reads(reads, output)
    return

if __name__ == '__main__':
  CleanNs().main(sys.argv[1], sys.argv[2])
