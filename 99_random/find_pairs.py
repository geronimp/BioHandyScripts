import sys
from Bio import SeqIO
import os

forward=sys.argv[1]
reverse=sys.argv[2]
forward_out=os.path.basename(sys.argv[1]).split('.')[0] + '_pair.fastq'
reverse_out=os.path.basename(sys.argv[2]).split('.')[0] + '_pair.fastq'


forward_dict = SeqIO.to_dict(SeqIO.parse(forward, "fastq"))
reverse_dict = SeqIO.to_dict(SeqIO.parse(reverse, "fastq"))

pairs=[x for x in forward_dict.keys() if x in reverse_dict.keys()]

if any(pairs):
  with open(forward_out, 'w') as f_out:
    for record in pairs:
      SeqIO.write(forward_dict[record], f_out, "fastq")
  with open(reverse_out, 'w') as r_out:
    for record in pairs:
      SeqIO.write(reverse_dict[record], r_out, "fastq")
else:
  exit(0)
