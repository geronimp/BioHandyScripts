from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='''Filter an alignment.''')
parser.add_argument('--id', type=float, help='ID cutoff to filter.', required=True)
parser.add_argument('--alignment', type=str, help='Alignment', required=True)

args = parser.parse_args()

alignment = SeqIO.parse(args.alignment, "fasta")

for sequence in alignment:
    s = sequence.seq
    c = 1-float(s.count('-'))/float(len(s))
    if c>args.id:
        print '>%s' % sequence.id
        print sequence.seq