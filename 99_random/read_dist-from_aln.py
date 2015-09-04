from Bio import SeqIO
import sys
alignment = []
for x in sys.argv[1].split(','):
  alignment += list(SeqIO.parse(x, 'fasta'))

length=len(str(alignment[0].seq))
profile=[0]*length
for record in alignment:
  for idx, nucl in enumerate(str(record.seq)):
    if nucl != '-':
      profile[idx]+=1
for idx in profile: print idx
