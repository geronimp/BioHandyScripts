import sys
from Bio import SeqIO
import code

heretics = sys.argv[1]
population = sys.argv[2]
cleansed_population = sys.argv[3]

heretic_list = []

for heretic in open(heretics, 'r'):
  heretic_list.append(heretic.rstrip())

heretic_set = set(heretic_list)

with open(cleansed_population, 'w') as cleansed:
  for member in list(SeqIO.parse(open(population, 'r'), 'fasta')):
    if member.id not in heretic_set:
      SeqIO.write(member, cleansed, "fasta")
    else:
      continue #PURGED

print 'Nobody suspects the Spanish Inquisition...'


exit()
code.interact(local=locals())
