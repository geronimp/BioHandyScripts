import argparse

parser = argparse.ArgumentParser(description='''Extraction from GG database using accession ids specified''')
parser.add_argument('--accessions', type=str, help='Accession numbers of interes', required=True)
parser.add_argument('--gg_taxonomy', type=str, help='GG database to extarct from', required=True)
parser.add_argument('--mode', type=str, help='Extarct negative or positive hits', choices=['n', 'p'], required=True)
args = parser.parse_args()

ids=[x.rstrip() for x in open(args.accessions, 'r')]
set_ids = set(ids)

for i in [x.rstrip() for x in open(args.gg_taxonomy, 'r')]:
  splt=i.split()
  if args.mode == 'n':
    if splt[0] not in set_ids:
      print i
  if args.mode == 'p':
    if splt[0] in set_ids:
      print i


exit(1)
