import argparse

parser = argparse.ArgumentParser(description='''Extraction from GG database using accession ids specified''')
parser.add_argument('--accessions', type=str, help='Accession numbers of interes', required=True)
parser.add_argument('--gg_taxonomy', type=str, help='GG database to extarct from', required=True)
parser.add_argument('-v', help='Extarct negative or positive hits', action="store_true", default=False)
args = parser.parse_args()

tax_dict={}

for i in [x.rstrip() for x in open(args.gg_taxonomy, 'r')]:
  splt=i.split()
  tax_dict[splt[0]]=' '.join(splt[1:]).strip()


if args.v:
  print "not right now"
else:
  for i in [x.strip() for x in open(args.accessions)]:
    if i in tax_dict:
      print "%s\t%s" % (i, tax_dict[i])

