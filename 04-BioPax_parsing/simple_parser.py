import sys
import xml.etree.ElementTree as ET
Pathway = 'Pathway'
standardName = 'standardName'

biopax = ET.parse(sys.argv[1])
root = biopax.getroot()

tag = "{http://www.biopax.org/release/biopax-level3.owl#}"
synt = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
# Get pathway name
pathway_name = root.find(tag + Pathway).find(tag + standardName).text

ecs = []
for line in open(sys.argv[1]):
    if "eCNumber" in line:
        ecs.append(line.strip().split('>')[1].split('<')[0])
    
print('\t'.join([pathway_name, ' '.join(ecs)]))
