import sys
import xml.etree.ElementTree as ET

direction_dict = {"RIGHT": 'right', "LEFT": 'left'}
BiochemicalPathwayStep = 'BiochemicalPathwayStep'
nextStep = 'nextStep'
stepConversion = "stepConversion"
stepDirection = "stepDirection"
BiochemicalReaction = "BiochemicalReaction"
eCNumber = "eCNumber"
left = 'left'
right = 'right'
SmallMolecule = 'SmallMolecule'
standardName = 'standardName'
protein = 'Protein'
Pathway = 'Pathway'
ID = 'ID'

pathway = list()
pathway_info = dict()
reaction_info = dict()
compound_info = dict()

biopax = ET.parse(sys.argv[1])
root = biopax.getroot()

tag = "{http://www.biopax.org/release/biopax-level3.owl#}"
synt = '{http://www.w3.org/1999/02/22-rdf-syntax-ns#}'
# Get pathway name
pathway_name = root.find(tag + Pathway).find(tag + standardName).text
# Iterate through steps 
for bps in root.findall(tag + BiochemicalPathwayStep):
    id = bps.attrib[synt + ID]
    
    if id not in pathway_info:
        pathway_info[id] = dict()
    
    if len(pathway) == 0:
        pathway.append(id)
    import IPython
    IPython.embed()

    if pr == nextStep:
        nextstep_value = va.replace("#", "")
