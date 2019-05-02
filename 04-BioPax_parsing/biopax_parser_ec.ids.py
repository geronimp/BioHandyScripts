import sys
import os
input_biopax = sys.argv[1]
sep = '\t'

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

pathway = list()
pathway_info = dict()
reaction_info = dict()
compound_info = dict()

def rename(compound_info, compound_list):
	output_list = list()
	
	for compound_id in compound_list:

		if not compound_id.startswith(protein):
			compound_name = compound_info[compound_id]			
			output_list.append(compound_name[0])

	return output_list

def merger(pathway_info, reaction_info, compound_info):
	left_list = list()
	right_list = list()

	for key, entry in pathway_info.items():
		conversion = entry[stepConversion]
		entry.update(reaction_info[conversion])
		entry[left] = rename(compound_info, entry[left])
		entry[right] = rename(compound_info, entry[right])

	return pathway_info

io_biopax = open(input_biopax)
header = io_biopax.readline().split(sep)

for line in io_biopax:
	sline = line.split(sep)
	cl, id, pr, at, va, de = sline[:6]
	if cl == Pathway:

		if pr == standardName:
		
			pathway_name = de.strip().replace(' ', '_')
	
	if cl == BiochemicalPathwayStep:
	
		if id not in pathway_info:
			pathway_info[id] = dict()
			
		if len(pathway) == 0:
			pathway.append(id)
			
		if pr == nextStep:
			nextstep_value = va.replace("#", "")
			
			if id in pathway:
				
				if nextstep_value not in pathway:
					id_index = pathway.index(id)
					pathway = pathway[:id_index+1] + [nextstep_value] + pathway[id_index+1:] 
			
			elif nextstep_value in pathway:
				
				if id not in pathway:
					nextstep_index = pathway.index(nextstep_value)
					pathway = pathway[:nextstep_index] + [id] + pathway[nextstep_index:]

			else:
				pathway += [id, nextstep_value]
	
		if pr == stepConversion:
			pathway_info[id][stepConversion] = va.replace("#", "")
	
		if pr == stepDirection:
			pathway_info[id][stepDirection] = de.replace("#", "").strip().split('-TO-')

	if cl == BiochemicalReaction:

		if id not in reaction_info:
			reaction_info[id] = {eCNumber: list(), left: list(), right: list()}

		if pr == eCNumber:
			reaction_info[id][eCNumber].append(de.replace("#", "").strip())

		if pr == left:
			reaction_info[id][left].append(va.replace("#", ""))
		
		if pr == right:
			reaction_info[id][right].append(va.replace("#", ""))

	if cl == SmallMolecule:

		if id not in compound_info:
			compound_info[id] = list()
		
		if pr == standardName:
			compound_info[id].append(de.strip())

merged_pathway_info = merger(pathway_info, reaction_info, compound_info)

pathway_definition = list()

for step in pathway:
	ec_numbers = merged_pathway_info[step][eCNumber]
	step_definition = list()

	for ec in ec_numbers:
		step_definition.append(ec)

	pathway_definition.append(','.join(step_definition))

print("%s\t%s" % (os.path.basename(input_biopax.split('.')[0]), ' '.join(pathway_definition)))