#!/usr/bin/env python
###############################################################################
#                                                                             #
# This program is free software: you can redistribute it and/or modify        #
# it under the terms of the GNU General Public License as published by        #
# the Free Software Foundation, either version 3 of the License, or           #
# (at your option) any later version.                                         #
#                                                                             #
# This program is distributed in the hope that it will be useful,             #
# but WITHOUT ANY WARRANTY; without even the implied warranty of              #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the                #
# GNU General Public License for more details.                                #
#                                                                             #
# You should have received a copy of the GNU General Public License           #
# along with this program. If not, see <http://www.gnu.org/licenses/>.        #
#                                                                             #
###############################################################################
# *py <GPGK> <output>
#
# Provide the seqinfo and the taxonomy file within a graftM package.
#

import sys
import os
import tempfile
import gzip
import json
from Bio import SeqIO
from collections import Counter
from graftm.getaxnseq import Getaxnseq
gpkg = sys.argv[1]
output_file = sys.argv[2]
contents=json.load(open(os.path.join(gpkg, "CONTENTS.json")))
refpkg=os.path.join(gpkg, contents['refpkg'])
refpkg_contents=json.load(open(os.path.join(refpkg, "CONTENTS.json")))
pfam_annotations = "/srv/projects/abisko/Joel/99_phd/01_Projects/08_gpkgs_for_ko_groups/data/PFAM_annotation_of_GTDB_hits"
genomes_file = '/srv/db/img/latest/genomes/'
seqinfo_path=os.path.join(refpkg, refpkg_contents['files']['seq_info'])
taxonomy_path=os.path.join(refpkg, refpkg_contents['files']['taxonomy'])
unaligned_sequences=os.path.join(gpkg, contents['unaligned_sequence_database'])
aligned_sequences=os.path.join(refpkg, refpkg_contents['files']['aln_fasta'])
prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
GTDB_FILE = "/srv/projects/share/for_joel_gtdbgenomes_prokka/ace_genomes_prokka"
KO_ID=seqinfo_path.split('_')[0]
IMG_TAX = os.path.join("/srv/projects/abisko/Joel/99_phd/01_Projects/08_gpkgs_for_ko_groups/data/IMG_taxonomy/",
                       KO_ID + '.txt.gz')
gene_id_to_tax = {x.strip().split('\t')[1]:x.split('\t')[0] for x in  gzip.open(IMG_TAX)}
sprot_ids = set([x.strip() for x in gzip.open("/srv/projects/abisko/Joel/99_phd/01_Projects/08_gpkgs_for_ko_groups/flat/uniprot-all.list.gz")])

parsed_img_file = {x.split()[1]:os.path.join(genomes_file,
                                             x.split()[0].strip(), x.split()[0].strip()+'.gff')
                   for x in gzip.open(IMG_TAX)}
cog_annotation_file = {x.split()[1]: os.path.join(genomes_file, x.split()[0].strip(),
                                                  x.split()[0].strip()+'.cog.tab.txt')
                       for x in gzip.open(IMG_TAX)}
pfam_annotation_file = {x.split()[1]: os.path.join(genomes_file, x.split()[0].strip(),
                                                   x.split()[0].strip()+'.pfam.tab.txt')
                        for x in gzip.open(IMG_TAX)}
tigrfam_annotation_file = {x.split()[1]: os.path.join(genomes_file, x.split()[0].strip(),
                                                      x.split()[0].strip()+'.tigrfam.tab.txt')
                           for x in gzip.open(IMG_TAX)}

xref_annotation_file = {x.split()[1]: os.path.join(genomes_file, x.split()[0].strip(),
                                                      x.split()[0].strip()+'.xref.tab.txt')
                           for x in gzip.open(IMG_TAX)}

def parse_pfam_annotation_hmmoutput(domtblout_file):
    domtblout = {}
    for line in gzip.open(domtblout_file):
        if line.startswith('#'): continue
        seqname = line.split()[0]
        pfam_hit = line.split()[3]
        pfam_id = line.split()[4]
        if seqname in domtblout:

            domtblout[seqname]['annotation_ids'].append(pfam_id)
            domtblout[seqname]['annotations'].append(pfam_hit)
        else:
            domtblout[seqname] = {"annotation_ids":[pfam_id],
                                  "annotations":[pfam_hit]}
    return domtblout

def parse_xref_file(xref_file):
    xref = {}
    for line in open(xref_file):
        if line.startswith("gene_oid"): continue
        gene_id, db_name, id = tuple(line.split('\t'))
        if db_name == "UniProtKB":
            xref[gene_id] = id.strip()
    return xref

def parse_pfam_file(pfam_file):
    pfam = {}
    for line in open(pfam_file):
        if line.startswith("gene_oid"): continue
        line = line.strip().split('\t')
        gene_name = line[0]
        pfam_annotation_id = line[8]
        pfam_annotation = line[9]
        if gene_name in pfam:
            pfam[gene_name]['annotation_ids'].append(pfam_annotation_id)
            pfam[gene_name]['annotations'].append(pfam_annotation)
        else:
            pfam[gene_name] = {'annotation_ids':[pfam_annotation_id],
                               'annotations': [pfam_annotation]}
    return pfam

def parse_tigrfam_file(tigrfam_file):
    tigrfam = {}
    for line in open(tigrfam_file):
        if line.startswith("gene_oid"): continue
        line = line.strip().split('\t')
        gene_name = line[0]
        tigrfam_annotation_id = line[6]
        tigrfam_annotation = line[7].replace(' ','_')
                                    #.translate(None, "/\():-,[]!?;")
        if gene_name in tigrfam:
            tigrfam[gene_name]['annotation_ids'].append(tigrfam_annotation_id)
            tigrfam[gene_name]['annotations'].append(tigrfam_annotation)
        else:
            tigrfam[gene_name] = {'annotation_ids':[tigrfam_annotation_id],
                                  'annotations': [tigrfam_annotation]}
    return tigrfam

def parse_cog_file(cog_file):
    cog = {}
    for line in open(cog_file):
        if line.startswith("gene_oid"): continue
        line = line.strip().split('\t')
        gene_name = line[0]
        cog_annotation_id = line[9]
        cog_annotation = line[10].replace(' ','_')
                                    #.translate(None, "/\():-,[]!?;")
        if gene_name in cog:
            cog[gene_name]['annotation_ids'].append(cog_annotation_id)
            cog[gene_name]['annotations'].append(cog_annotation)
        else:
            cog[gene_name] = {'annotation_ids':[cog_annotation_id],
                              'annotations': [cog_annotation]}
    return cog

def parse_GTDB_gff_file(gff_file):

    gff = {}

    for line in open(gff_file):
        if line.startswith("##FASTA"):break
        if line.startswith('#'): continue
        line = line.strip().split('\t')
        if len(line) == 9:
            c = line[8].split(';')
            id = c[0][3:]
            if c[-1].startswith("product"):
                classification = c[-1][8:]
            else:
                classification = line[2]
            gff[id] = classification
    return gff


arb_db = {}
def parse_seq_ids(seq_ids):
    dbs = {'pfams': {},
           'cogs': {},
           'tigrfams': {},
           'gffs': {}}

    for seq_id in seq_ids:

        if '~' in seq_id:
            seq, genome_id = tuple(seq_id.split('~'))
            split_seq_id = seq.split('_')
            contig_id = split_seq_id[0]
            seq_name = '_'.join(split_seq_id[-3:][1:])
            gff_path = os.path.join(GTDB_FILE,
                                  genome_id + '.gff')
            parsed_gff_file = parse_GTDB_gff_file(gff_path)

            if seq_id in parsed_gtdb_pfam_annotation:
                pfam_classifications = ','.join(parsed_gtdb_pfam_annotation[seq_id]['annotations'])
                pfam_ids = ','.join(parsed_gtdb_pfam_annotation[seq_id]['annotation_ids'])
            else:
                pfam_classifications = None
                pfam_ids = None

            arb_db[seq_id] = {}
            tmp={"pfam_ids":pfam_ids,
            "pfam_classifications":pfam_classifications,
            "cog_ids":None,
            "cog_classifications":None,
            "tigrfam_ids":None,
            "tigrfam_classifications":None,
            "source":"GTDB",
            "swiss_prot":"No"}
            arb_db[seq_id].update(tmp)

        else:
            parsed_gff_file = parse_GTDB_gff_file(parsed_img_file[seq_id])
            parsed_pfam_file = parse_pfam_file(pfam_annotation_file[seq_id])
            parsed_tigrfam_file = parse_tigrfam_file(tigrfam_annotation_file[seq_id])
            parsed_cog_file = parse_cog_file(cog_annotation_file[seq_id])

            if seq_id in parsed_pfam_file:
                pfam_ids = ','.join(parsed_pfam_file[seq_id]['annotation_ids'])
                pfam_classifications = ','.join(parsed_pfam_file[seq_id]['annotations'])
            else:
                pfam_ids = None
                pfam_classifications = None
            if seq_id in parsed_cog_file:
                cog_ids = ','.join(parsed_cog_file[seq_id]['annotation_ids'])
                cog_classifications = ','.join(parsed_cog_file[seq_id]['annotations'])
            else:
                cog_ids = None
                cog_classifications = None
            if seq_id in parsed_tigrfam_file:
                tigrfam_ids = ','.join(parsed_tigrfam_file[seq_id]['annotation_ids'])
                tigrfam_classifications = ','.join(parsed_tigrfam_file[seq_id]['annotations'])
            else:
                tigrfam_ids = None
                tigrfam_classifications = None

            if seq_id in xref_annotation_file:
                if os.path.isfile(xref_annotation_file[seq_id]):
                    parsed_xref_file = parse_xref_file(xref_annotation_file[seq_id])
                    if seq_id in parsed_xref_file:
                        uniprot_id = parsed_xref_file[seq_id]
                        swiss_prot = ("Yes" if uniprot_id in sprot_ids else "No")
                    else:
                        swiss_prot = "No"
                else:
                    swiss_prot = "No"
            else:
                swiss_prot = "No"



            seq_name = seq_id
            arb_db[seq_id] = {}
            tmp={"pfam_ids":pfam_ids,
            "pfam_classifications":pfam_classifications,
            "cog_ids":cog_ids,
            "cog_classifications":cog_classifications,
            "tigrfam_ids":tigrfam_ids,
            "tigrfam_classifications":tigrfam_classifications,
            "source":"IMG",
            "swiss_prot": swiss_prot}
            arb_db[seq_name].update(tmp)
        arb_db[seq_id]["Classification"] = '_'.join(parsed_gff_file[seq_name].split())



parsed_gtdb_pfam_annotation = parse_pfam_annotation_hmmoutput(os.path.join(pfam_annotations, KO_ID+'_gtdb_pfam_annotation.domtblout.txt.gz'))
parsed_unaligned_sequences = SeqIO.to_dict(SeqIO.parse(open(unaligned_sequences), "fasta"))
parsed_aligned_sequences=SeqIO.to_dict(SeqIO.parse(open(aligned_sequences), "fasta"))

gtns = Getaxnseq()

parsed_taxonomy = {}
ids =[]
for id, tax in gtns.\
            read_taxtastic_taxonomy_and_seqinfo(open(taxonomy_path),
                                                open(seqinfo_path)).iteritems():

    ids.append(id)
    if '~' in id:
        gene_name, genome_id = id, id.split('~')[1]
        
        if genome_id in gene_id_to_tax:
            raise Exception("Genome ID encountered twice: %s" % genome_id)
        else:
            gene_id_to_tax[id] = genome_id
    #
    parsed_taxonomy[id] = tax
sequence_ids = parse_seq_ids(ids)

for seq_id in arb_db:

  tax = [x.translate(None, "/\:-,[]!?;") for x in  parsed_taxonomy[seq_id]]
  tax = '; '.join([x.translate(None, "/\:-,[]!?;")
                            for x in  parsed_taxonomy[seq_id]] + prefixes[len(tax):])
  arb_db[seq_id].update({"tax_string":tax})

genome_list = [gene_id_to_tax[x] for x in arb_db.keys()]

para_db = {key:range(1, 1+item) for key,item in Counter(genome_list).iteritems() if item > 1}

for name in arb_db.keys():
    if gene_id_to_tax[name] in para_db:
        para='p'+str(para_db[gene_id_to_tax[name]].pop(0))
    else:
        para="0"
    genome_id = (name.split('~')[0] if '~' in name else gene_id_to_tax[name])
    arb_db[name]["paralog"]=para
    arb_db[name]["genome_id"]=genome_id

output=["name\tfull_name\tClassification\ttax_string\tpfam_ids\tpfam_classifications\
\tcog_ids\tcog_classifications\ttigrfam_ids\ttigrfam_classifications\tsource\tSwiss_Prot_annotated\
\tsequence_length\tdomain\tphylum\tclass\torder\tfamily\tgenus\tspecies\tparalog\tgenome_id\tAlignment"]
for name, entry in arb_db.iteritems():
    d,p,c,o,f,g,s=tuple(arb_db[name]["tax_string"].split('; '))
    if name in parsed_aligned_sequences:
        line = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % (name,
                                                                           parsed_unaligned_sequences[name].description,
                                                                           arb_db[name]["Classification"],
                                                                           arb_db[name]["tax_string"],
                                                                           arb_db[name]['pfam_ids'],
                                                                           arb_db[name]['pfam_classifications'],
                                                                           arb_db[name]['cog_ids'],
                                  arb_db[name]['cog_classifications'],
                                  arb_db[name]['tigrfam_ids'],
                                  arb_db[name]['tigrfam_classifications'],
                                  arb_db[name]["source"],
                                  arb_db[name]["swiss_prot"],
                                  str(len(parsed_unaligned_sequences[name].seq)),
                                  d,p,c,o,f,g,s,
                                  arb_db[name]['paralog'],
                                  arb_db[name]['genome_id'],
                                  str(parsed_aligned_sequences[name].seq))
        output.append(line)
with open(output_file, "w") as out:
  for l in output:
    out.write(l + '\n')
