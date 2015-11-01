#!/usr/bin/env python
###############################################################################
#
# fundec.py <tree> <database_file> <%id cutoff>
# I should probably iclude some reguular expression recognition to make the
# counts of each annotation more robust. I imagine that as is this
# code will only cluster the major branches. I should really set the
# %ID cutoff high to avoid over clustering. I dont want that to come
# back to bite.
#
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
 
__author__ = "Joel Boyd"
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
 
###############################################################################

import sys
import argparse
import logging
import os
import shutil

from skbio.tree import TreeNode
from collections import Counter

###############################################################################
############################### - Exceptions - ################################

class BadTreeFileException(Exception):
    pass

###############################################################################
################################ - Classes - ##################################

class FunDec:
    
    def __init__(self, db_name):
        self.db_name = db_name
        
    # ------------------------------- Parsing ------------------------------- #
    def _parse_db_file(self, database):
        db_dict = {}
        for line in open(database):
            if line.startswith('name'): continue
            line=line.strip().split('\t')
            db_dict[line[0]]= {'name':line[0],
                               'full_name':line[1],
                               'annotation':line[2],
                               'Domain':line[3],
                               'Phylum':line[4],
                               'Class':line[5],
                               'Order':line[6],
                               'Family':line[7],
                               'Genus':line[8],
                               'Species':line[9],
                               'sequence_length':line[10],
                               'Alignment':line[11],
                               'classificaton':['%s_homolog' % self.db_name]}
            
        return db_dict

    def _open_tree(self, tree_path):
        tree_obj=TreeNode.read(open(tree_path))
        return tree_obj
    
    # ---------------------------- Manipulating ----------------------------- #
    def _gather_annotations(self, tips, db):
        'Gather annotations into a list, return relative abundances of each'
        annotations = []
        for tip in tips:
            annotations.append(db[tip.name.replace(' ', '_')]['annotation'])
        counted_annotations = Counter(annotations)
        total_annotations=float(sum(counted_annotations.values()))
        counted_annotations_relab={}
        for annotation, count in counted_annotations.iteritems():
            counted_annotations_relab[annotation]=\
                                        float(count)/total_annotations
        return counted_annotations_relab

    def  _gather_ancestors(self, node):
        previous_clade = []
        for ancestor in node.ancestors():
            if ancestor.name:
                if ancestor.name.startswith('a__'):
                    previous_clade.append(ancestor.name)
        return previous_clade
    
    def _write_tree(self, tree, output_path):
        tree.write(output_path,
                   format = "newick")

    # =============================== MAIN ================================== #
    def main(self, tree_path, database_path, cutoff, output_path):
        database = self._parse_db_file(database_path)
        tree     = self._open_tree(tree_path)
        
        for node in tree.preorder():
            tip_annotations = self._gather_annotations(node.tips(),
                                                       database)
            if any([x for x in tip_annotations.values() if x > cutoff]):
                previous_clades = self._gather_ancestors(node)
                key, _ = max(tip_annotations.iteritems(), key=lambda x:x[1])
                new_classification = "a__%s" % key
                if new_classification not in previous_clades:
                    node.name = new_classification
                    
        self._write_tree(tree, output_path)

###############################################################################
############################### - Functions - #################################

def check_args(args):
    # Constants
    tree_extensions = (".tree", ".tre")
    # Checks
    if args.tree.endswith(tree_extensions):
        args.base = '.'.join(args.tree.split('.')[:-1])
        if not args.output:
            args.output = args.base + '_functional_decoration.tree'            
    else:
        raise BadTreeFileException("Tree file without one of the following \
Extensions was provided: %s" % ' '.join(tree_extensions))
    return args

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s %(levelname)s: %(message)s', 
                        datefmt='%m/%d/%Y %I:%M:%S %p')   
    parser = argparse.ArgumentParser(description='''Functionally decorate a  
tree as best as possible''')
    parser.add_argument('--tree', 
                        help='Path to newick tree file',
                        required=True)
    parser.add_argument('--database', 
                        help='Path to datasbase file', 
                        required=True)
    parser.add_argument('--cutoff', 
                        help='This fraction of tips from a given node must \
have consistent taxonomy in order to be clustered (default = 0.9)', 
                        type=float, 
                        default=0.9)
    parser.add_argument('--output', 
                        help='Output tree to this file path.')
    
    args = check_args(parser.parse_args())
        
    fd = FunDec(args.base)
    status = fd.main(args.tree,
                     args.database,
                     args.cutoff,
                     args.output)

    if status == 0:
        exit(0)
    else:
        exit(1)
    