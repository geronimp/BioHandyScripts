#!/usr/bin/env python
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Josh Daly"
__copyright__ = "Copyright 2015"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Development"

###############################################################################

# system imports
import argparse
import subprocess	
import re
import os
from Bio import SeqIO

def doWork( fasta_file ):

	genome_name = os.path.basename(fasta_file)
	genome_base = os.path.splitext(genome_name)[0]
	output_mhc_protein_path = genome_base + '.MHC.faa'
	
	for record in SeqIO.parse(fasta_file, "fasta"):
		count=0
		m2 = re.findall('C..CH', str(record.seq))
		m3 = re.findall('C...CH', str(record.seq))
		m4 = re.findall('C....CH', str(record.seq))
		m5 = re.findall('C.....CH', str(record.seq))
		print record.id +'\t' + str(len(m2)), str(len(m3)), str(len(m4)), str(len(m5))
		

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('fasta_file', help="")
	args = parser.parse_args()
	# do what we came here to do
	doWork(args.fasta_file)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
