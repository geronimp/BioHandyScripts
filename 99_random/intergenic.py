#!/usr/bin/env python
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2017"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
__version__ = "0.0.1"

###############################################################################

# Imports
import argparse
import logging
import subprocess
import os

###############################################################################

debug={1: logging.CRITICAL,
       2: logging.ERROR,
       3: logging.WARNING,
       4: logging.INFO,
       5: logging.DEBUG}

###############################################################################

class Intergenic:
	
	HASH	= '#'
	TAB 	= '\t'
	NEWLINE	= '\n'
	CDS 	= 'CDS'
	
	def _parse_gffs(self, gffs):
		'''
		Extract the coding regions of input GFF files
		
		Parameters
		----------
		gffs - List. list of paths to gff files to parse
		
		Output
		------
		Dictionary with each contig name as the keys, and a set of positions
		of genic regions as enties.
		'''
		
		genic_regions = {}

		for gff in gffs:
			logging.info('Parsing GFF file: %s' % (gff))
		
			for line in open(gff):
		
				if line.startswith(self.HASH): # Skip headers and comment lines
					continue
		
				contig, program, region_type, start, end, score, strand, phase, attributes = line.split(self.TAB)
		
				if region_type == self.CDS:
					range_limits = [int(start), int(end)]
					genic_region = range( min(range_limits), max(range_limits) )

					if contig not in genic_regions:
						genic_regions[contig] = set([])
						
					for position in genic_region:
						genic_regions[contig].add(position)

		return genic_regions

	def _filter_bam(self, contig_regions, bams):
		'''
		Extract reads that map to intergenic regions.
		
		Parameters
		----------
		contig_regions 	- Dict. dictionary with each contig name as the keys, and a set of positions
						  of genic regions as enties.
		bam 			- List. List of paths to bam files to parse.

		Output
		------
		List of lists, each containing the column entires for the line to be written to the output line.
		'''

		output_lines = []
		
		for bam in bams:
			contig_counts 	= {contig:[0,0] for contig in contig_regions.keys()}
			bam_name      	= os.path.basename(bam)
			cmd    			= 'samtools view %s' % bam	
			logging.info('Filtering BAM file: %s' % (bam_name))
			
			for line in subprocess.check_output(cmd, shell=True).strip().split(self.NEWLINE):
				start  			= int(line.split(self.TAB)[3])
				end    			= start + len(line.split(self.TAB)[9])
				mapped_region 	= range(start, end)
				contig 			= line.split(self.TAB)[2]
				
				try:
					contig_regions[contig]
				except:
					import IPython ; IPython.embed()

				if len(contig_regions[contig].intersection(mapped_region))==0:
					contig_counts[contig][0] += 1 # Intergenic
				else:
					contig_counts[contig][1] += 1 # Genic
			
			for contig, counts in contig_counts.items():
				output_lines.append([bam_name, contig, str(counts[0]), str(counts[1])])

		return output_lines

	def _write(self, output_lines, output_file_path):
		'''
		Write results to file		

		Parameters
		----------
		output_file_path - String. Write results to this file.
		'''

		header = ['bam', 'contig', 'intergenic', 'genic']
		
		with open(output_file_path, 'w') as out_io:
			out_io.write(header + self.NEWLINE)
		
			for line in output_lines:
				out_io.write('\t'.join(line) + self.NEWLINE)

	def do(self, gffs, bams, output):

		logging.info('Parsing GFFs')
		genic_regions = self._parse_gffs(gffs)
		
		logging.info('Filtering BAMs')
		filtered_reads = self._filter_bam(genic_regions, bams)
		
		logging.info('Writing results to output: %s' % (output))
		self._write(filtered_reads, output)

		logging.info('Done')

###############################################################################

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='''Count intergenic vs genic reads using GFF and BAM files''')
	
	parser.add_argument('--gff', type=str, help='space separated list of gff files', nargs='+', required=True)
	parser.add_argument('--bam', type=str, help='bam file', nargs = '+', required=True)
	parser.add_argument('--output', type=str, help='output file', required=True)
	parser.add_argument('--log', help='Output logging information to file', type=str, default=False)
	parser.add_argument('--verbosity', help='1 - 5, 1 being silent, 5 being noisy indeed. Default = 4', type=int, default=4)
	
	args = parser.parse_args()
	
	if args.log:
		if os.path.isfile(args.log): 
			raise Exception("File %s exists" % args.log)
		logging.basicConfig(filename=args.log, level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	else:
		logging.basicConfig(level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
	
	i = Intergenic()
	i.do(args.gff, args.bam, args.output)

	exit(0)
