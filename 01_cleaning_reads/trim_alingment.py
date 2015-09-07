#!/usr/bin/env python

__author__ = "Joel Boyd"
__copyright__ = "Copyright 2015"
__credits__ = ["Joel Boyd"]
__license__ = "GPL3"
__maintainer__ = "Joel Boyd"
__email__ = "joel.boyd near uq.net.au"
__status__ = "Development"
__version__ = "0.0.1"


################################################################################
################################# - Imports - ##################################

# System imports
import argparse
import logging
from Bio import SeqIO

################################################################################
################################# - Setup - ####################################

debug={1:logging.CRITICAL,
       2:logging.ERROR,
       3:logging.WARNING,
       4:logging.INFO,
       5:logging.DEBUG}

################################################################################
################################## - Code - ####################################

class AlignmentTrimmer:
    def __init__(self, alignment):
        
        logging.info("Reading alignment")
        self.alignment = list(SeqIO.parse(open(alignment), "fasta"))
        self.aln_len = len(self.alignment[0].seq)
        self.nseq = len(self.alignment)
        logging.debug("Alignment length: %i" % self.aln_len)
        logging.debug("Number of sequences: %i" % self.nseq)
    
    def _write(self, output_file_path, to_filter):
        logging.info("Writing trimmed sequences to file: %s" % output_file_path)
        sequences = [x for x in self.alignment if x not in to_filter]
        SeqIO.write(sequences, open(output_file_path, "w"), "fasta")
    
    def trim(self, output, cutoff):

        to_trim = []
        for idx in range(0, self.aln_len):
            pos_covered = []
            for record in self.alignment:
                if record.seq[idx] == '-':
                    continue
                else:
                    pos_covered.append(record)
            ncovered = len(pos_covered)
            coverage = float(ncovered)/self.nseq
            if coverage < cutoff:
                logging.debug("Position %i has coverage %s which is below %f, removing %i sequences aligned to this position" % (idx, coverage, cutoff, ncovered))
                for rec in pos_covered:
                    if rec not in to_trim:
                        to_trim.append(rec)
                
        logging.info("Excluded %i sequences from the alignment" % len(to_trim))
        self._write(output, to_trim)

################################################################################
################################ - argparser - #################################
          
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Remove infrequent positions from alignment''')
    parser.add_argument('--alignment', type=str, help='Input alignment', required=True)
    parser.add_argument('--output', type=str, help='Output alignment', required=True)
    parser.add_argument('--cutoff', type=int, help='cutoff (at each position in the hmm, any reads that have bases at this position will \
                                                    be cut out if there < this percent of reads align to that position. Default = 0.05)', default=0.0005)

    parser.add_argument('--log', help='Output logging information to file', type=str, default=False)
    parser.add_argument('--verbosity', help='1 - 5, 1 being silent, 5 being noisy indeed. Default = 4', type=int, default=4)

    args = parser.parse_args()
    
    if args.log:
        if os.path.isfile(args.log): raise Exception("File %s exists" % args.log)
        logging.basicConfig(filename=args.log, level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=debug[args.verbosity], format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    
    at=AlignmentTrimmer(args.alignment)
    
    at.trim(
            args.output,
            float(args.cutoff)
            )

    exit()

