#!/usr/bin/env python


## USAGE::
## clip_Nshift_primer.py <IN> <PRIMER> <OUT>
# Created by Joel Boyd, 23/06/2015
# Australian Centre for ecogenomics

import sys
import regex
import subprocess
import tempfile
from Bio import SeqIO

class primerClipper:
    def __init__(self): pass
    
    def _translateDegenerate(self, primer):
        ## Does the job
        degenerates={"R": ["A", "G"],
                     "Y": ["C", "T"],
                     "S": ["G", "C"],
                     "W": ["A", "T"],
                     "K": ["G", "T"],
                     "M": ["A", "C"],
                     "B": ["C", "G", "T"],
                     "D": ["A", "G", "T"],
                     "H": ["A", "C", "T"],
                     "V": ["A", "C", "G"],
                     "N": ["A", "C", "G", "T"]}
        all=[]
        for nucl in primer:
            if nucl not in degenerates:
                if len(all) > 0:
                    for idx,i in enumerate(all):
                        all[idx]=i+nucl
                else:
                    pass
                    all.append(nucl)
            else:
                p =  degenerates[nucl]
                all = all*len(p)
                x= len(all)/len(p)
                app_nucl= [p[i//x] for i in range(len(p)*x)]
                for idx, i in enumerate(all):
                    all[idx]=i+app_nucl[idx]
        return set(all)
    
    def _get_seqs(self, path):
        
        with tempfile.NamedTemporaryFile(suffix='.fastq') as fastq:
            cmd=['zcat', path, '>', fastq.name]
            subprocess.check_call(' '.join(cmd), shell=True)
            fastq.flush()
            records=list(SeqIO.parse(open(fastq.name, 'r'), 'fastq'))
        return records

    
    def _get_ind(self, read, primerset):
        read=str(read)
        for primer in primerset:
            if primer in read:
                return "exact_hits", read.index(primer)+len(primer)
        else: # no break
            for primer in primerset:
                match=regex.findall("(%s){e<=5}" % (primer), read)
                if any(match):
                    return "fuzzy_hits", read.index(max(match, key=len))+len(primer)
            else:
                return "undetected", None
    

    def _print_summary(self, counter):
        print "seqs total:                         %i" % counter["total"]
        print "seqs with exact primer hit:         %i\t(%s%% total)" % (counter["exact_hits"], round(100*(float(counter["exact_hits"])/float(counter["total"])),2))
        print "seqs with fuzzy match to primer:    %i\t(%s%% total)" % (counter["fuzzy_hits"], round(100*(float(counter["fuzzy_hits"])/float(counter["total"])),2))
        print "seqs with undetected primer:        %i\t(%s%% total)" % (counter["undetected"], round(100*(float(counter["undetected"])/float(counter["total"])),2))


    def main(self, path_to_reads, primer, path_to_output):
        records=self._get_seqs(path_to_reads)
        primer_seq=self._translateDegenerate(primer)
        
        counter={
                 "total":len(records),
                 "exact_hits":0,
                 "fuzzy_hits":0,
                 "undetected":0
                 }
        
        with open(path_to_output, 'w') as out:
            for record in records:
                cat,ind=self._get_ind(record.seq, primer_seq)
                counter[cat]+=1
                if ind:
                    record=record[ind:]
                    SeqIO.write(record, out, "fastq")
        
        self._print_summary(counter)
        
primerClipper().main(sys.argv[1], sys.argv[2], sys.argv[3])
