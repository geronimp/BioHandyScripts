#!/usr/bin/env python


## USAGE::
## clip_Nshift_primer.py <IN> <PRIMER> <OUT>
# Created by Joel Boyd, 23/06/2015
# Australian Centre for ecogenomics
import difflib
import sys
import regex
import subprocess
import tempfile
import argparse
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
            records=SeqIO.to_dict(SeqIO.parse(open(fastq.name, 'r'), 'fastq'))
        return records

    
    def _get_ind(self, records, primerset, fuzzy):
        n_no_ht = 0
        n_ht    = 0
        fz_ht   = 0
        no_ht   = []
        ht      = []
        for key, record in records.iteritems():
            read=str(record.seq)
            h=False
            for primer in primerset:
                if primer in read:
                    ht.append(record[read.index(primer)+len(primer):])
                    n_ht+=1
                    h=True
                else:
                    continue
            if h == False:
                no_ht.append(record)
        if fuzzy:
            for primer in primerset:
                a={record.name: regex.findall("(%s){e<=5}" % (primer), str(record.seq)) for record in no_ht}
                a={key:item for key, item in a.iteritems() if item}
                
                for key,item in a.iteritems():
                    if len(item)>1:    
                        best_match=[difflib.SequenceMatcher(None, x, primer).ratio() for x in item]
                        best=item[ best_match.index(max(best_match)) ] 
                        ind = str(records[key].seq).index(best) + len(best)
                    else:
                        best=item[0]
                        ind = str(records[key].seq).index(best) + len(best)
                        
                    record=records[key]
                    no_ht.remove(record)
                    ht.append(record[ind:])
                    fz_ht+=1
            
            n_no_ht=len(no_ht)
            return ht, n_ht, fz_ht, n_no_ht
           
        else:
            n_no_ht=len(no_ht)
            return ht, n_ht, fz_ht, n_no_ht

    def _print_summary(self, counter, fuzzy):
        print "seqs total:                         %i" % counter["total"]
        print "seqs with exact primer hit:         %i\t(%s%% total)" % (counter["exact_hits"], round(100*(float(counter["exact_hits"])/float(counter["total"])),2))
        if fuzzy:
            print "seqs with fuzzy match to primer:    %i\t(%s%% total)" % (counter["fuzzy_hits"], round(100*(float(counter["fuzzy_hits"])/float(counter["total"])),2))
        else:
            print "seqs with fuzzy match to primer:    N/A" 
        print "seqs with undetected primer:        %i\t(%s%% total)" % (counter["undetected"], round(100*(float(counter["undetected"])/float(counter["total"])),2))


    def main(self, path_to_reads, primer, path_to_output, fuzzy):
        records=self._get_seqs(path_to_reads)
        primer_seq=self._translateDegenerate(primer)
        
        counter={
                 "total"      :len(records),
                 "exact_hits" :0,
                 "fuzzy_hits" :0,
                 "undetected" :0
                 }
        
        with open(path_to_output, 'w') as out:
            
            hits, numhits, fuzzynum, miss =self._get_ind(records, primer_seq, fuzzy)
            counter["exact_hits"]=numhits
            counter["fuzzy_hits"]=fuzzynum
            counter["undetected"]=miss
            
            SeqIO.write(hits, out, "fastq")
        
        self._print_summary(counter, fuzzy)
        
parser = argparse.ArgumentParser(description='''clip N shift primers''')
parser.add_argument('--input', type=str, help='input sequences', required=True)
parser.add_argument('--output', type=str, help='output sequences', required=True)
parser.add_argument('--primer', type=str, help='primer sequence', required=True)
parser.add_argument('--accept_fuzzy', help='Accept fuzzy hits', action="store_true", default=False)
args = parser.parse_args()

primerClipper().main(args.input, args.primer, args.output, args.accept_fuzzy)
