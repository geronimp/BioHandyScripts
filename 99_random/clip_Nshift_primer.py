## USAGE::
## zcat <FILE> | python clip_Nshift_primer.py <PRIMER_SEQUENCE> <OUT_FILE> <STDIN>
# Created by Joel Boyd, 23/06/2015
# Australian Centre for ecogenomics


import sys
import regex
def translateDegenerate(primer):
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
    return set(all), len(primer)

primer_set, primer_length=translateDegenerate(sys.argv[1])
first_round_misses = {}
fastq_line = 0
misses     = 0
fuzzy_hits = 0
undetected = 0
exact      = 0
total      = 0
with open(sys.argv[2], 'w') as out:
    for line in sys.stdin:
        if fastq_line == 0:
            header = line
            fastq_line+=1
            total +=1
            continue
        elif fastq_line == 1:
            seq = line
            fastq_line+=1
            continue
        elif fastq_line == 2:
            spacer=line
            fastq_line+=1
            continue
        elif fastq_line == 3:
            qual=line
            hit=False
            for primer in primer_set:
                if primer in seq:
                    ind = seq.index(primer)
                    ind=primer_length+ind
                    out.write(header.strip())
                    out.write(seq[ind:].strip())
                    out.write(spacer.strip())
                    out.write(qual[ind:].strip())
                    hit=True
            if hit==True:
                fastq_line=0
                exact+=1
                continue
            else:
                first_round_misses[header]=[seq,spacer,qual]
                misses+=1
                fastq_line=0
                continue
    if first_round_misses:
        for header, item in first_round_misses.iteritems():
            for primer in primer_set:
                hit=False
                fuzzy_match=regex.findall("(%s){e<=5}" % (primer), item[0])
                if len(fuzzy_match) ==1:
                    ind = str(item[0]).index(fuzzy_match[0])
                    ind=primer_length+ind
                    out.write(header)
                    out.write(item[0][ind:])
                    out.write(item[1])
                    out.write(item[2][ind:])
                    hit=True
                if hit==True:
                    fuzzy_hits+=1
                    break
                else:
                    continue
            if hit==False:
                undetected+=1
                continue

print "seqs total:                         %i" % total
print "seqs with exact primer hit:         %i\t(%s%% total)" % (exact, round(100*(float(exact)/float(total)),2))
print "seqs missing exact hit of primer:   %i\t(%s%% total)" % (misses, round(100*(float(misses)/float(total)),2))
print "seqs with fuzzy match to primer:    %i\t(%s%% total)" % (fuzzy_hits, round(100*(float(fuzzy_hits)/float(total)),2))
print "seqs with undetected primer:        %i\t(%s%% total)" % (undetected, round(100*(float(undetected)/float(total)),2))
