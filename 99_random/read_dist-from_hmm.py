import os
import sys

hmm_output_lines = [x.strip().split() for x in open(sys.argv[1]) \
                    if x[0] != "#"]

hmm_length = int(hmm_output_lines[0][5])
bit_profile = [[] for x in range(0, hmm_length)]

for entry in hmm_output_lines:
  bit=float(entry[13])
  hmm_from = int(entry[15])
  hmm_to = int(entry[16])
  entry_span=range(min([hmm_from, hmm_to]),max([hmm_from, hmm_to]))

  for position in entry_span:
    bit_profile[position].append(bit)


print "HMM_position,bit_score_average,frequency"
for idx, pos in enumerate(bit_profile):
  if pos:
    pos_cov = len(pos)
    bit_avg = str(round(sum(pos)/pos_cov, 2))
    print "%i,%s,%i" %  (idx, bit_avg, pos_cov)
  else:
    print "%i,%i,%i" %  (idx, 0, 0)
