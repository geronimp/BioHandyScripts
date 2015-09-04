import os
import sys

hmm_output_files = [[x.strip().split() for x in open(y) if x[0] != "#"] for y in sys.argv[1].split(',')]
seen={}

for hmm_output_lines in hmm_output_files:
  for entry in hmm_output_lines:
    if entry[0] in seen:
      if float(seen[entry[0]][13]) < float(entry[13]):
        seen[entry[0]] = entry
    else:
      seen[entry[0]] = entry
hmm_output_lines = seen.values()
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
