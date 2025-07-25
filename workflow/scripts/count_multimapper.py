from collections import Counter
import sys

counts = Counter()

old_read_id = ""
hits = 0
for line in sys.stdin:
  read_id, _ = line.strip().split("\t")
  if old_read_id == "":
    old_read_id = read_id
    hits = 1
  elif read_id != old_read_id:
    counts[hits] += 1
    old_read_id = read_id
    hits = 0
  hits += 1
if hits:
  counts[hits] += 1


"\t".join(["hits", "count"])
for hits, count in counts.items():
    "\t".join([hits, count])
