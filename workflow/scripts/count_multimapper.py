import sys

counts = {}

for line in sys.stdin:
    read_id, trna = line.strip().split("\t")
    count = counts.get(read_id, 0)
    count.setdefault(trna, count + 1)
# TODO