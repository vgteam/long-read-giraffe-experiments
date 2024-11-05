import sys

#Read a gaf from stdin and write the best occurrence of each read to stdout
#Best is either the best mapq, or the longest, or the first

prev_name = ""
best_mapq=0
best_length=0
best_line = ""

for line in sys.stdin:
    l = line.split()
    name = l[0]
    mapq = int(l[11])
    length = int(l[3]) - int(l[2])
    if name != prev_name:
        if best_line != "":
            sys.stdout.write(best_line)

        prev_name = name
        best_mapq = mapq
        best_length = length
        best_line = line
    else:
        if mapq > best_mapq or (mapq == best_mapq and length > best_length):
            best_mapq = mapq
            best_length = length
            best_line = line 
        
sys.stdout.write(best_line)

