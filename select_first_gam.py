import json
import sys

#Read a gam json from stdin and write the first occurrence of each read to stdout

seen_readnames = set()

for line in sys.stdin:
    name = json.loads(line)["name"]
    if name not in seen_readnames: 
        sys.stdout.write(line)
        seen_readnames.add(name)
