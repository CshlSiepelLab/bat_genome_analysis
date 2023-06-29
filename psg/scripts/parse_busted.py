import sys
import json
with open(sys.argv[1]) as f:
    data = json.load(f)
codon = 0
for ER in data["Evidence Ratios"]["optimized null"][0]:
    #print(ER)
    codon += 1
    if int(ER) > 2:
        print(codon,ER)
    #if value["Uncorrected P-value"]:
    #    print(key,value["Uncorrected P-value"])
