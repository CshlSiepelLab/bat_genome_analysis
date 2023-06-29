import sys
import json
with open(sys.argv[1]) as f:
    data = json.load(f)
for key,value in data["branch attributes"]["0"].items():
    if value["Uncorrected P-value"]:
        print(key,value["Uncorrected P-value"])
