import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="case ID")
args = parser.parse_args()

qc = args.filename
qc = open(qc, 'r+')
lines = qc.readlines()
print(len(lines))
qc.close()

output = "ccfCombined.txt"
if os.path.isfile(output):
    lines = lines[1:]
    with open(output, 'a') as target:
        for line in lines:
            target.write(str(line))
else:
    with open(output, 'w') as target:
        for line in lines:
            target.write(str(line))
