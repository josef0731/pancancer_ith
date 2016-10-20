import os

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("case", help="filename (case)")
args = parser.parse_args()

snp = args.case + "parsed.maf"
file = open(snp, 'r+')
lines = file.readlines()
file.close()

newlines = []
newlines.append(lines[0])
lines = lines[1:]

for line in lines:
	line = line.split("\t")
	chrom = line[4]
	start = line[5]
	mafindex = chrom + "_" + start
	line = "\t".join(line)
	line = mafindex + "\t" + line
	newlines.append(line)

firehose = "/storage_disk1/home/jocfng/Protected_Mutations/COADsomaticMut/firehose/gdac.broadinstitute.org_COAD.Mutation_Packager_Oncotated_Raw_Calls.Level_3.2016012800.0.0/" + args.case + ".maf"
file = open(firehose, 'r+')
lines = file.readlines()
file.close()
lines = lines[4:]
firehoseindices = []

for line in lines:
	line = line.split("\t")
	chrom = line[4]
	start = line[5]
	firehoseindex = chrom + "_" + start
	firehoseindices.append(firehoseindex)

newnewline = []
newnewline.append(newlines[0])

for line in newlines:
	line = line.split("\t")
	mafindex = line[0]
	if mafindex in firehoseindices:
		line = "\t".join(line[1:])
		newnewline.append(line)

print(args.case)

output = args.case + "parsedTRIMMED.maf"
output = open(output, 'w')
for newline in newnewline:
	output.write(newline)
output.close()