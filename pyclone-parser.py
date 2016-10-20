import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("case", help="case ID i.e. pfgXXX or TCGA-xx-xxxx")
args = parser.parse_args()

vcf = "APPENDED_" + args.case + "T.vcf"
vcf = open(vcf, 'r+')
lines = vcf.readlines()[1:]

newlines = []

mutation_id = []
ref_counts = []
var_counts = []
chr = []
start = []
end = []
status = []

retain = []
for i in range(22):
        retain.append(str(i + 1))
        
for i in range(len(lines)):
        if lines[i][0:1] != "##" and any(retained in lines[i].split("\t")[0] for retained in retain):
                newlines.append(lines[i])

for i in range(len(newlines)):
        newlines[i] = newlines[i].split("\t")
        newlines[i][10] = newlines[i][10].split(":")
        var_counts.append(newlines[i][10][1].split(",")[1])
        ref_counts.append(newlines[i][10][1].split(",")[0])
        chr.append(newlines[i][0])
        start.append(int(newlines[i][1]))
        end.append(int(newlines[i][1]))
        mutation_id.append(str(newlines[i][0]) + "_" + str(newlines[i][1]))
vcf.close()
                
cnv = str(args.case) + "_oncosnp.txt"
cnv = open(cnv, 'r+')
cnvlines = cnv.readlines()

major_cn = [0]*len(mutation_id)
minor_cn = [0]*len(mutation_id)

for i in range(len(cnvlines)):
        cnvlines[i] = cnvlines[i].split('\t')
                
for j in range(len(mutation_id)):
        for k in range(len(cnvlines)):
                if chr[j] == str(cnvlines[k][0]): 
                        if start[j] >= int(cnvlines[k][1]) and end[j] <= int(cnvlines[k][2]):
                                major_cn[j] = cnvlines[k][3]
                                minor_cn[j] = cnvlines[k][4]

if len(mutation_id) > 10000:
	sample = np.random.choice(len(mutation_id), 10000)
else:
	sample = range(len(mutation_id))

target = "pyclone.mutation.tsv"
with open(target,'w') as writing:
        writing.write("mutation_id" + "\t" + "ref_counts" + "\t" + "var_counts" + "\t" + "normal_cn" + "\t" + "minor_cn" + "\t" + "major_cn" + "\n")
        for j in sample:
                if major_cn[j] != 0 and minor_cn[j] != 0:
                        writing.write(str(mutation_id[j]) + "\t" + str(ref_counts[j]) + "\t" + str(var_counts[j]) + "\t2\t" + str(minor_cn[j]) + "\t" + str(major_cn[j]) + "\n") 
