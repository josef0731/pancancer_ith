import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("filename", help="filename (TCGA level 3 methylation file)")
args = parser.parse_args()

methyl = args.filename
original = open(methyl, 'r')
lines = original.readlines()
lines = lines[1:]

print "Reading methylation data ..."

beta = []
IllumID = []	#col0

for line in lines:
	line = line.split("\t")
	IllumID.append(line[0])
	beta.append(line[1])

ID = []	#col0
infinium = []	#col6
strands = []	#col16
build = []	#col10
chr = []	#col11
posit = []	#col12
name = []	#col21
accession = []	#col22
group = []	#col23
islandName = []	#col24
relation = []	#col25

print "Methylation data read. Now read array master file..."

standard = "wgEncodeHaibMethyl450CpgIslandDetailsUPDATED.txt"
standard = open(standard, 'r')
standard = standard.readlines()
for entry in standard:
	entry = entry.split(",")
	ID.append(entry[0])
	infinium.append(entry[6])	#col6
	strands.append(entry[16])
	build.append(entry[10])	#col10
	chr.append(entry[11])	#col11
	posit.append(entry[12])	#col12
	name.append(entry[21])	#col21
	accession.append(entry[22])	#col22
	group.append(entry[23])	#col23
	islandName.append(entry[24])	#col24
	relation.append(entry[25])	#col25

illuminaID = []
infiniumType = []
strand = []
genomeBuild = []
chrom = []
position = []
UCSCgeneName = []
UCSCgeneAccession = []
UCSCgeneGroup = []
UCSCislandName = []
relationToUCSCisland = []
	
print "Now do the matching ..."
	
for count in range(len(lines)):
	if IllumID[count] in ID:
		index = ID.index(IllumID[count])
		illuminaID.append(IllumID[count])
		infiniumType.append(infinium[index])
		strand.append(strands[index])
		genomeBuild.append(build[index])
		chrom.append(chr[index])
		position.append(posit[index])
		UCSCgeneName.append(name[index])
		UCSCgeneAccession.append(accession[index])
		UCSCgeneGroup.append(group[index])
		UCSCislandName.append(islandName[index])
		relationToUCSCisland.append(relation[index])

if not os.path.exists("TCGA_COAD_methyl"):
	os.makedirs("TCGA_COAD_methyl")

print "Matching done. Writing to file..."
	
new = "TCGA_COAD_methyl/" + methyl
new = open(new, 'w')
new.write("IllumID\tinfiniumDesignType\tgenomeBuild\tchrom\tposition\tUCSCgeneName\tUCSCgeneAccession\tUCSCgeneGroup\tUCSCislandName\trelationToUCSCisland\n")
for count in range(len(illuminaID)):
	newline = str(illuminaID[count]) + "\t" + str(strand[count]) + "\t" + str(infiniumType[count]) + "\t" +  str(genomeBuild[count]) + "\t" + str(chrom[count]) + "\t" + str(position[count]) + "\t" + str(UCSCgeneName[count]) + "\t" + str(UCSCgeneAccession[count]) + "\t" + str(UCSCgeneGroup[count]) + "\t" + str(UCSCislandName[count]) + "\t" + str(relationToUCSCisland[count]) + "\n"
	new.write(newline)
new.close()

print "All done."	