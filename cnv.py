cnv = open('cnv_data.txt','r+')
lines = cnv.readlines()
lines = lines[1:]

for i in range(len(lines)):
	lines[i] = lines[i].split("\t")
	lines[i][3] = lines[i][3].split(";")
	for j in range(len(lines[i][3])):
		lines[i][3][j] = lines[i][3][j].split(",")

cnv_id = []
cnv_region = []		
		
for i in range(len(lines)):
	pair = []
	ranges = []
	pair.append(lines[i][3][0][0])
	cnvs = len(lines[i][3])
	for j in range(1,cnvs):
		if int(lines[i][3][j][0][1:]) != int(lines[i][3][j - 1][0][1:]) + 1:		# test for consecutive overlapping SSMs
			pair.append(lines[i][3][j - 1][0])
			ranges.append(pair)
			cnv_id.append(lines[i][0])
			pair = []
			pair.append(lines[i][3][j][0])
	if ranges == []:								# for CNVs specifying only 1 region
		pair.append(lines[i][3][len(lines[i][3]) - 1][0])
		ranges.append(pair)
		cnv_id.append(lines[i][0])
	cnv_region.append(ranges)

cnv_region = sum(cnv_region,[])

cnv.close()

ssm = open('ssm_data.txt','r+')
lines = ssm.readlines()
lines = lines[1:]

for i in range(len(lines)):
	lines[i] = lines[i].split("\t")
	lines[i][1] = lines[i][1].split("_")

coordinate = []	
	
for x in range(len(cnv_region)):
	start = cnv_region[x][0]
	end = cnv_region[x][1]
	triplet = []							#chromsoome number + start + end of each CNV region
	for i in range(len(lines)):
		if lines[i][0] == start:
			triplet.append(lines[i][1][0])			# chromosome number
			triplet.append(lines[i][1][1])			# start position
		if lines[i][0] == end:
			triplet.append(lines[i][1][1])			# end position	
			coordinate.append(triplet)	
ssm.close()

cytoband = "../../cytoBand.txt"
cytoband = open(cytoband,'r+')
cytoband = cytoband.readlines()

for w in range(len(cytoband)):
        cytoband[w] = cytoband[w].split("\t")
        cytoband[w][0] = cytoband[w][0][3:]

filename = "cnv_coordinates.txt"

with open(filename,'w') as target:
	target.write("cnv_id" + "\t" + "chrs" + "\t" + "start" + "\t" +"end" + "\t" + "cytoband" + "\n")
	for z in range(len(coordinate)):
		target.write(cnv_id[z] + "\t" + coordinate[z][0] + "\t" + coordinate[z][1] + "\t" + coordinate[z][2] + "\t")
		for w in range(len(cytoband)):
			band = []
			if coordinate[z][0] == cytoband[w][0]:
				if int(coordinate[z][1]) >= int(cytoband[w][1]) and int(coordinate[z][1]) <= int(cytoband[w][2]):
					band.append(coordinate[z][0] + cytoband[w][3])
					target.write(str(band))
				if int(coordinate[z][1]) <= int(cytoband[w][1]) and int(coordinate[z][2]) >= int(cytoband[w][2]):
					band.append(coordinate[z][0] + cytoband[w][3])
					target.write(str(band))
				if int(coordinate[z][2]) >= int(cytoband[w][1]) and int(coordinate[z][2]) <= int(cytoband[w][2]):
					target.write("\n")	
					
