methyl = "wgEncodeHaibMethyl450CpgIslandDetails.txt"
original = open(methyl, 'r')
lines = original.readlines()

newlines = []
for line in lines:
	if ",cg" in line:
		newlines.append(line)
		
new = "wgEncodeHaibMethyl450CpgIslandDetailsUPDATED.txt"
new = open(new, 'w')
new.write(lines[0])
count = 0
for newline in newlines:
	newline = newline.split(",")
	for j in range(len(newline)):
		if newline[j] == "":
			newline[j] = "-"
	if count == 0:
		print newline[4]
	newline = ",".join(newline)		
	new.write(newline)
	count = count + 1
new.close()	