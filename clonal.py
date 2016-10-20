import os
prefixed = [filename for filename in os.listdir('.') if filename.startswith("tree_0_") and filename.endswith(".tex")]

filename = str(prefixed[0])
target = open(filename, 'r+')
line = target.readlines()

j = 0
while(line[21 + j] != "\\hline\n"):
	j = j + 1
target_lines = [None]*j
freq = [None]*j

for y in range(j): 
	target_lines[y] = line[21 + y].split(" & ")
	target_lines[y][1] = target_lines[y][1].split(", ")
	freq[y] = target_lines[y][-1]
	freq[y] = freq[y][:-3]
target.close()

#freq = target_lines[-1]
#freq = freq[:-3]

newfile = "clone_member.txt"	
target = open(newfile, 'w+')
target.write("clone_num" + "\t" + "id" + "\t" + "frequency" + "\n")
for z in range(len(target_lines)):
        for w in range(len(target_lines[z][1])):
                target.write(target_lines[z][0] + "\t" + target_lines[z][1][w] + "\t" + freq[z] + "\n")
target.close()
