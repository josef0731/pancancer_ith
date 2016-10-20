from collections import Counter
import os

plist = []
for file in os.listdir(os.getcwd()):
	if file.endswith(".tex"):
		plist.append(file)

filename = "consensus.txt"
target = open(filename, 'w')
for j in range(len(plist)):
	with open(plist[j]) as f:
		line = f.readlines()
	repeat = plist[j].split("_")
	prob = repeat[2].split(".t")
	prob = prob[0]
	prob = float(prob) / 0.002
	for i in range(int(prob)):
		target.write(repeat[1] + "\t" + str(int(prob)) + "\t" + line[10])
target.close()

target = open(filename, 'r+')
structures = target.readlines()
for y in range(len(structures)):
	structures[y] = structures[y].split("\t")[2]
count = Counter(structures)
target.write(str(count))
top3 = count.most_common()
coverage = 0
n_tree = 0
if len(count) > 3:
	for k in range(3):
		coverage += top3[k][1]
	if coverage >= 350:
		n_tree = 3
	elif len(count) >= 6:
		for m in range(3,6):
			coverage += top3[m][1]
			if coverage >= 350:
				n_tree = m
		if n_tree == 0:
			n_tree = 6
			for n in range(6,len(count)):
				if top3[n][1] == top3[5][1]:
					n_tree += 1
	else:
		for m in range(3,len(count)):
			coverage += top3[m][1]
			if coverage >= 350:
				n_tree = m
		if n_tree == 0:
			n_tree = len(count)
else:
	n_tree = len(count)

merge_trees = [0]*n_tree
merge_count = [0]*n_tree
top_trees = [0]*n_tree
repeat_count = [0]*n_tree
parents = [0]*n_tree

target = open(filename, 'r+')
line = target.readlines()

#scan for top tree structure
top = str(count)
i = 11
j = 0
comma = 0
print(top)
for merge in range(n_tree):
	while(top[j] != ";"):
		j = j + 1
	while(comma < (len(top) - 2) and top[comma] != ","):
		comma = comma + 1
	merge_count[merge] = top[j+6 : comma] 
	merge_trees[merge] = top[i : j+1]

	#parsing tree structure
	start = 0
	previous = 0
	clone = 0
	parent = []
	tree = merge_trees[merge]
	child = 0
	while(child != -1):
		child = tree.find("child", start)
		if child != -1:
			count = previous + tree.count('{', start, child) - tree.count('}',start, child) #total balance = 0 then parent = 0
			parent.append(count)
			if previous != 0:
				previous = count
			else:
				previous = clone
			clone = clone + 1
		start = child + 1

	parents[merge] = parent

	merge_count[merge] = int(merge_count[merge])
	top_trees[merge] = [0]*merge_count[merge]
	repeat_count[merge] = [0]*merge_count[merge]
	i = comma + 4
	comma = comma + 1
	j = i
	cnt = 0
	for x in range(len(line) - 1):
		split = line[x].split("\t")
		if split[2][:-1] == merge_trees[merge]:
			top_trees[merge][cnt] = split[0]
			repeat_count[merge][cnt] = split[1]
			cnt = cnt + 1
target.close()

newfile = "top_tree_stat.txt"	
target = open(newfile, 'w+')
target.write("rank" + "\t" + "tree#" + "\t" + "clone#" + "\t" + "#ssm" + "\t" + "#cnv" + "\t" + "freq" + "\t" + "parent" + "\n")

for merge in range(n_tree):
	target_lines = [0]*merge_count[merge]
	
	for i in range(merge_count[merge]):
		k = 0
		with open("tree_" + top_trees[merge][i] + "_" + str(float(repeat_count[merge][i]) * 0.002) + ".tex") as f:
			line = f.readlines()
		f.close()
		while(line[21 + k] != "\\hline\n"):								# copy all subclones
			k = k + 1
		target_lines[i] = [0]*k
		for y in range(k): 
			target_lines[i][y] = line[21 + y].split(" &")
			target_lines[i][y][1] = target_lines[i][y][1][1:]
			target_lines[i][y][2] = target_lines[i][y][2][1:]
			target_lines[i][y][3] = target_lines[i][y][3][0:5]
			if target_lines[i][y][3][4] == "$": target_lines[i][y][3] = target_lines[i][y][3][:-2]
			target.write(str(merge) + "\t" + str(top_trees[merge][i]) + "\t" + target_lines[i][y][0] + "\t" + target_lines[i][y][1] + "\t" + target_lines[i][y][2] + "\t" + target_lines[i][y][3] + "\t" + str(parents[merge][y]) + "\n")
target.close()