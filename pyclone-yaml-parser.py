import os
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("case", help="case ID i.e. pfgXXX")
args = parser.parse_args()

case = args.case

yaml = open('pyclone-config.yaml','r+')
lines = yaml.readlines()
yaml.close()

current_dir = os.getcwd()
lines[1] = "working_dir: /pathosy01/disk1/home/jocfng/STADCOADUCEC/HKU100/pyclone/" + args.case + "\n"
lines[55] = "  " + str(case) + ":\n"
lines[56] = "    mutations_file: " + str(case) + ".yaml\n"

filename = 'pyclone-config.yaml'

with open(filename,'w') as target:
	for z in range(len(lines)):
		target.write(lines[z]) 
