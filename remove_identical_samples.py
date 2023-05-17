#!/usr/bin/env python 

####Remove positions where both parent samples are "1/1"
import argparse

ap=argparse.ArgumentParser()
ap.add_argument('-i',type=str,required=True,help="Input VCF")
parse=ap.parse_args()

outfile=[]
removed=[]

for x in open(parse.i):
	if x.startswith("#"):
		x=x.strip()
		outfile.append(x)
		removed.append(x)
		continue
	x=x.strip().split("\t")
	if x[9][0:3] == "1/1" and x[10][0:3] == "1/1":
		removed.append("\t".join(x))
	else:
		outfile.append("\t".join(x))

#print(outfile)

with open("removed_data.vcf", 'w') as file:
	for x in removed:
		file.write(f"{x}\n")

with open("kept_data.vcf", 'w') as file:
	for x in outfile:
		file.write(f"{x}\n")