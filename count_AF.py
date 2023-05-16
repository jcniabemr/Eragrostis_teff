#!/usr/bin/env python 

####Count allele frequency 

import argparse 
ap=argparse.ArgumentParser()
ap.add_argument('-i',required=True,type=str,help='Input variant file')
parse=ap.parse_args()

count_data=[]

for x in open(parse.i):
	if x.startswith("#"):
		continue 
	x=x.strip().split("\t")
	pos="_".join([x[0],x[1]])
	count0=0
	count1=0
####Number of samples 
	for y in range(9,300):
		v1=x[y].split("/")[0]
		v2=x[y].split("/")[1].split(":")[0]
		if int(v1)==0:
			count0+=1
		else:
			count1+=1
		if int(v2)==0:
			count0+=1
		else:
			count1+=1
	count_data.append("\t".join([pos,str(count0),str(count1)]))

#print(count_data)

with open("AF_results.txt", 'w') as file:
	for x in count_data:
		file.write(f"{x}\n")