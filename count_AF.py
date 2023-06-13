#!/usr/bin/env python 

####Count allele frequency 

import argparse 
ap=argparse.ArgumentParser()
ap.add_argument('-i',required=True,type=str,help='Input variant file')
parse=ap.parse_args()

count_data=[]
remove_data=[]
new_vcf=[]


with open(parse.i) as f:
	for x in f:
		if x.startswith("#"):
			continue 
		x=x.strip().split("\t")
		pos="_".join([x[0],x[1]])
		count0=0
		count1=0
		for y in range(9,len(x)):
			if x[y][0:3] != "./.":
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
		count_data.append("\t".join([pos,str(float(count0/(count0+count1)*100)),str(float(count1/(count0+count1)*100))]))

#print(count_data)

with open("AF_results.txt", 'w') as file:
	for x in count_data:
		file.write(f"{x}\n")

for x in count_data:
 	x=x.split("\t")
 	if float(x[1]) < 10 or float(x[2]) < 10:
 		remove_data.append(x[0])

with open("remove_data_list.txt", 'w') as file:
	for x in remove_data:
		file.write(f"{x}\n")

with open(parse.i) as f:
	for x in f:
		if x.startswith("#"):
			x=x.strip()
			new_vcf.append(x)
			continue
		x=x.strip().split("\t")
		if "_".join([x[0],x[1]]) not in remove_data:
			new_vcf.append("\t".join(x))

with open("filtered_variants.vcf",'w') as f:
	for x in new_vcf:
		f.write(f"{x}\n") 