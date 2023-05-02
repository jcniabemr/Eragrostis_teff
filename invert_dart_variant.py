#!/usr/bin/env python 

####Invert variants in dartseq file based on orientation in the genome sequence SNPs
####Run using only presorted files eg bcftools sort, else ref and alt will be wrong order 

####Import funcions 
import argparse

####Parse files
ap=argparse.ArgumentParser()
ap.add_argument('--dart',type=str,required=True,help="dartseq file")
ap.add_argument('--vcf',type=str,required=True,help="sequencing SNP file")
parse=ap.parse_args()

####Create lists
ref_vars=[]
dart_data=[]
zipped_data=[]
outfile=[]

####Open files 
with open(parse.vcf) as vcf:
	for x in vcf:
		if x.startswith("#"):
			x=x.strip()
			continue
		x=x.strip().split("\t")
		ref_vars.append(x[3]+"_"+x[4])

with open(parse.dart) as dart:
	for x in dart:
		if x.startswith("#"):
			x=x.strip()
			outfile.append(x)
			continue
		dart_data.append(x)

####Zip ref SNPs and dart data together as a list using map 
data=list(map(list, zip(ref_vars, dart_data)))
for x in data:
	zipped_data.append("\t".join(map(str,x)))

for x in zipped_data:
	x=x.strip().split()
	if x[0][0] == x[4]:
		outfile.append("\t".join(map(str,x[1:])))
	else:
		x[4]=x[4].replace(x[4],x[0][0])
		x[5]=x[5].replace(x[5],x[0][2])
		for i in range(10,254):
			if x[i] == "0/0":
				x[i]=x[i].replace(x[i],"1/1")
		outfile.append("\t".join(map(str,x[1:]))) 

with open("inverted_final_dart_data2.vcf",'w') as f:
	for x in outfile:
		f.write(f"{x}\n")