#!/usr/bin/env python 

####Import functions 
import argparse 

####Parse arguments 
ap=argparse.ArgumentParser()
ap.add_argument('--dart',required=True,type=str,help="dartseq SNP file")
ap.add_argument('--seq',required=True,type=str,help="SNP called from sequencing")
parse=ap.parse_args()

####Create sets 
dart_data=set()
vcf_data=set()

####Open and extract data 
with open(parse.dart) as f:
	for x in f:
		if x.startswith("#"):
			continue
		x=x.strip().split("\t")
		dart_data.add(str(x[0])+"_"+str(x[1]))
with open(parse.seq) as f:
	for x in f:
		if x.startswith("#"):
			continue
		x=x.strip().split("\t")
		vcf_data.add(str(x[0])+"_"+str(x[1]))

####Find common data 
common_data=set(dart_data.intersection(vcf_data))

#print(len(dart_data)) #10116
#print(len(vcf_data)) #74034
#print(len(common_data)) #2663

####Subset SNP files based on common variants 
def subset_variants (infile):
	data=[]
	with open(infile) as f:
		for x in f:
			if x.startswith("#"):
				data.append(x)
				continue
			col=x.split("\t")[0]+"_"+x.split("\t")[1]
			if col in common_data:
				data.append(x)
		return (data)

vcf_list=subset_variants(parse.seq) 
dart_list=subset_variants(parse.dart) 

####Write files 
with open("Kora_Asgori_common_SNP.vcf", 'w') as f:
	for x in vcf_list:
		f.write(f"{x}")

with open("Dart_common_SNP.vcf", 'w') as f:
	for x in dart_list:
		f.write(f"{x}")
