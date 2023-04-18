#!/usr/bin/env python

####Script to extract info from SNP call vcf

####Import fuctions 
import argparse 
import pandas as pd 

####Parse files 
ap=argparse.ArgumentParser()
ap.add_argument('-snp',type=str,required=True,help="snpinfile")
ap.add_argument('-indel',type=str,required=True,help="indel file")
args=ap.parse_args()

snp_infile=(args.snp)
indel_file=(args.indel)

####Create required lists 
depth_data=[]
indel_data=[]

####Sep out depth data 
with open (snp_infile) as file:
	for x in file:
		if x.startswith("#"):
			continue
		x=x.replace("\n","")
		x=x.split(":")
		depth_data.append(x[5] + "\t" + x[8])
	
df=pd.DataFrame([x.strip().split("\t") for x in depth_data])
df.to_csv("depth_data.txt",header=False, index=False, sep="\t")

####Create SNP exclusion ranges 
with open (indel_file) as file:
	for x in file:
		if x.startswith("#"):
			continue
		x=x.replace("\n","")
		x=x.split("\t")
		lower_value=int(x[1])-int(5)
		upper_val=len(x[3]) + int(x[1]) + int(5)
		indel_data.append(x[0]+ "\t" + str(lower_value) + "\t" + str(upper_val))

df=pd.DataFrame([x.strip().split("\t") for x in indel_data])
df.to_csv("snp_excluson_data.txt",header=False, index=False, sep="\t")
