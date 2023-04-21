#!/usr/bin/env python

####Script to extract depth data from VCF file and create 5bp window around INDELs for filtering SNP call VCF

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
S1_hom_depth=[]
S1_het_depth=[]
S2_hom_depth=[]
S2_het_depth=[]
qual_data=[]
indel_data=[]

####Sep out depth data 
with open (snp_infile) as file:
	for x in file:
		if x.startswith("#"):
			continue
		x=x.replace("\n","")
		x=x.split(":")
		depth_data.append(x[5] + "\t" + x[8])
		if (x[3][2:]).strip() == "1/1":
			S1_hom_depth.append(x[5])
		else:
			if (x[3][2:]).strip() == "0/1":
				S1_het_depth.append(x[5])
		if (x[6][-3:]).strip() == "1/1":
			S2_hom_depth.append(x[8])
		else:
			if (x[6][-3:]).strip() == "0/1":
				S2_het_depth.append(x[8])
c=0
for lst in depth_data, S1_hom_depth, S1_het_depth, S2_hom_depth, S2_het_depth:
	df=pd.DataFrame([x.strip().split("\t") for x in lst])
	c+=1
	if c == 1:
		df.to_csv("depth_data.txt",header=False,index=False,sep="\t")
	elif c==2:
		df.to_csv("S1_hom_depth_data.txt",header=False,index=False,sep="\t")
	elif c==3:
		df.to_csv("S1_het_depth_data.txt",header=False,index=False,sep="\t")
	elif c==4:
		df.to_csv("S2_hom_depth_data.txt",header=False,index=False,sep="\t")
	elif c==5:
		df.to_csv("S2_het_depth_data.txt",header=False,index=False,sep="\t")
   
####Sep out qual data 
with open (snp_infile) as file:
	for x in file:
		if x.startswith("#"):
			continue
		x=x.replace("\n","")
		x=x.split("\t")
		qual_data.append(x[5])

df=pd.DataFrame([x.strip().split("\t") for x in qual_data])
df.to_csv("qual_data.txt",header=False,index=False,sep="\t")

####Create SNP exclusion ranges 
with open (indel_file) as file:
	for x in file:
		if x.startswith("#"):
			continue
		x=x.replace("\n","")
		x=x.split("\t")
		lower_value=int(x[1]) - int(5)
		if "," in x[4]:
			y=x[4].split(",")
			upper_val=int(len(max(y, key=len))) + int(x[1]) + int(5) 
		else:
			upper_val=int(len(x)) + int(x[1]) + int(5)
		indel_data.append(x[0]+ "\t" + str(lower_value) + "\t" + str(upper_val))

df=pd.DataFrame([x.strip().split("\t") for x in indel_data])
df.columns=["chromosome","start","end"]
df.to_csv("snp_excluson_data.txt",index=False,sep="\t")