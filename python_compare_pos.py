#!/usr/bin/env python

####Compare the positions of variants in different files 

import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--dart1",required=True,type=str,help="Dart file 1")
ap.add_argument("--dart2",type=str,help="Dart file 2")
ap.add_argument("--vcf",type=str,help="vcf")
parse=ap.parse_args()

####Initialise sets
newData=set()
existing=set()

for x in open(parse.dart1):
	if x.startswith("*"):
		continue
	x=x.strip().split(",")
	newData.add("_".join([x[4],x[5],x[0].split(":")[-1].replace(">","_")]))

####DArTseq vs DArTseq
for x in open(parse.dart2):
	if x.startswith("*"):
		continue
	x=x.strip().split(",")
	existing.add("_".join([x[6],x[7],x[0].split(":")[-1].replace(">","_")]))

# ####DArTseq vs VCF
# # for x in open(parse.vcf):
# # 	if x.startswith("#"):
# # 		continue 
# # 	x=x.strip().split("\t")
# # 	existing.add("_".join([x[0],x[1]]))

#print(len(newData))
#print(len(existing))
#print(len(newData.intersection(existing)))
matching=newData.intersection(existing)

data1=open("file_1_results.txt", 'w')
for x in open(parse.dart1):
	if x.startswith("*"):
		continue 
	x=x.strip().split(",")
	if "_".join([x[4],x[5],x[0].split(":")[-1].replace(">","_")]) in matching:
		data1.write(",".join(x)+"\n")
data1.close()

data2=open("file_2_results.txt", 'w')
for x in open(parse.dart2):
	if x.startswith("*"):
		continue 
	x=x.strip().split(",")
	if "_".join([x[6],x[7],x[0].split(":")[-1].replace(">","_")]) in matching:
		data2.write(",".join(x)+"\n")
data2.close()