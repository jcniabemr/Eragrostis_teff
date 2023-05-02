#!/usr/bin/env python 

####Import funcions 
import argparse 

####Parse args 
ap=argparse.ArgumentParser()
ap.add_argument('--odart',required=True,type=str,help="dart vcf")
ap.add_argument('--dart',required=True,type=str,help="origional dart file")
parse=ap.parse_args()

####Create lists 
pos=[]
complemented_out=[]

####Possible complements 
complements={
"A":"T",
"T":"A",
"G":"C",
"C":"G"
}

####Open and extract
with open(parse.odart) as f:
	for x in f:
		if x.startswith("*") or x.startswith("AlleleID"):
			continue
		x=x.strip().split(",")
		if x[11] == "Minus":
			pos.append(x[6]+"_"+x[8])

with open(parse.dart) as f:
	for x in f:
		if x.startswith("#"):
			x=x.strip()
			complemented_out.append(x)
			continue
		x=x.strip().split("\t")
		if x[0]+"_"+x[1] in pos:
			REF=x[3]
			ALT=x[4]
			x[3]=x[3].replace(x[3], complements[REF])
			x[4]=x[4].replace(x[4], complements[ALT])
		else:
			pass
		complemented_out.append("\t".join(map(str,x)))


####Write file 
with open("complment_Dart_common_SNP.vcf", 'w') as f:
	for x in complemented_out:
		f.write(f"{x}\n")