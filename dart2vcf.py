#!/usr/bin/env python 

####Convert dartfile to VCF format 
####Import functions 
import argparse 
import pandas as pd

####Set arguments
ap=argparse.ArgumentParser()
ap.add_argument('-i',type=str,required=True,help="dartseq file")
parse=ap.parse_args()

####Create lists
out_file_info=[]
position=[]

####Read in data 
df=pd.read_csv(parse.i,skiprows=6,low_memory=False)

####Drop cols not found in VCF
df = df.drop(columns=['CloneID','SnpPosition','AlleleSequenceRef','AlleleSequenceSnp','TrimmedSequenceRef','TrimmedSequenceSnp','ChromPosTag_Eragrostis_CogeV3','AlnCnt_Eragrostis_CogeV3','Strand_Eragrostis_CogeV3','CallRate','OneRatioRef','OneRatioSnp','FreqHomRef','FreqHomSnp','FreqHets','PICRef','PICSnp','AvgPIC','AvgCountRef','AvgCountSnp','RepAvg'])

####Duplicate this col for later 
df['second_snp'] = df.loc[:, 'SNP']
df['FILTER']='PASS'
df['INFO']='.'
df['FORMAT']='GT'

####Re-order cols 
cols = list(df.columns)
cols = cols[-4:] + cols[:-4]
df = df[cols]

def swap_columns(df,a,b,c,d,e,f,g,h,i):
    col_list = list(df.columns)
    r,s,t,u,v,w,x,y,z = col_list.index(a),col_list.index(b),col_list.index(c),col_list.index(d),col_list.index(e),col_list.index(f),col_list.index(g),col_list.index(h),col_list.index(i)
    col_list[r],col_list[s],col_list[t],col_list[u],col_list[v],col_list[w],col_list[x],col_list[y],col_list[z] = col_list[w],col_list[x],col_list[v],col_list[z],col_list[r],col_list[y],col_list[s],col_list[t],col_list[u]
    df = df[col_list]
    return df
df = swap_columns(df,'second_snp','FILTER','INFO','FORMAT','AlleleID','Chrom_Eragrostis_CogeV3','ChromPosSnp_Eragrostis_CogeV3','AlnEvalue_Eragrostis_CogeV3','SNP')

####Convert back to list for further editing 
reformed=df.T.reset_index().values.T.tolist()

####Create set of duplicated positions 
for x in reformed:
	position.append(str(x[0])+"_"+str(x[1]))
duplicates = set(x for x in position if position.count(x) > 1)

####Edit new reformed data to be consistent with VCF
for x in reformed:
	a="\t".join(map(str,x))
	if a.startswith("nan"): 
		continue 
	true_id=x[2].split("|")[0]
	x[2]=x[2].replace(x[2],true_id)
	true_ref=x[3].split(">")[0][-1]
	x[3]=x[3].replace(x[3],true_ref)
	true_alt=x[4].split(">")[-1]
	x[4]=x[4].replace(x[4],true_alt)
	if x[0] == "Chrom_Eragrostis_CogeV3":
		x[0]=x[0].replace(x[0],"#CHROM")
	if x[1] == "ChromPosSnp_Eragrostis_CogeV3":
		x[1]=x[1].replace(x[1],"POS")
	if x[2] == "AlleleID":
		x[2]=x[2].replace(x[2],"ID")
	if x[3] == "P":
		x[3]=x[3].replace(x[3],"REF")
	if x[4] == "second_snp":
		x[4]=x[4].replace(x[4],"ALT")
	if x[5] == "AlnEvalue_Eragrostis_CogeV3":
		x[5]=x[5].replace(x[5],"QUAL")
	x[5]=str(x[5])
	if x[5] != "QUAL":
		x[5]=x[5].replace(x[5], "999")	
	for i in range(9,253):
		if x[i] == "0":
			x[i]=x[i].replace(x[i],"0/0")
		elif x[i] == "1":
			x[i]=x[i].replace(x[i],"0/1")
		elif x[i] == "2":
			x[i]=x[i].replace(x[i],"1/1")
		elif x[i] == "-":
			x[i]=x[i].replace(x[i],"./.")
	if str(x[0])+"_"+str(x[1]) not in duplicates:
		out_file_info.append("\t".join(map(str,x)))
	
####Write file 
with open("dart_mapping.vcf", 'w') as f:
	for x in out_file_info:
		f.write(f"{x}\n")



