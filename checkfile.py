#!/usr/bin/env python 

import sys
import itertools 

# f=open(sys.argv[1])

# for x in f:
# 	x=x.replace("\n","")
# 	if x[0] == "#":
# 		pass
# 	else:
# 		x=x.split("\t")
# 		print (x[3]+"_"+x[4])
# f.close()

total=0
matches=0

f1=open(sys.argv[1])
f2=open(sys.argv[2])

for x,y in zip(f1,f2):
	total+=1
	if x == y:
		matches+=1 

print (matches, total)

f1.close()
f2.close()