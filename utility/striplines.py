#!/usr/bin/python 

import sys

if len(sys.argv) > 1:
	fin=open(sys.argv[1],'r')
else:
	fin=sys.stdin

linee=[l.strip() for l in fin]
output="\n".join(linee)
fin.close()

if len(sys.argv) > 2:
	fout=open(sys.argv[2],'w')
else:
	fout=sys.stdout

fout.write(output)
fout.close()
