#!/usr/bin/python
import re
import sys
from itertools import izip

HA=open("usa_HA.fasta","w")
NA=open("usa_NA.fasta","w")
VAC=open("vaccine_data.txt","w");

data_americana=re.compile(r"(\d{2})/(\d{2})/(\d{4})");
data_giusta=re.compile(r"(\d{4})/(\d{2})/(\d{2})");
data_incompleta=re.compile(r"(\d{4})");

def correct_date(data):
	if re.match(data_giusta,data):
		corretta=data
	elif len(data)==4:
		corretta=data+"/01/01"
	else:
		m=re.match(data_americana,data)
		if m:
			corretta="%s/%s/%s" % (m.group(3),m.group(2),m.group(1))
		else:
			print "Very wrong date %s" % data
	return(corretta)


righe=[]
try:
	filename=sys.argv[1]
except:
	sys.stderr.write("Usage : python %s inputfile\n" % sys.argv[0])
	raise SystemExit(1)


#for line in open(filename):
#	stripped=line.strip()
#	if not stripped:
#		continue
#	if (stripped[0]=='>'):
#		righe.append([stripped])
#	else:
#		righe[-1].append(stripped)
#
#righe=map(lambda l: (l[0],''.join(l[1:])),righe)
f=open(filename);
pulite=[line.strip() for line in f]
righe = zip(pulite[::2],pulite[1::2])

print "Read %d sequences" % len(righe)

arch={}
regexp=re.compile(r">(?P<date>.+?)\|(?P<accession>.+?)" \
			r"\|(?P<strain>.+?)\|(?P<protein>.+?)" \
			r"\|(?P<country>.+?)(?P<other>\|.*)?$");

for header,sequence in righe:
	m=re.search(regexp,header).groupdict();
	strain=m["strain"];
	date=correct_date(m["date"]);
	if m["other"]:
		other=m["other"]
	else:
		other=""
	prot=m["protein"];
	new={"strain": strain, prot: sequence, "date": date, "other": other};
	if strain in arch:
		arch[strain].update(new)
	else:
		arch[strain]=new

strains=arch.keys()
strains.sort(key=lambda strain: arch[strain]["date"])

already_printed=set('');
printed=skipped=0;

for k in strains:
	if not (("NA" in arch[k]) and ("HA" in arch[k])):
		print "Missing a sequence of 2"
		continue
	seq_na=arch[k]["NA"]
	seq_ha=arch[k]["HA"]
	if (seq_na in already_printed) and (seq_ha in already_printed):
		skipped+=1
		continue
	already_printed.add(seq_na);
	already_printed.add(seq_ha);
	#print arch[k]['date']
	print >>NA, ">%s|%s|NA%s\n%s" % (arch[k]["date"],arch[k]["strain"], \
			arch[k]["other"],seq_na) ;
	print >>HA, ">%s|%s|HA%s\n%s" % (arch[k]["date"],arch[k]["strain"], \
			arch[k]["other"],seq_ha) ;
	printed+=1;

print "Printed %d, skipped %d/%d (%.1f%%)" % (printed,skipped, len(strains),skipped*100.0/len(strains))
