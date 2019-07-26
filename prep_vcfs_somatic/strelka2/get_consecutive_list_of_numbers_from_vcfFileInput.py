#!/usr/bin/env python3

from itertools import groupby
from operator import itemgetter
from os import path
from sys import exit,argv
import VCF
import cyvcf2

try:
	## MUST BE A VCF file
	finput = argv[1]
except IndexError as ie:
	exit("{}\nUSAGE: $0 $vcf_file ".format(ie))

if not path.exists(finput):
	msg = "ERROR: FNF {}".format(finput)
	raise IOError(msg)

d = {}

for v in VCF.lines(finput):
	if v['CHROM'] in d:
		d[v['CHROM']].append(v['POS'])
	else:
		d[v['CHROM']] = [v['POS']]

with open("{}.consPos.txt".format(finput), 'wt') as of:
	for key,val in d.items():
		## make sure all positions are integer ; if not raise error
		try:
			data = [int(i) for i in val]
		except ValueError as ve:
			exit("ERROR: {}".format(e))
		# https://stackoverflow.com/questions/2361945/detecting-consecutive-integers-in-a-list
		for k, g in groupby(enumerate(data), lambda ix: ix[0] - ix[1]):
			cn = list(map(itemgetter(1), g))
			if len(cn)>1:
				for p in cn:
					of.write("{}\t{}\n".format(str(key), p))
of.close()
exit()
