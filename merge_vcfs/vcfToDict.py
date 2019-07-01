#!/usr/bin/env python3

### vcfMerger2
###
### MIT License
###
### Copyright (c) 2018 Translational Genomics Research Institute
###
### Permission is hereby granted, free of charge, to any person obtaining a copy
### of this software and associated documentation files (the "Software"), to deal
### in the Software without restriction, including without limitation the rights
### to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
### copies of the Software, and to permit persons to whom the Software is
### furnished to do so, subject to the following conditions:
###
### The above copyright notice and this permission notice shall be included in all
### copies or substantial portions of the Software.
###
### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
### AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
### OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
### SOFTWARE.
###
### Major Contributors: Christophe Legendre'a0
### Minor Contributors:



import re
import cyvcf2
from collections import defaultdict
import logging as log


class vcfToDict:
	""" manage input VCF ; read the VCF into Dictionary and process the data to make it ready for the merging"""

	def __init__(self, fvcf, toolname):

		self.fvcf = fvcf
		self.toolname = toolname

		self.headers,self.header_minus_chrom_line,self.header_chrom_line,self.header_other_info,self.contigs,\
		self.header_filters,self.header_info,self.header_format = self.getHeaderOnly()

		self.samplenames = self.getSampleNames()
		## we init the dictionary that will contain the loci to be merged with the others dictionaries from the other vcfs
		self.ds = {} ;

	def getHeaderChromLine(self):
		log.info("capturing header LINE for " + self.toolname)
		with(open(self.fvcf,'r')) as f:
			for line in f:
				if line.startswith("#CHROM"):
					return line

	def getHeaderOnly(self):
		"""
		only get the header lines from the vcf to check first the presence of the FLAGS that should be common
		between the VCF files ; All the VCFs files should be up to sVCF specifications ;
		if not we exit the program and ask user to update the VCF to specs compatible with this program ;
		program will return a value (not decided yet on which one) that tells us if compatible
		"""
		log.info("capturing headers for " + self.toolname)
		headers = ""
		header_other_info = []
		header_filters = []
		header_info = []
		header_format = []
		contigs = [] ;
		header_minus_chrom_line = "" ;

		regex = re.compile(r",assembly.*>$", re.IGNORECASE)
		with(open(self.fvcf,'r')) as f:
			for line in f:
				line = line.strip()
				headers = ''.join([ headers,line ])
				# log.info(headers.splitlines()[0:5])
				if not line.startswith("#CHROM") and line.startswith("##"): ## not really useful and not use later ;
					header_minus_chrom_line = ''.join([ header_minus_chrom_line,line ])
				if line.startswith("#CHROM"):
					header_chrom_line = line
					continue
				if line.startswith("##contig"):
					line = regex.sub(">", line) ## we remove any character after the string ",assembly" including assembly ; due to lancet ; and to avoid contigs duplicates
					contigs.append(line)
					continue
				if line.startswith("##FILTER"): header_filters.append(line) ; continue
				if line.startswith("##INFO"):
					header_info.append(line)
					continue
				if line.startswith("##FORMAT"): header_format.append(line) ; continue
				if line.startswith("##"): header_other_info.append(line) ; continue
				## we assume that the VCF is up to specs, and reaching a line without ## at the beginning means we have reached either '##CHROM' line  or  a variant record
				if not line.startswith("##"): break ;

		return headers,\
		       header_minus_chrom_line,\
		       header_chrom_line,\
		       header_other_info,\
		       contigs,\
		       header_filters,\
		       header_info,\
		       header_format

	def getSampleNames(self):
		''' By default in the VCF specs, if a FORMAT column exists, it has to be the 9th column; 
		Then starting from the 10th column, it is sample information;
		As we deal with 0-base info in python we use [9:] here
		'''
		return tuple(self.header_chrom_line.strip().split('\t')[9:])
		

	def readVCF(self):
		"""
		read the whole vcf file into memory and create the dictionnary that contains the mutations loci
		normally somatic calls are smaller fiels compared to germline calls and therefore can be hold in memory
		even if you have up to 10 vcfs with 1 million lines each. This will of course require up to 16GB
		"""
		return(cyvcf2.VCF(self.fvcf))

	def dictOfLoci(self, ovcf):
		"""
			we populate the dictionary with every variant record using current tool.
		"""
		log.info("creating variant dictionary for " + self.toolname)
		for variant in ovcf:
			self.ds[str(variant.CHROM)+"__"+str(variant.POS)+"_"+str(variant.REF)+"_"+str("".join(x for x in variant.ALT)) ] = [self.toolname, variant ]
		return(self.ds)

	def workflow(self):
		log.info(self.getHeaderOnly())


class Dict(defaultdict):
	''' Class that defines MultiDimensional Dictionary such as: myDict[A][B][C][D] = value  '''
	def __init__(self):
		defaultdict.__init__(self, Dict)
	def __repr__(self):
		return dict.__repr__(self)

