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


from sys import exit
from sys import argv
from os import path
import getopt
from cyvcf2 import VCF, Writer
import numpy as np
import logging as log
import warnings
from math import isnan

global AR_threshold_for_GT
AR_threshold_for_GT = 0.90  ## value HARDCODED but dynamically modify it with --threshold_AR option

class Genotype(object):
    '''
    # genotypes = [Genotype(li) for li in variant.genotypes ]
	# print genotypes
	# which shows: [./., ./., ./., 1/1, 0/0, 0/0, 0/0]
    '''
    __slots__ = ('alleles', 'phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123456."[a] for a in self.alleles)
    __repr__ = __str__

class GenotypeInv(object):

	def __init__(self, li):

		self.allele1 = li[0]
		self.phased = bool(0) if li[1] == "/" else bool(1)
		self.allele2 = li[2]
		self.GT = []

	def get_gt_numpy_compatible(self):
		self.GT = [] ## we need to reinit the GT list here otherwise shared by all instances. Weird because we reinitiated it already in the _init_
		if self.phased:
			if self.allele1 != "." :
				self.GT.append((2*int(self.allele1)) + 3)
			else:
				self.GT.append(1)
			if self.allele2 != ".":
				self.GT.append((2*int(self.allele2)) + 3)
			else:
				self.GT.append(1)
		else:
			if self.allele1 != ".":
				self.GT.append((2*int(self.allele1)) + 2)
			else:
				self.GT.append(0)
			if self.allele2 != ".":
				self.GT.append((2*int(self.allele2)) + 2)
			else:
				self.GT.append(0)
		return self.GT



def usage(scriptname, opts):
	print("\nUSAGE: \npython3 " + scriptname + '  --help')
	print("python3 " + scriptname + " -i vardictJava.somatic.vcf  --tumor_column 11 --normal_column 10  -o "
	                                "anyUserGivenName_for_prepped_updated_vcf.vcf\n")
	print("#" * 40 + "\nWARNING WARNING: This script is to be used only with somatic calls in vcf having NORMAL "
	                 "sample in column 10 and TUMOR sample in column 11; if not like this, update your vcf file to these "
	                 "specs; and the input vcf has to be decomposed as well if needed.\n" + "#" * 40)

def parseArgs(scriptname, argv):

	new_vcf_name = None
	column_tumor, column_normal = None, None

	warnings.simplefilter(action='ignore', category=FutureWarning)
	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)


	try:
		opts, args = getopt.getopt(argv,"hi:t:o:",[ "vcf=", \
		                                            "tumor_column=", "normal_column=", \
		                                            "debug", "outfilename=", \
		                                            "threshold_AR=", \
		                                            "help"] )
		log.info(opts)

	except getopt.GetoptError:
		usage(scriptname, opts)
		exit(2)
	for opt, arg in opts:
		if opt == '-h' or opt == "--help":
			usage(scriptname, opts)
			exit()
		elif opt in ("", "--threshold_AR"):
			try:
				global AR_threshold_for_GT
				AR_threshold_for_GT = float(arg)
			except TypeError:
				log.info("ERROR: threshold values MUST be integer or float")
				exit(2)
		elif opt in ("", "--tumor_column"):
			column_tumor = arg
		elif opt in ("", "--normal_column"):
			column_normal = arg
		elif opt in ("-o","--outfilename"):
			new_vcf_name = arg.strip()
		elif opt in ("-i", "--vcf"):
			fvcf = arg
			if not path.exists(fvcf):
				exit("ERROR: FNF --> " +fvcf)
		else:
			exit("Unknown Option --> " + opt )

	if column_tumor is None or column_normal is None:
		exit(
			"Please Provide column number for tumor and normal Samples; should be 10 and 11  - or -  11 and 10;\n"
			"options are: --tumor_column and --normal_column;\nAborting. ")
	if column_normal == column_tumor:
		exit("ERROR: number for the columns Tumor and Normal MUST be different")


	return(fvcf, column_tumor, column_normal, new_vcf_name)

def update_header(vcf):
	'''
	We modify the current header in the vcf object with the new fields or modifying old fields
	:param v: cyvcf2 VCF object
	:return v: cyvcf2 VCF object
	'''

	## if Adding Fields to INFO field
	vcf.add_info_to_header(
		{'ID': 'OGT', 'Description': ''.join([
			                                      'Original VaradictJava GT fields for each sample before '
			                                      'reassigning the GT value based on AR threshold (GT_0/1 < AR_',
			                                      str(AR_threshold_for_GT),
			                                      ' and GT_1/1 >= AR_', str(AR_threshold_for_GT), ' )']),
		 'Type': 'String', 'Number': '.'})

	## if Adding AR to FORMAT field
	vcf.add_format_to_header({'ID': 'AR', 'Description': 'Alt Allelic Ratio', 'Type': 'Float', 'Number': '1'})

	return vcf

def is_obj_nan(obj):
	if isnan(obj):
		return True
	return False

def get_GT_value_from_AR(AR_value, GT_value):
	'''
		return the GT value according to AR threshold values
		This is based off TGen's current thresholds of assigning the genotype
		1/1 if AR>=0.90 and 0/1 if AR<0.90

		Genotype representation for cyvcf2
		0 --> Unknown ; 1 --> Unknown_phased
		2 --> 0     ; 3 --> 0_PHASED_if_secondValue
		4 --> 1     ; 5 --> 1_PHASED_if_secondValue
		6 --> 2     ; 7 --> 2_PHASED_if_secondValue
		[0,0] == ./. ; [1,1] == .|.
		[1,0] == ./. ; [0,1] == .|.
		[2,2] == 0/0 ; [2,2] == 0/0
		[2,3] == 0|0 ; [3,2] == 0/0
		[2,4] == 0/1 ; [4,2] == 1/0
		[2,5] == 0|1 ; [5,2] == 1/0
		[2,6] == 0/2 ; [6,2] == 2/0
		[3,6] == 0/2 ; [6,3] == 2|0
		[4,6] == 1/2 ; [6,4] == 2/1
		[5,6] == 1/2 ; [6,5] == 2|1
		[9,8] == 3/3 ; [8,9] == 3|3

		[2,2] == 0/0 ; [4,4] == 1/1
		[3,3] == 0|0 ; [5,5] == 1|1
		[6,6] == 2/2 ; [7,7] == 2|2
		[8,8] == 3/3 ; [9,9] == 3|3

		return [int(2),int(4)] ; --> 0/1
		return [int(4),int(4)] ; --> 1/1

		'''
	log.debug("AR =" + str(AR_value) + " ---  AR_threshold = " + str(AR_threshold_for_GT))



	try:
		if AR_value < AR_threshold_for_GT:
			if "|" in GT_value:
				return [2, 5]
			return [2, 4]
		if AR_value >= AR_threshold_for_GT:
			if "|" in GT_value:
				return [5, 5]
			return [4, 4]
	except ValueError:
		print("ERROR: AR value not a number")
	except TypeError:
		print("ERROR: AR value not of type Float")

def get_GT_value_from_GT_value(GT_value):
	'''
		return the numpy array compatible GT value according to string GT value
		See also Genotype representation for cyvcf2 in Class Genotype
	'''

	dico_mapping_GT = {
		"./.": [0, 0],
		"0/0": [2, 2],
		"0/1": [2, 4], "1/0": [4, 2], "1/1": [4, 4],
		"0/2": [2, 6], "2/0": [6, 2], "2/2": [6, 6],
		"0/3": [2, 8], "3/0": [8, 2], "3/3": [8, 8],
		".|.": [1, 1],
		"0|0": [3, 3],
		"0|1": [2, 5], "1|0": [4, 3],
		"0|2": [2, 7], "2|0": [6, 3],
		"0|3": [2, 9], "3|0": [8, 3], "3|3": [8, 9],
		"1|1": [4, 5], "1|2": [4, 7], "1|3": [4, 9], "1|4": [4, 11],
		"2|2": [6, 7], "2|3": [6, 9], "2|4": [6, 11], "2|5": [6, 13],

	}  ## unused value ; kept only for the mapping informatino

	x = GenotypeInv(list(GT_value))
	try:
		return x.get_gt_numpy_compatible()
	except ValueError:
		print("ERROR: GT value ")
	except TypeError:
		print("ERROR: GT value not of right type ")
	except Exception as e:
		print("ERROR: Unknown Error ; Check with the Author :-( ; "+str(e))

def process_GTs(tot_number_samples, v, col_tumor, col_normal):
	'''
	Reassign GT value based on ala TGen threshold for AR value using _th_AR_for_GT_ CONSTANT

	:param tot_number_samples:
	:param v: variant record
	:param col_tumor:
	:param col_normal:
	:return: updated variant record
	'''

	if tot_number_samples != 2:
		msg = "Expected 2 Samples in VCF found {}. We are suppose to treat the vcf as a SOMATIC vcf and expect two samples;  Aborting.".format(tot_number_samples)
		raise Exception(msg)
	## capturing original GTs and adding them to INFO field
	v.INFO["OGT"] = ','.join([ str(Genotype(li)) for li in v.genotypes ])
	## ReAssiging GT with value based on AR thresholds comparison to CONSTANT threshold value
	GTs = [[0], [0]]  ; # need to init  list as we used index later for the list to replace values
	GTOs = [ str(Genotype(li)) for li in v.genotypes ]
	ARs = v.format('AR')
	idxN = 0 if col_normal == 10 else 1
	idxT = 1 if col_tumor == 11 else 0
	## we need to keep the order of the information based on the index; so the list GTs MUST be ordered;
	GTs[idxN] = get_GT_value_from_GT_value(GTOs[idxN]) ## we do not modify the GT field for the Normal sample
	GTs[idxT] = get_GT_value_from_AR(ARs[idxT][0], GTOs[idxT]) ## we do modify the GT field for the Tumor Sample based on defined threshold
	v.set_format('GT', np.array(GTs))
	log.debug("v after reassigning GT: " + str(v))
	return v


def process_GTs_Deprecated(tot_number_samples, v, col_tumor, col_normal):
	'''
	Reassign GT value based on ala TGen threshold for AR value using _th_AR_for_GT_ CONSTANT

	:param tot_number_samples:
	:param v:
	:param col_tumor:
	:param col_normal:
	:return: v
	'''

	## capturing original GTs and adding them to INFO field
	v.INFO["OGT"] = ','.join([ str(Genotype(li)) for li in v.genotypes ])
	## ReAssiging GT with value based on AR thresholds comparison to CONSTANT threshold value
	GTs = []
	ARs = v.format('AR')
	for sidx in range(tot_number_samples):
		GTs.append(get_GT_value_from_AR(ARs[sidx][0]))  ## cyvcf returns an array when using v.format('AR')
	v.set_format('GT', np.array(GTs))
	return v

def add_new_flags(v, column_tumor, column_normal, tot_number_samples):
	'''
	Calculate the AR for each sample in the variant record v
	The Total number of Sample in the VCF file is given by the variable tot_number_samples
	We assume that the numberof sample does not vary form one record to another as recommended in VCF specs

	:param tot_number_samples: number of samples in vcf (should be 2)
	:param v: cyvcf2 object which is going ot be updated
	:param col_tumor: column number of the Tumor sample in the vcf file (should be 11 as the Tgen's convention
	decided on column 11 for the tumor sample)
	:param col_normal: column number of the Constitutional (or Normal) sample in the vcf file (should be 10 as the
	Tgen's convention decided on column 10 for the tumor sample)
	:return: v
	'''

	idxT = 0 if int(column_tumor) == 10 else 1
	idxN = 1 if int(column_normal) == 11 else 0
	log.debug("___".join(str(x) for x in [ idxT, idxN ]) )
	if 'AF' in v.FORMAT:
		AR_tumor = v.format('AF')[idxT][0]
		AR_normal = v.format('AF')[idxN][0]
	if 'DP' in v.FORMAT:
		DP_tumor = v.format('DP')[idxT][0]
		DP_normal = v.format('DP')[idxN][0]

	if is_obj_nan(AR_tumor): AR_tumor = 0.00 ## takes care of possible abscence of the AF flag for unknown reason
	if is_obj_nan(AR_normal): AR_normal = 0.00
	if is_obj_nan(DP_tumor): DP_tumor = 0.00
	if is_obj_nan(DP_normal): DP_normal = 0.00

	if idxT == 0: ## as we do not loop over samples automatically, we need to decide the order of ARs value based on
		#  the sample column numbers given by user.
		ARs = [float(AR_tumor), float(AR_normal)]
	else:
		ARs = [ float(AR_normal), float(AR_tumor) ]

	log.debug("\t".join([ str(x) for x  in [ idxT, idxN , AR_tumor, AR_normal, DP_tumor, DP_normal ] ] ))
	v.set_format('AR', np.array(ARs))

	return process_GTs(tot_number_samples, v, column_tumor, column_normal)

#@#########
#@ MAIN  ##
#@#########
if __name__ == "__main__":

	vcf_path, column_tumor, column_normal, new_vcf_name = parseArgs(argv[0], argv[1:])  ; ## tth means tuple of thresholds
	vcf = VCF(vcf_path)

	print(' '.join([ "Constant AR threshold is: ", str(AR_threshold_for_GT) ]) )

	if new_vcf_name is None:
		new_vcf = '.'.join([str(vcf_path), "uts.vcf"])
	else:
		new_vcf = new_vcf_name


	vcf = update_header(vcf)

	# create a new vcf Writer using the input vcf as a template.
	w = Writer(new_vcf, vcf)

	tot_number_samples = len(vcf.samples)
	if tot_number_samples != 2:
		exit("ERROR: Number of Sample greater than 2; Expected 2 samples only TUMOR and NORMAL")

	log.info("looping over records ...")
	for v in vcf: ## v for variant which represents one "variant record"
		v = add_new_flags(v, column_tumor, column_normal, tot_number_samples)
		if v is not None:
			w.write_record(v)

	w.close()
	vcf.close()
	log.info("work completed")
	log.info('new vcf is << {} >>'.format(new_vcf))
	exit()
