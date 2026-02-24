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

import sys, os
import getopt
from sys import argv  # Used to bring in the feature argv, variables or arguments
from cyvcf2 import VCF, Writer
import numpy as np
import logging as log
import warnings
from math import isnan

global AR_threshold_for_GT
AR_threshold_for_GT: float = 0.90  ## value HARDCODED but dynamically modify it with --threshold_AR option


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
		return sep.join("0123456789."[a] for a in self.alleles)
	
	__repr__ = __str__


class GenotypeInv(object):
	
	def __init__(self, li):
		
		self.allele1 = li[0]
		self.phased = bool(0) if li[1] == "/" else bool(1)
		self.allele2 = li[2]
		self.GT = []
	
	def get_gt_numpy_compatible(self):
		self.GT = []  ## we need to reinit the GT list here otherwise shared by all instances. Weird because we reinitiated it already in the _init_
		if self.phased:
			if self.allele1 != ".":
				self.GT.append((2 * int(self.allele1)) + 3)
			else:
				self.GT.append(1)
			if self.allele2 != ".":
				self.GT.append((2 * int(self.allele2)) + 3)
			else:
				self.GT.append(1)
		else:
			if self.allele1 != ".":
				self.GT.append((2 * int(self.allele1)) + 2)
			else:
				self.GT.append(0)
			if self.allele2 != ".":
				self.GT.append((2 * int(self.allele2)) + 2)
			else:
				self.GT.append(0)
		return self.GT


def usage(scriptname):
	print("USAGE: \n" + scriptname + ' -i <VCF file [M]>  --tumor_column INT --normal_column INT ')
	print("")
	print("options available:")
	print(" -i|--fvcf  [ Mandatory, no default value, String Filename full or relative path expected ]\n", \
	      "--tumor_column  [ Mandatory, no default value, Integer Expected ]\n", \
	      "--normal_column [ Mandatory, no default value, Integer Expected ]\n", \
	      "-o|--outfilename  [ Optional, no default value, String Expected ]\n", \
	      "--threshold_AR [ Optional; default value:0.9 ; float expected ]\n", \
	      "--debug [Optional, Flag, for debug only; increase verbosity ]\n", \
	      )
	print("")
	
	print("#" * 40 + "\nWARNING WARNING\n" + "#" * 40)
	print(" This script is to be used only with somatic snvs vcf having NORMAL sample in column 10 and TUMOR sample "
	      "in column 11;")
	print("if not like this, update manually your vcf file to these specs; and the vcf has to be decomposed as "
	      "well ... ")
	print("or you may use the script 'prep_vcf_somatic.sh' to do it for you\n")


# print("List of INDEL Field already in Original Mutect2's VCF FORMAT columns: "
#       "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB ")
# print("List of SNV Field already in Original Mutect2's VCF FORMAT columns: "
#       "GT:AD:AF:F1R2:F2R1:MBQ:MFRL:MMQ:MPOS:SA_MAP_AF:SA_POST_PROB ")
# print("we will add DP to that column and << recalculate >> GT ala TGen; the original GTs are going to be "
#       "transfered to INFO field;\n We also add AR to Genotype fields")
# print("NOTE: the threshold_AR value is designed to assign the genotype to GT; if AR<=0.9, GT=0/1; above, "
#       "GT=1/1 ; you have to be aware that 0/0 does not exist in that instance.")


def parseArgs(scriptname, argv):
	new_vcf_name = None
	column_tumor, column_normal = None, None
	warnings.simplefilter(action='ignore', category=FutureWarning)
	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)
	
	try:
		opts, args = getopt.getopt(argv, "hi:o:", ["fvcf=", \
		                                           "debug", "outfilename=", \
		                                           "tumor_column=", "normal_column=", \
		                                           "threshold_AR=", "help"
		                                           ])
		log.info(opts)
	except getopt.GetoptError:
		usage(scriptname)
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h', "--help"):
			usage(scriptname)
			sys.exit()
		elif opt in ("", "--threshold_AR"):
			try:
				global AR_threshold_for_GT
				AR_threshold_for_GT = float(arg)
			except TypeError:
				log.info("ERROR: threshold values MUST be integer or float")
				sys.exit(2)
		elif opt in ("", "--tumor_column"):
			column_tumor = int(arg)
		elif opt in ("", "--normal_column"):
			column_normal = int(arg)
		elif opt in ("", "--debug"):
			print("DEBUG MODE")
			log.basicConfig(format=FORMAT_LOGGING, level=log.DEBUG)  ## not working due to incorrect implementation ##TODO
		elif opt in ("-o", "--outfilename"):
			new_vcf_name = arg.strip()
		elif opt in ("-i", "--vcfs"):
			fvcf = arg
			if not os.path.exists(fvcf):
				sys.exit("ERROR: FNF --> " + fvcf)
	
	log.debug(AR_threshold_for_GT)
	log.debug("normal_column = {} and tumor_column = {}".format(str(column_normal), str(column_tumor)))
	
	if column_tumor is None or column_normal is None:
		usage(scriptname, opts)
		sys.exit(
			"Please Provide column number for tumor and normal Samples; should be 10 and 11  - or -  11 and 10 respectively; Aborting. ")
	if column_normal == column_tumor:
		sys.exit("ERROR: number for the columns Tumor and Normal MUST be different")
	
	values_to_return = (fvcf, new_vcf_name, column_tumor, column_normal)
	
	return values_to_return


def update_header(vcf):
	'''
	We modify the current header in the vcf object with the new fields or modifying old fields
	:param v: cyvcf2 VCF object
	:return v: cyvcf2 VCF object
	'''
	## if Adding Fields to INFO field
	vcf.add_info_to_header(
		{'ID': 'OGT', 'Description': ''.join(['Original DeepSomatic GT fields for each sample before reassigning the GT value based on AR threshold (GT_0/1 < AR_', str(AR_threshold_for_GT),
		                                      ' and GT_1/1 >= AR_', str(AR_threshold_for_GT), ' )']), 'Type': 'String', 'Number': '.'})
	
	vcf.add_info_to_header(
		{'ID': 'OGQ', 'Description': ''.join(['Original DeepSomatic GQ fields for each sample. Conditional genotype quality']), 'Type': 'Integer', 'Number': '1'})
	
	vcf.add_info_to_header(
		{'ID': 'OPL', 'Description': ''.join(['Original DeepSomatic PL fields for each sample. Phred-scaled genotype likelihoods rounded to the closest integer']), 'Type': 'Integer', 'Number': 'G'})
	
	vcf.add_info_to_header(
		{'ID': 'OVAF', 'Description': ''.join(['Original DeepSomatic VAF fields for each sample. Variant allele fractions.']), 'Type': 'Float', 'Number': 'A'})
	
	## if Adding Fields to FORMAT field
	vcf.add_format_to_header({'ID': 'AR',
	                          'Description': 'Alt tier1 Allelic Ratios for each sample in same order as list of samples found in VCF beyond column FORMAT',
	                          'Type': 'Float', 'Number': '1'})
	
	return vcf


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
		print("ERROR: Unknown Error ; Check with the Author :-( ; " + str(e))


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
		msg = "Expected 2 Samples in VCF found {}. We are supposed to treat the vcf as a SOMATIC vcf and expect two samples; Sorry for the inconvenience.  Aborting.".format(tot_number_samples)
		raise Exception(msg)
	## capturing original GTs and adding them to INFO field
	v.INFO["OGT"] = ','.join([str(Genotype(li)) for li in v.genotypes])
	# log.info(f'v.format("GQ") == {v.format("GQ")}')
	# log.info(f'type of v.format("GQ") == {type(v.format("GQ"))}')
	#
	# log.info(f'v.format("PL") == {v.format("PL")}')
	# log.info(f'type of v.format("PL") == {type(v.format("PL"))}')
	#
	# log.info(f'v.format("VAF") == {v.format("VAF")}')
	# log.info(f'type of v.format("VAF") == {type(v.format("VAF"))}')
	
	# v.INFO["OGQ"] = ','.join(v.format('GQ'))  ## NOT WORKING --> RAISE ERROR numpy.ndarray
	# v.INFO["OPL"] = ','.join(v.format('PL'))
	# v.INFO["OVAF"] = ','.join(v.format('VAF'))
	return v


def process_records(tot_number_samples, var_rec, col_tumor, col_normal):
	"""
	:param tot_number_samples: number of samples in the VCF ; can be extracted from len(v.samples)
	:type tot_number_samples: int
	:param var_rec: variant_record object from cyvcf2.VCF
	:type var_rec: cyvcf2.VCF
	:param col_tumor: integer value representing the column number of the Tumor sample
	:type col_tumor: int
	:param col_normal: integer value representing the column number of the Tumor sample
	:type col_normal: int
	:return: variant record
	"""
	
	v = add_new_flags_AR_to_FORMAT(tot_number_samples, var_rec, col_tumor, col_normal)
	if v is None:
		return v
	return process_GTs(tot_number_samples, v, col_tumor, col_normal)


def is_obj_nan(obj):
	"""
        check if the object has the value NaN
    """
	log.debug("obj test for nan: {}".format(str(obj)))
	if isnan(obj):
		return True
	return False


def if_dot_assign_value_zero(obj, idx):
	"""
    check if the value is a dot and reassign a zero
    """
	try:
		if obj == "." or is_obj_nan(float(obj)):
			return 0
		elif str(obj) != "." and str(obj) != "-2147483648" and str(obj) != "./.":
			return obj
		else:
			# Negative value was causing issue and adding more bytes to a file; So we switch to Zero number
			return 0
	except TypeError:
		return 0


def add_new_flags_AR_to_FORMAT(tot_number_samples, v, column_tumor, column_normal, ):
	"""
    Calculate the AR for each sample in the variant record v
    The Total number of Sample in the VCF file is given by the variable tot_number_samples
    We assume that the number of sample does not vary form one record to another as recommended in VCF specs
    """
	
	idxT = 0 if int(column_tumor) == 10 else 1
	idxN = 1 if int(column_normal) == 11 else 0
	log.debug("___".join(str(x) for x in [idxT, idxN]))
	
	# in that version 0.7.4 of Octopus, the fields AD and AF and ADP are already present and AD, AF and ADP now follow the specs.
	# For each allele, there is a number now except for ADP with only one value which is supposed to represent Depth for an assigned allele but actually represent the depth for
	# the current sample in the FORMAT field whereas DP in INFO represent the total depth across all the samples (normally)
	# ##FORMAT=<ID=ADP,Number=1,Type=Integer,Description="Number of reads overlapping the position that could be assigned to an allele">
	# ##FORMAT=<ID=AF,Number=R,Type=Float,Description="Empirical allele frequency (AD / ADP)">
	# ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Empirical allele depth">
	# So we do not need to check if AD, AF or ADP are present, we assume they are there in our use case from the TGen pipeline

	if 'VAF' in v.FORMAT:
		AR_tumor = v.format('VAF')[idxT]
		AR_normal = v.format('VAF')[idxN]
		if idxT == 0:
			ARs = [AR_tumor, AR_normal]
		else:
			ARs = [AR_normal, AR_tumor]
	else:  # the else is entered only if there is no 'AD' flag in the octopus vcf
		log.error("VAF not Found in the VCF; Are you sure you are using a vcf generated by a version of DeepSomatic v1.9.0-gpu or up?; Aborting.")
		exit(2)
	
	# checking the values after processing and before adding them to the variant object v
	log.debug("ARs  are  : " + str(ARs))
	v.set_format('AR', np.array(ARs))
	log.debug("updated v object with new ARs: " + str(v))
	
	return process_GTs(tot_number_samples, v, column_tumor, column_normal)


# @#########
# @ MAIN  ##
# @#########
if __name__ == "__main__":
	
	vcf_path, new_vcf_name, column_tumor, column_normal = parseArgs(argv[0], argv[1:])
	log.info(f'recap input data after parsing arguments:')
	log.info(f'column_number normal and tumor: {column_normal} & {column_tumor}')
	log.info(' '.join(["Constant AR threshold is: ", str(AR_threshold_for_GT)]))
	
	if new_vcf_name is None:
		new_vcf = '.'.join([str(os.path.splitext(vcf_path)[0]), "uts.vcf"])
	else:
		new_vcf = new_vcf_name
	
	vcf = VCF(vcf_path)
	
	# we first Add/Modify/Update Fields to the Header
	update_header(vcf)
	
	# create a new vcf Writer using the input vcf as a template.
	w = Writer(new_vcf, vcf)
	
	tot_number_samples = len(vcf.samples)
	
	log.info("looping on records ...")
	for v in vcf:  ## v for variant which represents one "variant record"
		log.debug(f'v: {v}')
		rec = process_records(tot_number_samples, v, column_tumor, column_normal)
		if rec is not None:
			w.write_record(rec)
	
	w.close()
	vcf.close()
	log.info("work completed")
	log.info('new vcf is << {} >>'.format(new_vcf))
	sys.exit()
