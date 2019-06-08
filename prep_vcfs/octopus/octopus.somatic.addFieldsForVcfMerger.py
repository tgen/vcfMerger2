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

## WHAT DOES THIS SCRIPT DO?
## In order to <<harmonized>> the GENOTYPE flags, we add, modifed, remove or replace the Genotype Octopus' Flags
## We want some common following flags in the Genotype fields for each of our vcf that need to be merged:
## GT:DP:AR:AD are the flags we want to be common;
## We are not removing





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
AR_threshold_for_GT = 0.90  ## value HARDCODED

class Genotype(object):
    __slots__ = ('alleles', 'phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    ## genotypes = [Genotype(li) for li in variant.genotypes ]
    ## print genotypes
    ## which shows: [./., ./., ./., 1/1, 0/0, 0/0, 0/0]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)
    __repr__ = __str__

def usage(scriptname, opts):
	print("\nUSAGE: \npython3 " + scriptname + '  --help')
	print("python3 " + scriptname + " -i octopus.somatic.snvs.pass.vcf  --tumor_column 11 --normal_column 10  -o updated_vcf.vcf\n")


	print("")
	print("options available:")
	print(" -i|--fvcf  [ Mandatory, no default value, String Filename full or relative path expected ]\n", \
	      "--tumor_column  [ Mandatory, no default value, Integer Expected ]\n", \
	      "--normal_column [ Mandatory, no default value, Integer Expected ]\n", \
	      "-o|--outfilename  [ Optional, no default value, String Expected ]\n", \
	      "--threshold_AR [ Optional; default value:0.9 ; float expected ]\n", \
	      )
	print("")

	print("#"*40 + "\nWARNING WARNING\n" + "#"*40 )
	print("1) This script is to be used only with somatic snvs vcf having NORMAL sample in "
	                 "column 10 and TUMOR sample in column 11; if not like this, update your vcf file to these "
	                 "specs;\n2) and the vcf has to be decomposed as well.\n")


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
			except TypeValueError:
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
		usage(scriptname, opts)
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
		{'ID': 'OCGT', 'Description': ''.join([
			                                      'Original Mutect2 GT fields for each sample before reassigning the GT value based on AR threshold (GT_0/1 < AR_',
			                                      str(AR_threshold_for_GT),
			                                      ' and GT_1/1 >= AR_', str(AR_threshold_for_GT), ' )']),
		 'Type': 'String', 'Number': '.'})

	## Adding AR and new AD to FORMAT field
	vcf.add_format_to_header({'ID': 'AR', 'Description': 'Alt Allelic Ratios for each sample in same order as list of samples found in VCF beyond column FORMAT', 'Type': 'Float', 'Number': '1'})
	vcf.add_format_to_header({'ID': 'AD',
	                          'Description': 'Reformatted Allele Depth according to specs (AD=ADP-ADO,ADO)',
	                          'Type': 'Integer', 'Number': '2'})
	vcf.add_format_to_header({'ID': 'ADO',
	                          'Description': 'Original Octopus Allele Depth Value',
	                          'Type': 'Integer', 'Number': '1'})

	return vcf

def is_obj_nan(obj):
	if isnan(obj):
		return True
	return False

def get_GT_value_from_AR(AR_value):
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
			return [2,5]
		if AR_value >= AR_threshold_for_GT:
			return [5,5]
	except ValueError:
		print("ERROR: AR value not a number")
	except TypeError:
		print("ERROR: AR value not of type Float")
	else:
		return [1,1]

def process_GTs(tot_number_samples, v, col_tumor, col_normal):
	'''
	Reassign GT value based on ala TGen threshold for AR value using _th_AR_for_GT_ CONSTANT

	:param tot_number_samples:
	:param v:
	:param col_tumor:
	:param col_normal:
	:return:
	'''

	## capturing original GTs and adding them to INFO field
	v.INFO["OCGT"] = ','.join([ str(Genotype(li)) for li in v.genotypes ])
	## ReAssiging GT with value based on AR thresholds comparison to CONSTANT threshold value
	GTs = []
	ARs = v.format('AR')
	for sidx in range(tot_number_samples):
		GTs.append(get_GT_value_from_AR(ARs[sidx][0]))
	v.set_format('GT', np.array(GTs))
	log.debug("v after reassigning GT: "+str(v))
	return v

def check_if_PS_in_FORMAT_field(vcf_cyobj, input_vcf_path, new_vcf_name, list_of_fields_to_check):
	v1 = next(iter(vcf_cyobj))
	log.info("Checking PS flag presence in FORMAT ...")
	# Minimum_Expected_Fields_in_FORMAT_not_MANAGE_by_the_CODE
	ExpectedFlags="GT:DP"
	for FIELD in list_of_fields_to_check:
		if not FIELD in v1.FORMAT:
			log.warning(
				FIELD+" tag is ABSENT from the FORMAT field of OCTOPUS\nPlease Check you have run Octopus with appropriate options (i.e, with random forest option and no other filtering option; random forest filtering add ALL the adequate fields to FORMAT columns )")
			if ':'.join(v1.FORMAT) == ExpectedFlags:
				log.info("FORMAT field is equivalent to {}, but we request to have these flags at least: {}; And we try to manage the rest, such as AD and AR; Would be better if you could run random_Forest_Filtering; ".format(v1.FORMAT, ExpectedFlags))
				log.warning(
					"We assume the vcf has already been prepared for vcfMerger2 and therefore just copy the vcf by assigning the decomposed expected filename output")
				from shutil import copyfile
				copyfile(input_vcf_path, new_vcf_name)
				exit()
			else:
				log.error(FIELD+" flag NOT found in FORMAT field; Aborting VCF preparation.")
				exit(FIELD+" flag Absent")
		else:
			log.info(FIELD+" flag Found")

def if_dot_assign_negative_value(obj):
	if obj != ".":
		return obj
	else:
		return -2


def add_new_flags(v, column_tumor, column_normal, filter, tot_number_samples):
	'''
	Calculate the AR for each sample in the variant record v
	The Total number of Sample in the VCF file is given by the variable tot_number_samples
	We assume that the numberof sample does not vary form one record to another as recommended in VCF specs
	'''

	idxT = 0 if int(column_tumor) == 10 else 1
	idxN = 1 if int(column_normal) == 11 else 0
	log.debug("___".join(str(x) for x in [ idxT, idxN ]) )

	if 'AD' in v.FORMAT:
		log.debug(str(v))
		log.debug("AD____TUMOR is " + str(v.format('AD')[idxT]))
		log.debug("AD____TUMOR is "+str(v.format('AD')[idxT][0]))
		## capturing Original AD and ADP
		AD_tumor = v.format('AD')[idxT]
		ADP_tumor = v.format('ADP')[idxT]
		AD_normal = v.format('AD')[idxN]
		ADP_normal = v.format('ADP')[idxN]


		# AD_tumor = -1
		# ADP_tumor = 600
		#AD_normal = -1
		#ADP_normal = 9999

		log.debug(str(v))
		log.debug("ADP --->  " + str(ADP_normal) + " -----  " + str(ADP_tumor))
		## Re-Allocationg ADs to ADOs, new tag for Original Octopus AD flags and values
		AD_tumor = int(if_dot_assign_negative_value(AD_tumor))
		AD_normal = int(if_dot_assign_negative_value(AD_normal))
		ADP_tumor = int(if_dot_assign_negative_value(ADP_tumor))
		ADP_normal = int(if_dot_assign_negative_value(ADP_normal))
		ADOs = [AD_tumor,AD_normal] if idxT == 0 else [AD_normal,AD_tumor]
		log.debug("ADOs --->>>  " + str(ADOs) + " <<<<<<<<-----  ")
		v.set_format('ADO', np.array(ADOs))

		## Calculate AR (allele ration using AD and ADP)
		try:
			log.debug("before calcul AR, AD = ---->>    " + str(AD_normal) + " -----  " + str(AD_tumor))
			log.debug("before calcul AR, ADP = ---->>    " + str(ADP_normal) + " -----  " + str(ADP_tumor))

			if AD_tumor != "." and AD_tumor >=0:
				AR_tumor = round(float(AD_tumor/ADP_tumor), 2)
				## Reformmating AD to expected VCF specs for that Reserved AD field, using the original AD and ADP values
				AD_tumor = [ADP_tumor - AD_tumor, ADP_tumor]
			else:
				AR_tumor = int(-2)
				AD_tumor = [0, ADP_tumor]
		except ZeroDivisionError:
			log.debug("division by zero!")
			AD_tumor = [0, 0]
			AR_tumor = float(0.00)
		try:
			if AD_normal != "." and AD_normal >=0:
				AR_normal = round(float(AD_normal / ADP_normal), 2)
				AD_normal = [ADP_normal - AD_normal, ADP_normal]
			else:
				AR_normal = int(-2)
				log.debug("AD in the else of the try --> " + str(ADP_normal))
				AD_normal = [0, ADP_normal]

		except ZeroDivisionError:
			log.debug("division by zero!")
			AD_normal = [0, 0]
			AR_normal = float(0.00)

		log.debug("AR --->> " + str(AR_normal) + " -----  " + str(AR_tumor))
		log.debug("AD = ---->>    " + str(AD_normal) + " -----  " + str(AD_tumor))
		log.debug("DP tumor is  : "+str(v.format('DP')[idxT][0]))
		log.debug("DP normal is  : " + str(v.format('DP')[idxN][0]))
		DP_tumor = v.format('DP')[idxT][0] #if ',' in v.format('DP')[idxT] else v.format('DP')[idxT]		## returns numpy.str_
		DP_normal = v.format('DP')[idxN][0] #if ',' in v.format('DP')[idxN][0] else v.format('DP')[idxN]		## returns numpy.str_
		#DP_tumor = v.format('DP')[idxT]
		#DP_normal = v.format('DP')[idxN]
		log.debug(str(DP_normal) + " -----  " + str(DP_tumor))
		if is_obj_nan(float(AR_tumor)): AR_tumor = 0.00
		if is_obj_nan(float(AR_normal)): AR_normal = 0.00
		if is_obj_nan(int(DP_tumor)): DP_tumor = 0.00
		if is_obj_nan(int(DP_normal)): DP_normal = 0.00

		if idxT == 0:
			ARs = [AR_tumor, AR_normal]
			ADs = [AD_tumor, AD_normal]
		else:
			ARs = [AR_normal, AR_tumor]
			ADs = [AD_normal, AD_tumor]

	else:
		# Because Octopus does not provide enough information to calculate AD, we assign default
		# values of 0,0 ## can be discussed and modify if users think differently
		dummy_value = int(-2)  ## set dummy value for the AD when AD is absent or not capturable from octopus' vcf.
		AR_tumor = [dummy_value, dummy_value]
		AR_normal = [dummy_value, dummy_value]
		DP_tumor = v.format('DP')[idxT][0]
		DP_normal = v.format('DP')[idxN][0]
		AD_tumor = [dummy_value, dummy_value]
		AD_normal = [dummy_value, dummy_value]

		ADs = [AD_normal, AD_tumor]
		ARs = [AR_normal, AR_tumor]

	log.debug("ADs  is  : " + str(ADs))
	log.debug("ARs  is  : " + str(ARs))

	log.debug("\t".join([ str(x) for x  in [ idxT, idxN , AR_tumor, AR_normal, DP_tumor, DP_normal ] ] ))
	v.set_format('AR', np.array(ARs))
	v.set_format('AD', np.array(ADs))

	return process_GTs(tot_number_samples, v, column_tumor, column_normal)

#@#########
#@ MAIN  ##
#@#########
if __name__ == "__main__":

	vcf_path, column_tumor, column_normal, new_vcf_name = parseArgs(argv[0], argv[1:])  ; ## tth means tuple of thresholds
	vcf = VCF(vcf_path)

	log.info(' '.join([ "Constant AR threshold is: ", str(AR_threshold_for_GT) ]) )

	if new_vcf_name is None:
		new_vcf = '.'.join([str(vcf_path), "uts.vcf"])
	else:
		new_vcf = new_vcf_name

	## checking if PS flag is still present in the VCF genotype fields
	check_if_PS_in_FORMAT_field(vcf, vcf_path, new_vcf_name, ["PS"])

	vcf = update_header(vcf)

	# create a new vcf Writer using the input vcf as a template.
	w = Writer(new_vcf, vcf)

	tot_number_samples = len(vcf.samples)
	if tot_number_samples != 2:
		exit("ERROR: Number of Sample greater than 2; Expected 2 samples only TUMOR and NORMAL")

	log.info("looping over records ...")
	for v in vcf: ## v for variant which represents one "variant record"
		v = add_new_flags(v, column_tumor, column_normal, filter, tot_number_samples)
		if v is not None:
			w.write_record(v)

	w.close()
	vcf.close()
	log.info("work completed")
	log.info('new vcf is << {} >>'.format(new_vcf))
	exit()
