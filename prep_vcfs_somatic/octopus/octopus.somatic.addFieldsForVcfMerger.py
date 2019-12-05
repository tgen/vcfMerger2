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

		try:
			## we added the if statements because Octopus may put only one letter in GT when dealing with chrY
			## so we decided to rebuilt the Genotype to be consistent here
			if len(li) != 3:
				if li[0] == "." or li[0] == "0":
					li = ['0', '|', '0']
				else:
					li = [0, "|", li[0]]

			self.allele1 = li[0]
			self.phased = bool(0) if li[1] == "/" else bool(1)
			self.allele2 = li[2]
			self.GT = []
		except Exception:
			print("LI==" + li)
			exit(1)

	def get_gt_numpy_compatible(self):
		self.GT = []  ## we need to reinit the GT list here otherwise shared by all instances. Weird because we reinitiated it already in the _init_ ; I am probably missing knowledge in some python features behaviour for classes.
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

	print("#" * 40 + "\nWARNING WARNING\n" + "#" * 40)
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
		opts, args = getopt.getopt(argv, "hi:t:o:", ["vcf=", \
		                                             "tumor_column=", "normal_column=", \
		                                             "debug", "outfilename=", \
		                                             "threshold_AR=", \
		                                             "help"])
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
				AR_threshold_for_GT = round(float(arg), 4)
			except TypeValueError:
				log.info("ERROR: threshold values MUST be integer or float")
				exit(2)
		elif opt in ("", "--tumor_column"):
			column_tumor = int(arg)
		elif opt in ("", "--normal_column"):
			column_normal = int(arg)
		elif opt in ("-o", "--outfilename"):
			new_vcf_name = arg.strip()
		elif opt in ("-i", "--vcf"):
			fvcf = arg
			if not path.exists(fvcf):
				exit("ERROR: FNF --> " + fvcf)
		else:
			exit("Unknown Option --> " + opt)

	if column_tumor is None or column_normal is None:
		usage(scriptname, opts)
		exit(
			"Please Provide column number for tumor and normal Samples; should be 10 and 11  - or -  11 and 10;\n"
			"options are: --tumor_column and --normal_column;\nAborting. ")
	if column_normal == column_tumor:
		exit("ERROR: number for the columns Tumor and Normal MUST be different")

	return (fvcf, column_tumor, column_normal, new_vcf_name)


def update_header(vcf):
	'''
	We modify the current header in the vcf object with the new fields or modifying old fields
	:param v: cyvcf2 VCF object
	:return v: cyvcf2 VCF object
	'''

	## if Adding Fields to INFO field
	vcf.add_info_to_header(
		{'ID': 'OGT', 'Description': ''.join([
			'Original Octopus GT fields for each sample before reassigning the GT value based on AR threshold (GT_0/1 < AR_',
			str(AR_threshold_for_GT),
			' and GT_1/1 >= AR_', str(AR_threshold_for_GT), ' )']),
		 'Type': 'String', 'Number': '.'})

	## Adding AR and new AD to FORMAT field
	## NOTE: if the field already exist in the Header, it will not be replaced or update; You must rename the Field apready present in the VCF to add specifically the following fields to the vcf HEADER
	vcf.add_format_to_header({'ID': 'AR', 'Description': 'Alt Allelic Ratios for each sample in same order as list of samples found in VCF beyond column FORMAT', 'Type': 'Float', 'Number': '1'})
	vcf.add_format_to_header({'ID': 'AD',
	                          'Description': 'Reformatted Allele Depth according to specs (AD=ADP-ADO,ADO)',
	                          'Type': 'Integer', 'Number': '.'})
	vcf.add_format_to_header({'ID': 'ADO',
	                          'Description': 'Original Octopus Allele Depth Value',
	                          'Type': 'Integer', 'Number': '1'})

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
		log.error("ERROR: AR value not a number")
	except TypeError:
		log.error("ERROR: AR value not of type Float")


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
		log.error("ERROR: GT value ")
	except TypeError:
		log.error("ERROR: GT value not of right type ")
	except Exception as e:
		log.error("ERROR: Unknown Error ; Check with the Author :-( ; " + str(e))


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
	v.INFO["OGT"] = ','.join([str(Genotype(li)) for li in v.genotypes])
	## ReAssiging GT with value based on AR thresholds comparison to CONSTANT threshold value
	GTs = [[0], [0]];  # need to init  list as we used index later for the list to replace values
	GTOs = [str(Genotype(li)) for li in v.genotypes]
	ARs = v.format('AR')
	idxN = 0 if col_normal == 10 else 1
	idxT = 1 if col_tumor == 11 else 0
	## we need to keep the order of the information based on the index; so the list GTs MUST be ordered;
	GTs[idxN] = get_GT_value_from_GT_value(GTOs[idxN])  ## we do not modify the GT field for the Normal sample
	GTs[idxT] = get_GT_value_from_AR(ARs[idxT], GTOs[idxT])  ## we do modify the GT field for the Tumor Sample based on defined threshold
	v.set_format('GT', np.array(GTs))
	log.debug("v after reassigning GT: " + str(v))
	return v


def check_if_PS_in_FORMAT_field(vcf_cyobj, input_vcf_path, new_vcf_name, list_of_fields_to_check):

	try:
		# v1 = next(iter(vcf_cyobj))
		v1 = vcf_cyobj

		# Minimum_Expected_Fields_in_FORMAT_not_MANAGE_by_the_CODE
		ExpectedFlags = "GT:DP"
		for FIELD in list_of_fields_to_check:
			if not FIELD in v1.FORMAT:
				log.warning(
					FIELD + " tag is ABSENT from the FORMAT field of OCTOPUS\nPlease Check you have run Octopus with appropriate options (i.e, with random forest option and no other filtering option; random forest filtering add ALL the adequate fields to FORMAT columns )")
				if ':'.join(v1.FORMAT) == ExpectedFlags:
					log.info("FORMAT field is equivalent to {}, but we request to have these flags at least: {}; And we try to manage the rest, such as AD and AR; Would be better if you could run random_Forest_Filtering; ".format(v1.FORMAT, ExpectedFlags))
					log.warning(
						"We assume the vcf has already been prepared for vcfMerger2 and therefore just copy the vcf by assigning the decomposed expected filename output")
					from shutil import copyfile
					copyfile(input_vcf_path, new_vcf_name)
					exit()
				else:

					log.error(FIELD + " flag NOT found in FORMAT field; Aborting VCF preparation.")
					exit(FIELD + " flag Absent")
			else:
				log.info(FIELD + " flag Found")
	except StopIteration as si:
		log.info(si)
		log.warning("no records")
	log.info("Checking PS flag presence in FORMAT ...")


def if_dot_assign_negative_value(obj):
	if obj != ".":
		return obj
	else:
		return 0



def add_new_flags(v, column_tumor, column_normal, filter, tot_number_samples):
	'''
	Calculate the AR for each sample in the variant record v
	The Total number of Sample in the VCF file is given by the variable tot_number_samples
	We assume that the numberof sample does not vary form one record to another as recommended in VCF specs
	'''

	idxT = 0 if int(column_tumor) == 10 else 1
	idxN = 1 if int(column_normal) == 11 else 0
	log.debug("___".join(str(x) for x in [idxT, idxN]))

	if 'AD' in v.FORMAT:
		log.debug(str(v))
		log.debug("AD____TUMOR is " + str(v.format('AD')[idxT]))
		log.debug("AD____TUMOR is " + str(v.format('AD')[idxT][0]))
		log.debug("AD____NORMAL is " + str(v.format('AD')[idxN]))
		log.debug("AD____NORMAL is " + str(v.format('AD')[idxN][0]))
		## capturing Original AD and ADP
		AD_tumor = v.format('AD')[idxT]
		ADP_tumor = v.format('ADP')[idxT]
		AD_normal = v.format('AD')[idxN]
		ADP_normal = v.format('ADP')[idxN]

		log.debug("ADP____TUMOR is " + str(v.format('ADP')[idxT]))
		log.debug("ADP____NORMAL is " + str(v.format('ADP')[idxN]))

		DP_tumor = v.format('DP')[idxT][0]  # if ',' in v.format('DP')[idxT] else v.format('DP')[idxT]		## returns numpy.str_
		DP_normal = v.format('DP')[idxN][0]  # if ',' in v.format('DP')[idxN][0] else v.format('DP')[idxN]		## returns numpy.str_
		# DP_tumor = v.format('DP')[idxT]
		# DP_normal = v.format('DP')[idxN]
		log.debug(str(DP_normal) + " -----  " + str(DP_tumor))


		log.debug(str(v))
		log.debug("ADP --->  " + str(ADP_normal) + " -----  " + str(ADP_tumor))

		## Re-Allocationg ADs to ADOs, new tag for Original Octopus AD flags and values
		AD_tumor = int(if_dot_assign_negative_value(AD_tumor))
		AD_normal = int(if_dot_assign_negative_value(AD_normal))
		ADP_tumor = int(if_dot_assign_negative_value(ADP_tumor))
		ADP_normal = int(if_dot_assign_negative_value(ADP_normal))
		ADOs = [AD_normal, AD_tumor] if idxT == 1 else [AD_tumor, AD_normal]
		log.debug("ADOs --->>>  " + str(ADOs) + " <<<<<<<<-----  ")  ## ADO stands for AD Old values
		v.set_format('ADO', np.array(ADOs))

		## Calculate AR (allele ratio using AD and ADP) for TUMOR SAMPLE
		try:
			log.debug("before calcul AR, AD = ---->>  AD_N, AD_T  " + str(AD_normal) + " -----  " + str(AD_tumor))
			log.debug("before calcul AR, ADP = ---->> ADP_N, ADP_T   " + str(ADP_normal) + " -----  " + str(ADP_tumor))

			if AD_tumor != "." and AD_tumor >= 0:
				AR_tumor = round(float(AD_tumor / ADP_tumor), 4)
				## Reformmating AD to expected VCF specs for that Reserved AD field, using the original AD and ADP values
				AD_tumor = [ADP_tumor - AD_tumor, AD_tumor]
			else:
				AR_tumor = round(float(0.0), 4)
				AD_tumor = [ADP_tumor - AD_tumor, AD_tumor]
		except ZeroDivisionError:
			log.debug("division by zero!")
			AD_tumor = [0, 0]
			AR_tumor = round(float(0.0), 4)
		## Calculate AR (allele ratio using AD and ADP) for NORMAL SAMPLE
		try:
			if AD_normal != "." and AD_normal >= 0:
				AR_normal = round(float(AD_normal / ADP_normal), 4)
				AD_normal = [ADP_normal - AD_normal, 0]
			else:
				AR_normal = int(0)
				log.debug("AD in the else of the try --> where AD is either a dot or equals to zero ; AD ==" + str(AD_normal))
				log.debug("ADP in the else of the try --> where AD is either a dot or equals to zero ; ADP == " + str(ADP_normal))
				AD_normal = [ADP_normal, 0]

		except ZeroDivisionError:
			log.debug("division by zero!")
			AD_normal = [0, 0]
			AR_normal = round(float(0.0), 4)

		log.debug("AR --->> " + str(AR_normal) + " -----  " + str(AR_tumor))
		log.debug("AD = ---->>    " + str(AD_normal) + " -----  " + str(AD_tumor))
		log.debug("DP tumor is  : " + str(v.format('DP')[idxT][0]))
		log.debug("DP normal is  : " + str(v.format('DP')[idxN][0]))

		if is_obj_nan(round(float(AR_tumor), 4)): AR_tumor = round(float(0.0), 4)
		if is_obj_nan(round(float(AR_normal), 4)): AR_normal = round(float(0.0), 4)
		if is_obj_nan(int(DP_tumor)): DP_tumor = round(float(0.0), 4)
		if is_obj_nan(int(DP_normal)): DP_normal = round(float(0.0), 4)

		if idxT == 0:
			ARs = [AR_tumor, AR_normal]
			ADs = [AD_tumor, AD_normal]
		else:
			log.debug("idxT == " + str(idxT))
			ARs = [AR_normal, AR_tumor]
			ADs = [AD_normal, AD_tumor]

		log.debug("ARs --->> (N,T) " + str(ARs))
		log.debug("ADs --->> (N,T)" + str(ADs))


	else:  ## the else is entered only if there is no 'AD' flag in the octopus vcf
		# Because Octopus does not provide enough information to calculate AD, we assign default if AD or ADP or any other information needed to get AD is not present in vcf
		# values of 0,0 ## can be discussed and modify if users think differently

		dummy_value = int(0)  ## set dummy value for the AD when AD is absent or not capturable from octopus' vcf.
		# AR_tumor = [dummy_value, dummy_value] ; old line, Why did I assign two values to AR where only one was enough???
		# AR_normal = [dummy_value, dummy_value]

		log.debug("In Else because NO 'AD' flag found in vcf; Flag NOT FOUND")

		if 'MAP_VAF' in v.FORMAT:
			AR_tumor = v.format('MAP_VAF')[idxT]
			AR_normal = v.format('MAP_VAF')[idxN] if not isnan(v.format('MAP_VAF')[idxN]) else [round(float(0.0), 4)]
		else:
			AR_tumor = [dummy_value]
			AR_normal = [dummy_value]

		DP_tumor = v.format('DP')[idxT][0]
		DP_normal = v.format('DP')[idxN][0]
		AD_tumor = [dummy_value, dummy_value]
		AD_normal = [dummy_value, dummy_value]

		ADs = [AD_normal, AD_tumor]
		ARs = [AR_normal, AR_tumor]

	## checking the values after processing and before adding them to the variant object v
	log.debug("ADs  are   : " + str(ADs))
	log.debug("ARs  are  : " + str(ARs))
	log.debug("\t".join([ str(x) for x in [ "idxT", "idxN", "AR_tumor", "AR_normal", "DP_tumor", "DP_normal" ] ]))
	log.debug("\t".join([str(x) for x in [idxT, idxN, AR_tumor, AR_normal, DP_tumor, DP_normal]]))
	v.set_format('AR', np.array(ARs))
	v.set_format('AD', np.array(ADs))

	log.debug("updated v object with new ARs and ADs:  " + str(v))

	return process_GTs(tot_number_samples, v, column_tumor, column_normal)


# @#########
# @ MAIN  ##
# @#########
if __name__ == "__main__":

	vcf_path, column_tumor, column_normal, new_vcf_name = parseArgs(argv[0], argv[1:]);  ## tth means tuple of thresholds
	vcf = VCF(vcf_path)

	log.info(' '.join(["Constant AR threshold is: ", str(AR_threshold_for_GT)]))

	if new_vcf_name is None:
		new_vcf = '.'.join([str(vcf_path), "uts.vcf"])
	else:
		new_vcf = new_vcf_name

	## test if no variants in vcf:
	try:
		vtest = next(iter(vcf))  ## we just try to check if we have no variant at all in VCF
	except StopIteration as si:
		## we bring the header upto specs and write down the output file
		log.warning("No Variants found in VCF; Creating Final Empty VCF now ..." + str(si))
		vcf = update_header(vcf)
		w = Writer(new_vcf, vcf)
		w.close()
		exit()

	## checking if PS flag is still present in the VCF genotype fields
	# vcf = VCF(vcf_path)  ## as we have already consumed once the generator; if in case there is only one variant (edge case encountered already), we need to re-init vcf here; need to think about another way of doing these tests without consuming the vcf object
	# vcf = update_header(vcf)
	check_if_PS_in_FORMAT_field(vtest, vcf_path, new_vcf_name, ["PS"])

	vcf = VCF(vcf_path)  ## as we have already consumed twice the generator; we do not want to lose any variant
	vcf = update_header(vcf)

	# create a new vcf Writer using the input vcf as a template.
	w = Writer(new_vcf, vcf)

	tot_number_samples = len(vcf.samples)
	if tot_number_samples != 2:
		exit("ERROR: Number of Sample greater than 2; Expected 2 samples only TUMOR and NORMAL")

	log.info("looping over records ...")
	for v in vcf:  ## v for variant which represents one "variant record"
		v = add_new_flags(v, column_tumor, column_normal, filter, tot_number_samples)
		log.debug("v object right before writing it to the VCF ...")
		log.debug(str(v))
		log.debug("---")
		if v is not None:
			w.write_record(v)

	w.close()
	vcf.close()
	log.info("work completed")
	log.info('new vcf is << {} >>'.format(new_vcf))
	exit()
