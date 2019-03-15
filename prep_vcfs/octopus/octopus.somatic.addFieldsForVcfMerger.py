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

	## if Adding AR to FORMAT field
	vcf.add_format_to_header({'ID': 'AR', 'Description': 'Alt tier1 Allelic Ratios for each sample in same order as list of samples found in VCF beyond column FORMAT', 'Type': 'Float', 'Number': 'A'})
	vcf.add_format_to_header({'ID': 'AD',
	                          'Description': 'Allele Depth',
	                          'Type': 'String', 'Number': '.'})

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
			return [2,4]
		if AR_value >= AR_threshold_for_GT:
			return [4,4]
	except ValueError:
		print("ERROR: AR value not a number")
	except TypeError:
		print("ERROR: AR value not of type Float")

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
	return v

def add_new_flags(v, column_tumor, column_normal, filter, tot_number_samples):
	'''
	Calculate the AR for each sample in the variant record v
	The Total number of Sample in the VCF file is given by the variable tot_number_samples
	We assume that the numberof sample does not vary form one record to another as recommended in VCF specs
	'''

	idxT = 0 if int(column_tumor) == 10 else 1
	idxN = 1 if int(column_normal) == 11 else 0
	log.debug("___".join(str(x) for x in [ idxT, idxN ]) )

	AR_tumor = v.format('MAP_VAF')[idxT][0]		## returns numpy.float32
	AR_normal = v.format('MAP_VAF')[idxN][0]	## returns numpy.float32
	#DP_tumor = v.format('DP')[idxT][0] if ',' in v.format('DP')[idxT] else v.format('DP')[idxT]		## returns numpy.str_
	#DP_normal = v.format('DP')[idxN][0] if ',' in v.format('DP')[idxN][0] else v.format('DP')[idxN]		## returns numpy.str_
	DP_tumor = v.format('DP')[idxT]
	DP_normal = v.format('DP')[idxN]
	if is_obj_nan(float(AR_tumor)): AR_tumor = 0.00
	if is_obj_nan(float(AR_normal)): AR_normal = 0.00
	if is_obj_nan(int(DP_tumor)): DP_tumor = 0.00
	if is_obj_nan(int(DP_normal)): DP_normal = 0.00

	if idxT == 0:
		ARs = [float(AR_tumor), float(AR_normal)]
	else:
		ARs = [ float(AR_normal), float(AR_tumor) ]
	ADs = [ (0,0), (0,0) ] ## HARDCODED information ;
	# Because Octopus does not provide enough information to calculate AD, we assign default
	# values of 0,0 ## can be discussed and modify if users think differently

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
