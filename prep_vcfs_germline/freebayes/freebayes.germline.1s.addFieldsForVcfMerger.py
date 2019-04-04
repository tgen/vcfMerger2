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
import argparse
from cyvcf2 import VCF, Writer, VCFReader
import numpy as np
import logging as log
#from myGenotype import Genotype


global AR_threshold_for_GT
AR_threshold_for_GT = 0.90  ;  ##  value HARDCODED het_0/1 < 90%  and   homalt_1/1 >= 90%

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

class UniqueStore(argparse.Action):
	"""
	class grabbed from stackOverflow (2018-11-08)
	https://stackoverflow.com/questions/23032514/argparse-disable-same-argument-occurences
	Thanks To the Community
	We override the function __call__ from argparse to check if one option is given more than once
	"""

	def __call__(self, parser, namespace, values, option_string):
		if getattr(namespace, self.dest, self.default) is not None:
			parser.error(option_string + " appears several times.  Please modify your options.")
		setattr(namespace, self.dest, values)

def update_header(vcf):
	'''
	We modify the current header in the vcf object with the new fields or modifying old fields
	:param v: cyvcf2 VCF object
	:return v: cyvcf2 VCF object
	'''

	## Adding New Fields to FORMAT field
	vcf.add_format_to_header({'ID': 'AR', 'Description': 'Alternate Allelic Ratios for each sample in same order as list of samples found in VCF beyond column FORMAT ; AR calculus is (AO/(AO+RO))',
	                          'Type': 'Float', 'Number': '1'})
	vcf.add_format_to_header({'ID': 'AD',
	                          'Description': 'Allele Depth Ref,Alt',
	                          'Type': 'Float', 'Number': '1'})
	vcf.add_format_to_header({'ID': 'GT',
	                          'Description': ''.join([ 'GT fields for each sample after reassigning the GT value based on AR threshold (GT_0/1 < AR_', str(AR_threshold_for_GT) ,
	                                                   ' and GT_1/1 >= AR_', str(AR_threshold_for_GT),  ' )' ]) ,
	                          'Type': 'String', 'Number': 'A'})
	vcf.add_info_to_header({'ID': 'OGT',
	                          'Description': ''.join([
		                                                 'Old GT (OGT) field for each sample after reassigning the GT value based on AR threshold (GT_0/1 < AR_',
		                                                 str(AR_threshold_for_GT),
		                                                 ' and GT_1/1 >= AR_', str(AR_threshold_for_GT), ' )']),
	                          'Type': 'String', 'Number': 'A'})

	# vcf.add_format_to_header({'ID': 'AD',
	#                           'Description': 'Allele Depth capture from ref_tier1 and alt_tier1 ',
	#                           'Type': 'Integer', 'Number': 'R'})

	return vcf

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
	log.debug("AR = " + str(AR_value) + " ---  AR_threshold = " + str(AR_threshold_for_GT))

	try:

		if AR_value < AR_threshold_for_GT:
			return [2,4]
		if AR_value >= AR_threshold_for_GT:
			return [4,4]
	except ValueError:
		print("ERROR: AR value not a number")
	except TypeError:
		print("ERROR: AR value not of type Float")

def process_GTs(tot_number_samples, v):
	'''
	Reassign GT value based on ala TGen threshold for AR value using _th_AR_for_GT_ CONSTANT

	:param tot_number_samples:
	:param v:
	:return: v
	'''

	## capturing original GTs and adding them to INFO field
	v.INFO["OGT"] = ','.join([ str(Genotype(li)) for li in v.genotypes ])
	## ReAssiging GT with value based on AR thresholds comparison to CONSTANT threshold value
	GTs = []
	ARs = v.format('AR')
	for sidx in range(tot_number_samples):
		GTs.append(get_GT_value_from_AR(ARs[sidx][0]))
	v.set_format('GT', np.array(GTs))
	return v

def update_flags(tot_number_samples, v):

	if v.INFO['DP'] == 0:
		return None
	list_flag_required_for_calculating_AR = ["RO", "AO"]  ## HARDCODED b/c specific to Strelka2 for INDELS calls but NOT SNVs.
	isRequired_Flag_Present = []
	for flag in list_flag_required_for_calculating_AR:
		isRequired_Flag_Present.append(bool(flag in v.FORMAT))

	if len(isRequired_Flag_Present) == len(list_flag_required_for_calculating_AR):
		# '''We need to deal here with INDEL processing and not just returning the record '''
		return process_snvs_records(tot_number_samples, v)
		# if v is None: return v
		# return process_GTs(tot_number_samples, v)
	else:
		return None

def process_indels_records(tot_number_samples, v):
	'''
	Calculate the AR for indel in each sample in the variant record v
	The Total number of Samples in the VCF file is given by the variable tot_number_samples
	We assume that the number of sample does not vary form one record to another as recommended in VCF specs.
	Within that function we then can calculate GTs because we have AR available
	'''
	## in Haplotype caller the Flags are the same for Indels and SNVs (so far)

	ARs, ADs, GTs, DPs = [], [], [], []

	## loop through samples to calculate AR for each one
	for sidx in range(tot_number_samples):
		if v.format('DP') is not None:
			DP = v.format('DP')[sidx][0]
			DPs.append(DP)
		else:
			DPs.append('0')
		if v.format('AD') is not None:
			AD = v.format('AD')[sidx]
			ADs.append(AD)
			try:
				AR = float(AD[0]/(AD[0]+AD[1])) if AR[0] != 0 else 0.00
				ARs.append(AR)
			except ZeroDivisionError:
				log.debug("You can't divide by zero!")
				AR = float(0.00)  ## so we make it zero manually
		else:
			ADs.append('0,0')
			ARs.append(float(0.00))

		log.error("AR={} ; AD={} ; locus= {}".format(str(AR),str(AD),str(str(v.CHROM)+":"+str(v.POS))))
		ARs.append(AR)
		ADs.append(AD)
		GTs.append(get_GT_value_from_AR(AR))
		log.error("GT={} ".format(get_GT_value_from_AR(AR)))
	v.set_format('GT', np.array(GTs))
	v.set_format('AR', np.array(ARs))
	v.set_format('AD', np.array(ADs))

	## returning the updated variant; if v is None, this means the variant got filtered out based on rules
	return v

def process_snvs_records(tot_number_samples, v):
	'''
	Calculate the AR for each sample in the variant record v
	The Total number of Sample in the VCF file is given by the variable tot_number_samples
	We assume that the number of sample does not vary form one record to another as recommended in VCF specs.
	Within that function we then can calculate GTs because we have AR available
	'''
	ARs, ADs, DPs, GTs, OGTs = [], [], [], [], []

	## get REF and ALT bases ; note: ALT is a list in cyvcf2, not a character as multiple ALT can exist
	## here we only deal with the first ALT ## TODO implement AR for each ALT
	## loop through samples to calculate AR for each one
	for sidx in range(tot_number_samples):
		refCounts = v.format('RO')[sidx][0]
		altCounts = v.format('AO')[sidx][0]
		try:
			AR = float(altCounts/(refCounts + altCounts)) if altCounts != 0 else 0.00
		except ZeroDivisionError:
			log.info("refcounts {}, altCounts {} locus {}:{}".format( refCounts, altCounts, str(v.CHROM), str(v.POS) ) )
			log.info("You can't divide by zero!")
			AR = float(0.00)   ## so we make it zero manually
		except Exception:
			AR = float(0.00)  ## so we force it zero

		AD = (int(altCounts), int(refCounts))
		ADs.append(AD)
		ARs.append(AR)
		GTs.append(get_GT_value_from_AR(AR))
	#	OGTs.append(v.format('GT')[sidx])

	# log.info(str(OGTs))
	log.debug(str(GTs))
	v.INFO["OGT"] = ','.join([str(Genotype(li)) for li in v.genotypes])
	v.set_format('GT', np.array(GTs))
	v.set_format('AD', np.array(ADs))
	v.set_format('AR', np.array(ARs))

	## returning the updated variant; if v is None, this means the variant got filtered out based on rules
	return v


def check_inputs(vcf, threshold_AR):
	"""
	:param lvcfs:
	:param ltoolnames:
	:param ltpo:
	:param acronyms:
	:param delim:
	:param lprepped_vcf_outfilenames:
	:return:
	"""
	if vcf is None :
		log.info("ERROR: Found list of input vcfs to be prepped or to be merged empty")
		sys.exit("ERROR: list of vcfs empty")
	if not os.path.exists(vcf):
		log.error("ERROR:  FILE NOT FOUND  --->  Check your input for vcf:" + vcf)
		sys.exit(-1)
	try:
		threshold_AR = float(threshold_AR)
	except ValueError:
		log.error("threshold_AR value not a Float")
	if threshold_AR < 0.00 or threshold_AR > 1.00:
		raise ValueError("threshold_AR Value MUST be interval [0.00, 1.00] ")


	log.info(' '.join(["Constant AR threshold is: ", str(AR_threshold_for_GT)]))
	log.info("input vcf = " + str(vcf))

def make_parser_args():
	parser = argparse.ArgumentParser(description='Processing arguments and options ...')
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	isRequired = True


	required.add_argument('-i','--vcf',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='Input VCF to be updated to specs for vcfMerger2')

	optional.add_argument('-o','--prep-outfilename',
	                      required=False,
	                      action=UniqueStore,
	                      help='outfilename for the prepped vcf')

	optional.add_argument('--threshold-AR',
	                      required=False,
	                      action=UniqueStore,
	                      help='AlleRatio threshold value to assign genotype; 0/1 if less than threshold, 1/1 if equal or above threshold; default is 0.90 ; range ]0,1] ')

	print(str(parser.prog) + "   " + str(parser.description))
	return parser

def main(args):

	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)

	vcf = None
	if args["vcf"]:
		vcf = str(args["vcf"])
		vcf_path = os.path.abspath(vcf)
		log.info("vcf given:\t\t{}".format(str(vcf)))
		log.info("vcf_path captured:\t{}".format(str(vcf_path)))

	prep_outfilename = None
	if args["prep_outfilename"]:
		prep_outfilename = str(args["prep_outfilename"])
		log.info("tool-specific filename for prepped upto vcfMerger2-specs vcf: " + str(prep_outfilename))

	TH_AR = 0.90
	if args['threshold_AR']:
		TH_AR = args['threshold_AR']
		if isinstance(TH_AR, (int, float, complex)):
			raise Exception("Threshold-AR must be a float or integer value between 0 and 1 (range ]0,1]). Check your inputs.")
		log.info("user given threshold for AR: " + str(TH_AR))

	check_inputs(vcf, TH_AR)

	vcf = VCFReader(vcf_path)
	filebasename = str(os.path.splitext(vcf_path)[0])
	if prep_outfilename is None:
		new_vcf = '.'.join([filebasename, "uts.vcf"])
	else:
		new_vcf = prep_outfilename

	# we first Add/Modify/Update Fields to the Header
	update_header(vcf)

	# create a new vcf Writer using the updated 'vcf' object above as a template for the header (mostly).
	w = Writer(new_vcf, vcf)
	tot_number_samples = len(vcf.samples)

	c=0
	log.info("looping over records ...")
	for v in vcf: ## v for variant which represents one "variant record"
		c = c+1
		if c % 100000 == 0:
			log.info("processed {}".format(str(c)))
		v = update_flags(tot_number_samples, v)
		if v is not None:
			w.write_record(v)

	w.close()
	vcf.close()
	log.info("work completed")
	log.info('new vcf is << {} >>'.format(new_vcf))
	sys.exit()

#@#########
#@ MAIN  ##
#@#########
if __name__ == "__main__":

	parser = make_parser_args()
	args = vars(parser.parse_args())  # vars() function returns a dictionary of key-value arguments
	print(str(args))
	main(args)



