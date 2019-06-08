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
from sys import argv	;# Used to bring in the feature argv, variables or arguments
from cyvcf2 import VCF, Writer, VCFReader
import numpy as np
import logging as log
import warnings

global AR_threshold_for_GT
AR_threshold_for_GT = 0.90  ;  ##  value HARDCODED het_0/1 < 90%  and   homalt_1/1 >= 90%

def usage(scriptname):
	print("USAGE: \n")
	print("Example minimum requirements:\n python3 " + scriptname + " -i somatic.snvs.pass.vcf "
	      "--tumor_column 11  "
	      "--normal_column 10 ")
	print("")
	print("options available:")
	print(" -i|--fvcf  [ Mandatory, no default value, String Filename full or relative path expected ]\n", \
	      "-o|--outfilename  [ Optional, no default value, String Expected ]\n", \
	      "--tumor_column  [ Mandatory, no default value, Integer Expected ]\n", \
	      "--normal_column [ Mandatory, no default value, Integer Expected ]\n", \
	      "--threshold_AR [ Optional; default value:0.9 ; float expected ]\n", \
	      "--debug [Optional, Flag, for debug only; increase verbosity ]\n", \
	      "--pass [Optional, Flag, will create a file wil only the Strelka's-specific PASS variants]\n", \
	      )

def parseArgs(scriptname, argv):

	## init some variables
	column_tumor, column_normal = None, None
	new_vcf_name = None
	generate_vcf_pass_calls_only = False

	warnings.simplefilter(action='ignore', category=FutureWarning)
	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)


	try:
		opts, args = getopt.getopt(argv,"hi:o:",[ "fvcf=", "help",\
		                                            "debug", \
		                                            "pass", \
		                                            "outfilename=", \
		                                            "tumor_column=", \
		                                            "normal_column=",\
		                                            "threshold_AR=" \
		                                            ] )
		log.info(opts)
	except getopt.GetoptError:
		usage(scriptname)
		sys.exit(2)

	for opt, arg in opts:
		if opt == '-h':
			usage(scriptname)
			sys.exit()
		if opt == '--help':
			usage(scriptname)
			sys.exit()
		elif opt in ("", "--pass"):
			generate_vcf_pass_calls_only = True
		elif opt in ("", "--threshold_AR"):
			try:
				global AR_threshold_for_GT
				AR_threshold_for_GT = float(arg)
				log.info("threshold AR reassigned by user: "+str(AR_threshold_for_GT))
			except TypeValueError:
				log.info("ERROR: threshold values MUST be integer or float")
				sys.exit(2)
		elif opt in ("", "--tumor_column"):
			column_tumor = int(arg)
		elif opt in ("", "--normal_column"):
			column_normal = int(arg)
		elif opt in ("","--debug"):
			log.info("in DEBUG options test")
			newloglevel = "DEBUG"
			newLevel = getattr(log, newloglevel.upper(), None)
			log.basicConfig(level=newLevel, format=FORMAT_LOGGING)
		elif opt in ("-o","--outfilename"):
			new_vcf_name = arg.strip()
		elif opt in ("-i", "--vcfs"):
			fvcf = arg
			if not os.path.exists(fvcf):
				sys.exit("ERROR: FNF --> " +fvcf)



	log.debug(AR_threshold_for_GT)
	log.debug("normal_column = {} and tumor_column = {}".format( str(column_normal), str(column_tumor) ))

	if column_tumor is None or column_normal is None:
		usage(scriptname, opts)
		sys.exit(
			"Please Provide column number for tumor and normal Samples; should be 10 and 11  - or -  11 and 10 respectively; Aborting. ")
	if column_normal == column_tumor:
		sys.exit("ERROR: number for the columns Tumor and Normal MUST be different")

	return(fvcf, new_vcf_name, column_tumor, column_normal, generate_vcf_pass_calls_only)

def update_header(vcf):
	'''
	We modify the current header in the vcf object with the new fields or modifying old fields
	:param v: cyvcf2 VCF object
	:return v: cyvcf2 VCF object
	'''

	## Adding New Fields to FORMAT field
	vcf.add_format_to_header({'ID': 'AR', 'Description': 'Alt tier1 Allelic Ratios for each sample in same order as list of samples found in VCF beyond column FORMAT ; AR calculus is alt_tier1/(alt_tier1 + ref_tier1)',
	                          'Type': 'Float', 'Number': '1'})
	vcf.add_format_to_header({'ID': 'GT',
	                          'Description': ''.join([ 'GT fields for each sample after reassigning the GT value based on AR threshold (GT_0/1 < AR_', str(AR_threshold_for_GT) ,
		                                       ' and GT_1/1 >= AR_', str(AR_threshold_for_GT),  ' )' ]) ,
	                          'Type': 'String', 'Number': 'A'})
	vcf.add_format_to_header({'ID': 'AD',
	                          'Description': 'Allele Depth capture from ref_tier1 and alt_tier1 ',
	                          'Type': 'Integer', 'Number': 'R'})

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
	v.INFO["OMGT"] = ','.join([ str(Genotype(li)) for li in v.genotypes ])
	## ReAssiging GT with value based on AR thresholds comparison to CONSTANT threshold value
	GTs = []
	ARs = v.format('AR')
	for sidx in range(tot_number_samples):
		GTs.append(get_GT_value_from_AR(ARs[sidx][0]))
	v.set_format('GT', np.array(GTs))
	return v

def update_flags(tot_number_samples, v, column_tumor, column_normal):

	## We check first if we deal with SNV or INDELS as the FLAGs in FORMAT are different:
	list_flag_indels = [ "TAR", "TIR", "TOR", "DP50", "FDP50", "SUBDP50", "BCN50" ]  ## HARDCODED b/c specific to Strelka2 for INDELS calls but NOT SNVs.
	isAnyOfIndelFlag = []
	for F in list_flag_indels:
		isAnyOfIndelFlag.append(bool(F in v.FORMAT))

	if any(isAnyOfIndelFlag):
		#'''We need to deal here with INDEL processing and not just returning the record '''
		return process_indels_records(tot_number_samples, v, column_tumor, column_normal)
	else:
		return process_snvs_records(tot_number_samples, v, column_tumor, column_normal)

def process_indels_records(tot_number_samples, v, column_tumor, column_normal):
	'''
	Calculate the AR for indel in each sample in the variant record v
	The Total number of Samples in the VCF file is given by the variable tot_number_samples
	We assume that the number of sample does not vary form one record to another as recommended in VCF specs.
	Within that function we then can calculate GTs because we have AR available
	'''
	'''
	excerpt from: https://www.biostars.org/p/51965/
	Per the Strelka paper, tier1 counts are simply more stringent than tier2. I'd recommend using tier1 counts...
	So, as in a normal VCF, DP would be your total depth across the site. In a normal VCF, the AD tag is a comma-delimited
	list of REF and ALT counts. There is a caveat mentioned at this link, but a workaround is to use the first entry
	under TIR as your ALT count, and the remainder in DP as your REF count.
    '''
	'''
	The closest number to what you are looking for is TAR. Strictly this is reads strongly support the reference 
	allele or an alternate overlapping indel. So TIR/(TIR+TAR) will give you approximately fraction of reads supporting
	the variant allele among all reads strongly supporting one local haplotype.
	TAR might not look like the right value in some contexts compared to the way reads appear in IGV because of 
	local repeat structures.Many reads provide very similar support to two or more alleles even though 
	they are mapped/appear in IGV in such a way so as to apparently support the reference.
	'''
	'''	
	TGen NOTE: we found out that sometimes DP value is less that TIR value (if 100% alternate indel) or 
	DP value is less then TAR value when no indels present like in normal sample, or DP value is less
	than TAR+TIR even though we could expect DP=TAR+TIR+TOR ; this might be due to what was explained in the 
	BioStar thread that due to error in strelka DP reads count may be less thatn TAR+TIR+TOR ; so calulating the REF 
	becomes tricky.
	'''

	ARs, ADs, GTs = [], [], []

	try:
		## loop through samples to calculate AR for each one
		for sidx in range(tot_number_samples):
			DP = v.format('DP')[sidx][0]
			TAR= v.format('TAR')[sidx][0] ## we only manage tier1 values
			TIR = v.format('TIR')[sidx][0]  ## we only manage tier1 values
			try:
				log.debug("TIR={}, DP={}, locus {}:{}".format(str(TIR), str(DP),  str(v.CHROM), str(v.POS)))
				AR = float(TIR/(TAR+TIR)) if TIR>0 else 0
				#AD = (int(max(DP-TIR, 0)), int(TIR))
				AD = (int(TAR), int(TIR))
			except ZeroDivisionError:
				log.debug("TAR= {}, TIR={}, DP={}, AR={} ; locus {}:{}".format(str(TAR), str(TIR),str(DP),str(AR),str(v.CHROM),str(v.POS)))
				log.debug("You can't divide by zero!")
				AR = float(0.0)  ## so we make it zero manually

			log.debug("AR={} ; AD={} ; locus= {}".format(str(AR),str(AD),str(str(v.CHROM)+":"+str(v.POS))))
			ARs.append(AR)
			ADs.append(AD)
			GTs.append(get_GT_value_from_AR(AR))
			log.debug("GT={} ".format(get_GT_value_from_AR(AR)))
		v.set_format('GT', np.array(GTs))
		v.set_format('AR', np.array(ARs))
		v.set_format('AD', np.array(ADs))

	except Exception as e:
		log.error("record raising ERROR: {}".format(str(v)))
		raise(e)
	## returning the updated variant; if v is None, this means the variant got filtered out based on rules
	return v

def process_snvs_records(tot_number_samples, v, column_tumor, column_normal):
	'''
	Calculate the AR for each sample in the variant record v
	The Total number of Sample in the VCF file is given by the variable tot_number_samples
	We assume that the number of sample does not vary form one record to another as recommended in VCF specs.
	Within that function we then can calculate GTs because we have AR available
	'''
	ARs, ADs, DPs, GTs = [], [], [], []

	## need to know the column number to be sure that we filter on the correct sample
	# col_tumor = 11 ; col_normal = 10 or vice-versa ;
	# idxT = 0 if column_tumor == 10 else 1
	# idxN = 0 if column_normal == 10 else 1

	try:
		## get REF and ALT bases ; note: ALT is a list nt a character as multiple ALT can exist
		## here we only deal with the first ALT ## TODO implement AR for each ALT
		ref_tag = ''.join([v.REF, "U"])
		alt_tag = ''.join([v.ALT[0], "U"])
		## loop through samples to calculate AR for each one
		for sidx in range(tot_number_samples):
			refCounts = v.format(ref_tag)[sidx]
			altCounts = v.format(alt_tag)[sidx]
			## we only consider tier1, therefore only the index 0 is needed
			ref_tier1 = int(refCounts[0])
			alt_tier1 = int(altCounts[0])
			try:
				AR = float(alt_tier1/(alt_tier1 + ref_tier1))
			except ZeroDivisionError:
				log.debug("refcount {}, altCounts {}, refTier1 {}, altTier1 {}, locus {}:{}".format( refCounts, altCounts, str(ref_tier1), str(alt_tier1), str(v.CHROM), str(v.POS) ) )
				log.debug("You can't divide by zero!")
				AR = float(0.0)  ## so we make it zero manually

			ARs.append(AR)
			AD = (ref_tier1, alt_tier1)
			ADs.append(AD)
			GTs.append((get_GT_value_from_AR(AR)))


		v.set_format('GT', np.array(GTs))
		v.set_format('AR', np.array(ARs))
		v.set_format('AD', np.array(ADs))

	except Exception as e:
		log.error("record raising ERROR: {}".format(str(v)))
		raise(e)
	## returning the updated variant; if v is None, this means the variant got filtered out based on rules
	return v


#@#########
#@ MAIN  ##
#@#########
if __name__ == "__main__":

	vcf_path, new_vcf_name, column_tumor, column_normal, generate_vcf_pass_calls = parseArgs(argv[0], argv[1:])  ; ## tth means tuple of thresholds

	log.info(' '.join(["Constant AR threshold is: ", str(AR_threshold_for_GT)]))
	log.info("vcf_path = " + str(vcf_path))

	vcf = VCFReader(vcf_path)
	filebasename = str(os.path.splitext(vcf_path)[0])
	if new_vcf_name is None:
		new_vcf = '.'.join([filebasename, "uts.vcf"])
	else:
		new_vcf = new_vcf_name

	if generate_vcf_pass_calls:
		## whether we filter or not, we may output a vcf with the Strelka's PASS variant only
		new_vcf_pass_only = ".".join([ new_vcf[:-4], "pass.vcf" ])
		print("filename for vcf with PASS calls only: " + str(new_vcf_pass_only))

	# we first Add/Modify/Update Fields to the Header
	update_header(vcf)

	# create a new vcf Writer using the updated 'vcf' object above as a template for the header (mostly).
	w = Writer(new_vcf, vcf)
	if generate_vcf_pass_calls: wpo = Writer(new_vcf_pass_only, vcf)

	tot_number_samples = len(vcf.samples)

	log.info("looping over records ...")
	for v in vcf: ## v for variant which represents one "variant record"
		v = update_flags(tot_number_samples, v, column_tumor, column_normal)
		if v is not None:
			w.write_record(v)
			if generate_vcf_pass_calls and v.FILTER is None: ## with cyvcf2, if FILTER has value PASS, it returns None
				wpo.write_record(v)

	w.close()
	if generate_vcf_pass_calls: wpo.close()
	vcf.close()
	log.info("work completed")
	log.info('new vcf is << {} >>'.format(new_vcf))
	sys.exit()

