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
from sys import argv	;# Used to bring in the feature argv, variables or arguments
from os import path
import getopt
from cyvcf2 import VCF
import numpy as np
import logging as log
import warnings
from collections import defaultdict


class Genotype(object):
    __slots__ = ('alleles', 'phased')

    def __init__(self, li):
        self.alleles = li[:-1]
        self.phased = li[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123456789."[a] for a in self.alleles)
    __repr__ = __str__
        # genotypes = [Genotype(li) for li in variant.genotypes ]
        # print genotypes
        # which shows: [./., ./., ./., 1/1, 0/0, 0/0, 0/0]

def usage(scriptname, opts):
	print("\nUSAGE: \npython3 " + scriptname + '  --help')
	print("python3 " + scriptname + " -i octopus.somatic.snvs.pass.vcf  --normal_column 10 --tumor_column 11 -o updated_vcf.vcf\n")
	print("#" * 40 + "\nWARNING WARNING: This script is to be used only with somatic snvs vcf having NORMAL sample in "
	                 "column 10 and TUMOR sample in column 11; if not like this, update your vcf file to these "
	                 "specs; and the vcf has to contain Phasing Information from Phaser with the PI and PW and PS fields in FORMAT.\n" + "#" * 40)

def parseArgs(scriptname, argv):
	new_vcf_name = None
	column_tumor, column_normal = None, None


	warnings.simplefilter(action='ignore', category=FutureWarning)
	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.WARNING)


	try:
		opts, args = getopt.getopt(argv,"hi:o:",[ "vcf=",
		                                            "tumor-col=", "normal-col=", \
		                                            "debug", "outfilename=", \
		                                            "help"] )
		log.info(opts)
		print("%" * 20 + "\n" + str(opts) + "\n" + "%" * 20)
	except getopt.GetoptError as go:
		log.error(go)
		usage(scriptname, opts)
		exit(2)
	for opt, arg in opts:
		if opt == '-h' or opt == "--help":
			usage(scriptname, opts)
			exit()
		elif opt in ("", "--tumor-col"):
			column_tumor = arg
		elif opt in ("", "--normal-col"):
			column_normal = arg
		elif opt in ("","--debug"):
			print("DEBUG")
			log.basicConfig(format=FORMAT_LOGGING, level=log.DEBUG)
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
			"Please Provide column number for tumor and normal Samples; should be 10 and 11  - or -  11 and 10 respectively;\n"
			"options are: --tumor_column and --normal_column;\nAborting. ")
	if column_normal == column_tumor:
		exit("ERROR: number for the columns Tumor and Normal MUST be different")


	return(fvcf, column_tumor, column_normal, new_vcf_name)

def collapse_variants(LOV):
	'''
	After going through several tests, it looks like the variants in the LOV object can be considered as
	a block-substitution variant, and therefore need to be collapsed into one line;
	To keep track of which variants record were collapsed we add the INFO field  BLOCSUBSFROM  which gathers
	a comma-separated list of CHROM_POS  variant record information found in the original octopus output file.
	The First record is taken as a template and the ALT is going to be replaced with the concatenated string
	of all the ALTs in the LOV records

	:param LOV:  List Variant Object (a.k.a variant record in vcf file
	:return: variant as a list of String ready to be written into the output vcf
	'''
	newREF = ""
	newALT = ""
	BLOCSUBSFROM= ""
	delim = ","
	for rec in LOV:
		newREF = newREF + rec.REF
		newALT=newALT+rec.ALT[0]
		locus = str(rec.CHROM)+"_"+str(rec.POS)
		BLOCSUBSFROM = delim.join([ BLOCSUBSFROM, locus ]) if BLOCSUBSFROM != "" else locus
	newV = str(LOV[0]).split("\t")
	newV[4-1] = newREF
	newV[5-1] = newALT
	newV[8-1] = ';'.join([ newV[8-1], "BLOCSUBSFROM="+BLOCSUBSFROM ])
	return str('\t'.join(newV))

def processing_variants_as_block_substitution(LOV, w):
		'''
		We consider that ALL the variants within the list_Of_Variant_Object in argument belong to the same BLOCK
		We need to "merge" the records to create only one.
		cyvcf2 does NOT allow to modify the ALT (attribute 'ALT' of 'cyvcf2.cyvcf2.Variant' objects is not writable.)
		The only thing we could think of so far with our basic knowledge of Python, is to convert the Variant object
		into a list of Strings (variant record is tab-separated) and modify the column 5 (index 4 in python) with
		the value of the ALT alleles from the variant in the list_Of_Variant_Object argument.
		The LAST test we need to perform before considering the variants as a BLOCK is to look at the AR
		If their AR is too FAR apart should we consider it as a block or not? ; we Hardcoded values; DP, MAP_VAF and is_indel
		:param LOV: List Of Variant object (a.k.a variant record in vcf file)
		:param w: writer object to write processed record into output vcf file
		'''

		## we will need to check if:
		## 1) we only deal with SNV in all ALTs ##TODO: should we tests for DEL in all ALTs or INS in all ALTs?? or a mix of SNV and Indels?
		## 2) need to check if AR within x% from each other calls otherwise might mean it is two UN-related events
		## 3) should we check for a minimum DEPTH to consider worth combining records or accurate and specs?
		## 4) create a new updated record using the first variant in the block
		##      -- update ALT field
		##      -- update DP field or ADD a field in info field that might help keeping track of data

		tests = []
		for rec in LOV: ## HARDCODED Here Below for TESTS  ; QUESTION: Why did we choose MAP_VAF>=0.2 ?? Oups
			if int(rec.format('DP')[1]) >= 10 \
					and float(rec.format('MAP_VAF')[1]) >= 0.20 \
					and not rec.is_indel:
				tests.append(True)
			else:
				tests.append(False)

		if all(tests):
			log.debug('\n'.join([ str(rec).strip() for rec in LOV ] ))
			w.write(collapse_variants(LOV))
		else:
			log.debug("BLOCK FAILED DP_TUMOR>=10 and/or MAP_VAF_TUMOR >= 0.20 or is_not_indel; So we did not collapse the variants in LOV obj: "+'\n'.join([str(rec) for rec in LOV]))
			for rec in LOV:
				w.write(str(rec))


def check_if_PS_in_FORMAT_field(vcf_cyobj, input_vcf_path, new_vcf_name):
	print("Checking PS flag presence in FORMAT ...")
	try:
		psid = vcf_cyobj.get_header_type('PS')
		print("psid == "+str(psid))
		if not psid['Description'] == '"Phase Set"':
			raise KeyError("PS Phase Set flag Absent in VCF; Aborting Recomposition of Records")
	except KeyError as ke:
		log.error("PS tag is not present in the FORMAT field of the VCF; We assume that the VCF has not been processed for phasing.\nvcf_in = {} \nvcf_out = {}".format(input_vcf_path, new_vcf_name))
		log.error(ke)
	print("PS flag FOUND in Header ...")


def check_for_block_substitution(vcf, column_tumor, w):
	print("looping over records to capture and concatenate Block Substitutions Variants ...")

	idxT = 0 if int(column_tumor) == 10 else 1

	dico_PS = defaultdict(list)
	k = -1 ## init the key value to -1 for the first next if statement

	for v in vcf:  ## v for variant which represents one "variant record" ; !!! WARNING: We can consume the iterator only once. Once that loop is done, the VCF object will be empty and we cannot loop over vcf object anymore
		## we gather the variant with the same PhaseSet Value
		new_k = '_'.join([ v.CHROM ,v.format('PS')[idxT] ])
		if k!=-1 and k != new_k and len(dico_PS[k])==1:
			## dico has only one item, so we write it to output vcf right away and re-init dico for olk key
			print("dico has only one item, so we write it to output vcf right away and re-init dico for olk key BEGIN:")
			print(str(v))
			print(str(k))
			print(str(dico_PS[k]))
			print("dico has only one item, so we write it to output vcf right away and re-init dico for olk key END")
			for rec in dico_PS[k]:
				w.write(str(rec))
			k = new_k
			dico_PS[k] = []
		elif k!=-1 and k != new_k and len(dico_PS[k])>1:
			## dico has only one item, so we write it to output vcf right away and re-init dico for olk key
			## this mean we already processed value for the previous PS value and we reached a new PS value, so the key k is different
			## we also check if the dico has more than one records if not we print record
			print("processing_variants_as_block_substitution BEGIN:")
			print(str(v))
			print(str(k))
			print(str(dico_PS[k]))
			print("processing_variants_as_block_substitution END")
			processing_variants_as_block_substitution(dico_PS[k], w)
			k = new_k
			dico_PS[k] = []  ## we reinit the dico for the current k value


		if len(dico_PS[k]) == 0:  ## if dico[k] empty we add the current variant
			print("dico[k] empty we add the current variant:")
			print(str(v))
			dico_PS[k].append(v)

		elif len(dico_PS[k]) == 1: ## this mean we already had one variant with same PS
			## we compare the variant POSITION(implies that the VCF was sorted by variant records)
			## if position is shifted by one, this means the vairant calls are right next to each other and we keep both in the dictionnary
			## if not, we only keep the current variant in the dico for that PS and continue to next variant in the loop

			## Let's get 'GT' for the current variant and the variant in position [-1] in dico
			## the variant of the Tumor should only be considered here
			## same genotype means similar to same PHASE
			lgenotype_current  = [str(Genotype(li)) for li in v.genotypes][1]
			lgenotype_previous = [str(Genotype(li)) for li in dico_PS[k][-1].genotypes][1]
			## we then compare the current locus position (v.POS) with the last locus added to the dictionary
			## if consecutives loci and genotype is identical, good chances we deal with variants that can
			## constitute a block-substitution
			if int(v.POS) == int(dico_PS[k][-1].POS)+1 and lgenotype_current == lgenotype_previous:
				dico_PS[k].append(v)  ## we gather the variant with the same PhaseSet Value, Same genotype in Tumor and consecutive POS
			else:
				## the variants in the dico either do not have the same genotype or were not consecutive
				## so we have to process what we have, and as it is only one variant in the dico,
				## we write it here, directly to the output file;
				print("variants in the dico either do not have the same genotype or were not consecutive BEGIN:")
				print(str(v))
				print(str(k))
				print(str(dico_PS[k]))
				print("variants in the dico either do not have the same genotype or were not consecutive END")
				for rec in dico_PS[k]:
					w.write(str(rec))
				dico_PS[k] = [ v ] ## we re-init the Value to the current variant as the previous variant is definitely not right next to the current one even though in the same PhaseSet

		elif len(dico_PS[k]) > 1:
			## this mean we already have at least two variants right next to each other in the dico
			## Either the 3rd or nth variant is right after the position of the latest variant in the list
			## or if not this implies to statements:
			## 1) the current variant does not belong to the same block-substitution as the last variant already in the Dictionary
			## 2) we have to process the variants currently in the dictionary has "potential" block-substitution and the current variant will be added to the dico_PS[k] replacing the old variants
				## 2a) if Genotype of each consecutive position is the same; it should be considered as a block
				## 2b) if Genotype of each consecutive position is NOT the same; we have to Process the variants that
				##      are already in the dictionary has belonging to the same Block-Substittion Variant
			lgenotype_current = [str(Genotype(li)) for li in v.genotypes][1]
			lgenotype_previous = [str(Genotype(li)) for li in dico_PS[k][-1].genotypes][1]
			if int(v.POS) == int(dico_PS[k][-1].POS) + 1 and lgenotype_current == lgenotype_previous:
				dico_PS[k].append(v)  ## we gather the variant with the same PhaseSet Value
			print("len(dico_PS[k]) > 1 BEGIN:")
			print(str(v))
			print(str(k))
			print(str(dico_PS[k]))
			print("len(dico_PS[k]) > 1 END")
			else:
				processing_variants_as_block_substitution(dico_PS[k], w)  ##dicoPS[k] is a list of Variants
				dico_PS[k] = [ v ]  ## we re-init the Value to the current variant as the previous variant is definitely not right next to the current one even though in the same PhaseSet

	## we need to write the data from the last phase set captured
	if len(dico_PS[k]) == 1:
		for rec in dico_PS[k]:
			w.write(str(rec))
	else:
		processing_variants_as_block_substitution(dico_PS[k], w)

def collapse_variants_with_same_PS(LOV):
	'''
	After going through several tests, it looks like the variants in the LOV object can be considered as
	a block-substitution variant, and therefore need to be collapsed into one line;
	To keep track of which variants record were collapsed we add the INFO field  BLOCSUBSFROM  which gathers
	a comma-separated list of CHROM_POS  variant record information found in the original octopus output file.
	The First record is taken as a template and the ALT is going to be replaced with the concatenated string
	of all the ALTs in the LOV records

	:param LOV:  List Variant Object (a.k.a variant record in vcf file
	:return: variant as a list of String ready to be written into the output vcf
	'''
	newREF = ""
	newALT = ""
	BLOCSUBSFROM= ""
	newAR = float(0)
	delim = ","
	# print("*"*200) ; 	print(str(LOV)) ; print("*"*200) ;
	for rec in LOV:
		newREF = newREF + rec.REF
		newALT=newALT+rec.ALT[0]
		# locus = str(rec.CHROM)+"_"+str(rec.POS)
	#	BLOCSUBSFROM = delim.join([ BLOCSUBSFROM, locus ]) if BLOCSUBSFROM != "" else locus
		# newARs = process_newAR() ## TODO: WE WILL NEED TO process AR for each sample
	newV = str(LOV[0]).split("\t")
	newV[4-1] = newREF
	newV[5-1] = newALT
	#newV[8-1] = ';'.join([ newV[8-1], "BLOCSUBSFROM="+BLOCSUBSFROM ])
	return str('\t'.join(newV))

def recompose_consecutive_blocks(vcf, column_tumor, w):
	"""

	:param vcf:
	:param column_tumor:
	:param w:
	:return: None
	"""
	log.info("looping over records to capture and concatenate Block Substitutions Variants ...")

	idxT = 1 if int(column_tumor) == 11 else 0
	dico_PS = defaultdict(list)
	# count = 0
	previous_pos = 0
	same_phased_but_not_consecutive = False ## spnc

	for v in vcf:
		if v.format('PS') is not None:
			curpos = v.POS
			k = int(v.format('PS')[idxT]) if not None else 0
			# print(str(k))
			if k in dico_PS and (curpos == previous_pos+1):
				if same_phased_but_not_consecutive:
					k = str(k) + "spnc"
				dico_PS[k].append(v)
			elif k not in dico_PS:
				dico_PS[k] = [v]
				same_phased_but_not_consecutive = False
			else: ## this mean we have the same phase but the position are not consecutives ... so we cannot make a block with previous position ## issue MMRF_2490_2_BM_Exome
				same_phased_but_not_consecutive = True
				k = str(k) + "spnc"
				dico_PS[k] = [ v ]
			# count += 1
			# if count >100:
			# 	break
			# print(dico_PS)
			previous_pos = curpos

		else:
			# print(100*"***"+"\nonly unphased variants here ...")
			#dico_PS[0].append(v)
			w.write(str(v))
			# print(dico_PS[0])
	# print(str(dico_PS[1][0]))
	# print(str(dico_PS[1][1]))
	for k in dico_PS:
		##dico_PS[k] is a list of variants object
		# print("k= "+str(k)+" <---> value = "+str(dico_PS[k]))
		w.write(collapse_variants_with_same_PS(dico_PS[k]))


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@##

#@#########
#@ MAIN  ##
#@#########
if __name__ == "__main__":

	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)


	vcf_path, column_tumor, column_normal, new_vcf_name = parseArgs(argv[0], argv[1:])  ; ## tth means tuple of thresholds

	vcf = VCF(vcf_path)
	tot_number_samples = len(vcf.samples)
	if tot_number_samples > 2 or tot_number_samples < 2:
		exit("ERROR: Number of Samples is different from 2; Expected 2 samples only; NORMAL in column 10 and TUMOR in column 11 (in that order in the VCF) ; Aborting Here.")

	if new_vcf_name is None:
		new_vcf = '.'.join([str(vcf_path), "recompose.vcf"])
	else:
		new_vcf = new_vcf_name


#	## checking if PS flag is still present in the VCF genotype fields

	check_if_PS_in_FORMAT_field(vcf, vcf_path, new_vcf_name)
	### WARNING WARNING: Because here we read at least once the vcf object which is a generator opbject, we LOSE the first variant;
	## we need to recreate the generator here
	vcf = VCF(vcf_path)
	## Adding Fields to INFO field
	# vcf.add_info_to_header({'ID': 'BLOCSUBSFROM', 'Description': 'List of Original Position in VCF That were used to make current block substitution variant', 'Type': 'String', 'Number': '.'})

	# create a new vcf Writer using the input vcf as a template.
	w = open(new_vcf, 'w')
	## writing the updated header immediately to the new output vcf file
	w.write(str(vcf.raw_header))


	log.info("looping over records ...")
	#check_for_block_substitution(vcf, column_tumor, w)
	# try:
	recompose_consecutive_blocks(vcf, column_tumor, w)
	# except KeyError as ke:
	# 	log.error("MISSING FIELD in VCF: \n"+ke)
	# 	exit(-1)

	w.close()
	vcf.close()

	exit()

