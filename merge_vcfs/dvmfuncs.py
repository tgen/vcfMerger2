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
### Major Contributors: Christophe Legendre
### Minor Contributors:


import logging as log
import re
import sys
import os
import csv
import subprocess
#from os import path (see line775)

from myGenotype import Genotype
from natsort import natsorted


def get_dictOfLoci(vcfToDict_instance):
	"""
	process the VCF file to capture ALL the loci (aka variants) in a dico and return a tuple containing one dictionary.
	:param vcfToDict_instance: an instance of vcfToDict Class
	:return: tuple with one dictionary of loci [[ KEY as CHROM_POS_REF_ALT, and VALUE as the associated line of the variant call
	"""
	return tuple(vcfToDict_instance.dictOfLoci(vcfToDict_instance.readVCF()),)

def compareTuples(LoT, msg):
	""" Compare Tuples Contents ; Here we compare the list of Samplenames from each VCF ; if they do not match, we exit; CASE SENSITIVE
	t1 = ('NORMAL', 'TUMOR') ; t2 = ('NORMAL', 'TUMOR')  ; LOT= [ t1 , t2 ] ; compareTuples(LOT)
	:param: LoT stands for List Of Tuple
	:param: msg is a string giving information about what we compare here, such as "toolnames" or samplenames
	:return: None
	"""
	x = 0
	for y in range(1,len(LoT)):
		if not LoT[x] == LoT[y]:
			log.info("\n" + "-"*41 + "\n!! VCF's "+ msg + " Differ; Aborting !!\n" + "-"*41 )
			log.info(str(x) + " vs. " + str(y))
			log.info(" <-- vs --> ".join([str(LoT[x]),str(LoT[y])]))
			if msg == "SampleNames":
				log.info("ERROR: SAMPLES should be in the EXACT same order and exact same CASE (case sensitive comparison)")
			if msg == "CONTIGS":
				log.info("ERROR: CONTIGS number or/and CONTIGS names are not identical from one vcf to another")
			log.info("Correct and modify the inputs VCF files to get them up to specs for the current merging tool")
			sys.exit(2)
	log.info("Comparison of << " + msg + " >> between VCFs: PASSED")

def get_list_contig_from_fasta_dict_file(dicofilename):
	'''
	Extract Conti Names from fasta.dict file and return the list of contigs ordered the same way it is in the fasta.dict file
	:param dicofilename:
	:return:  list of contig names
	'''
	lcontigs = []
	with open(dicofilename, newline='') as f:
		data = csv.reader(f, delimiter='\t')
		for line in data:
			if line[0] == '@SQ':
				lcontigs.append(str(line[1].replace("SN:", "")))
	return lcontigs

def addHeaderToOutVcf(t_obj, f):
		log.info("in function")

def prefix_headers_information_line_with_toolname(myHeaderString, toolname):
	"""
	Prefixing the Flags in the INFO fields with the name of the tool they come from
	:param myHeaderString: a string representing a Header line from a vcf header so should strart with 2 pound signs
	:param toolname: string name of the tool that called the variant
	:return: String
	"""
	return re.sub(r"<ID=", ''.join(["<ID=", toolname.upper(), "_"]), str(myHeaderString))

def prefix_headers_other_information_line_with_toolname(myHeaderString, toolname):
	"""
	Prefixing the Flags in the INFO fields with the name of the tool they come from
	:param myHeaderString: a string representing a Header line from a vcf header so should strart with 2 pound signs
	:param toolname: string name of the tool that called the variant
	:return: String
	"""
	return re.sub(r"^##", ''.join(["##", toolname.upper(), "_"]), str(myHeaderString))

def create_new_header_for_merged_vcf(tuple_objs, command_line, vcfMerger_Format_Fields_Specific, vcfMerger_Info_Fields_Specific, dico_map_tool_acronym, list_contig_from_fastadict_captured_as_is):
	'''
	Manage the Headers from all the input VCFs and recreate a brand new updated Header
	:param: tuple of vcf2dict objects which containg ALL the information about each input vcfs
	:param: String of the Command line captured by sys.argv value
	:param: vcfMerger_Format_Fields_Specific is a list of strings with the strings being ready to be added to header lines
	:param: vcfMerger_Info_Fields_Specific is a list of strings with the strings being ready to be added to header lines
	:return: new Fully Formated Header as String
	'''

	## here we will parse the object and capture header from each tool and ...
	## updating the header as necessary such as
	## prefixing the INFO IDS with toolname ;
	## we will also need to add the new header such as the command line that geenrated the out vcf file.
	## we will need to ad only the FORMAT field from the list of common field found in FORMAT

	log.info("creating new header")
	lh = [] ## list headers
	l_contigs = []
	## capture infos and assign values
	fileformat="##fileformat=VCFv4.2" ## harcoded
	from time import gmtime, strftime
	filedate = "##fileDate="+str(strftime("%Y%m%d", gmtime()))
	command_line='##cmdLine="'+command_line+'"\n'

	lh.append(fileformat)
	lh.append(filedate)
	## process contigs separately to remove duplicates
	for vtdo in tuple_objs:  ## list of vcfToDict objects
		# print("vtdo.contigs is of type : "+str(type(vtdo.contigs)))
		for contig in vtdo.contigs:
			l_contigs.append(contig)
	## removing duplicates with the set function
	l_contigs = set(l_contigs)
	## Manipulate l_contigs to have a sortable object by key and values
	dtemp = {} ## dico with key as contig names and values thetail of the string
	for item in l_contigs:
		strip_item = item.replace('##contig=<ID=', '').replace(">", '')   ## need to strip off the prefix and suffix
		if not "," in strip_item:
			strip_item = strip_item+","
		k, v = strip_item.split(',', 1)
		v = v + ">"
		if k in dtemp:
			dtemp[k].append(v)
		else:
			dtemp[k] = [v]
	## The Contigs are not well managed here; Need to Improve ##TODO
	## Here below we test if the values are more than one (should be one) and contains the keyword "length" as expected ;
	## If not, we should capture expection ##TODO
	for k, v in dtemp.items():
		if len(v)>1:
			for litem in v:
				if "length" in litem:
					newval = [ litem ]
					break
			dtemp[k] = newval


	## performing a sort of a dictionary with a list of contigs
	index_map = {v: i for i, v in enumerate(list_contig_from_fastadict_captured_as_is)}

	try: ## if an error is raised here, it is mostly because the a contig present in the input vcfs is absent from the fasta dictionnary file
		d3 = sorted(dtemp.items(), key=lambda pair: index_map[pair[0]])
	except KeyError as e:
		log.error("KeyError: ({0})".format(e))
		log.info("ERROR raised because a contig present in the input vcfs is actually absent from the given fasta dictionary file")
		exit()

	## rebuilding the contigs header lines after the correct sorting
	nlc = [] ## new list contig
	for pair in d3:
		nlc.append(''.join(['##contig=<ID=', pair[0], ",", str(pair[1][0])]))

	## adding the contigs to the list of strings called "lh" ; We DO NOT SORT or touch the list of contigs to keep the order defined in the fasta dictionary above
	for contig in nlc:
		lh.append(contig)
	## prefixing the header with the toolname, the same way the INFO Fields Flag are prefixed
	reference=""
	log.info("tuple_objs is length : {}".format(str(len(tuple_objs))))

	for vtdo in tuple_objs:  ## list of vcfToDict objects

		## capturing the ##reference informatino from the tool which has precedence
		if reference == "":
			indices = [i for i, s in enumerate(vtdo.header_other_info) if '##reference=' in s]
			if indices == None or len(indices) == 0:
				reference = ""
				# log.error("ERROR: Line ##reference is missing in your input vcf file for tool {}".format(vtdo.toolname) )
				# sys.exit(-1)
			else:
				reference = vtdo.header_other_info[indices[0]]
			log.info("reference is: {}".format(reference if reference != "" else "Reference Line Not Defined In {} Vcf ".format(vtdo.toolname)))

		toolname_or_acronym = get_acronym_for_current_tool(vtdo.toolname, dico_map_tool_acronym)
		for s in vtdo.header_filters:
			lh.append(prefix_headers_information_line_with_toolname(s, toolname_or_acronym))
		for s in vtdo.header_info:
			lh.append(prefix_headers_information_line_with_toolname(s, toolname_or_acronym))
		for s in vtdo.header_format:
			lh.append(prefix_headers_information_line_with_toolname(s, toolname_or_acronym))
		for s in vtdo.header_other_info:
			lh.append(prefix_headers_other_information_line_with_toolname(s, toolname_or_acronym))
		## if LOSSLESS, the column QUAL, FILTER, ID, and some others are ADDED to the variant record
		## this creates NEW fields prefixed with the toolname
		for COLUMN in ["FILTER", "QUAL", "ID"]:
			## ##INFO=<ID=SEURAT_AR1,Number=1,Type=Float,Description="Allele frequency of ALT allele in normal">
			stringline = ''.join(["##INFO=<ID=", toolname_or_acronym, "_", COLUMN,
			                      ',Number=.,Type=String,Description=',
			                      '"Represents lossless data from tool ',vtdo.toolname,' or (if given acronym: aka ', toolname_or_acronym,
			                      'for column ', COLUMN,'">'] )
			lh.append(stringline)
		## Here when LOSSLESS is enabled, fields that were in format of the secondary tools, are added to
		## the INFO field with the following format: TOOLNAME_Sx_FIELDofINTEREST
		## where x represents an indice of the Sample starting at 1 up to n.
		## if dealing with TUMOR_NORMAL, we should have S1 and S2 (respectively from column 10 and 11 in vcf)
		## Now we need to implement this here TODO: add the test if lossless enabled;
		##1) we capture the Format column, aka column number 9 for the current tool and prefixed it with tool names
		## and Sample number

		numberOfSamples = len(vtdo.samplenames)
		for S in vtdo.header_format:
			## return the first indice where the pattern is in the string
			idx1 = S.find(',')
			idx2 = S[:idx1].rfind("=")
			FIELD = (S[idx2+1:idx1])
			for i in range(1,numberOfSamples+1):
				newField = '_'.join([ toolname_or_acronym, "S"+str(i), FIELD ])
				# print(newField)
				stringline = ''.join(["##INFO=<ID=", newField,
			                      ',Number=.,Type=String,Description=',
			                      '"lossless data from defined tool">'])
				lh.append(stringline)


	for item in vcfMerger_Format_Fields_Specific:
		lh.append(item)
	for item in vcfMerger_Info_Fields_Specific:
		lh.append(item)
	if reference is not None or reference != "":
		lh.append(reference)

	lh.append(command_line)
	return lh  ## returns a list

def add_new_header_to_merged_file(f, list_lines_header, header_chrom_line):
	"""
	:param list_lines_header:  list of string objects ; each string is a line that has to be added to the merged vcf file ;
	 but the CHROM line should not be in it.
	:param f: string representing the filename of the output merged vcf file
	:param header_chrom_line: CHROM line to be added to the merged vcf file
	:return: None
	"""
	log.info("adding new_header to out vcf file")
	with open(f, "w") as outf:
		outf.write('\n'.join([str(item) for item in
		                    list_lines_header]));  ## header lines already have \n embedded in it so no need for extra
		outf.write(str(header_chrom_line))

def output_list_variant_sorted_by_contigs_as_same_order_as_in_fastdict_file(dd, ordered_l_contigs_ref_genome_fasta_dict):
	"""
	takes keys from an ordered dictionary which has our special keys formated as follow CHR__POS_REF_ALT
	and sort the keys according to the ordered list of contigs given and returns the ordered dictionary;
	Warning: sorting the contigs and the position in a different way that sorting by natural sort will not be compatible with bcftools sort
	:param ordered_l_contigs_ref_genome_fasta_dict:
	:return: ordered list of dictionary's keys
	"""
	dtemp = {}
	##1) we extract the contigs from the keys
	for item in dd.keys():  ## we capture the keys from dd dico and slipt the contig names from the pos_ref_alt string using the keyword '__'
		k, v = item.split('__', 1)
		if k in dtemp:
			dtemp[k].append(v)
		else:
			dtemp[k] = [v]
	##2) we sort the contigs using the ordered contig list
	index_map = {v: i for i, v in enumerate(ordered_l_contigs_ref_genome_fasta_dict)}
	try:  ## if an error is raised here, it is mostly because the a contig present in the input vcfs is absent from the fasta dictionnary file
		ordered_list_of_list = sorted(dtemp.items(), key=lambda pair: index_map[pair[0]])
	except KeyError as e:
		log.error("KeyError: ({0})".format(e))
		log.info(
			"ERROR raised because a contig present in the input vcfs is actually absent from the given fasta dictionary file")
		exit()
	##3) we rebuilt the keys from the ordered list of list
	d4 = []
	for key, value in ordered_list_of_list:
		## normally, the list representing the values are ordered in ascending order but if we want them reordered, please uncomment next line
		value = natsorted(value)
		for pos in value:
			d4.append(str(key) + "__" + str(pos))

	##4) again we have an ordered list of string which maps the list of keys found in dd, but probably in a different order
	index_map = {v: i for i, v in enumerate(d4)}
	try:  ## if an error is raised here, it is mostly because the a contig present in the input vcfs is absent from the fasta dictionnary file
		ordered_list_of_list = sorted(dd.items(), key=lambda pair: index_map[pair[0]])
	except KeyError as e:
		log.error("KeyError: ({0})".format(e))
		log.info(
			"ERROR raised WHY?? because a contig present in the input vcfs is actually absent from the given fasta dictionary file")
		exit()
	##5) final rebuilt of the ordered list of keys
	sod = []
	for key, value in ordered_list_of_list:
		sod.append(key)

	return sod

def isOfTYPE(ref, alt, v):  # we need to elaborate here ALL the possibilities to Define an INDEL, SNV or Complex
	# Variant ;
	''' Check what type of variant we are dealing with
	:param ref: REF column from the VCF (of type STRING)
	:param alt: ALT column from the VCF (as our prerequisite for the vcf is to be decomposed format, only one ALT
	should be available ; ALT is of type LIST  --->  according to cyvcf2 module
	:return: string value predefined (hardcoded; should be one of snv, indel, complex or unknown)
	'''
	### NOTE:
	### In our case we should have only ONE item within alt column as we have decomposed ALL the input vcf files (using
	#  vt). At least all the VCF that needed to be decomposed :-)
	### so alt should have only one item of variable length depending if this is indel of snv;
	### Remember, WE HAVE DECOMPOSED the VCF, therefore, only ONE Variant per line
	### If not decomposed, ALT or REF may have more than ONE variant per REF or ALT and therefore having comma
	### in it ; as  ALT is of type LIST, we use alt[0] to capture the item in the list.

	try:
		if alt == []:
			alt = "."
		log.debug("len(ref_==" + str(len(ref)) + " and len(alt)==" + str(len(alt[0])))
		log.debug("ref_==" + str(ref) + " and alt == " + str(alt[0]))
		if len(ref) > 1:
			log.debug("len(ref)>1")
			if alt[0] == ".":
				return "del"
			if len(ref) > len(alt[0]):
				log.debug("len(ref) == 1 and alt == . " + " --> del3")
				return "del"
			if len(ref) < len(alt[0]):
				return "ins"
			if len(ref) == len(alt[0]):
				return "blocksub"
		if len(ref) == 1 and alt == ".":
			log.debug("len(ref) == 1 and alt == . " + " --> del3")
			return "del"
		if len(ref) == 1 and len(alt[0]) > 1:
			log.debug("len(ref) == 1 and len(alt[0]) > 1" + " --> ins2")
			return "ins"
		if len(ref) == 1 and len(alt[0]) == 1:
			log.debug("len(ref_==1 and len(alt)==1")
			if len(alt[0]) == 1 and alt[0] is not "." \
					and (ref.capitalize() in ['A', 'T', 'C', 'G']) \
					and (alt[0].capitalize() in ['A', 'T', 'C', 'G']):
				return "snv"
			if alt[0] is "." and ref.capitalize() in ['A', 'T', 'C', 'G']:
				return "del"
			if len(alt[0]) > 1 and ref.capitalize() in ['A', 'T', 'C', 'G']:
				return "ins"
			return "unknown"
		else:
			log.debug("len(ref_==" + str(len(ref)) + " and len(alt)==" + str(len(alt[0])))
			return "complex"

	except Exception as e:
		log.info("ERROR in function << isOfType >> ({0}): {1}".format(e.errno, e.strerror))
		return "unknown"

def renameINFO_IDS(obj_variant, toolname):
	"""
	Prefixing the Flags in the INFO fields with the name of the tool they come from
	:param obj_variant: cyvcf2 object repreent the variant object
	:param toolname: string name of the tool that called the variant
	:return: String
	"""
	return re.sub(r"^", ''.join([toolname.upper(), "_"]), str(obj_variant).split("\t")[8-1].replace(
		";", ''.join([";", toolname.upper(), "_"]))).strip(';')

def addINFO_FromSecondaryToolsToNewRebuiltINFO_Field(currentNewRebuiltINFO, obj_variant, toolname):
	"""
	Get the Format field information from the remaining tools and add them to the newINFO column
	:param obj_variant: cyvcf2 object represent the variant object
	:param toolname: string name of the tool that called the variant
	:param currentNewRebuiltINFO:
	:return: NewCurrentRebuiltINFO_String
	"""
	#column_number=8 ; indice = column_number-1 ; ## python is 0-based ; column 8 is INFO column in VCF
	delim = "" if len(currentNewRebuiltINFO) == 0 else ";" ;
	newInfoToAppend = re.sub(r"^", ''.join([ toolname.upper() , "_"]),
	              str(obj_variant).split("\t")[(8-1)].replace(";",''.join([";", toolname.upper(), "_"])))
	return delim.join([currentNewRebuiltINFO, newInfoToAppend])

def addFORMAT_FromSecondaryToolsToNewRebuiltINFO_Field(currentNewRebuiltINFO, tv, toolname, totnum_samples):
	"""
	Get the Format field information from the remaining tools and add them to the newINFO column
	:param: currentNewRebuiltINFO:
	:param: tv or obj_variant: cyvcf2 object represent the variant object for current tool
	:param: toolname: string name of the tool that called the variant
	:param: number of samples in the VCF
	:return: NewCurrentRebuiltINFO_String
	"""
	#column_number=9 ; indice = column_number-1 ; ## python is 0-based ; column 9 is FORMAT column in VCF
	#Columns 10 and beyond are the format for each sample in the VCF
	delim = ";" # if len(currentNewRebuiltINFO) == 0 else ";" ;
	list_format_fields = str(tv).split("\t")[(9 - 1)].split(":")
	for idx_col_sample in range(1, totnum_samples+1): ## loop over the column that represent the SAMPLES
		list_value_current_field = str(tv).split("\t")[((9 + idx_col_sample) - 1)].split(":") ## column 10 is the
		# first sample, column 11, the second sample, etc.
		#log.info("VFF --> " + str(list_value_current_field))
		log.debug(toolname + "  VFF --> " + str(list_value_current_field))

		for idx_field in range(0, len(list_format_fields)-1):
			prefix_name="_".join( [ toolname, "S"+str(idx_col_sample), list_format_fields[idx_field].strip("\n")] )
			value_associated_to_prefix_name = list_value_current_field[idx_field].strip("\n")
			currentNewRebuiltINFO = delim.join([currentNewRebuiltINFO,
			                                    "=".join([str(prefix_name), str(value_associated_to_prefix_name)])
			                                      ])

	return currentNewRebuiltINFO.strip(';')

def add_ID_FILTER_QUAL_FromSecondaryToolsToNewRebuiltINFO_Field(currentNewRebuiltINFO, tv, toolname):
	"""
	Get the Format field information from the remaining tools and add them to the newINFO column
	:param: currentNewRebuiltINFO:
	:param: tv or obj_variant: cyvcf2 object represent the variant object for current tool
	:param: toolname: string name of the tool that called the variant
	:return: RebuiltINFO_String
	"""
	#column_number=9 ; indice = column_number-1 ; ## python is 0-indexed ; column 9 is FORMAT column in VCF
	#Columns 10 and beyond represent the format values for each sample in the VCF
	delim = ";"
	prefix_name_ID = "_".join([toolname, "ID" ])
	ID = str(tv).split("\t")[(3 - 1)][0]
	prefix_name_QUAL = "_".join([toolname, "QUAL"])
	QUAL = str(tv).split("\t")[(6 - 1)][0]
	prefix_name_FILTER = "_".join([toolname, "FILTER"])
	FILTER = str(tv).split("\t")[(7 - 1)].replace(";",",")
	#local_rebuilt = ";".join([ "=".join([ prefix_name_ID, ID ]), "=".join([ prefix_name_QUAL, QUAL ] ) ])
	local_rebuilt = str(prefix_name_ID) + "=" + str(ID) + ";" + \
	                str(prefix_name_QUAL) + "=" + str(QUAL) + ";" + \
	                str(prefix_name_FILTER) + "=" + str(FILTER)
	return delim.join([currentNewRebuiltINFO, local_rebuilt]).strip(';')

def process_GT_field(field, totnum_samples, dicField, v):
	'''
	:param: field to Process
	:param: Total number of sample in the VCFs (number should be identical between the vcfs)
	:param: dictField is a dictionary of RECONSTRUCTED values for the ALL the field in FORMAT ; We therefore update
	this	dictionnary everytime we encounter a new Field in the v.FORMAT value
	:param: v stands for variant and is the object variant created by cyvcf2 module
	:return: the updated dictField
	'''
	if field not in v.FORMAT:
		for i in range(totnum_samples):
			if not len(dicField[i]) == 0:
				dicField[i] = ":".join([dicField[i], str("./.")])
			else:
				dicField[i] = str("./.")
	else:
		log.debug("%" * 10 + " IN ELSE " + "%" * 10)
		log.debug(v)
		log.debug(v.format('GT'))
		log.debug(v.FORMAT)
		log.debug(v.genotypes)

		gts = v.genotypes
		log.debug("GENOTYPE gts is:" + str(gts))
		my_genotypes = [Genotype(li) for li in gts]

		log.debug("TYPE of my_genotypes :::: " + str(type(my_genotypes)))
		log.debug("MY GENOTYPES ARE: ")
		log.debug(my_genotypes)
		log.debug(str(my_genotypes[0]))

		log.debug(str(range(totnum_samples)))

		for i in range(totnum_samples):
			log.debug(" index i == "  + str(i))
			if not len(dicField[i]) == 0:
				# dicField[i] = ":".join([dicField[i], str('/'.join(["x", "x"]))])
				dicField[i] = ":".join([dicField[i], str(my_genotypes[i])])
			else:
				# dicField[i] = str('/'.join(str(e) for e in ["y","y"]))
				dicField[i] = str(my_genotypes[i])

def process_Unknown_field(field, totnum_samples, dicField, v):
	'''
		:param: field to Process
		:param: Total number of sample in the VCFs (number should be identical between the vcfs)
		:param: dictField is a dictionary of RECONSTRUCTED values for the ALL the field in FORMAT ; We therefore update
		this	dictionnary everytime we encounter a new Field in the v.FORMAT value
		:param: v stands for variant and is the object variant created by cyvcf2 module
		:return: the updated dictField
	'''
	if field not in v.FORMAT:
		log.debug("%" * 10 + " UNKOWN not FOUND in FORMAT " + "%" * 10)
		for i in range(totnum_samples):
			if not len(dicField[i]) == 0:
				dicField[i] = ":".join([dicField[i], str(".")])
			else:
				dicField[i] = str(".")

def process_known_field(field, totnum_samples, dicField, v):
	'''
	:param: field to Process
	:param: Total number of sample in the VCFs (number should be identical between the vcfs)
	:param: dictField is a dictionary of RECONSTRUCTED values for the ALL the field in FORMAT ; We therefore update
	this	dictionnary everytime we encounter a new Field in the v.FORMAT value
	:param: v stands for variant and is the object variant created by cyvcf2 module
	:return: the updated dictField
	'''
	nfor = v.format(field).tolist() ## return a numpy array  ; we need to manage this array for each recaptured field
	log.debug("In process known for flag << "+ field +" >> :  "+ str(nfor))
	for i in range(len(nfor)):
		if not len(dicField[i]) == 0:
			dicField[i] = ":".join([dicField[i], str(','.join(str(e) for e in nfor[i] if e is not ","))])
		else:
			dicField[i] = str(','.join(str(e) for e in nfor[i] if e is not "," ))
		log.debug(str(dicField[i]))

def process_extra_format_fields_from_winner_tool(currentNewRebuiltINFO, field, totnum_samples, tv, toolname):
	'''
		:param: field to Process
		:param: Total number of sample in the VCFs (number should be identical between the vcfs)
		:param: currentNewRebuiltINFO: string of the new rebuilding INFO field
		:param: v stands for variant and is the object variant created by cyvcf2 module
		:return: the updated currentNewRebuiltINFO variable (aka String)
		'''

	delim = ";" # if len(currentNewRebuiltINFO) == 0 else ";" ;
	for idx_col_sample in range(1, totnum_samples+1): ## loop over the column that represent the SAMPLES
		list_values_current_field=tv.format(field).tolist()## return a numpy array  ; we need to manage this array
		# for each recaptured values

		for idx_field in range(0, totnum_samples):
			prefix_name="_".join( [ toolname, "S"+str(idx_col_sample), field] )
			if isinstance(list_values_current_field[idx_field], list):
				value_associated_to_prefix_name = str(list_values_current_field[idx_field][0])
			else:
				value_associated_to_prefix_name = str(list_values_current_field[idx_field])

			if value_associated_to_prefix_name == "nan":
				value_associated_to_prefix_name="."
			currentNewRebuiltINFO = delim.join([currentNewRebuiltINFO,
			                                    "=".join([str(prefix_name), str(value_associated_to_prefix_name)])
			                                    ])
	return currentNewRebuiltINFO.strip(';')

def get_acronym_for_current_tool(tool, dico_map_tool_acronym):
	"""
	return the given acronym of toolname if acronym for that tool got defined; otherwise return toolname
	:param tool:
	:param map:
	:return: string_of_interest [ acronym or toolname ]
	"""
	d = dico_map_tool_acronym
	if tool in d.keys() and ( d[tool] is not None or d[tool] is not ""):
		log.debug("RETURN ACRONYM: " + d[tool])
		return d[tool]
	else:
		log.debug("RETURN TN: " + tool )
		return tool

def rebuiltVariantLine(LV, dico_map_tool_acronym, lossless, list_Fields_Format, totnum_samples = 1):

	"""
	built a lossless variant line using all the variant INFO from each tool.
	:param LV: variant information object for the Winner Tool; v represents an entire line of a vcf
	:param dico_map_tool_acronym is a dico with keys as toolname and values as acronyms
	:param: lossless, boolean; if True, will keep ALL the information from teh secondary tools and add that
	information to the INFO field in the  merged VCF
	:param  As this version is dedicated to Somatic, we expect two Sample by default in the FORMAT Columns.
	:param: list_Fields_Format, Format Fields that wll be present in the Merged VCF ; if no information about hte
	field exist in any of the tools' FORMAT columns, the field is set with dots to indicat absence of information or value
	:return: String of the new whole variant line
	"""
	'''
	wtv is the variant information object for the Winner Tool [aka tool with highest precedence , wtv ];
	newinfo is the new info field for that variant ;
	The format field should be captured from the tool having precedence because normally ALL the input vcfs
	should be up to specs or should have at least be updated to be up to specs for using this current merging tool.
	So the FORMAT field should already contains the FLAGS that are common to ALL the inputs VCF; FORMAT column should be identical between the input VCFs
	'''

	# LV = List Variant from each Tool @ the Same Position;  calls are from minimum 1 to maximum N tools at the same
	# position;
	wtn = LV[0][0]  ; ## wtv stands for Winner Tool Name
	wtv = LV[0][1]  ; ## wtv stands for Winner Tool Variant (here Variant refers to the Object Variant Created by
	# cyvcf2

	## 1) we first deal with rebuildign INFO field
	## Lets' deal with the winner tool first

	newINFO = renameINFO_IDS(wtv, get_acronym_for_current_tool(wtn, dico_map_tool_acronym))  ## LV[0][0] is the toolname which has precedence for that current call

	INFOfromOtherTools = ""

	# print("\n" + "#" * 50 + "\nnew ------> ALT is :" + str(wtv.ALT))

	## we process the Winner Tool first

	# for field in wtv.FORMAT:
	# 	if field not in list_Fields_Format:  ## we add the data that are not going to be in the new FORMAT
	# 		INFOfromOtherTools=process_extra_format_fields_from_winner_tool(INFOfromOtherTools, field, totnum_samples, wtv, wtn)

	## Here we will capture the INFO field from the other tools because they came second in the list;
	## we process the calls from the other tools
	## we add their data from ALL the fields and column if lossless ; or just from their INFO field otherwise
	for i in range(1, len(LV)):  # i represent the index of the variant in the list LV for current
		# position ; If we have only one tool which called that variant, we never enter this for loop;
		tn = get_acronym_for_current_tool(LV[i][0], dico_map_tool_acronym)  # In the list of 2 elements, the first is the toolname to remember the file of origin
		tv = LV[i][1]  # the second element is the variant object (tv stands for tool variant)
		INFOfromOtherTools = addINFO_FromSecondaryToolsToNewRebuiltINFO_Field(INFOfromOtherTools, tv, tn)

		if lossless:
			INFOfromOtherToolsForLossless = ""
			INFOfromOtherToolsForLossless = addFORMAT_FromSecondaryToolsToNewRebuiltINFO_Field(
				INFOfromOtherToolsForLossless, tv, tn, totnum_samples)
			INFOfromOtherToolsForLossless = add_ID_FILTER_QUAL_FromSecondaryToolsToNewRebuiltINFO_Field(
				INFOfromOtherToolsForLossless, tv, tn)
			INFOfromOtherTools=";".join([INFOfromOtherTools, INFOfromOtherToolsForLossless])

	# DELETE orints HERE was for debugging the GENOTYPEs management
	# gts = wtv.genotypes
	# print(str(type(gts)))
	# print("GENOTYPE gts is:")
	# print(str(gts))
	# print(str(gts[0]))
	# print(str(type(gts[0])))


	#ListFieldsToProcessForOurFORMATColumn=[ "GT" , "DP" , "AR", "AD" ] #, "AU", "CU", "TU", "GU" ]
	ListFieldsToProcessForOurFORMATColumn=list_Fields_Format

	if ListFieldsToProcessForOurFORMATColumn[0] is not "GT":
		raise("ERROR: GT field MUST be the first field in the FORMAT column according to VCF Specs; Aborting Merging")

	### same comment as earlier, the list of FIELDS or FLAGS in that FORMAT column may be set up by USER in tool wrappers
	### or we lock it down to specific flags such as GT, DP, AF, and AR.

	#format(self, field, vtype=None)
	#format returns a numpy array for the requested field.
	#The numpy array shape will match the requested field. E.g. if the fields has number=3, then the shape will be (n_samples, 3).

	## we need the number of total samples managed here ...  and the number of values defined in the header for that field in order to manage the number
	## ; This is shit ; they should have given us the option of capturing only the format field itself in this cyvcf2.
	## otherwise we may consider that the number of samples captured by getSampleNames in vcfToDict class can be captured here.
	##For now let's assign it manually:

	lFieldsForColumn9 = ""
	dicField = {} ## to capture and represent the value in current processed field for each Sample in VCF

	for i in range(totnum_samples):
		dicField[i] = ""

	## processing the FORMAT using the Winner tool
	for field in ListFieldsToProcessForOurFORMATColumn:
		## we are just rebuilding the column FORMAT with FLAGS names, i.e. column #9 in the VCF given in variable
		# above ListFieldsToProcessForOurFORMATColumn=[ "GT", "DP", "AR", "AD" ]
		lFieldsForColumn9 = ":".join([lFieldsForColumn9, field]) if not len(lFieldsForColumn9) == 0 else "".join([lFieldsForColumn9, field])

		## We First Process the Tool that gets the Precedence (aka, winner tool according to user given precedence
		## order
		## wtv represents the variant line of VCF object captured using the cyvcf2 module
		## wtv.format captures the format fields and flags for current variant (see cyvcf2 module)
		if field == "GT":
			process_GT_field(field, totnum_samples, dicField, wtv)
			continue

		if field not in wtv.FORMAT:
			process_Unknown_field(field, totnum_samples, dicField, wtv)
			continue

		## we process the field that exists in the vcf of the tool which got precedence for that variant
		process_known_field(field, totnum_samples, dicField, wtv)

	## CAPTURING ALL THE INFORMATION FOR ALL THE COLUMNS IN the FINAL MERGED VCF
	## get the new column number 9
	newFORMATcolumn=lFieldsForColumn9
	## get the STRING for the column 10 and beyond ; Sorting MUST be performed to have the Keys and therefore the
	# samples in the correct order ; We used an integer as key sort easy to sort;
	newFORMATusingWinnerTool="\t".join([ str(dicField[k]) for k in sorted(dicField.keys()) ])

	## The primary tool is taken care in the parent loop and after that loop
	CC = '='.join(["CC", str(len(LV))])  ## CC is the flag that stands for CALLERS_COUNT
	TPCE = '='.join(["TPCE", wtn])  ## TPCE stands for TOOL PRECEDENCE ;
	# TOOLS_AS_FLAG=";".join([wtn, ";".join([ LV[i][0] for i in range(1, len(LV)) ] ) ])
	TOOLS_AS_FLAG = ";".join([LV[i][0] for i in range(0, len(LV))])


	### We also could keep the original user's input that tell us about the order of precedence chosen by the user;
	### we will add this to the VCF header as a comment: ##TOOL_PRECEDENCE='>'.join([str(item) for item in ltoolnames])

	## wtv.ALT = [ '.' ] ## strangely, if the ALT in vcf is a dot, cyvcf2 retrun an empty list. I am sure I am
	## missing something in the concept of returning an empty list instead of the dot itself.
	## attribute 'ALT' of 'cyvcf2.cyvcf2.Variant' objects is not writable

	TYPE = str("TYPE=" + isOfTYPE(wtv.REF, wtv.ALT, wtv))
	newINFO = ';'.join([x for x in [TYPE, CC, TPCE, TOOLS_AS_FLAG, newINFO, INFOfromOtherTools] if (x != "")])


	CHROM = wtv.CHROM if wtv.CHROM != None else sys.exit("CHROMOSOME Unknown for position " + str(wtv.POS) + "; "
	                                                                                                       "Aborting")
	POS = wtv.POS if wtv.POS != None else sys.exit("POSITION Unknown for position " + str(CHROM) + "; Aborting")
	ID = wtv.ID if wtv.ID != None else "."
	REF = wtv.REF if wtv.REF != None else sys.exit("REFERENCE Nucleotide Unknown for position " + str(CHROM) + ":" + str(POS) + "; Aborting")

	if wtv.ALT == []:
		ALT =  "."
	elif len(wtv.ALT)==1:
		ALT = str(wtv.ALT[0])
	else:
		ALT = ','.join([ x for x in wtv.ALT ])
	log.debug(str(ALT))
	if not ALT:
		sys.exit("ALT Nucleotide Unknown for position " + str(CHROM) + ":" + str(POS) + "; Aborting")
	QUAL = wtv.QUAL if wtv.QUAL != None else "."

	#prefixing each FILTER information with toolname or acronym given by the user
	log.debug(str(wtv.FILTER))
	FILT = ';'.join([get_acronym_for_current_tool(wtn, dico_map_tool_acronym) +"_"+FLTV for FLTV in wtv.FILTER.split(';') ]) if wtv.FILTER != None else "."
	#FILT = wtv.FILTER if wtv.FILTER != None else "."  ## note: cyvcf2 consider the keyword PASS as 'no filter' and
	# therefore, when PASS is present, wtv.FILTER returns  None, which in our case here is a dot or will be force to
	#  be a dot.

	elements = [str(e) for e in [CHROM, POS, ID, REF, ALT, QUAL, FILT, newINFO, newFORMATcolumn,
	                             newFORMATusingWinnerTool]]

	return '\t'.join(elements)

def get_colors_for_venns(number):
	dictColors = {
		2: ("#d8b365", "#5ab4ac"),
		3: ("#dfc27d", "#80cdc1", "#018571"),
		4: ("#a6611a", "#dfc27d", "#80cdc1", "#018571"),
		5: ("#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e"),
		6: ("#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", "#5ab4ac", "#01665e"),
		7: ("#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e"),
		8: ("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e"),
		9: ("#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30"),
		10: ("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")
	}
	try:
		return dictColors[number]
	except KeyError:
		log.info("ERROR: Invalid Number; expected Integer between 2 and 10; Aborting Venn Diagram Creation")
		sys.exit("ERROR: Invalid Number; expected Integer between 2 and 10; Aborting Venn Diagram Creation")



def make_venn(ltoolnames, lbeds, variantType="Snvs_and_Indels", venn_title="", saveOverlapsBool=True,
              upsetBool=False, dirout=None, prefixPngFilenames="vcfMerger2"):

	names = ','.join([name for name in ltoolnames])
	## TODO we could check if any of the tools or any of the vcfs filenames already contains a comma; if so raise error


	numberOfTools = len(ltoolnames)
	if len(lbeds) != numberOfTools:
		log.info("WARNING: Number of Tools and number of BED files do NOT match; we skip the creation of Venn")
		return 0

	type = "genomic"

	## as intervene only uses the bed coordinates to get the sets, if two same positions are listed but one is a SNV and
	## the second is an INDEL, intervene considers that case as only ONE item; so we lose info here; we need to
	## have the REF and ALT taken into account as well;
	## so we convert our 5-columns bed files  into a simple list of strings and reassigned the new files to the lbeds object
	bedTolist = True;
	if bedTolist:
		newFileList = []
		for bed in lbeds:
			# mycmd = [ "sed ", "", " 's/\\t/_/g' ", bed, " > ", bed+".asList " ]
			pattern = "\\t"
			replacement = "_"
			full_string_sed = ["s/" + pattern + "/" + replacement + "/g"]
			args = ["sed"] + full_string_sed + [bed]
			log.info(str(args))
			log.info(" ".join([x for x in args]))
			# check_output will run the command and store to result
			log.info("*" * 50)
			log.info("full command to convert bed to StringList")
			fout = open(bed + ".asList", "w")
			process = subprocess.Popen(args, shell=False, universal_newlines=False, stdout=fout)
			process.wait()
			log.info(str(process.returncode))
			if process.returncode is not 0:
				sys.exit("Conversion BED to LIST FAILED; Aborting.")
			newFileList.append(bed + ".asList")
			fout.close()
		lbeds = newFileList
		type = "list";  ## need to update the type here as now we have a list of string in the file for intervene
		log.info("new list of inputs for intervene:")
		log.info(",".join(lbeds))

	colors = list(get_colors_for_venns(numberOfTools))

	title = venn_title + " [ " + variantType + "]\""
	figtype = "png"
	dpi = 300
	bordercolors = ["black"] * numberOfTools
	fontsize = 20
	project = prefixPngFilenames + "." + str(numberOfTools) + "_tools_" + variantType.replace(" ", "_") + "." + str(figtype) ;  ## this is actually the name of the png image file while the output_name is the folder where the intervene results are going into

	# Define the type of venn
	if numberOfTools >= 5:
		upsetBool = True

	#output_name = "upsetPlot_" + str(numberOfTools) + "_tools_"+ variantType.replace(" ", "_")  if upsetBool else "venn_" + str(numberOfTools) + "_tools_" + variantType.replace(" ", "_") ## this is the name of the directory created by intervene where the filename (aka project) above will be in.
	output_name = "SummaryPlots_" + str(numberOfTools) + "_tools_"+ variantType.replace(" ", "_") ## this is the name of the directory created by intervene where the filename (aka project) above will be in.
	import os
	if dirout is None:
		dirout = os.path.basename(os.path.curdir)
	output_name = os.path.sep.join([dirout, output_name])

	# Define command and arguments
	command = 'intervene'

	## common arguments
	common_args = ["--type", type,
	        "--names", names,
	        "--figtype", figtype,
	        "--dpi", str(dpi),
	        "--project", project,
	        "--output", output_name
	        ]


	if upsetBool:
		vtype = "upset"
		type_specific_additional_args = ["--ninter", "50",
		                                 "--sbcolor", "#d8b365",
		                                 "--mbcolor", "#5ab4ac",
		                                 "--showzero",
		                                 "--order", "freq", "--scriptonly" ]
		# type_specific_additional_args = type_specific_additional_args + ["nsets=2"]
		if saveOverlapsBool:
			type_specific_additional_args = type_specific_additional_args + ["--save-overlaps"]

	else:
		vtype = "venn"
		type_specific_additional_args = ["--bordercolors",  ",".join([color for color in bordercolors]),
		                                  "--colors", ','.join([color for color in colors]),
		                                 "--fontsize", str(fontsize),
		                                 "--title", title]
		if saveOverlapsBool:
			type_specific_additional_args = type_specific_additional_args + ["--save-overlaps"]




	## Build subprocess command
	mycmd = [command, vtype]
	mycmd = mycmd + common_args + type_specific_additional_args
	mycmd = mycmd + ["--input"] + lbeds

	list_commands = []
	# for C in [mycmd, mycmdsnvs, mycmdindels]:
	for C in [mycmd]:
		if C is not None:
			list_commands.append(C)

	for mycmd in list_commands:
		log.info(str(mycmd))
		log.info(" ".join([x for x in mycmd]))
		# check_output will run the command and store to result
		#import subprocess
		log.info("*" * 50)
		log.info("full command run intervene")
		process = subprocess.Popen(mycmd, shell=False, universal_newlines=False)
		process.wait()
		log.info(str(process.returncode))
		if process.returncode is not 0:
			sys.exit("Venn or Upset Creation FAILED")


	## update Rscript to colorize the intersection of all tools
	if upsetBool:
		## before doing any further steps, checking here if at least one of the bed file has NO variants; If so, UpSetR package has a bug when dealing with sets with zeros.
		## in order to avoid getting an error and breaking any pipeline, we will skip Running the Rscript;
		## but we need to fake the Folders and Files in order to keep consistency (mostly for pipelines)
		for bedfile in lbeds:
			## fast way of counting lines without loading file in memory
			number_lines = sum(1 for line in open(bedfile))
			if number_lines == 0:
				cmd_touch = "touch " + output_name + os.path.sep + project + "_PlotNotGeneratedBecauseFoundToolWithoutVariants.png"
				log.info(str(cmd_touch))
				log.warning("WARNING: Upset plots NOT generated because one of the given tool returned no-variant; Please Check the Rscript to generate manually the plot for "+output_name)
				os.system(cmd_touch)
				return None


		## Processing with the Rscript created by Intervene for UpSet plots
		list_tools = ",".join(["\""+tool+"\"" for tool in ltoolnames ])
		pattern = "nsets"
		replacement = "queries=list(list(query=intersects, params=list("+list_tools+"),color=\"red\", active=T)), matrix.dot.alpha=0.5, point.size = 3, nsets"

		full_string_sed = [ "s/" + pattern + "/" + replacement + "/" ]
		scriptname = project + "_upset.R"
		file_path = [ output_name + "/" + scriptname ]

		log.info(replacement)
		log.info(file_path)
		args = ["sed", "-i" ] + full_string_sed +  file_path
		log.info(str(args))
		log.info("Running Sed command")
		process = subprocess.Popen(args,  shell=False, universal_newlines=False)
		process.wait()
		if process.returncode is not 0:
			sys.exit("Upset Creation FAILED")
		log.info("Running Rscript Command")
		import os
		log.info(os.path.abspath("."))
		process = subprocess.Popen(str("Rscript " + file_path[0]),  shell=True, universal_newlines=False)
		process.wait()
		log.info("process return code: " + str(process.returncode))
		if process.returncode is not 0:
			log.info("return code not zero in << if upsetBool>>")
			sys.exit("Upset Creation FAILED")

	if not upsetBool:
		## this means with deal with Venns and therefore teh sets directory got created with too limited access permissions we need to change
		modify_acces_permissions_to_files_recursively(output_name)

	## annotate the images created by make_venn function
	## We expect at least three files, snvs+indels, snvs_only and indels_only

	project = os.path.splitext(project)[0]+"."+figtype+"_"+vtype+os.path.splitext(project)[1]
	log.info("output_name = "+ output_name +"/"+ project)
	path_to_image_file=os.path.realpath(os.path.join(output_name,project))
	log.info("Venn Annotation in progress for ... "+path_to_image_file)
	add_annotation_to_image(path_to_image_file, ltoolnames, lbeds)


def modify_acces_permissions_to_files_recursively(path):
	import os
	for root, dirs, files in os.walk(path):
		for d in dirs:
			os.chmod(os.path.join(root, d), 0o755)
		for f in files:
			os.chmod(os.path.join(root, f), 0o755)

def get_os_fonts():

	if sys.platform.startswith('linux'):
		fonts_sys_dir = "/usr/share/fonts/gnu-free"
		res = os.walk(fonts_sys_dir)
	elif sys.platform.startswith('win'):
		fonts_sys_dir = os.path.join(os.environ['WINDIR'], 'Fonts')
		res = os.walk(fonts_sys_dir)
	elif sys.platform.startswith('darwin'):
		fonts_sys_dir = "/Library/Fonts/"
		res = os.walk(fonts_sys_dir)
	log.info("\n\nAvailable Fonts on Your System:")
	log.info("\n".join(sorted([ str(f) for f in next(res)[2]])))
	log.info("\n")

def get_os_specific_system_font(my_os, user_font=None):
	fpf = None
	try:
		if user_font is not None and os.path.exists(user_font):  ## relative or full path
			return user_font
		if my_os.startswith("linux"):
			default_dir_font_lin = "/usr/share/fonts/gnu-free"
			default_font_lin = "FreeSansBold.ttf"
			fpf = os.path.join(default_dir_font_lin, default_font_lin)
		elif my_os.startswith("darwin"):
			default_dir_font_mac = "/Library/Fonts/"
			default_font_mac = "Arial Bold.ttf"
			fpf = os.path.join(default_dir_font_mac, default_font_mac)
		elif my_os.startswith("win"):
			default_dir_font_win = os.path.join(os.environ['WINDIR'], 'Fonts')
			default_font_win = "ARLRDBD.TTF"
			fpf = os.path.join(default_dir_font_win, default_font_win)

		if fpf is not None and os.path.exists(fpf):
			log.info("Found Font: " + fpf + " and PATH to it EXISTS")
			return fpf
		else:
			try:
				get_os_fonts()
			except Exception as e:
				log.error(e)
			log.error(
				"ERROR: vcfMerger2 could not determine automatically either the correct system, "
				"could not find the correct font directory or could not find the default preset system font;\nAborting;\nPlease "
				"install < " + fpf + " > equivalent on your system ")
			raise FileNotFoundError
	except (TypeError, AttributeError) as n:
		log.error(n)
	except FileNotFoundError as fnf:
		log.error(fnf)

def check_if_files_in_list_exist(list_file):
	"""
	Check if the file in the given list of file exit and return only the files that exist or empty list if none exit
	:param list_file: object list of files
	:return: list_of_file or empty list
	"""
	lf = []
	for f in list_file:
		try:
			if os.path.exists(f):
				lf.append(f)
		except FileExistsError as fee:
			log.warning("FILE NOT EXIST: "+f+"  with error traceback: "+fee)
			return None
		except FileNotFoundError as fnf:
			log.warning("FILE NOT FOUND: "+f+"  with error traceback: "+fnf)
			return None
	return lf

def add_annotation_to_image(finput_image, ltoolnames, list_of_files_with_variants):
	'''
	Adds text to image; Here the Text will be the toolname and the number of variants;
	ltoolnames and list_of_files_with_variants MUST be paired [ (t1, t2, t3) to (lv1, lv2, lv3) ]
	:param image: image to update (i.e. filename of the image)
	:param ltoolnames is the list of toolnames
	:param list_of_files_with_variants: Expecting Bedfiles without Any header or comment lines; just data of interest.
	:return: None
	'''
	lvarfiles = check_if_files_in_list_exist(list_of_files_with_variants)

	if lvarfiles is None or lvarfiles == []:
		log.error("None of the expected png files for annotation were found; Skipping Image file annotation;")
		## we just do not do any annotation;
		return None

	if len(ltoolnames) != len(lvarfiles):
		msg = "ERROR: number of toolnames MUST match the number of Given files that contain the variants"
		log.error(msg)
		raise(msg)
	lanno = []
	try:
		for pair in zip(ltoolnames, lvarfiles):
			tn = pair[0]
			N = sum(1 for i in open(pair[1], 'rb'))
			lanno.append(" :  ".join([tn, str(N)]))
			log.info(" -- ".join([str(x) for x in [tn, N]]))
		log.info(str(lanno))

		from PIL import Image, ImageDraw, ImageFont
		import os
		# create Image object with the input image
		image = Image.open(finput_image)

		# initialise the drawing context with the image object as background
		draw = ImageDraw.Draw(image)

		# create font object with the font file and specify desired size
		font = ImageFont.truetype(get_os_specific_system_font(sys.platform), size=40)
		#font = ImageFont.load_default(size=40)

		# starting position of the message
		(x, y) = (150, 200)
		message = "\n".join(lanno)
		color = 'rgb(0, 0, 0)'  # black color
		# draw the message on the background
		draw.text((x, y), message, fill=color, font=font)
		# save the edited image
		anno_image_name = os.path.splitext(os.path.realpath(finput_image))[0]+ ".anno" + os.path.splitext(os.path.realpath(
			finput_image))[1]
		image.save( anno_image_name )
		## uncomment line below if we decide to keep only the annotated image file
		os.rename(anno_image_name, finput_image)
	except ImportError as ie:
		raise(ie)
	except FileNotFoundError as fnf:
		raise(fnf)
	except Exception as e:
		raise(e)
