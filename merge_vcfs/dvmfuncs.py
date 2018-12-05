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


import re, sys
from myGenotype import Genotype
import logging as log

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
			log.info("Samples should be in the EXACT same order and exact same CASE (case sensitive comparison)")
			log.info("Correct and modify the inputs VCF files to get them up to specs for the current merging tool")
			sys.exit(2)
	log.info("Comparison of << " + msg + " >> between VCFs: PASSED")

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

def create_new_header_for_merged_vcf(tuple_objs, command_line, vcfMerger_Format_Fields_Specific, vcfMerger_Info_Fields_Specific ):
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
	lh = []
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
	l_contigs = set(l_contigs)
	## adding the contigs to the list of strings called "lh"
	for contig in sorted(l_contigs):
		lh.append(contig)
	## prefixing the header with the toolname, the same way the INFO Fields Flag are prefixed
	reference=""
	for vtdo in tuple_objs:  ## list of vcfToDict objects
		if reference == "":
			indices = [i for i, s in enumerate(vtdo.header_other_info) if '##reference=' in s]
			reference = vtdo.header_other_info[indices[0]]
		for s in vtdo.header_filters:
			lh.append(prefix_headers_information_line_with_toolname(s, vtdo.toolname))
		for s in vtdo.header_info:
			lh.append(prefix_headers_information_line_with_toolname(s, vtdo.toolname))
		for s in vtdo.header_format:
			lh.append(prefix_headers_information_line_with_toolname(s, vtdo.toolname))
		for s in vtdo.header_other_info:
			lh.append(prefix_headers_other_information_line_with_toolname(s, vtdo.toolname))
		## if LOSSLESS, the column QUAL, FILTER, ID, and some others are ADDED to the variant record
		## this creates NEW fields prefixed with the toolname
		for COLUMN in ["FILTER", "QUAL", "ID"]:
			## ##INFO=<ID=SEURAT_AR1,Number=1,Type=Float,Description="Allele frequency of ALT allele in normal">
			stringline = ''.join(["##INFO=<ID=", vtdo.toolname, "_", COLUMN,
			                      ',Number=.,Type=String,Description=',
			                      '"Represents lossless data from tool ',vtdo.toolname,
			                      'for column ', COLUMN,'">'] )
			lh.append(stringline)
		## Here when LOSSLESS is enabled, fields that were in format of the secondary tools, are added to
		## the INFO field with the following format: TOOLNAME_Sx_FIELDofINTEREST
		## where x represents an indice of the Sample starting at 1 up to n.
		## if dealing with TUMOR_NORMAL, we should have S1 and S2 (respectively from column 10 and 11 in vcf)
		## Now we need to implement this here TODO: add the test if lossless enabled;
		##1) we capture the Format column, aka column number 9 for the current tool and prefixed it with tool names
		## and Sample number
		toolname = vtdo.toolname
		numberOfSamples = len(vtdo.samplenames)
		for S in vtdo.header_format:
			## return the first indice where the pattern is in the string
			idx1 = S.find(',')
			idx2 = S[:idx1].rfind("=")
			FIELD = (S[idx2+1:idx1])
			for i in range(1,numberOfSamples+1):
				newField = '_'.join([ toolname, "S"+str(i), FIELD ])
				# print(newField)
				stringline = ''.join(["##INFO=<ID=", newField,
			                      ',Number=.,Type=String,Description=',
			                      '"lossless data from defined tool">'])
				lh.append(stringline)


	for item in vcfMerger_Format_Fields_Specific:
		lh.append(item)
	for item in vcfMerger_Info_Fields_Specific:
		lh.append(item)
	lh.append(reference)

	lh.append(command_line)
	return lh  ## returns a list

def add_new_header_to_merged_file(f, list_lines_header,header_chrom_line):
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
		# log.info("%" * 10 + " IN ELSE " + "%" * 10)
		# log.info(v)
		# log.info(v.format('GT'))
		# log.info(v.FORMAT)
		# log.info(v.genotypes)
		gts = v.genotypes
		log.debug("GENOTYPE gts is:" + str(gts))
		my_genotypes = [Genotype(li) for li in gts]
		# log.info("TYPE of my_genotypes :::: " + str(type(my_genotypes)))
		# log.info("MY GENOTYPES ARE: ")
		# log.info(my_genotypes)
		for i in range(totnum_samples):
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

def rebuiltVariantLine(LV, lossless, list_Fields_Format, totnum_samples = 2 ):

	"""
	built a lossless variant line using all the variant INFO from each tool.
	:param LV: variant information object for the Winner Tool; v represents an entire line of a vcf
	:param: lossless, boolean; if True, will keep ALL the information from teh secondary tools and add that
	information to the INFO field in the  merged VCF
	:param  As this versio is dedicated to Somatic, we expect two Sample by default in the FORMAT Columns.
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
	newINFO = renameINFO_IDS(wtv, wtn)  ## LV[0][0] is the toolname which has precedence for that current call

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
		# position ; If we have only one tool that called that variant, we never enter this for loop;
		tn = LV[i][0]  # In the list of 2 elements, the first is the toolname to remember the file of origin
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
	TOOLS_AS_FLAG = ";".join([LV[i][0] for i in range(1, len(LV))])


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
	FILT = ';'.join([ wtn+"_"+FLTV for FLTV in wtv.FILTER.split(';') ]) if wtv.FILTER != None else "."
	#FILT = wtv.FILTER if wtv.FILTER != None else "."  ## note: cyvcf2 consider the keyword PASS as no filter and
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

def make_venn(ltoolnames, lbeds, delim, saveOverlapsBool=False, upsetBool=False):
	#TODO we could check if any of the tools or any of the vcfs filenames already contains a comma; if so raise error
	names = "\""+','.join([name for name in ltoolnames])+"\""
	beds = ' '.join([name for name in lbeds])
	numberOfTools = len(ltoolnames)
	type = "\"genomic\""
	colors = get_colors_for_venns(numberOfTools)
	title = "\"Venn using " + str(numberOfTools) + " variant callers\""
	figtype = "png"
	dpi = 300
	bordercolors = ','.join(["\"black\""]*numberOfTools)
	fontsize = 20
	project = "vcfMerger2_"+str(numberOfTools)+"_tools."+str(figtype) ;  ## this is actually the name of the png image file while the output_name is the folder where the intervene results are going into
	output_name = "upsetPlot_"+str(numberOfTools)+"_tools" if upsetBool else "venn_"+str(numberOfTools)+"_tools"


	# Define command and arguments
	command = 'intervene'
	# Define the type of venn
	if numberOfTools >= 5:
		upsetBool = True
	if upsetBool:
		vtype = "upset"
		additional_args = " ".join(["--ninter", "150",
		                            "--sbcolor", "\"#d8b365\"",
		                            "--mbcolor", "\"#5ab4ac\"",
		                            "--showzero",
		                            "--showsize",
		                            "--order", "freq"
		                            ])
	else:
		vtype = "venn"
		saveOverlaps = "--save-overlaps" if saveOverlapsBool else ""
		additional_args = ' '.join(["--bordercolors", bordercolors,
		                            "--colors", ','.join([ "\""+color+"\"" for color in colors ]),
		                            saveOverlaps
		                            ])

	log.info(str(additional_args))
	log.info(str(beds))
	print(str(type(beds)))

	args = ' '.join(["--input", beds,
	                 "--type", type,
	                 "--names", names,
	                 "--title", title,
	                 "--figtype", figtype,
	                 "--dpi", str(dpi),
	                 "--fontsize", str(fontsize),
	                 "--project", project,
	                 "--output", output_name
	                 ])

	# Build subprocess command
	cmd = [command, vtype, args, additional_args]
	cmd = [command, vtype, args]
	#cmd = [ "intervene", "venn", "--input", "SLK.prepped.intervene.bed", "MUT.prepped.intervene.bed", "--type", "genomic", "--names", "STRELKA2,MUTECT2", "--title", "\"Venn using 2 variant callers\"", "--dpi", "300" ] ## THAT LINE WORKED
	log.info(str(cmd))
	log.info(" ".join([x for x in cmd]))
	# check_output will run the command and store to result
	import subprocess
	print("*"*50)
	print("full command run intervene")
#	subprocess.call(cmd, shell=True, universal_newlines=True)
	process = subprocess.run(cmd, shell=False, universal_newlines=False)
	#process.wait()
	print(str(process.returncode))
	if process.returncode is not 0:
		sys.exit("Venn Creation FAILED")


