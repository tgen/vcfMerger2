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


import argparse
import json
import subprocess
import os
from os import path
import sys
import logging as log
import time

# CAPTURED VARIABLES AUTOMATICALLY
## capturing the current path of the current script
scriptDirectory = os.path.dirname(os.path.realpath(__file__))
## as the project should be installed by the user and not modified by the user, we know where the prep_vcf.sh script is
prep_script_path = os.path.join(os.path.dirname(scriptDirectory), "prep_vcfs", "prep_vcf.sh")
prep_vcf_functions_script_path = os.path.join(os.path.dirname(scriptDirectory), "prep_vcfs", "prep_vcf_functions.sh")
vcfmerger_tool_path = os.path.join(os.path.dirname(scriptDirectory), "merge_vcfs", "vcfMerger.py")
snpsift_filter_script_path = os.path.join(os.path.dirname(scriptDirectory), "prep_vcfs", "utils", "filter_vcf.sh")


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

def quote_str(s):
	'''
	This is an important function to get the empty string with appropriate quotes in order to
	avoid errors with the bash command in the subprocesses functions ; empty string MUST show in the bash command

	:param s: input is a string (empty or not)
	:return: a string (with new quotes around or as is)
	'''
	if s == "":
		return ''.join(["\"\'", s, "\'\""])
	return s

def double_quote_str(s):
	'''
	This is an important function to get the empty string with appropriate quotes in order to avoid
	error with the bash or python commands in the subprocess ; empty string MUST show in the bash command

	:param s: input is a string (empty or not)
	:return: a double-quoted string
	'''
	return ''.join(["\"", s, "\""])

def make_json(data, f):
	"""
	create the json file for further processing ... I know, I know, I could avoid this intermediate step
	but I started this script with json file manually made as input of this script. So read_json function
	got created first along with other functions; so I am recycling :-)
	:param data: dictionary or other object that can be converted to json format
	:param f: filename
	:return: filename f
	"""
	with open(f, "w") as write_file:
		json.dump(data, write_file, indent=4)
	return f

def read_json(f):
	"""
	read the json file with specific format for the vcfMerger

	:param f: json file as input
	:return: dictionary of converted json data into a dictionary

	"""
	with open(f, "r") as read_file:
		data = json.load(read_file)
	log.debug(str(data) + '\n')
	return data


def make_data_for_json(lvcfs,
                       ltoolnames,
                       normal_sname,
                       tumor_sname,
                       ref_genome_fasta,
                       lossy,
                       ltpo=None,
                       lacronyms=None,
                       lprepped_vcf_outfilenames=None,
                       lbams=None,
                       lcontigs=None,
                       filter_string_for_snpsift=None,
                       TH_AR=0.9,
                       do_venn=False):
	# TODO : Check if tool precedence is different from order of toolnames
	# if different, reorder the list;
	# otherwise, currently the order of precedence is the same as the toolnames given list
	data = {}
	for tool_idx in range(len(ltoolnames)):
		data[ltoolnames[tool_idx]] = {}

		if lacronyms is None:
			data[ltoolnames[tool_idx]]['tool_acronym'] = ""
		else:
			data[ltoolnames[tool_idx]]['tool_acronym'] = lacronyms[tool_idx]

		data[ltoolnames[tool_idx]]['dir_work'] = "."

		# if os.path.dirname(lvcfs[tool_idx]) is None or os.path.dirname(lvcfs[tool_idx]) == "":
		# 	data[ltoolnames[tool_idx]]['dir_work'] = "."
		# elif "/" in lvcfs[tool_idx] and (lvcfs[tool_idx].startswith("./") or lvcfs[tool_idx].startswith("/")):
		# 	if os.path.exists(lvcfs[tool_idx]):
		# 		data[ltoolnames[tool_idx]]['dir_work'] = os.path.dirname(lvcfs[tool_idx])
		# 	else:
		# 		sys.exit("ERROR: vcf path NOT FOUND from current directory '{}'\npath must be the basename, relative path or full path to the vcf file".format(lvcfs[tool_idx]))
		# else:
		# 	sys.exit("ERROR: vcf pathname is INVALID '{}'\npathname must be the basename, relative path or full path to the vcf file".format(lvcfs[tool_idx]))
		data[ltoolnames[tool_idx]]['normal'] = normal_sname
		data[ltoolnames[tool_idx]]['tumor'] = tumor_sname
		data[ltoolnames[tool_idx]]['vcf'] = lvcfs[tool_idx] ; # manages wherever is the vcf (relative of full path to the current directory)
		data[ltoolnames[tool_idx]]['ref_genome_fasta'] = ref_genome_fasta

		if lprepped_vcf_outfilenames is not None:
			data[ltoolnames[tool_idx]]['prepped_vcf_outfilename'] = lprepped_vcf_outfilenames[tool_idx]
		else:
			data[ltoolnames[tool_idx]]['prepped_vcf_outfilename'] = ""

		data[ltoolnames[tool_idx]]['vcf_indels'] = ""
		data[ltoolnames[tool_idx]]['vcf_snvs'] = ""

		if lbams is not None:
			data[ltoolnames[tool_idx]]['bam'] = lbams[tool_idx]
		else:
			data[ltoolnames[tool_idx]]['bam'] = ""

		if lcontigs is not None:
			data[ltoolnames[tool_idx]]['contigs_file'] = lcontigs[tool_idx]
		else:
			data[ltoolnames[tool_idx]]['contigs_file'] = ""

		# lossy do not need to be added to the json file because it applies to vcfMerger not prep
		# we keep it here just in case we use the json file as a reminder of what was run
		data[ltoolnames[tool_idx]]['lossy'] = lossy
		data[ltoolnames[tool_idx]]['filter_string_snpsift'] = filter_string_for_snpsift
		data[ltoolnames[tool_idx]]['threshold_AR'] = TH_AR
		data[ltoolnames[tool_idx]]['do_venn'] = do_venn

	return data

def filter_vcf(data, path_jar_snpsift):
	"""
	1) parse the data object
	2) prepare command for running the filter bash script
	3) call and run bash script
	4) update the data object with new vcf file

	:param data: dictionary of converted json data
	:return: updated data object
	"""
	for tool in data.keys():

		vcf = data[tool]['vcf']
		log.info("%" * 10 + "  " + str(tool).upper() + "  " + "%" * 10)
		log.info("Filtering vcf ...")
		log.info("input: \ttool\t==\t{}".format(str(tool)))
		log.info("input: \tvcf\t==\t{}".format(str(vcf)))
		mycmd = ["bash", snpsift_filter_script_path, path_jar_snpsift, str("\""+data[tool]['filter_string_snpsift']+"\""), vcf]
		log.info(str(mycmd))
		log.info(" ".join([x for x in mycmd]))
		log.info(("Running filter stage for vcf: {}".format(vcf)))
		subprocess_cmd(' '.join([str(x) for x in mycmd]))
		new_vcf_name = os.path.basename(os.path.splitext(vcf)[0]+".filt.vcf")
		log.info("Expected new filename for the input vcfs for the next stage is: ".format(str(new_vcf_name)))
		data[tool]['vcf'] = os.path.abspath(new_vcf_name)

	return data


def parse_json_data_and_run_prep_vcf(data, dryrun):
	"""
	1) parse the data from the json file
	2) run prep_vcf program for each tool's vcf

	:param data: dictionary of converted json data
	:return: list of expected prep_vcfs from each tool unless error
	"""

	for tool in data.keys():

		log.info("%" * 10 + "  " + str(tool).upper() + "  " + "%" * 10)
		log.info("recap inputs captured from json file")
		log.info("input: \ttool\t==\t{}".format(str(tool)))

		for var, val in data[tool].items():
			log.info("input:\t{} \t==\t{}".format(var, quote_str(val)))
		print()

		# check if outfilename for prep_vcf is Null or None
		if data[tool]['prepped_vcf_outfilename'] is None or data[tool]['prepped_vcf_outfilename'] == "":
			msg = "prepped_vcf_outfilename has not been defined injson file for tool {} ; Aborting.".format(tool)
			sys.exit(str(msg))
		# check if outfilename will be outputted in current working folder or not; important for vcfMerger2.0 to
		# know where the prep vcfs are located (users can provide full path for prep vcf outfilename otherwise ; no relative path in json )
		if os.path.curdir != data[tool]['dir_work']:
			if not str(data[tool]['prepped_vcf_outfilename']).startswith("/"):
				data[tool]['prepped_vcf_outfilename'] = os.path.join(os.path.abspath(os.path.curdir),
				                                                     os.path.basename(
					                                                     data[tool]['prepped_vcf_outfilename']))

		# make string from all options for the future bash command line
		cmdLine = ' '.join(
			[
				"-d", quote_str(data[tool]['dir_work']),
				'--toolname', quote_str(tool),
				'--normal-sname', quote_str(data[tool]['normal']),
				'--tumor-sname', quote_str(data[tool]['tumor']),
				'--vcf', quote_str(data[tool]['vcf']),
				'-g', quote_str(data[tool]['ref_genome_fasta']),
				'-o', quote_str(data[tool]['prepped_vcf_outfilename']),
				'--vcf-indels', quote_str(data[tool]['vcf_indels']),
				'--vcf-snvs', quote_str(data[tool]['vcf_snvs']),
				'--bam', quote_str(data[tool]['bam']),
				'--contigs-file', quote_str(data[tool]['contigs_file']),

			]
		)

		## capture if user wants to make the Venn/upset plots later on before merging vcfs
		if data[tool]['do_venn']:
			cmdLine = ' '.join([cmdLine, "--make-bed-for-venn"])

		# capture threshold AR found in json
		TH_AR = data[tool]['threshold_AR']
		if TH_AR is not None and TH_AR != "" and TH_AR != 0.9:
			cmdLine = ' '.join([cmdLine, "--threshold-AR", TH_AR])

		# display the command line for log purposes
		log.info(prep_script_path + " " + cmdLine)
		my_command = ' '.join(["bash", prep_script_path, cmdLine])
		logFilename = "log_prep_vcf_{}.logs".format(tool)
		subp_logfile = open(logFilename, "w")
		if not dryrun:
			log.info("")
			log.info("%" * 10 + " prep {} vcf ".format(tool).upper() + "%" * 10)
			log.info("logging prep steps to file: " + str(logFilename))
			process = subprocess.Popen(my_command,
			                           shell=True,
			                           universal_newlines=True,
			                           stdout=subp_logfile,
			                           stderr=subp_logfile)
			process.wait()
			print(str(process.returncode))
			subp_logfile.close()
			if process.returncode is not 0:
				sys.exit("{} FAILED for tool {} ".format(prep_script_path, tool))



def subprocess_cmd(command):
	os.system(command)

def prepare_bed_for_venn(vcf):
	'''
	if no beds have been provided to vcfMerge2.py using --beds option and --do-venn has been enabled, and ...
	the skip-prep-vcf is used, we can make the beds from the vcf file(s) provided. As we already have the
	function << prepare_input_file_for_Venn >> in the bash script named << prep_vcf.sh >>, we will source the function
	and run it using system.command()
	:param vcf:
	:return: none
	'''

	# Build subprocess command
	mycmd = ["source ", prep_vcf_functions_script_path, " && ", "prepare_input_file_for_Venn ", vcf]
	log.info(str(mycmd))
	log.info(" ".join([x for x in mycmd]))
	log.info(("Running bash function prep_input_file_for_venn command"))
	subprocess_cmd(''.join([ str(x) for x in mycmd]))


def merging_prepped_vcfs(data, merged_vcf_outfilename, delim, lossy, dryrun, do_venn, lbeds, skip_prep_vcfs, ):
	"""

	:param data:
	:param merged_vcf_outfilename:
	:return: None
	"""

	list_vcfs = ""
	list_tools = ""
	list_tools_acronyms = ""
	for tool in data.keys():
		vcf = data[tool]['vcf']
		prepped_vcf = data[tool]['prepped_vcf_outfilename']
		acronym = data[tool]['tool_acronym'] if data[tool]['tool_acronym'] != "" else str(tool).upper()
		list_tools = delim.join([list_tools, tool]) if list_tools != "" else tool
		vcf_to_add = prepped_vcf if prepped_vcf != "" else vcf
		list_vcfs = delim.join([list_vcfs, vcf_to_add]) if list_vcfs != "" else vcf_to_add
		list_tools_acronyms = delim.join([list_tools_acronyms, acronym]) if list_tools_acronyms != "" else acronym

	log.info(str(list_tools_acronyms))
	my_command = ' '.join(["python", vcfmerger_tool_path,
	                       "--toolnames", double_quote_str(list_tools),
	                       "--vcfs", double_quote_str(list_vcfs),
	                       "-o", merged_vcf_outfilename,
	                       "-a", double_quote_str(list_tools_acronyms)
	                       ])

	if lossy:
		my_command = my_command + " --lossy"

	if do_venn:
		for tool in data.keys():
			if not skip_prep_vcfs:
				list_beds = delim.join([str(os.path.splitext(vcf)[0]+".intervene.bed") for vcf in list_vcfs.split(delim)]) ## extension intervene.bed defined in prep_vcf.sh
			elif lbeds == "":
				log.info("processing vcf to make bed for intervene tool: "+str(tool))
				log.info("trying to create on the fly the bed file using function in prep_vcf.sh script")
				print(tool + " __ prepare_bed_for_venn __  " + data[tool]['vcf'])
				prepare_bed_for_venn(data[tool]['vcf'])
				list_beds = delim.join([str(os.path.splitext(vcf)[0] + ".intervene.bed") for vcf in list_vcfs.split(delim)])
			else:
				list_beds = lbeds


		if len(list_beds) == 0:
			sys.exit("ERROR: --do-venn provided, but list_beds file EMPTY ; check you inputs. Aborting.")
		log.info(double_quote_str(list_beds))
		if len(list_beds.split(delim)) != len(list_tools.split(delim)):
			sys.exit("ERROR: --do-venn provided, but number of bed files ({}) DO NOT matched number of tools ({})  ; check you input or contact your IT/HPC admin".format(list_beds,list_tools))
		my_command = my_command + " --do-venn --beds " + double_quote_str(list_beds)

	log.info(double_quote_str(list_tools))
	log.info(double_quote_str(list_vcfs))
	log.info("merging vcfs ...")
	log.info(my_command)

	if not dryrun:
		log.info("")
		log.info("%" * 10 + " starting vcfMerger ... ".upper() + "%" * 10)
		process = subprocess.Popen(my_command, shell=True)
		process.wait()
		log.info("vcfMerger2.0 exit value: " + str(process.returncode))
		if process.returncode is not 0:
			sys.exit("{} FAILED with vcfs files: {} ".format(prep_script_path, list_vcfs))

def check_inputs(lvcfs, ltoolnames, ltpo=None, lacronyms=None, lprepped_vcf_outfilenames=None, lbeds=None ):
	"""

	:param lvcfs:
	:param ltoolnames:
	:param ltpo:
	:param acronyms:
	:param delim:
	:param lprepped_vcf_outfilenames:
	:return:
	"""
	if lvcfs is None:
		log.info("ERROR: Found list of input vcfs to be prepped or to be merged empty")
		sys.exit("ERROR: list of vcfs empty")
	if len(lvcfs) == 1:
		log.info("ERROR: Found list of input vcfs with ONE vcfs only; Minimumn number of vcfs must be TWO;")
		sys.exit("ERROR: list of vcfs length of 1 vcf only")
	if len(lvcfs) != len(ltoolnames):
		log.info("ERROR: Found {} vcfs provided and {} tools given ".format(len(lvcfs), len(ltoolnames)))
		sys.exit(
			"ERROR: list vcfs files MUST be equal to the number of tools given ; check if delimiter is adequate and do not interfere")
	if ltpo is not None and len(ltoolnames) != len(ltpo):
		log.info("ERROR: Found {} tools in precedence list and {} toolnames given".format(len(ltpo), len(ltoolnames)))
		sys.exit(
			"ERROR: Number of toolnames MUST be equal to the number of tools given in the list of precedence ; "
			"check if delimiter is adequate and do not interfere with splitting the given lists of tools")
	if lacronyms is not None and len(ltoolnames) != len(lacronyms):
		log.info("ERROR: Found {} in acronyms list and {} toolnames given".format(len(lacronyms), len(ltoolnames)))
		sys.exit(
			"ERROR: Number of toolnames MUST be equal to the number of acronyms given ;\n"
			"ERROR: check if delimiter is adequate and do not interfere with splitting the given lists of tools")
	if lprepped_vcf_outfilenames is not None and len(ltoolnames) != len(lprepped_vcf_outfilenames):
		log.info(
			"ERROR: Found {} in list of given prep-filenames and {} toolnames given".format(len(lprepped_vcf_outfilenames),
			                                                                         len(ltoolnames)))
		sys.exit(
			"ERROR: Number of toolnames MUST be equal to the number of intermediate prep-outfilenames given ;\ncheck if delimiter is adequate and do not interfere with splitting the given lists of tools")

def main(args, cmdline):

	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)

	delim = "|"
	if args["delim"]:
		delim = args["delim"]
		if len(delim) != 1:
			exit("the list delimiter given with --delim should be one character only and NOT a space; found {}".format(str(len(delim))))
		log.info("delimiter is:\t" + delim)

	lvcfs = []
	if args["vcfs"]:
		lvcfs = str(args["vcfs"]).split(delim)
		log.info("ordered list of vcfs given:\t\t{}".format(str(lvcfs)))

	ltoolnames = []
	if args["toolnames"]:
		ltoolnames = str(args["toolnames"]).split(delim)
		ltoolnames = [x.lower() for x in ltoolnames]
		log.info("ordered list of toolnames given:\t{}".format(str(ltoolnames)))

	ref_genome_fasta = None
	if args["refgenome"]:
		ref_genome_fasta = args["refgenome"]
		if not os.path.exists(ref_genome_fasta):
			sys.exit("reference genome {} NOT FOUND; Aborting.".format(ref_genome_fasta))

	list_tool_precedence_order = None
	if args["precedence"]:
		list_tool_precedence_order = str(args["precedence"]).split(delim)
		list_tool_precedence_order = [x.upper() for x in list_tool_precedence_order]
		log.info("given precedence: {} ".format(str(list_tool_precedence_order)))

	lacronyms = None
	if args["toolacronyms"]:
		lacronyms = str(args["toolacronyms"]).split(delim)
		lacronyms = [x.upper() for x in lacronyms]
		log.info("tool acronyms given:\t\t\t{}".format(lacronyms))

	lossy = False
	if args["lossy"]:
		log.info("lossy enabled")
		lossy = True

	normal_sname = None
	if args['normal_sname']:
		normal_sname = args['normal_sname']

	tumor_sname = None
	if args['tumor_sname']:
		tumor_sname = args['tumor_sname']

	lprepped_vcf_outfilenames = None
	if args["prep_outfilenames"]:
		lprepped_vcf_outfilenames = str(args["prep_outfilenames"]).split(delim)
		log.info("tool-specific filenames for intermediate upto vcfMerger2-specs: " + str(lprepped_vcf_outfilenames))
		if args["skip_prep_vcfs"] == True:
			lprepped_vcf_outfilenames = None

	lbams = None
	if args["bams"]:
		lbams = str(args["bams"]).split(delim)
		log.info("ordered list of bams given:\t{}".format(str(lbams)))

	lcontigs = None
	if args["contigs_file_for_vcf_header"]:
		lcontigs = str(args["contigs_file_for_vcf_header"]).split(delim)
		log.info("ordered list of contigs files given:\t{}".format(str(lcontigs)))

	merged_vcf_outfilename = None
	if args["merged_vcf_outfilename"]:
		merged_vcf_outfilename = str(args["merged_vcf_outfilename"])
		log.info("filename for the merged output vcf will be: " + merged_vcf_outfilename)

	skip_prep_vcfs = False
	if args["skip_prep_vcfs"]:
		skip_prep_vcfs = args["skip_prep_vcfs"]
		log.info("skip_prep_vcfs:" + str(skip_prep_vcfs))

	filter_string_for_snpsift = None
	if args["filter"] is not None:
		filter_string_for_snpsift = args["filter"]
		log.info("filter string to be used with snpSift: \"" + str(filter_string_for_snpsift) +"\"")

	path_jar_snpsift = None
	if args["filter"] is not None and args['path_jar_snpsift'] is not None:
		path_jar_snpsift = args['path_jar_snpsift']
		log.info("Path to provided snpSift.jar file:" + str(filter_string_for_snpsift))
		if not os.path.exists(path_jar_snpsift):
			raise Exception("ERROR: snpSift.jar FILE NOT FOUND. Aborting!")
	elif args["filter"] is not None and args['path_jar_snpsift'] is None:
		raise Exception("Please provide the Full PATH to a snpSift.jar file using the option --path-jar-snpsift. Aborting!")
	else:
		log.info("Well, you provided the path to snpSift for nothing as you have not set the filter option. :-) ")


	skip_merge = False
	if args["skip_merge"]:
		skip_merge = args["skip_merge"]
		log.info("skip_merge:" + str(skip_merge))

	TH_AR = 0.90
	if args['threshold_AR']:
		TH_AR = args['threshold_AR']
		if isinstance(TH_AR, (int, float, complex)):
			raise Exception("Threshold-AR must be a float or integer value between 0 and 1 (range ]0,1]). Check your inputs.")
		log.info("user given threshold for AR: " + str(TH_AR))

	lbeds = ""
	if args["beds"]:
		lbeds = str(args["beds"]).split(delim)
		log.info("ordered list of beds given:\t\t{}".format(str(lbeds)))

	do_venn = False
	if args["do_venn"]:
		do_venn = True
		log.info("make venn enabled")

	dryrun = False
	if args["dry_run"]:
		dryrun = args["dry_run"]
		log.info("dry-run?:" + str(dryrun))

	# if inFileJson is None:
	##@@@@@@@@@
	## MAIN  ##
	##@@@@@@@@@
	check_inputs(lvcfs, ltoolnames, ltpo=list_tool_precedence_order, lacronyms=lacronyms,
	             lprepped_vcf_outfilenames=lprepped_vcf_outfilenames, lbeds=lbeds)

	data = make_data_for_json(lvcfs,
	                          ltoolnames,
	                          normal_sname,
	                          tumor_sname,
	                          ref_genome_fasta,
	                          lossy,
	                          ltpo=list_tool_precedence_order,
	                          lacronyms=lacronyms,
	                          lprepped_vcf_outfilenames=lprepped_vcf_outfilenames,
	                          lbams=lbams,
	                          lcontigs=lcontigs,
	                          filter_string_for_snpsift=filter_string_for_snpsift,
	                          TH_AR=TH_AR,
	                          do_venn=do_venn)
	json_filename = "vcfMerger.json"
	make_json(data, json_filename)
	#inFileJson = make_json(data, json_filename)
	#data = read_json(inFileJson) ## uncomment for debugging if necessary ; data is already created above

	# if filter_string_for_snpsift is not None:
	# 	data = filter_vcf(data, path_jar_snpsift)
	# 	print(str(data))

	if not skip_prep_vcfs:
		log.info("**** prep vcf steps section ***".upper())
		parse_json_data_and_run_prep_vcf(data, dryrun)
		log.info("**** merge process section  ****".upper())
	else:
		log.info("**** SKIPPED prep vcfs step SKIPPED ****")

	## FILTERING STEP for PREPPED VCFS (note: Filtering must not be applied ot un-prepped vcfs unless users would like to
	if filter_string_for_snpsift is not None:
		log.info("Performing vcf filtering for ALL the provided prepped vcfs ...")
		data = filter_vcf(data, path_jar_snpsift)
		log.info(str(data))


	if not skip_merge:
		merging_prepped_vcfs(data, merged_vcf_outfilename, delim, lossy, dryrun, do_venn, lbeds, skip_prep_vcfs)
	else:
		log.info("**** SKIPPED merge step SKIPPED ****")

	if not skip_prep_vcfs and not skip_merge:
		log.info("prep and merge vcfs Elapsed time in seconds:  {}".format(str(int(round((time.time() - start_time))))))
		log.info("prep and merge vcfs completed successfully")
	elif skip_prep_vcfs:
		log.info("merge vcfs Elapsed time in seconds:  {}".format(str(int(round((time.time() - start_time))))))
		log.info("merge vcfs completed successfully")
	elif skip_merge:
		log.info("prep Elapsed time in seconds:  {}".format(str(int(round((time.time() - start_time))))))
		log.info("prep vcfs completed successfully")

def make_parser_args():
	parser = argparse.ArgumentParser(description='Processing options ...')
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')
	isRequired = True


	required.add_argument('--vcfs',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='List of vcfs file delimited by DELIM character; default DELIM is "|" (piped character) ; if needed re-assign DELIM '
	                           'using --delim option')
	required.add_argument('--toolnames',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='List of vcfs file delimited by DELIM character; default DELIM is pipe unless --delim '
	                           'option is used using --delim option')
	required.add_argument('-g','--refgenome',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='reference genome used with bcftools norm ; must match reference used for alignment')

	required.add_argument('--normal-sname',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='expected name of normal sample in vcf file')

	required.add_argument('--tumor-sname',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='expected name of tumor sample in vcf file')

	optional.add_argument('--prep-outfilenames',
	                      required=False,
	                      action=UniqueStore,
	                      help='delim-separated names to the tool-specific prepared vcf files')
	optional.add_argument('-c', '--precedence',
	                      required=False,
	                      action=UniqueStore,
	                      help=' sorted delim-separated list of the toolnames as listed in --toolnames ; '
	                           'This list stipulates an order of precedence for the tools different from the '
	                           'default order given by the --toolnames list')

	optional.add_argument('--bams',
	                      required=False,
	                      action=UniqueStore,
	                      help='List of bams necessary for capturing contigs if not present in input vcf; otherwise put empty_string as value for each tool ')

	optional.add_argument('--contigs-file-for-vcf-header',
	                      required=False,
	                      action=UniqueStore,
	                      help='List of contigs necessary for capturing adding them to tool vcf header if needed; otherwise put empty_string as value for each tool  ;do not provide if bam file is given instead ')

	optional.add_argument('-a', '--toolacronyms',
	                      required=False,
	                      action=UniqueStore,
	                      help='List of Acronyms for toolnames to be used as PREFIXES in INFO field ; same DELIM as --vcfs ')

	required.add_argument('-o', '--merged-vcf-outfilename',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='outfilename for the merge vcf (can be relative or full path)')

	optional.add_argument('--delim',
	                      required=False,
	                      action=UniqueStore,
	                      help='delimiter which will be use to create the arguments value for the vcfMerger2.0 tool ; default is "|" (a.k.a pipe character)')

	optional.add_argument('--threshold-AR',
	                      required=False,
	                      action=UniqueStore,
	                      help='AlleRatio threshold value to assign genotype; 0/1 if less than threshold, 1/1 if equal or above threshold; default is 0.90 ; range ]0,1] ')

	optional.add_argument('--lossy',
	                      required=False,
	                      help='This will create a lossy merged vcf by only keeping the information from the tool with first hand precedence',
	                      action='store_true')

	optional.add_argument('--skip-prep-vcfs',
	                      required=False,
	                      action='store_true',
	                      help=' skip the step for preparing vcfs up to specs and only run the merge step; implies all prep-vcfs are ready already ; same options and inputs required as if prep step was run ')

	optional.add_argument('--filter',
	                      required=False,
	                      action=UniqueStore,
	                      help='enable vcf filtering using snpSift tool; A string argument is passed to this filter option; \
	                      This string MUST be formatted as if you were using it for snpSift, i.e. we used "as_is" the string provided; \
	                      Therefore a valid format is mandatory; Check snpSift manual ; IMPORTANT NOTE: the filtering MUST use FLAG and TAGs \
	                      that are COMMON to ALL the tools involved ; The Prepped vcfs have some common flags in the GENOTYPE fields (so far GT, DP, AR, AD), \
	                      and the CC flag is common to ALL records. The filtering takes place after the prep_vcf step and before merging the vcfs; Example of String: \
	                      \"( GEN[TUMOR].AR >= 0.10 ) & ( GEN[NORMAL].AR <= 0.02 ) & ( CC >= 2 ) & ( GEN[TUMOR].DP >= 10 & GEN[NORMAL].DP>=10 )\", where NORMAL or TUMOR can be replaced with appropriate indices or other given names')

	optional.add_argument('--path-jar-snpsift',
	                      required=False,
	                      action=UniqueStore,
	                      help='Provide full Path of the snpSift.jar you want to use for filtering the vcf before merging them')

	optional.add_argument('--skip-merge',
	                      required=False,
	                      action='store_true',
	                      help='enabling this flag prevents doing the merging step [useful if only the prep step needs to be done ]')

	optional.add_argument('--beds',
	                      help='list of bed files to be used to make Venns or Upset plots; requires to enable --do-venn as well to validate making Venn/upset plots ; list MUST be delimited by DELIM character (--delim or default delim)',
	                      action=UniqueStore)
	optional.add_argument('--do-venn',
	                      help='using the bed files listed in --beds option, Venns or Upset plot will be created ; need to match the number of tools listed in --toolnames ',
	                      action='store_true')

	optional.add_argument('-n', '--dry-run',
	                      required=False,
	                      action='store_true',
	                      help=' perform a dry run to just see the commands lines that will be run behind the scenes ')

	print(str(parser.prog) + "   " + str(parser.description))
	return parser


if __name__ == '__main__':
	# capture time for calculating vcfMerger's runtime
	start_time = time.time()
	# capture commandline for adding it to vcf header
	cmdline = ' '.join(sys.argv)
	parser = make_parser_args()
	args = vars(parser.parse_args())  # vars() function returns a dictionary of key-value arguments
	print(str(args))
	main(args, cmdline)
	sys.exit()
