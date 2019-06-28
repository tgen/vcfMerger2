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


import argparse
import gzip
import json
import logging as log
import os
import re
import shutil
import subprocess
import sys
import time

# CAPTURED VARIABLES AUTOMATICALLY
## capturing the current path of the current script
scriptDirectory = os.path.dirname(os.path.realpath(__file__))
## as the project should be installed by the user and not modified by the user, we know where the prep_vcf.sh script is
prep_vcf_script_path = os.path.join(os.path.dirname(scriptDirectory), "prep_vcfs", "prep_vcf.sh")
prep_vcf_functions_script_path = os.path.join(os.path.dirname(scriptDirectory), "prep_vcfs", "prep_vcf_functions.sh")
prep_germline_vcf_script_path = os.path.join(os.path.dirname(scriptDirectory), "prep_vcfs_germline",
                                             "prep_vcf_germline.sh")
prep_germline_vcf_functions_script_path = os.path.join(os.path.dirname(scriptDirectory), "prep_vcfs_germline",
                                                       "prep_vcf_functions_germline.sh")
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


def is_gzip(path, magic_number=b'\x1f\x8b'):
	"""Returns True if the path is gzipped. Fucntion taken from JetStream and created by Ryan Richolt"""
	if os.path.exists(path) and not os.path.isfile(path):
		err = 'This should only be used with regular files because otherwise ' \
		      'it will lose some data.'
		raise ValueError(err)

	with open(path, 'rb') as fp:
		if fp.read(2) == magic_number:
			return True
		else:
			return False


def check_if_vcf_is_compressed(lvcfs, dirout):
	'''check if vcfs are comporessed and if so, uncompress the vcf in current working directory
	:param lvcfs
	:return updated lvcfs
	'''

	for i in range(len(lvcfs)):
		vcf = lvcfs[i]
		if is_gzip(vcf):
			uvcf = os.path.join(dirout, os.path.basename(os.path.splitext(vcf)[0]))
			log.info("vcf file after decompression: " + uvcf)
			with gzip.open(vcf, 'r') as f_in, open(uvcf, 'wb') as f_out:
				shutil.copyfileobj(f_in, f_out)
			lvcfs[i] = uvcf
	return lvcfs


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


def make_data_for_json(lvcfs, ltoolnames, normal_sname, tumor_sname,
                       ref_genome_fasta, lossy, germline_snames=None,
                       ltpo=None, lacronyms=None, lprepped_vcf_outfilenames=None,
                       lbams=None, lcontigs=None, filter_string_for_snpsift=None,
                       TH_AR=0.9, do_venn=False, venn_title="", skip_prep_vcfs=False,
                       dirout=None, delete_temps=False):
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
		data[ltoolnames[tool_idx]]['tool_precedence_order'] = "" if ltpo is None else ltpo

		if dirout is not None:
			data[ltoolnames[tool_idx]]['dir_work'] = dirout
		else:
			dirout = os.path.realpath(os.curdir)

		# if os.path.dirname(lvcfs[tool_idx]) is None or os.path.dirname(lvcfs[tool_idx]) == "":
		# 	data[ltoolnames[tool_idx]]['dir_work'] = "."
		# elif "/" in lvcfs[tool_idx] and (lvcfs[tool_idx].startswith("./") or lvcfs[tool_idx].startswith("/")):
		# 	if os.path.exists(lvcfs[tool_idx]):
		# 		data[ltoolnames[tool_idx]]['dir_work'] = os.path.dirname(lvcfs[tool_idx])
		# 	else:
		# 		sys.exit("ERROR: vcf path NOT FOUND from current directory '{}'\npath must be the basename, relative path or full path to the vcf file".format(lvcfs[tool_idx]))
		# else:
		# 	sys.exit("ERROR: vcf pathname is INVALID '{}'\npathname must be the basename, relative path or full path to the vcf file".format(lvcfs[tool_idx]))

		if germline_snames is not None:
			data[ltoolnames[tool_idx]]['germline_snames'] = germline_snames
		else:
			data[ltoolnames[tool_idx]]['germline_snames'] = ""
		data[ltoolnames[tool_idx]]['normal'] = normal_sname
		data[ltoolnames[tool_idx]]['tumor'] = tumor_sname

		data[ltoolnames[tool_idx]]['vcf'] = lvcfs[
			tool_idx];  # manages wherever is the vcf (relative of full path to the current directory)
		data[ltoolnames[tool_idx]]['ref_genome_fasta'] = ref_genome_fasta
		data[ltoolnames[tool_idx]]['ref_genome_fasta_dict'] = ref_genome_fasta_dict

		if lprepped_vcf_outfilenames is not None:
			data[ltoolnames[tool_idx]]['prepped_vcf_outfilename'] =  os.path.sep.join([str(dirout),  str(lprepped_vcf_outfilenames[tool_idx]) ])
		elif skip_prep_vcfs:
			data[ltoolnames[tool_idx]]['prepped_vcf_outfilename'] = lvcfs[tool_idx]
		else:
			raise("ERROR: No Prep-vcfs assign")

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
		data[ltoolnames[tool_idx]]['lossy'] = lossy;

		if filter_string_for_snpsift is not None:
			if len(re.findall("###", filter_string_for_snpsift)) == 0:  ## HARDCODED DELIMITER
				data[ltoolnames[tool_idx]]['filter_string_snpsift'] = filter_string_for_snpsift
			else:
				data[ltoolnames[tool_idx]]['filter_string_snpsift'] = filter_string_for_snpsift.split("###")[
					tool_idx];  ## HARDCODED DELIMITER
		else:
			data[ltoolnames[tool_idx]]['filter_string_snpsift'] = None

		data[ltoolnames[tool_idx]]['threshold_AR'] = TH_AR
		data[ltoolnames[tool_idx]]['do_venn'] = do_venn

		data[ltoolnames[tool_idx]]['venn_title'] = venn_title
		data[ltoolnames[tool_idx]]['delete_temps'] = delete_temps
	return data

def filter_unprepped_vcf(data, path_jar_snpsift):
	"""
	1) parse the data object
	2) prepare command for running the filter bash script
	3) call and run bash script
	4) update the data object with new vcf file

	:param data: dictionary of converted json data
	:return: updated data object
	"""
	CONSTANT_STRING_FOR_PASS_RECORDS = "( FILTER == 'PASS' )"  ## HARDCODED
	log.info("Filtering vcf ... " + CONSTANT_STRING_FOR_PASS_RECORDS)

	for tool in data.keys():
		dir_temp = data[tool]['dir_work']
		vcf = data[tool]['vcf']
		log.info("%" * 10 + "  " + str(tool).upper() + "  " + "%" * 10)

		log.info("input: \ttool\t==\t{}".format(str(tool)))
		log.info("input: \tvcf\t==\t{}".format(str(vcf)))
		mycmd = ["bash", snpsift_filter_script_path, path_jar_snpsift, "pass",
		         str("\"" + CONSTANT_STRING_FOR_PASS_RECORDS + "\""), dir_temp, vcf]
		log.info(str(mycmd))
		log.info(" ".join([x for x in mycmd]))
		log.info(("Running filter stage for vcf: {}".format(vcf)))
		subprocess_cmd(' '.join([str(x) for x in mycmd]))
		new_vcf_name = os.path.join(dir_temp, os.path.basename(os.path.splitext(vcf)[0] + ".pass.vcf"))
		log.info("Expected new filename for the input vcfs for the next stage is: ".format(str(new_vcf_name)))
		data[tool]['vcf'] = new_vcf_name
	return data

def filter_prepped_vcf(data, path_jar_snpsift):
	"""
	1) parse the data object
	2) prepare command for running the filter bash script
	3) call and run bash script
	4) update the data object with new vcf file

	:param data: dictionary of converted json data
	:return: updated data object
	"""
	for tool in data.keys():

		dir_temp = data[tool]['dir_work']
		vcf = data[tool]['prepped_vcf_outfilename']
		log.info("%" * 10 + "  " + str(tool).upper() + "  " + "%" * 10)
		log.info("Filtering vcf ...")
		log.info("input: \ttool\t==\t{}".format(str(tool)))
		log.info("input: \tvcf\t==\t{}".format(str(vcf)))
		mycmd = ["bash", snpsift_filter_script_path, path_jar_snpsift, "filt",
		         str("\"" + data[tool]['filter_string_snpsift'] + "\""), dir_temp, vcf]
		log.info(str(mycmd))
		log.info(" ".join([x for x in mycmd]))
		log.info(("Running filter stage for vcf: {}".format(vcf)))
		subprocess_cmd(' '.join([str(x) for x in mycmd]))
		new_vcf_name = os.path.join(dir_temp, os.path.basename(os.path.splitext(vcf)[0] + ".filt.vcf") )
		log.info("Expected new filename for the input vcfs for the next stage is: ".format(str(new_vcf_name)))
		data[tool]['prepped_vcf_outfilename'] = new_vcf_name
		data[tool]['vcf'] = new_vcf_name ## we consider that the input vcf is now the filtered vcf; which is also the prepped vcf ## TRICK here

	return data

def parse_json_data_and_run_prep_vcf_germline_parallel(tool, data, dryrun=False):
	"""
	1) parse the data from the json file
	2) run prep_vcf program for each tool's vcf

	:param tool: tool name for current data processing
	:param data: dictionary of converted json data
	:return: list of expected prep_vcfs from each tool unless error
	"""
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
			'--germline-snames', quote_str(data[tool]['germline_snames']),
			'--vcf', quote_str(data[tool]['vcf']),
			'-g', quote_str(data[tool]['ref_genome_fasta']),
			'-o', quote_str(data[tool]['prepped_vcf_outfilename']),
			'--bam', quote_str(data[tool]['bam']),
			'--contigs-file', quote_str(data[tool]['contigs_file']),

		]
	)

	## capture if user wants to make the Venn/upset plots later on before merging vcfs
	if data[tool]['do_venn']:
		cmdLine = ' '.join([cmdLine, "--make-bed-for-venn"])

	if data[0]['delete_temps']:
		cmdLine = ' '.join([cmdLine, "--delete-temps"])

	# capture threshold AR found in json
	TH_AR = data[tool]['threshold_AR']
	if TH_AR is not None and TH_AR != "" and TH_AR != 0.9:
		cmdLine = ' '.join([cmdLine, "--threshold-AR", TH_AR])

	# display the command line for log purposes
	log.info(prep_germline_vcf_script_path + " " + cmdLine)
	my_command = ' '.join(["bash", prep_germline_vcf_script_path, cmdLine])
	logFilename = "log_prep_vcf_{}.logs".format(tool)
	subp_logfile = open(logFilename, "w")
	if not dryrun:
		log.info("")
		log.info("%" * 10 + " prep {} vcf ".format(tool).upper() + "%" * 10)
		log.info("logging prep steps to file: " + str(logFilename))
		log.info("")
		process = subprocess.Popen(my_command,
		                           shell=True,
		                           universal_newlines=True,
		                           stdout=subp_logfile,
		                           stderr=subp_logfile)
		process.wait()

		log.info("return code value prep step for tool {} is: {}".format(tool, str(process.returncode)))
		log.info("prep step for tool {}: {} seconds".format(tool, str(int(round((time.time() - start_time))))))
		subp_logfile.close()
		if process.returncode is not 0:
			sys.exit("{} FAILED for tool {} ".format(prep_germline_vcf_script_path, tool))

def parse_json_data_and_run_prep_vcf_parallel(tool, data, dryrun=False):
	"""
	1) parse the data from the json file
	2) run prep_vcf program for each tool's vcf

	:param data: dictionary of converted json data
	:return: list of expected prep_vcfs from each tool unless error
	"""
	log.info("%" * 10 + "  " + str(tool).upper() + "  " + "%" * 10)
	log.info("recap inputs captured from json file")
	log.info("input: \ttool\t==\t{}".format(str(tool)))

	for var, val in data[tool].items():
		log.info("input:\t{} \t==\t{}".format(var, quote_str(val)))
	print()

	# check if outfilename for prep_vcf is Null or None
	if data[tool]['prepped_vcf_outfilename'] is None or data[tool]['prepped_vcf_outfilename'] == "":
		msg = "prepped_vcf_outfilename has not been defined in json file for tool {} ; Aborting.".format(tool)
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

	if data[tool]['delete_temps']:
		cmdLine = ' '.join([cmdLine, "--delete-temps"])

	# capture threshold AR found in json
	TH_AR = data[tool]['threshold_AR']
	if TH_AR is not None and TH_AR != "" and TH_AR != 0.9:
		cmdLine = ' '.join([cmdLine, "--threshold-AR", TH_AR])

	# display the command line for log purposes
	log.info(prep_vcf_script_path + " " + cmdLine)
	my_command = ' '.join(["bash", prep_vcf_script_path, cmdLine])
	logFilename = os.path.sep.join([ data[tool]['dir_work'] ,"log_prep_vcf_{}.logs".format(tool) ])
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
		log.info("return code value prep step for tool {} is {}".format(tool, str(process.returncode)))
		log.info("prep step for tool {}: {} seconds".format(tool, str(int(round((time.time() - start_time))))))
		subp_logfile.close()
		if process.returncode is not 0:
			sys.exit("{} FAILED for tool {} ".format(prep_vcf_script_path, tool))

def parse_json_data_and_run_prep_vcf_germline__DEPRECATED(data, dryrun=False):
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
				'--germline-snames', quote_str(data[tool]['germline_snames']),
				'--vcf', quote_str(data[tool]['vcf']),
				'-g', quote_str(data[tool]['ref_genome_fasta']),
				'-o', quote_str(data[tool]['prepped_vcf_outfilename']),
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
		log.info(prep_germline_vcf_script_path + " " + cmdLine)
		my_command = ' '.join(["bash", prep_germline_vcf_script_path, cmdLine])
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
			log.info(str(process.returncode))
			subp_logfile.close()
			if process.returncode is not 0:
				sys.exit("{} FAILED for tool {} ".format(prep_germline_vcf_script_path, tool))

def parse_json_data_and_run_prep_vcf__DEPRECATED(data, dryrun=False):
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
		log.info(prep_vcf_script_path + " " + cmdLine)
		my_command = ' '.join(["bash", prep_vcf_script_path, cmdLine])
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
			log.info(str(process.returncode))
			subp_logfile.close()
			if process.returncode is not 0:
				sys.exit("{} FAILED for tool {} ".format(prep_vcf_script_path, tool))


def subprocess_cmd(command):
	ev = os.system(command)
	log.info("return code from os.system is: "+ str(ev))
	if ev != 0:
		log.error("ERROR: subprocess command for command <<\"{}\">> FAILED".format(command))
		exit(ev)

def prepare_bed_for_venn(vcf, dirout):
	'''
	if no beds have been provided to vcfMerge2.py using --beds option and --do-venn has been enabled, and ...
	the skip-prep-vcf is used, we can make the beds from the vcf file(s) provided. As we already have the
	function << prepare_input_file_for_Venn >> in the bash script named << prep_vcf.sh >>, we will source the function
	and run it using system.command()
	:param vcf:
	:return: none
	'''

	# Build subprocess command;  ## I know this is defintely not the best implementation ever; I need to think about another strategy; Shame on me.
	# We prepare the bed files for snvs+indels venn, for snvs_only venns and for indels-only venns
	for FUNC in ["prepare_input_file_for_Venn", "prepare_input_file_for_Venn_SplitbyVariantType"]:
		mycmd = ["source", prep_vcf_functions_script_path, " && ", FUNC, " " , vcf, dirout]
		log.info(str(mycmd))
		log.info(" ".join([x for x in mycmd]))
		log.info(("Running bash function prep_input_file_for_venn command"))
		subprocess_cmd(' '.join([str(x) for x in mycmd]))


def merging_prepped_vcfs(data, merged_vcf_outfilename, delim, lossy, dryrun, do_venn, lbeds, skip_prep_vcfs, dirout, cmdline=None):
	"""

	:param data:
	:param merged_vcf_outfilename:
	:return: None
	"""

	list_vcfs = ""
	list_tools = ""
	list_tools_acronyms = ""
	list_precedence_order = ""
	for tool in data.keys():
		vcf = data[tool]['vcf']


		prepped_vcf = data[tool]['prepped_vcf_outfilename'] if not skip_prep_vcfs else  data[tool]['vcf']
		acronym = data[tool]['tool_acronym'] if data[tool]['tool_acronym'] != "" else str(tool).upper()
		list_tools = delim.join([list_tools, tool]) if list_tools != "" else tool
		vcf_to_add = prepped_vcf if prepped_vcf != "" else vcf
		list_vcfs = delim.join([list_vcfs, vcf_to_add]) if list_vcfs != "" else vcf_to_add
		list_tools_acronyms = delim.join([list_tools_acronyms, acronym]) if list_tools_acronyms != "" else acronym
		list_precedence_order = data[tool]['tool_precedence_order'] if list_precedence_order != "" else None
		if list_precedence_order is not None and ( list_precedence_order != "" and len(list_precedence_order) != len(data.keys())):
			list_precedence_order = ""
		ref_genome_fasta_dict = data[tool]['ref_genome_fasta_dict']

	log.info(str(list_tools_acronyms))
	my_command = ' '.join(["python", vcfmerger_tool_path,
	                       "--toolnames", double_quote_str(list_tools),
	                       "--vcfs", double_quote_str(list_vcfs),
	                       "-o", merged_vcf_outfilename,
	                       "-a", double_quote_str(list_tools_acronyms),
	                       "--dict", ref_genome_fasta_dict,
	                       "-d", dirout
	                       ])

	if list_precedence_order is not None and list_precedence_order != "":
		my_command = my_command + " --precedence " + double_quote_str(delim.join(list_precedence_order))

	if lossy:
		my_command = my_command + " --lossy "

	if data[tool]['delete_temps']:
		my_command = my_command + " --delete-temps "

	if do_venn:

		if data[tool]['venn_title'] is not None or data[tool]['venn_title'] != "":
			my_command = my_command + " --venn-title " + double_quote_str(data[tool]['venn_title']) + "  "

		for tool in data.keys():
			if not skip_prep_vcfs:
				list_beds = delim.join([str(os.path.splitext(vcf)[0] + ".intervene.bed") for vcf in
				                        list_vcfs.split(delim)])  ## extension intervene.bed defined in prep_vcf.sh
				log.info("list_bed for venn is: " + str(list_beds))
			elif lbeds == "":
				## as we skipped the prparation of vcfs, we already assigned in code before the vcfs to the prepped_vcf_outfilename field; so they should be equivalent
				log.info("for intervene tool, making bed from vcf intended to be merged : " + str(tool))
				log.info("trying to create on the fly the bed file using function in prep_vcf.sh script")
				log.info(tool + " __ prepare_bed_for_venn __  " + data[tool]['prepped_vcf_outfilename'] + " " + dirout)
				prepare_bed_for_venn(data[tool]['prepped_vcf_outfilename'], dirout)
				list_beds = delim.join(
					[os.path.join(dirout,str(os.path.splitext(os.path.basename(vcf))[0] + ".intervene.bed")) for vcf in list_vcfs.split(delim)])
			else:
				list_beds = lbeds

		if len(list_beds) == 0:
			sys.exit("ERROR: --do-venn provided, but list_beds file EMPTY ; check you inputs. Aborting.")
		log.info(double_quote_str(list_beds))
		if len(list_beds.split(delim)) != len(list_tools.split(delim)):
			sys.exit(
				"ERROR: --do-venn provided, but number of bed files ({}) DO NOT matched number of tools ({})  ; check you input or contact your IT/HPC admin".format(
					list_beds, list_tools))
		my_command = my_command + " --do-venn --beds " + double_quote_str(list_beds)

	if cmdline is not None or cmdline != "":
		my_command = my_command + " --cmdline " + double_quote_str(str(cmdline))

	log.info(double_quote_str(list_tools))
	log.info(double_quote_str(list_vcfs))
	log.info("merging vcfs ...")
	log.info(my_command)

	if not dryrun:
		log.info("")
		log.info("%" * 10 + " starting vcfMerger command... ".upper() + "%" * 10)
		process = subprocess.Popen(my_command, shell=True)
		process.wait()
		log.info("vcfMerger2.0 exit value: " + str(process.returncode))
		if process.returncode is not 0:
			sys.exit("{} FAILED with VCF files: {} ".format(vcfmerger_tool_path, list_vcfs))
		else:
			cpus = len(data.keys()) if len(data.keys()) > 2 else 2
			zvcf = str(merged_vcf_outfilename + ".gz")
			log.info("compressing vcf file using bcftools; final merged vcf name : " + zvcf)
			mycmd = ["bcftools view --threads", cpus, "-O z -o ", zvcf, merged_vcf_outfilename, ";",
			         "bcftools index --threads", cpus, "--tbi ", zvcf]
			subprocess_cmd(" ".join(str(x) for x in mycmd))


def check_path_to_vcfs(lvcfs):
	iterator = iter(lvcfs)
	while True:
		try:
			vcf = next(iterator)
			if not os.path.exists(vcf):
				log.error("ERROR:  FILE NOT FOUND  --->  Check your input for vcf:" + vcf)
				sys.exit(-1)
		except StopIteration:
			break  # Iterator exhausted: stop the loop
		else:
			log.info("VCF found: " + str(vcf))

def check_inputs(lvcfs, ltoolnames, ltpo=None, lacronyms=None, lprepped_vcf_outfilenames=None, lbeds=None,
                 germline=False, tumor_sname=None, normal_sname=None, germline_snames=None, merged_vcf_outfilename=None,
                 filter_by_pass=False, filter_string_for_snpsift=None,
                 path_jar_snpsift=None, fastaDicoFile=None):
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

	if fastaDicoFile is None:
		log.info("ERROR: fasta dictionary file MUST be provided to sort the contigs in the header correctly")
		sys.exit("ERROR: fasta dico file missing ; use --dict option and provide full path to the  dict file ")
	## Checking snpsift path to jar
	if (filter_by_pass or filter_string_for_snpsift is not None) and path_jar_snpsift is None:
		log.error("ERROR: You enabled a filter option but did not provide any path to snpSift.jar ; please provide the FULL PATH to SnpSift.jar file; Aborting!")
		sys.exit(-1)
	elif (filter_string_for_snpsift is not None or filter_by_pass) and path_jar_snpsift is not None:
		log.info("Path to provided snpSift.jar file:" + str(path_jar_snpsift))
		if not os.path.exists(path_jar_snpsift):
			raise Exception("ERROR: snpSift.jar FILE NOT FOUND. Aborting!")
	elif filter_string_for_snpsift is not None and path_jar_snpsift is None:
		raise Exception(
			"Please provide the Full PATH to a snpSift.jar file using the option --path-jar-snpsift. Aborting!")
	elif path_jar_snpsift is not None:
		if not os.path.exists(path_jar_snpsift):
			raise Exception("ERROR: snpSift.jar FILE NOT FOUND. Aborting!")
		log.info("Well, you provided the path to snpSift probably before the options for filtering... that is ok. Otherwise, well you have not set the filter option. and provided the path to snpSift.jar for nothing :-) ")
	else:
		log.info("No Path given for SnpSift")


	if filter_string_for_snpsift is not None and (len(filter_string_for_snpsift.split("###")) != len(ltoolnames) or len(filter_string_for_snpsift.split("###")) == 0 ):
		log.error("ERROR: Number of triple-pound separated Values in --filter option does NOT match the number of given toolnames or number of given vcfs; "
		          ". Check your inputs; Check if triple pouns are well used to separate out the values for that option.")
		sys.exit(-1)
	elif filter_string_for_snpsift is not None and ( len(filter_string_for_snpsift.split("###")) == 0 or filter_string_for_snpsift != ""):
		log.info(
			"FYI: Using the same filtering for all the vcfs; This implies all the VCFs contains the correct flags and fields.")


	#checking if merged_vcf_outfilename has relative or full path included
	if merged_vcf_outfilename is None:
		msg = "ERROR: Missing filename; Please provide a name for the merged final VCF; (fullpath, relative path or just basename)"
		log.error(msg)
		sys.exit(-1)
	dn = os.path.dirname(merged_vcf_outfilename)
	if dn == '':
		log.info("final vcf will be written in current directory: " + os.path.realpath(os.curdir))
	elif (dn != '' and dn != "." and dn is not None) and not os.path.exists(dn):
		msg = "ERROR: path to the merged vcf outfilename NOT FOUND relative to current path; Check your inputs ; folder << {} >> NOT FOUND. Please create that directory first.".format(
			dn)
		log.error(msg);
		sys.exit(-1)
	elif os.path.exists(dn):
		log.info("merged_vcf will be written into folder: " + dn)
	else:
		log.error("dn object has a weird assignment; Why am I here?? ; dn is <" + dn + ">")
		sys.exit(-1)
	log.info("filename for the uncompressed merged output vcf will be: " + merged_vcf_outfilename)
	log.info("filename for the bgzip-compressed merged output vcf will be: " + merged_vcf_outfilename + ".gz")

	if lprepped_vcf_outfilenames is not None and len(ltoolnames) != len(lprepped_vcf_outfilenames):
		log.info(
			"ERROR: Found {} in list of given prep-filenames and {} toolnames given".format(
				len(lprepped_vcf_outfilenames),
				len(ltoolnames)))
		sys.exit(
			"ERROR: Number of toolnames MUST be equal to the number of intermediate prep-outfilenames given ;\ncheck if delimiter is adequate and do not interfere with splitting the given lists of tools")
	check_path_to_vcfs(lvcfs)
	if germline and germline_snames is None:
		log.error(
			"ERROR: germline was enabled but no germline sample names given; please use option: --germline-snames and provide a DELIM-list of sample according to vcf content ")
		sys.exit(-1)
	if not germline:
		if tumor_sname is None or normal_sname is None or (tumor_sname is None and normal_sname is None):
			log.error(
				"ERROR: Somatic was enabled but no either/or/both tumor or/and normal sample names were not given; Provide options: --tumor-sname and --normal-sname with expected sample names")
			sys.exit(-1)
	if germline_snames is not None and germline is False and (tumor_sname is not None or normal_sname is not None):
		log.error(
			"ERROR: Ambiguous inputs and options; germline and soamtic analyses are EXCLUSIVE ; use --germline option to stipulate processing germlien calls and provide --germline-snames as well ; if only somatic, provide only tumor-sname and normal-sname")
		sys.exit(-1)


def main(args, cmdline):
	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)

	delim = "|"
	if args["delim"]:
		delim = args["delim"]
		if len(delim) != 1:
			exit("the list delimiter given with --delim should be one character only and NOT a space; found {}".format(
				str(len(delim))))
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

	ref_genome_fasta_dict = None
	if args["dict"]:
		ref_genome_fasta_dict = args["dict"]
		if not os.path.exists(ref_genome_fasta_dict):
			sys.exit("dictionnary file of reference genome {} NOT FOUND; Aborting.".format(ref_genome_fasta_dict))

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
		log.info(
			"tool-specific filenames for intermediate vcfMerger2-up-to-specs vcf: " + str(lprepped_vcf_outfilenames))
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

	skip_prep_vcfs = False
	if args["skip_prep_vcfs"]:
		skip_prep_vcfs = args["skip_prep_vcfs"]
		log.info("skip_prep_vcfs:" + str(skip_prep_vcfs))

	filter_by_pass = False
	if args['filter_by_pass']:
		filter_by_pass = True
		log.info("filtering by PASS variant enabled; This filter will be applied before the PREP_VCF step")

	filter_string_for_snpsift = None
	if args["filter"] is not None:
		filter_string_for_snpsift = args["filter"]
		log.info("filtering variant enabled and filter string to be used with snpSift: \"" + str(
			filter_string_for_snpsift) + "\"")

	path_jar_snpsift = None
	if args['path_jar_snpsift'] is not None:
		path_jar_snpsift = args['path_jar_snpsift']



	germline = False
	if args["germline"]:
		germline = True
		log.info("GERMLINE analysis enabled")

	germline_snames = None
	if args['germline_snames']:
		germline_snames = args['germline_snames']
		log.info("germline_snames == " + str(germline_snames))

	skip_merge = False
	if args["skip_merge"]:
		skip_merge = args["skip_merge"]
		log.info("skip_merge:" + str(skip_merge))

	TH_AR = 0.90 ; ## HARDCODED
	if args['threshold_AR']:
		TH_AR = args['threshold_AR']
		if isinstance(TH_AR, (int, float, complex)):
			raise Exception(
				"Threshold-AR must be a float or integer value between 0 and 1 (range ]0,1]). Check your inputs.")
		log.info("user given threshold for AR: " + str(TH_AR))

	delete_temps=False
	if args['delete_temps']:
		delete_temps = True
		log.info("Temps Files from prep_vcf stage will be DELETED")

	lbeds = ""
	if args["beds"]:
		lbeds = str(args["beds"]).split(delim)
		log.info("ordered list of beds given:\t\t{}".format(str(lbeds)))

	do_venn = False
	if args["do_venn"]:
		do_venn = True
		log.info("make venn enabled")
		list_executables = ['Rscript', 'intervene']
		check_if_executable_in_path(list_executables)


	venn_title = ""
	if args["venn_title"]:
		venn_title = args['venn_title']
		log.info("venn title will be: "+venn_title)
    
	dirout = os.curdir
	if args["dir_out"]:
		dirout = args["dir_out"]
		try:
			if not os.path.exists(dirout):
				log.info("creating temp directory recursively, present or not")
				os.makedirs(dirout, exist_ok=True)
			dirout = os.path.realpath(dirout)
			log.info(dirout + " created")
		except Exception as e:
			log.info(e)
			sys.exit(-1)
		log.info("TEMP FOLDER is:" + dirout )

	dryrun = False
	if args["dry_run"]:
		dryrun = args["dry_run"]
		log.info("dry-run?:" + str(dryrun))

	# if inFileJson is None:
	##@@@@@@@@@
	## MAIN  ##
	##@@@@@@@@@
	check_inputs(lvcfs, ltoolnames, ltpo=list_tool_precedence_order, lacronyms=lacronyms,
	             lprepped_vcf_outfilenames=lprepped_vcf_outfilenames, lbeds=lbeds,
	             germline=germline, tumor_sname=tumor_sname, normal_sname=normal_sname,
	             germline_snames=germline_snames, merged_vcf_outfilename=merged_vcf_outfilename,
	             filter_by_pass=filter_by_pass, filter_string_for_snpsift=filter_string_for_snpsift,
	             path_jar_snpsift=path_jar_snpsift, ref_genome_fasta_dict=ref_genome_fasta_dict)


	lvcfs = check_if_vcf_is_compressed(lvcfs, dirout)
	log.info(str(lvcfs))

	data = make_data_for_json(lvcfs,
	                          ltoolnames,
	                          normal_sname,
	                          tumor_sname,
	                          ref_genome_fasta,
	                          ref_genome_fasta_dict,
	                          lossy,
	                          germline_snames=germline_snames,
	                          ltpo=list_tool_precedence_order,
	                          lacronyms=lacronyms,
	                          lprepped_vcf_outfilenames=lprepped_vcf_outfilenames,
	                          lbams=lbams,
	                          lcontigs=lcontigs,
	                          filter_string_for_snpsift=filter_string_for_snpsift,
	                          TH_AR=TH_AR,
	                          do_venn=do_venn,
	                          venn_title=venn_title,
	                          skip_prep_vcfs=skip_prep_vcfs,
	                          dirout=dirout,
	                          delete_temps=delete_temps)
	json_filename = os.path.join(dirout, "vcfMerger2_somatic.json") if not germline else os.path.join(dirout, "vcfMerger2_germline.json")
	make_json(data, json_filename)
	# inFileJson = make_json(data, json_filename)
	# data = read_json(inFileJson) ## uncomment for debugging if necessary ; data is already created above

	## filter the PASS records before preparing the vcf for vcfMerger step
	if filter_by_pass:
		log.info("Performing PASS only filtering for ALL the provided input vcfs ...")
		data = filter_unprepped_vcf(data, path_jar_snpsift)
		log.info(str(data))

	if not skip_prep_vcfs:  ## PREP_VCF step Enabled
		log.info("**** prep vcf steps section ***".upper())
		# parse_json_data_and_run_prep_vcf(data, dryrun) if not germline else parse_json_data_and_run_prep_vcf_germline(
		# 	data, dryrun)
		time.sleep(1)
		## PARALLELIZING MAKING PREP FILES
		import multiprocessing
		funcToUse = parse_json_data_and_run_prep_vcf_parallel if not germline else parse_json_data_and_run_prep_vcf_germline_parallel
		procs = []
		for tool in data.keys():
			p = multiprocessing.Process(target=funcToUse, name=tool, args=(tool, data, dryrun))
			procs.append(p)

		for p in procs:
			p.start()
			time.sleep(1)

		procs_exit_codes = {}
		for p in procs:
			log.info("task " + str(p) + " started in background ...")
			p.join()

		for p in procs:
			procs_exit_codes[p.name] = p.exitcode

		STOP = bool(False)
		for name, ev in procs_exit_codes.items():
			msg = " prep ev value of {} for tool {} ".format(ev, name)
			if ev == 0:
				log.info(msg)
			else:
				log.error(msg + " --> FAILED")
				STOP = bool(True)

		if STOP:
			log.info("prep step Elapsed time before failure:  {} secs".format(str(int(round((time.time() - start_time))))))
			sys.exit("Aborting vcMerger2 because preps FAILED")

		log.info("prep step Total Elapsed time in seconds:  {}".format(str(int(round((time.time() - start_time))))))
		log.info("**** end PREP process section  ****".upper())

	else:
		log.info("**** SKIPPED prep vcfs step SKIPPED ****")

	## FILTERING STEP Enabled for Already PREPPED VCFS (note: Filtering must not be applied on un-prepped vcfs unless user knows what he/she is doing
	if filter_string_for_snpsift is not None:
		log.info("Performing vcf filtering for ALL the provided prepped vcfs ...")
		data = filter_prepped_vcf(data, path_jar_snpsift)
		log.info(str(data))

	if not skip_merge:  ## MERGING step Enabled
		merging_prepped_vcfs(data, merged_vcf_outfilename, delim, lossy, dryrun, do_venn, lbeds, skip_prep_vcfs, dirout, cmdline=cmdline)
	else:
		log.info("**** SKIPPED merge step SKIPPED ****")

	## MESSAGES END of WORK
	if not skip_prep_vcfs and not skip_merge:  ## we perform the ALL in one [ PREP_VCF + MERGE ]
		log.info("prep and merge vcfs Elapsed time in seconds:  {}".format(str(int(round((time.time() - start_time))))))
		log.info("prep and merge vcfs completed successfully")
	elif skip_prep_vcfs:  ## Here we only performed MERGING
		log.info("merge vcfs Elapsed time in seconds:  {}".format(str(int(round((time.time() - start_time))))))
		log.info("merge vcfs completed successfully")
	elif skip_merge:  ## here we only performed PREP_VCF
		log.info("prep Elapsed time in seconds:  {}".format(str(int(round((time.time() - start_time))))))
		log.info("prep vcfs completed successfully")


def make_parser_args():
	parser = argparse.ArgumentParser()
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
	required.add_argument('-g', '--refgenome',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='reference genome used with bcftools norm ; must match reference used for alignment')

	required.add_argument('--dict',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='dictionary file of reference genome; required to get correct order of contig names; this should be a .dict file created by picard or samtools sequence dictionary module')

	required.add_argument('-o', '--merged-vcf-outfilename',
	                      required=isRequired,
	                      action=UniqueStore,
	                      help='outfilename for the merge vcf (can be relative or full path)')

	optional.add_argument('--germline',
	                      required=False,
	                      action='store_true',
	                      help='option required if dealing with GERMLINE VCFs, otherwise data will be considered as Somatic calls')

	optional.add_argument('--germline-snames',
	                      required=False,
	                      action=UniqueStore,
	                      help='expected name of germline sample(s) in vcf file if option --germline is in use')

	optional.add_argument('--normal-sname',
	                      required=False,
	                      action=UniqueStore,
	                      help='expected name of normal sample in vcf file')

	optional.add_argument('--tumor-sname',
	                      required=False,
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

	optional.add_argument('--delete-temps',
	                      required=False,
	                      action='store_true',
	                      help="if set, temporary files created during the prep_vcf stage will be deleted")

	optional.add_argument('--filter-by-pass',
	                      required=False,
	                      action='store_true',
	                      help='enable vcf filtering of PASS variant using snpSift tool; String to snpSift hardcoded as (FILTER == "PASS") ')

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

	optional.add_argument('--venn-title',
	                      required=False,
	                      action=UniqueStore,
	                      help="Default is empty string")


	optional.add_argument('-d','--dir-out', '--dir-temp',
	                      help=' direcgtory where the outputs of vcfMerger2 will be written ',
	                      action=UniqueStore)

	optional.add_argument('-n', '--dry-run',
	                      required=False,
	                      action='store_true',
	                      help=' perform a dry run to just see the commands lines that will be run behind the scenes ')

	log.info(str(parser.prog) + "   " + str(parser.description))
	return parser


def check_if_executable_in_path(list_executables):
	for executable in list_executables:
		if shutil.which(executable) is None:
			sys.exit(str(executable) + "  NOT IN PATH ; Aborting;")


if __name__ == '__main__':
	print("   ")
	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)
	list_executables = ['bedtools', 'bcftools', 'samtools', 'python3']  ## list of executables required to be in path
	check_if_executable_in_path(list_executables)

	##TODO check if ALL intended python packages are present in python3

	# capture time for calculating vcfMerger's runtime
	start_time = time.time()
	# capture commandline for adding it to vcf header
	cmdline = ' '.join(sys.argv)
	parser = make_parser_args()
	args = vars(parser.parse_args())  # vars() function returns a dictionary of key-value arguments
	log.info(str(args))
	log.info(' '.join(["Command Line captured: ", cmdline]))
	# try:
	main(args, cmdline)
	# except Exception as e:
	# 	log.error("ERROR: Exception got Raised; Check you inputs and/or options".format(e))
	# 	sys.exit(1)
	sys.exit()
