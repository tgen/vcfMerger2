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


# from multiprocessing import Process, Queue, cpu_count
import argparse
from sys import exit
from sys import argv
from os import linesep
from os import path
from re import search
import dvmfuncs as dvm
import vcfToDict
import logging as log
from collections import OrderedDict
from natsort import natsorted
import time
import subprocess

'''
## STEPS so far:
## 1) we check if list of vcfs is correct and number of files >=2  ---> DONE
## 2) we check if vcfs are found --> DONE
## 3) we check if number of toolnames matched the number of vcf filenames. --> DONE
## 4) we create object that manages vcfs ; this object or "class" will be able to :
##	a) read the vcf by capturing the header only
##	b) read the entire vcf and put the variants into a dictionnary ; the Keys are going to be based off CHROM_POS_REF_ALT whether the ALT contains multiple alleles or not ; or JUST CHR_POS_REF ; we will use the module cyvcf2 from brent (https://github.com/brentp/cyvcf2) which is in part based off the pysam package to read vcfs
## 5) once all the dictionaries have been created for each vcf, we MERGE them using my trick
## 6) we need to check all the cases and the matching between the prevalence
## 7) how the prevalence is going to be managed? well, simple: the VCFs and toolnames have been given in a special ORDER; this order is kept in the order of the dictionaries merged one after the other;
##    because of that order, the data within the dictionnary is an ordrered list and therefore keep the order of the tool for each common variant; This implies that the user gives the toolnames and vcf files
##    in the order of the prevalence already. example: strelka, viper, mutect2. is different from viper, strelka , mutect2; 
## 
#######################################################################################################################
## cyvcf2 is a cython wrapper around htslib built for fast parsing of Variant Call Format (VCF) files.
## WARNING and IMPORTANT for cyvcf2 usage
## Attributes like variant.gt_ref_depths return a numpy array directly so they are immediately ready for downstream use.
## note that the array is backed by the underlying C data, so, once variant goes out of scope. The array will contain nonsense.
## To persist a copy, use: cpy = np.array(variant.gt_ref_depths) instead of just arr = variant.gt_ref_depths.
#######################################################################################################################
'''


class UniqueStore(argparse.Action):
	"""
	class grabbed from stackOverflow (2018-11-08)
	https://stackoverflow.com/questions/23032514/argparse-disable-same-argument-occurences
	Thanks To the Community
	We override the function __call__ from argparse to check if one option is given more than once
	"""

	def __call__(self, parser, namespace, values, option_string):
		if getattr(namespace, self.dest, self.default) is not None:
			parser.error(option_string + " appears several times.")
		setattr(namespace, self.dest, values)


def process_merging(lvcfs, ltoolnames, list_tool_precedence_order, lossless, merge_vcf_outfilename, cmdline):
	"""
	process the merging after having checked all the inputs
	:param lvcfs: list of vcf tools delimited by the delim value
	:param ltoolnames: list of names, toolnames, given in same order as the vcf
	:param list_tool_precedence_order: if specified, this list of toolnames given in a specific order (different or
	same as order inltoolnames, will be used for precedence order; NOTE: the precedence is only a term here that
	means the value in FORMAT fields will come from the first tool which call the variant.
	:param lossless: the INFO column will be FILLED with ALL the information coming from ALL the TOOLS and not only
	the first tool which got precedence
	:param merge_vcf_outfilename:  name of the vcf with the merged info
	:param cmdline: command line run to start the current tool
	:return: None (actually the function creates and populates the merge_vcf_outfilename)
	"""

	outputFilename = merge_vcf_outfilename
	tuple_objs = ()
	l_snames = []
	l_contigs = []

	ListFieldsToProcessForOurFORMATColumn = ["GT", "DP", "AR", "AD"]

	vcfMerger_Format_Fields_Specific = [
		'##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
		'##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth at locus in Sample">',
		'##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed from chosen prevalent tool">',
		'##FORMAT=<ID=AR,Number=1,Type=Float,Description="Allele frequency of ALT allele from chosen prevalent tool">'
	]

	TN_FLAGS = []
	for T in ltoolnames:
		TN_FLAG = str(''.join([
			'##INFO=<ID=' + T + ',Number=0,Type=Flag,Description="Toolname Flag means that position got '
			                    'called by this tool">']))
		TN_FLAGS.append(TN_FLAG)
	Additional_FLAGS = [
		'##INFO=<ID=CC,Number=1,Type=Integer,Description="CALLERS_COUNT,Number of tools calling this variant event '
		'out of a total of ' + str(len(ltoolnames)) + ' tools">',
		''.join(['##INFO=<ID=TPCE,Number=1,Type=String,Description="Tool that got precedence for called position; '
		         'user gave the following order for tool precedence: ', ', '.join([str(t) for t in
		                                                                           ltoolnames]),
		         '">']),
		'##INFO=<ID=TYPE,Number=1,Type=String,Description="Type of Variant">'
	]

	vcfMerger_Info_Fields_Specific = TN_FLAGS + Additional_FLAGS

	if list_tool_precedence_order is not None:
		'''here we sort and reassigned ltoolnames and lvcfs based on list_tool_precedence_order ; names of the 
		tools have to match 100%
		'''
		if len(list_tool_precedence_order) != len(ltoolnames):
			exit("ERROR: Tool Names in list precedence do not match 100% names in list toolnames ; check your "
			     "input\n" + "sorted_list_tool_precedence -> " + str(sorted(list_tool_precedence_order)) +
			     "\nsorted_list_tool_names ------> "
			     + str(sorted(ltoolnames)))
		indices = []
		for toolname in list_tool_precedence_order:
			indices.append(ltoolnames.index(toolname))
		## we reallocate/reorder the vcfs files the same order of the list_tool_precedence_order
		lvcfs = [lvcfs[i] for i in indices]
		ltoolnames = list_tool_precedence_order;  ## we re-assigned the list

	# the trick is here for the Tool Precedence!!! The user has given us an ordered list of
	# vcfs and toolnames in order of precedence or a specific PRECEDENCE order was given via --precedence
	# and we sort the vcf and add them to the tuple accordingly
	for i in range(len(lvcfs)):
		o = vcfToDict.vcfToDict(lvcfs[i], ltoolnames[i])  ## here we map the toolname and the vcf associated
		tuple_objs = tuple_objs + (o,)  ## we add instances of object vcfToDict to the tuple ; order FIFO is
		# equivalent to the order of precedence
		l_snames.append(o.samplenames)  ## we add tuples of samplenames to the list l_snames as a list of tuples
		l_contigs.append(sorted(o.contigs))

	# performing checks before processing data further
	dvm.compareTuples(l_snames,
	                  "SampleNames")  ## we cannot skip that one. If not matching, then modify vcf to get samples in
	# correct columns or with the same names across ALL the vcf files ;

	## UNCOMMENT NEXT LINE TO PUT THE CONTIGS CHECK BACK ON
#########	dvm.compareTuples(l_contigs, "CONTIGS")  ## we may add an option to skip that check ; even though we do not know
	# what could be the consequences of having different contigs ; we cannot think any so far.

	"""
	### we check here the presence of the expected MANDATORY fields in the FORMAT columns ;
	### Unfortunately as we do not read the entire VCFs file and therefore we do not have the object VCF created yet,
	### we cannot use the cyvcf2 API to check if an ID is defined in the VCF header or not, or in the variant or not;
	### So for now, we rely on our own vcf header capture as string; we therefore check the string;
	### BUT: this does not mean that the ID fields are necessary present in each variant;
	### If we want to check that presence, we will have to read the vcf files entirely see below "tuple_dicts = () loop" ;
	### and check every variant.
	### Or, second option, we will check while we merge and raise ERROR and either stop merging or skip that variant, or put NULL value for that field ;
	### for example: if AR does not exist, we set AR=.
	"""

	check_fields_definitions_in_header = True
	if check_fields_definitions_in_header:
		for flagid in ListFieldsToProcessForOurFORMATColumn:
			log.info(flagid)
			for tpo in tuple_objs:
				'''Check if flag we want to put in the format field have been defined in the VCF header'''
				res_search = search("".join(["ID=", flagid]), tpo.headers)
				if res_search is None:
					exit(
						"Id Flag " + flagid + " not Defined in header of vcf file " + tpo.fvcf + ".\nPlease bring the VCF up to specs before running this merging tool. Use a wrapper specific to your tool which has already been created by the Author of the current tool. Aborting!")


	# we process the files entirely after all the checks have PASSED successfully
	# we may make parallel this step But If we do, we lose the precedence order in the tuple_dicts variable and
	# this defies the purpose of that script
	tuple_dicts = ()
	for tpo in tuple_objs:
		tuple_dicts = tuple_dicts + (tpo.dictOfLoci(tpo.readVCF()),)

	# we merge the Loci from all the VCFs [Key + Value, where Key is defined as CHROM_POS_REF_ALT as assigned in the function "dictOfLoci" of class vcfToDict ]
	from collections import defaultdict
	dd = defaultdict(list)

	log.debug("-" * 41);
	log.debug(str(type(tuple_dicts)))

	for d in tuple_dicts:
		for key, value in d.items():
			try:
				dd[key].append(value)
			except KeyError:  ## I do not see why we should have an error here because we just list the Keys
				# from d dicts we created ; I put it probably because it happened?
				log.warning("KEY ERROR Detected - Skipping this values ; It should not have happened; please "
				            "report that to the Author")
	# NOTE: in the loop above, to get your .attrib, just change append(value) to append(value.attrib)
	# You may then want to make a normal dict out of the defaultdict so you have normal dict behavior for non-existent keys etc: dd = dict(dd)

	# 1) first we managed the Headers from all the tools
	log.info("processing headers of all the vcf files ...")
	list_lines_header = dvm.create_new_header_for_merged_vcf(tuple_objs,
	                                                         cmdline,
	                                                         vcfMerger_Format_Fields_Specific,
	                                                         vcfMerger_Info_Fields_Specific
	                                                         )
	# 2) we add the modified header lines to the output merger file
	log.info("adding the header to the out vcf file ...")
	dvm.add_new_header_to_merged_file(outputFilename, list_lines_header, tuple_objs[0].header_chrom_line + "\n")

	# 3) we process all the variants
	log.info("looping over variant calls, merging and writing back to file ... ")

	try:

		of = open(outputFilename, 'a')  # we open the output file with merged information here
		# sort dico by keys before iterating over it ...  ## normally the Keys are not sorted because we deal with a dictionary which do not keep the order
		dd = OrderedDict(sorted(dd.items()))
		natural_sorted_keys = natsorted(dd.keys())
		# dd.keys --> they are the KEYS that are represented by the PATTERN --> CHROM_POS_REF_ALT
		# dd.values --> represents the calls and their information from each tool having call the variant at position CHROM_POS
		# (the number of list in values may go from 1 to len(lvcfs); where len(lvcfs) represents the total number
		# of inputs vcfs and therefore ALL the tools would have called that variant )
		# wtv stands for Winning Tool Variant ; It always is the first one, as the tools have been sorted by
		# precedence given by the user
		# 3a) get the total number variants to process in order to calculate on the fly the value for the counter
		# steps
		tot_variants_count = len(dd)
		counter = 0
		# step is ~10% of tot_variants and round to the nearest nth value
		step = int(round(tot_variants_count / 10, -(len(str(round(tot_variants_count / 10))) - 1)))
		for K in [k for k in natural_sorted_keys]:  # sub is list__list__o.ovcf_variant ;
			counter += 1;
			if counter % step == 0:
				log.info("processed {} variants ...".format(counter))
			rebuilt_variant = dvm.rebuiltVariantLine(dd[K], lossless,
			                                         ListFieldsToProcessForOurFORMATColumn);  ## dd[K} represent a List of Variants (LV)
			of.write(rebuilt_variant + linesep)
		log.info("total processed variants: {}".format(counter))


	except IOError as e:
		log.info("Error I/O({0}): {1}".format(e.errno, e.strerror))
		of.close()
	else:
		of.close()


def main(args, cmdline):
	FORMAT_LOGGING = '%(levelname)s %(asctime)-15s %(module)s %(lineno)d\t %(message)s'
	log.basicConfig(format=FORMAT_LOGGING, level=log.INFO)

	# we may give the user the choice of the delimiter or the choice within several pre-determined list
	# of delimiters ; if we give the choice to the user, this variable HAS to be the first to be tested and assigned
	# as other variables depend upon it
	delim = "|"
	if args["delim"]:
		delim = args["delim"]
		if len(delim) != 1:
			exit("the list delimiter given with --delim should be one character only and NOT a space; found {"
			     "}".format(str(len(
				delim))))
		log.info("delimiter is:\t" + delim)

	lvcfs = []
	if args["vcfs"]:
		lvcfs = str(args["vcfs"]).split(delim)
		log.info("ordered list of vcfs given:\t\t{}".format(str(lvcfs)))

	ltoolnames = []
	if args["toolnames"]:
		ltoolnames = str(args["toolnames"]).split(delim)
		ltoolnames = [x.upper() for x in ltoolnames]
		log.info("ordered list of toolnames given:\t{}".format(str(ltoolnames)))

	list_tool_precedence_order = None
	if args["precedence"]:
		list_tool_precedence_order = str(args["precedence"]).split(delim)
		list_tool_precedence_order = [x.upper() for x in list_tool_precedence_order]
		log.info("precedence: {} ".format(str(list_tool_precedence_order)))

	lossless = True
	lossy = False
	if args["lossless"]:
		log.info("lossless enabled")
		lossless = True

	if args["lossy"]:
		log.info("lossy enabled")
		lossy = True
		lossless = bool(not lossy)

	if args["toolacronyms"]:
		lacronyms = str(args["toolacronyms"]).split(delim)
		lacronyms = [x.upper() for x in lacronyms]
		log.info("tool acronyms:\t\t\t{}".format(lacronyms))

	if args["outfilename"]:
		merge_vcf_outfilename = str(args["outfilename"])
		log.info("name of merged output vcf will be: " + merge_vcf_outfilename)

	lbeds = ""
	log.info("before --- ordered list of beds given:\t\t{}".format(str(lbeds)))
	if args["beds"]:
		lbeds = str(args["beds"]).split(delim)
		log.info("after --- ordered list of beds given:\t\t{}".format(str(lbeds)))

	do_venn = False
	if args["do_venn"]:
		do_venn = True
		log.info("make venn enabled")

	# CHECK POINTS
	if lossy and lossless:
		sys.exit("lossy and lossless are mutually exclusive options, please use one or the other but not both.")
	if lvcfs is None:
		exit("ERROR: please provide the list of vcf you want to prep or merged using --lvcfs option with value delimited with the DELIMITER (--delim)")
	log.debug("n(vcfs) = {} --- n(tools) = {}".format(str(len(lvcfs)),str(len(ltoolnames))))
	if len(lvcfs) != len(ltoolnames):
		exit("ERROR: number of  vcfs  MUST match number of  toolnames ; order sensitive; case insensitive; Aborting.")
	if lacronyms is not None and len(lacronyms) != len(ltoolnames):
		exit("ERROR: number of acronyms should match number of toolnames; order sensitive ; Aborting.")
	if list_tool_precedence_order is not None and len(list_tool_precedence_order) != len(ltoolnames):
		exit("ERROR: number of precedence toolname should match number of toolnames; Aborting.")
	if do_venn:
		if lbeds == "":
			exit("ERROR: list of bed files for making Venn/Upset plots MUST be provided while using --do-venn option")
		dvm.make_venn(ltoolnames, lbeds, saveOverlapsBool=False, upsetBool=False)
		exit()

	process_merging(lvcfs, ltoolnames, list_tool_precedence_order, lossless, merge_vcf_outfilename, cmdline)

def make_parser_args():
	parser = argparse.ArgumentParser(description='Processing options ...')
	parser._action_groups.pop()
	required = parser.add_argument_group('required arguments')
	optional = parser.add_argument_group('optional arguments')

	required.add_argument('--vcfs',
	                      required=True,
	                      action=UniqueStore,
	                      help='List of vcfs file delimited by DELIM character; default DELIM is space ; or assign it '
	                           'using --delim option')
	required.add_argument('--toolnames',
	                      required=True,
	                      action=UniqueStore,
	                      help='List of vcfs file delimited by DELIM character; default DELIM is pipe unless --delim '
	                           'option is used '
	                           'using --delim option ; same DELIM as --vcfs is required')
	required.add_argument('-o', '--outfilename',
	                      required=True,
	                      action=UniqueStore,
	                      help='given name to the merge vcf output file')
	optional.add_argument('-c', '--precedence',
	                      required=False,
	                      action=UniqueStore,
	                      help=' sorted delim-separated list of the toolnames defined also with --toolnames ')
	optional.add_argument('--delim',
	                      required=False,
	                      default="|",
	                      help=' sorted delim-separated list of the toolnames defined also with --toolnames ')
	optional.add_argument('-a', '--toolacronyms',
	                      required=False,
	                      action=UniqueStore,
	                      help='List of Acronyms for toolnames to be used as PREFIXES in INFO field ; same DELIM as --vcfs ')
	optional.add_argument('--lossless',
	                      help='This will create a lossless merged vcf by adding FORMAT columns of secondaries tools into the INFO field',
	                      action='store_true')
	optional.add_argument('--lossy',
	                      help='This will create a merged vcf with FORMAT information kept from only the tool having precedence at the dealt position. lossless and lossy are mutually exclusive',
	                      action='store_true')
	optional.add_argument('--beds',
	                      help='list of bed files to be used to make Venns or Upset plots; requires to enable --do-venn as well to validate making Venn/upset plots ; list MUST be delimited by DELIM character (--delim or default delim)',
	                      action=UniqueStore)
	optional.add_argument('--do-venn',
	                      help='using the bed files listed in --beds option, Venns or Upset plot will be created ; need to match the number and order of tools listed in --toolnames ',
	                      action='store_true')
	optional.add_argument('--verbose',
	                      help='Output verbose information',
	                      action='store_true')
	print(str(parser.prog) + "   " + str(parser.description))
	return parser


if __name__ == '__main__':
	# capture time for calculating vcfMerger's runtime
	start_time = time.time()
	# capture commandline for adding it to vcf header
	cmdline = ' '.join(argv);
	del argv

	# parsing arguments
	parser = make_parser_args()
	args = vars(parser.parse_args())  # vars() function returns a dictionary of key-value arguments
	main(args, cmdline)

	log.info("vcfMerger Elapsed time in seconds:  %s" % int(round((time.time() - start_time))))
	log.info("vcfMerger completed successfully")
	exit()
