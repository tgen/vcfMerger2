#!/bin/bash

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


trap "exit 1" TERM
export TOP_PID=$$


## CONSTANT VARIABLE (modified accordingly)
DIR_PATH_TO_SCRIPTS="$( dirname $0 )"


## CONSTANT VARIABLE : add of modify toolnames accordingly
VALID_TOOLNAMES="lancet,  mutect2,  octopus,  strelka2, vardict, or these abbreviations, LAN|LCT, MUT, OCT, SLK, VDJ|VDT, respectively"  ## if tools are later added, we will update this variable along with the function run_tool(), where the case statement will need to be updated.


function usage(){
	echo -e "\nUSAGE:"
	echo -e "`basename $0` \\
-d|--dir-work    	DIR_WORK directory where outputs will be written (needs to exist))  \\
-g|--ref-genome   	REFERENCE GENOME FASTA FILE  [Required] \\
-t|--toolname		Provide the toolname associated to the input vcf [REQUIRED]; see valid toolnames in prep_vcf_defaults.ini file, or use --print-defaults in in-line command \\
-o|--prepped-vcf-outfilename	Provide the name for the uptospecs vcf file that will be use as input for the vcfMerger2.0 tool \\
--print-defaults	Print default valid toolnames accepted so far \\
--vcf			vcf having all types of variants already (no need to concatenate) \\
--vcf-indels		tool's VCF with indels (.vcf or .vcf.gz) ; note: toDate, concerns strelka2 only \\
--vcf-snvs		tool's VCF with snvs (.vcf or .vcf.gz) ; note: toDate, concerns strelka2 only \\
--tumor-sname    	TUMOR SAMPLE NAME [Required]\\
--normal-sname   	NORMAL SAMPLE NAME [Required] \\
--contigs-file    	FILE_WITH_CONTIGS FORMATTED AS IT IS IN VCF HEADERS (Optional; depend on tool ; use the script 'convert_contig_list_from_bam_to_vcf_format.sh' located in utils directory to create the appropriate file ) [default is null ] \\
--do-not-normalize	disable normalization [default is enable] \\
--bam			BAM file to provide to generate intermediate contig file in case --contigs-file option is not provided but needed for current tool's vcf in process
--th-AR|--threshold-AR      AR value (float from 0.000001 to 1 ). Based on that value, the GT flag (genotype) will be assigned to 0/1 if below that threshold or 1/1 if equal or above that threshold [default value is 0.90 ] ;

"
}

function fexit(){
	kill -s TERM ${TOP_PID}
	exit -1
}

function check_ev(){
	if [[ $1 -ne 0 ]] ; then echo -e "ERROR: ${2} FAILED ;\nexit_value:${1} ; Aborting! " ; fexit ; fi
}

function checkDir(){
	local D=$1
	if [[ ! -e ${D} ]] ; then echo -e "DIR NOT FOUND << ${D} >>; Aborting!" ; fexit; fi ;
}

function checkFile(){
	local F=$1
	if [[ ! -e ${F} ]] ; then echo -e "FILE NOT FOUND << ${F} >>; Aborting!" ; fexit ; fi ;
}

function checkFileName(){
	local F=$1
	if [[ ${F} == "" ]] ;
	then
		echo -e "ERROR: Missing filename for one of the variable which is REQUIRED for that tool; Aborting."
		fexit
	fi
}

function checkIntegerOrFloat(){
    ## Check if a number is a Float or Integer ;
    local N=$1
    if [[ ! ${N} =~ ^[0-9]*$ && ! ${N} =~ ^[0-9]*\.?[0-9]*$  ]] ;
    then
        echo -e "ERROR: AR threshold MUST be a Float or an Integer value; Check your input ; Aborting."
        fexit
    fi
}

function checkIfRangeZeroOne(){
    local N=$1
    test_gtZero=$(bc -l <<< "${N} > 0.000")
    test_leOne=$(bc -l <<< "${N} <= 1.000")
    if [[ ( ${test_gtZero} -eq 1 && ${test_leOne} -eq 1 ) ]] ;
    then
        echo -e "AR interval ]0,1] check: PASSED"
    else
        echo  -e "ERROR: AR threshold MUST be in interval ]0,1] (0 exclude and 1 included) ; Check your input ; Aborting."
        fexit
    fi
}

function check_contig_file(){
	local TN="$1"
	local F="$2"
	echo -e "Checking if Contig Files has been provided for tool $( echo ${TN} | tr '[a-z]' '[A-Z]' ) ; provided contig file in check is: << ${F} >> or if not, BAM file should have been provided; if so, BAM_file == << ${F} >>. if not error is generated. please use options --contigs-file or --bam " 1>&2
	checkFileName "${F}" 1>&2
	checkFile "${F}" 1>&2
}

function init_some_vars(){
	LI="RECAP_INPUTS_USED:"
	TOOLNAME=""
	VCF_ALL_CALLS=""
	NORMAL_SNAME=""
	TUMOR_SNAME=""
	VCF_INDELS_FILE=""
	VCF_SNVS_FILE=""
	BAM_FILE=""
	NORMALIZE="yes"
	TARGETS_BED_FILE_FTT=""
	CONTIGS_FILE=""
	VCF_FINAL_USER_GIVEN_NAME=""
	TH_AR=""
}

function getOptions(){
# options may be followed by one colon to indicate they have a required argument
if ! options=`getopt -o hd:b:g:o: -l help,dir-work:,ref-genome:,tumor-sname:,normal-sname:,vcf-indels:,vcf-snvs:,vcf:,toolname:,prepped-vcf-outfilename:,bam:,contigs-file:,print-default-toolnames,do-not-normalize -- "$@" `
	then
	# something went wrong, getopt will put out an error message for us
		echo "ERROR in Arguments" ; usage
		fexit
	fi
	eval set -- "$options"
	while [[ $# -gt 0 ]]
	do
		# for options with required arguments, an additional shift is required
		case "$1" in
		-d|--dir-work) export DIR_WORK=$2 ; LI="${LI}\nDIR_WORK==\"${DIR_WORK}\"" ;  shift ;;
		--tumor-sname) export TUMOR_SNAME=$2 ; LI="${LI}\nTUMOR_SNAME==\"${TUMOR_SNAME}\"" ; shift ;; ## TUMOR SAMPLE NAME as represented in SM tag in BAM
		--normal-sname) export NORMAL_SNAME=$2 ; LI="${LI}\nNORMAL_SNAME==\"${NORMAL_SNAME}\"";  shift ;; ## NORMAL SAMPLE NAME as represented in SM tag in BAM
		--vcf) export VCF_ALL_CALLS="$2" ; LI="${LI}\nVCF_ALL_CALLS==\"${VCF_ALL_CALLS}\"";  shift ;;
		--vcf-indels)	export VCF_INDELS_FILE="$2"  ; LI="${LI}\nVCF_INDELS_FILE==\"${VCF_INDELS_FILE}\"";  shift ;;
		--vcf-snvs) export VCF_SNVS_FILE="$2"  ; LI="${LI}\nVCF_SNVS_FILE==\"${VCF_SNVS_FILE}\"";  shift ;;
		--bam) export BAM_FILE=$2 ; LI="${LI}\nBAM_FILE==\"${BAM_FILE}\"";  shift ;;
		-g|--ref-genome) export REF_GENOME_FASTA="${2}" ;  LI="${LI}\nREF_GENOME_FASTA==\"${REF_GENOME_FASTA}\"";  shift ;; ## FASTA FILE ; .fai file is needed
		--toolname) export TOOLNAME=$2 ; LI="${LI}\nTOOLNAME==\"${TOOLNAME}\"";  shift ;;
		--do-not-normalize) export NORMALIZE="no" ; LI="${LI}\nNORMALIZE==\"${NORMALIZE}\"" ;;
		--contigs-file) export CONTIGS_FILE="$2" ; LI="${LI}\nCONTIGS_FILE==\"${CONTIGS_FILE}\"";  shift ;; ## File containing the contigs in the same format as expected within a VCF file
		--print-default-toolnames) echo ${VALID_TOOLNAMES} ; exit ;;
		-o|--prepped-vcf-outfilename) export VCF_FINAL_USER_GIVEN_NAME="$2" ; LI="${LI}\nVCF_FINAL_USER_GIVEN_NAME==\"${VCF_FINAL_USER_GIVEN_NAME}\"";  shift ;;
		-h|--help) usage ; exit ;;
		(--) shift ;;
		(-*) echo -e "$0: error - unrecognized option $1\n\n" 1>&2   ; usage;  exit 1  ;;
		(*) break ; echo "$0: error --- unrecognized option $1" 1>&2 ; usage;  exit 1  ;;
		esac
		shift
	done
}

function recap_input(){
	#input recap
	LI="${LI}\nCURR_DIR==\"${PWD}\""
	echo -e "\n\n+------------------------------------------------+\n${LI[@]}\n+------------------------------------------------+\n\n"
}

function concatenate_snvs_indels(){
	local TOOLNAME=$1
	local SNVS=$2
	local INDELS=$3

	echo -e "## Concatenating ${TOOLNAME}'s somatic files ..." 1>&2
	fout_name=${TOOLNAME}.somatic.concat.vcf ## name of the concatenated vcf with snvs and indels from ${TOOLNAME}
	( zcat -f ${VCF_SNVS_FILE} | grep -E "^##" ; zcat -f ${VCF_INDELS_FILE} | grep -E "^##INFO|^##FORMAT" ; zcat -f ${VCF_INDELS_FILE} | grep -E "^#CHROM" ) > header_for_concat_${TOOLNAME}_vcf.txt
	( cat header_for_concat_${TOOLNAME}_vcf.txt ; zcat -f ${VCF_SNVS_FILE} | grep -vE "^#" ; zcat -f ${VCF_INDELS_FILE} | grep -vE "^#"   ) > ${fout_name}
	check_ev $? "concat vcfs indels and snvs" 1>&2
	rm header_for_concat_${TOOLNAME}_vcf.txt
	VCF=${fout_name}
	echo "${VCF}"
}

function check_and_update_sample_names(){
	##@@@@@@@@@@@@@@@@@@@@@@@
	## CHECK SAMPLES NAMES
	##@@@@@@@@@@@@@@@@@@@@@@@
	## so far the sample name should be NORMAL, TUMOR or the SM tag from BAM file ;
	## if not ##TODO revise the way we handle the sample names to make it much more generic;
	## will need discussion to do so, and will need to capture all the use cases possible ;

	local VCF=$1
	echo " in ${funcname} :  ${VCF} ${TOOLNAME} ${NORMAL_SNAME} ${TUMOR_SNAME} " 1>&2

	echo -e "## Checking the Sample names columns and swapping them if necessary (we assume that the VCF is a SOMATIC calls vcf )" 1>&2
	COL10_VALUE=`grep -m1 -E "^#CHROM" ${VCF} | cut -f10`
	COL11_VALUE=`grep -m1 -E "^#CHROM" ${VCF} | cut -f11`

	if [[ "${COL10_VALUE}" == "NORMAL" && "${COL11_VALUE}" == "TUMOR"  ]] ;
	then
		echo -e "\tsample name in column 10 is  NORMAL and column 11 is named TUMOR" 1>&2
		echo -e "\twe are updating the sample names appropriately here with the ones given by the user" 1>&2
		sed -i "/^#CHROM/ s/NORMAL/${NORMAL_SNAME}/ ; /^#CHROM/ s/TUMOR/${TUMOR_SNAME}/" ${VCF}

	elif [[ "${COL11_VALUE}" == "NORMAL" && "${COL10_VALUE}" == "TUMOR"  ]] ;
	then
		echo -e "\tsample name in column 10 is  TUMOR and column 11 is named NORMAL" 1>&2
		echo -e "\tas we decided that column 10 should be NORMAL sample and column 11 the TUMOR one" 1>&2
		echo -e "\twe RENAME and SWAPPED the sample names." 1>&2
		cat ${VCF} | awk -v TUMORSNAME=${TUMOR_SNAME} -v NORMALSNAME=${NORMAL_SNAME} -F"\t" '{OFS="\t" ; if($1~/^##/){print ; continue} ; if($1~/^CHROM/){ sub("TUMOR",TUMORSNAME,$10) ; sub("NORMAL",NORMALSNAME,$11) ;tempCol=$10 ; $10=$11; $11=tempCol ; print }  }' > temp_${TOOLNAME}.renamed_swapped_samples.vcf
		check_ev $? "swap column 10 and 11"  1>&2
		mv temp_${TOOLNAME}.renamed_swapped_samples.vcf ${VCF}

	elif [[ ( "${COL11_VALUE}" == "${TUMOR_SNAME}" && "${COL10_VALUE}" == "${NORMAL_SNAME}" ) || ( "${COL10_VALUE}" == "${TUMOR_SNAME}" && "${COL11_VALUE}" == "${NORMAL_SNAME}" )   ]] ;
	then
		echo -e "At least the sample name in column 10 is either Tumor or Normal; As we want normal sample in Column 10, we check if it is the case ..." 1>&2
		if [[ "${COL10_VALUE}" == "${TUMOR_SNAME}" ]]
		then
			echo "## we swap column 10 and 11 as we found that sample name in column 10 is the Tumor sample" 1>&2
			(grep -E "^##" ${VCF} ; grep -vE "^##" ${VCF} | awk -F"\t" '{OFS="\t" ; TEMP=$10 ; $10=$11 ; $11=TEMP ; print }' ) > temp_swap_samples_column.${TOOLNAME}.vcf
			if [[ $? -ne 0 ]] ; then echo -e "ERROR in swapping column; please check logs and your inputs to understand why it failed; Aborting LANCET post-processing; " ; fexit ; fi
			mv temp_swap_samples_column.${TOOLNAME}.vcf ${VCF}
		else
			echo -e "## We found that sample in column 10 is the NORMAL sample; we do not swap the columns" 1>&2
		fi
	else
		echo -e "ERROR: ${TOOLNAME}'s Sample name in VCF do NOT match 'NORMAL' or 'TUMOR' names OR any expected names already present in the VCF, sample name that was normally captured from SM tag in the BAM file;
		\nplease check your inputs; Aborting!" 1>&2
		fexit
	fi
	## return value
	echo "${VCF}"
}

function get_contigs_file(){
	local TN=$1
	local BAMF="$2"
	local CONTIGSF="$3"

	if [[ ${BAMF} == "" ]] ;
	then
		check_contig_file "${TN}" "${CONTIGSF}"
	elif [[ ${BAMF} != "" ]] ;
	then
		if [[ ! -e "${BAMF}" ]]
		then
			if [[ "${CONTIGSF}" == "" || ! -e "${CONTIGSF}"  ]]
			then
				echo -e "ERROR: at least one of the option --bam or --contigs-file has to be provided; Check your inputs; Aborting" ;
				fexit
			else
				export CONTIGS_FILE="${CONTIGSF}"
			fi
		else
			CONTIGSF=$( basename ${BAMF}).contigs.txt
			bash $(dirname $0 )/utils/convert_contig_list_from_bam_to_vcf_format.sh -i "${BAMF}" -o "${CONTIGSF}"
			check_ev $? "convert_contig_list_from_bam_to_vcf_format "
			export CONTIGS_FILE="${CONTIGSF}"
		fi
	fi

}

function add_Contigs(){
	local VCF=$1
	## Add missing Contigs lines to VCF Header
	fout_name=${VCF%.*}.contigs.vcf
	echo -e "## adding contigs to ${TOOLNAME}'s vcf header ..." 1>&2
	checkFile "${CONTIGS_FILE}" 1>&2
	( grep -E "^##" ${VCF} ; cat "${CONTIGS_FILE}" ; grep -m 1 -E "#CHROM" ${VCF} ; grep -vE "^#" ${VCF} ) > ${fout_name}
	check_ev $? "addContigs " 1>&2
	VCF=${fout_name}
	echo "${VCF}"
}

function look_for_block_substitution_in_octopus(){
	## Octopus Specific function ; collapse variants if they form a block-substitution
	local VCF=$1
	echo -e "## look for block substitution in ${TOOLNAME}'s vcf ... and collapse the variants "  1>&2
	fout_name=${VCF%.*}.blocs.vcf
	mycmd="python ${PYTHON_SCRIPT_OCTOPUS_BLOCK_SUBSTITUTION} --normal_column 10 --tumor_column 11 -i ${VCF} -o ${fout_name}"
	echo ${mycmd} 1>&2 ;
	eval ${mycmd} 1>&2 ;
	check_ev $? " python look_for_bloc_substitution_in_${TOOLNAME} " 2>&1
	VCF=${fout_name}
	echo ${VCF}
}

function decompose(){
	## decompose the vcf to get one ALT per line
	local VCF=$1
	echo -e "## vt decompose for somatic calls of ${TOOLNAME}'s vcf ..."  1>&2
	fout_name=${VCF%.*}.decomp.vcf
	vt decompose -s ${VCF} -o ${fout_name}
	check_ev $? "vt decompose " 2>&1
	VCF=${fout_name}
	echo "${VCF}"
}

function make_vcf_upto_specs_for_VcfMerger(){
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## updating ${TOOLNAME} VCF to specs for vcfMerger2
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	local VCF=$1
	echo -e "## prep vcf for vcfMerger2 ..." 1>&2
	fout_name=${VCF%.*}.prep.vcf
	mycmd="python ${PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER} -i ${VCF} --normal_column 10 --tumor_column 11 --outfilename ${fout_name}"
	if [[ ${TH_AR} != "" ]] ;
	then
	    mycmd="${mycmd} --threshold_AR ${TH_AR}"
	fi
	echo ${mycmd} 1>&2 ;
	eval ${mycmd} 1>&2 ;
	check_ev $? "${PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER} " 1>&2
	VCF=${fout_name}
	echo "${VCF}"
}

function normalize_vcf(){
	## normalize the VCF using bcftools
	local VCF=$1
	echo "in ${funcname}:  ${VCF} and refgenome is ${REF_GENOME_FASTA}" 1>&2
	echo "## Normalizing vcf ..."  1>&2
	fout_name=${VCF%.*}.norm.vcf
	echo -e "## bcftools sort -O v ${VCF} | bcftools norm -c x -f ${REF_GENOME_FASTA} -O v - > ${fout_name}" 1>&2
	bcftools sort -O v ${VCF} | bcftools norm -c x -f ${REF_GENOME_FASTA} -O v - > ${fout_name}
	check_ev $? "bcftools with ${VCF} " 1>&2
	VCF=${fout_name}
	echo "${VCF}"
}

function final_msg(){
	local VCF=$1
	if [[ ${VCF_FINAL_USER_GIVEN_NAME} != "" ]] ;
	then
		VCF_FINAL=${VCF_FINAL_USER_GIVEN_NAME}
	else
		VCF_FINAL=${TOOLNAME}.somatic.uts.vcf ; ## uts stands for up-to-specs for vcfMerger2
	fi
	cp ${VCF} ${VCF_FINAL}
	echo -e "\n##---------------------------------------------------------------------------##" 1>&2
	echo -e "Extension should represent the steps performed on the given input VCF" 1>&2
	echo -e "vcfMerger-compatible vcf file : << ${VCF} >>" 1>&2
	echo -e "copy/renamed to: ${VCF_FINAL} " 1>&2
	exit
}

function process_strelka2_vcf(){
	local VCF=$1
	VCF=$( check_and_update_sample_names ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}

function process_mutect2_vcf(){
	local VCF=$1
	VCF=$( check_and_update_sample_names ${VCF} )
	VCF=$( decompose ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}

function process_lancet_vcf(){
	local VCF=$1
	VCF=$( check_and_update_sample_names ${VCF} )
	VCF=$( add_Contigs ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}

function process_octopus_vcf(){
	local VCF=$1
	VCF=$( check_and_update_sample_names ${VCF} )
	VCF=$( add_Contigs ${VCF} )
	VCF=$( look_for_block_substitution_in_octopus ${VCF}) ## why do we put look for blocs before decompose? b/c we only use the first allele for collapsing block
	VCF=$( decompose ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}

function process_vardictjava_vcf(){
	local VCF=$1
	VCF=$( check_and_update_sample_names ${VCF} )
	VCF=$( add_Contigs ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}



function run_tool(){
	local TOOLNAME=$( echo $1 | tr '[A-Z]' '[a-z]')
	local VCF=$2
	local DIR_PATH_TO_SCRIPTS="${DIR_PATH_TO_SCRIPTS}"

	case $TOOLNAME in
		strelka2|slk)
			PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/strelka2/strelka2.somatic.addFieldsForVcfMerger.py"
			process_strelka2_vcf ${VCF}
			;;
		mutect2|mtc)
			PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/mutect2/mutect2.somatic.addFieldsForVcfMerger.py"
			process_mutect2_vcf ${VCF}
			;;
		lancet|lct|lan)
			PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/lancet/lancet.somatic.addFieldsForVcfMerger.py"
			get_contigs_file "${TOOLNAME}" "${BAM_FILE}" "${CONTIGS_FILE}"
			process_lancet_vcf ${VCF}
			;;
		octopus|oct)
			PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/octopus/octopus.somatic.addFieldsForVcfMerger.py"
			PYTHON_SCRIPT_OCTOPUS_BLOCK_SUBSTITUTION="${DIR_PATH_TO_SCRIPTS}/octopus/octopus.somatic.checkForBlockSubstitutionVariants.py"
			get_contigs_file "${TOOLNAME}" "${BAM_FILE}" "${CONTIGS_FILE}"
			process_octopus_vcf ${VCF}
			;;
		vardict|vdt|vdj)
			PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/vardictjava/vardictjava.somatic.addFieldsForVcfMerger.py"
			get_contigs_file "${TOOLNAME}" "${BAM_FILE}" "${CONTIGS_FILE}"
			process_vardictjava_vcf ${VCF}
			;;
		(*) echo -e "\nERROR: unrecognized toolname:  $1\nERROR: valid toolnames are << ${VALID_TOOLNAMES} >>" 1>&2 ; fexit  ;;
	esac
	echo "${SCRIPT_FP}" ## we return the full path of the script
}

function main(){

	##@@@@@@@@@@@@@@@@@##
	##  check inputs   ##
	##@@@@@@@@@@@@@@@@@##
	echo -e "## Checking inputs ..."
	if [[ ${TOOLNAME} == ""  ]] ; then echo -e "ERROR: --toolname has to be provided ; Aborting." ; fexit ; fi
	if [[ ${REF_GENOME_FASTA} == ""  ]] ; then echo -e "ERROR: reference genome option is required ; here we have a missing value; provide --ref-genome ; Aborting." ; fexit ; fi
	if [[ ${TUMOR_SNAME} == ""  ]] ; then echo -e "ERROR: TUMOR Sample Name MUST be provided ; Missing Values Found ; Check your inputs;  Aborting." ; fexit ; fi
	if [[ ${NORMAL_SNAME} == ""  ]] ; then echo -e "ERROR: NORMAL Sample Name MUST be provided ; Missing Values Found ; Check your inputs;  Aborting." ; fexit ; fi


	## check files and folders if exist
	checkDir ${DIR_WORK}
	checkFile ${REF_GENOME_FASTA}

	cd ${DIR_WORK}

	## check if we deal with indels and snvs in separated vcf or in all-in-one vcf
	if [[ ${VCF_ALL_CALLS} != "" ]] ;
	then
		checkFile ${VCF_ALL_CALLS}
		echo -e "CREATING SYMLINK in CURR DIR ${PWD}"
		ln -sf ${VCF_ALL_CALLS} $(basename ${VCF_ALL_CALLS}) &>/dev/null ## we create a symlimk in current working directory (in case original vcf folder is not writable)
		VCF=$(basename ${VCF_ALL_CALLS}) ## make basename vcf the new VCF name
		echo "processing vcf:  ${VCF}" 1>&2
		run_tool ${TOOLNAME} ${VCF}
	elif [[ ( ${VCF_SNVS_FILE} != "" && ${VCF_INDELS_FILE} != "" ) &&  ( -e ${VCF_SNVS_FILE} && -e ${VCF_INDELS_FILE} )  ]]
	then
		checkFile ${VCF_SNVS_FILE} ; checkFile ${VCF_INDELS_FILE}
		VCF=$( concatenate_snvs_indels ${TOOLNAME} ${VCF_SNVS_FILE} ${VCF_INDELS_FILE} )
		run_tool ${TOOLNAME} ${VCF}
	else
		echo -e "ERROR: Check your inputs ; VCF files information is missing or erroneous; Aborting!" ; 1>&2
		fexit
	fi

}



###@@@@@@@@@@@@@@
### START HERE
###@@@@@@@@@@@@@@
#if [[ `uname -s` == "Darwin" ]] ;
#then
#    type python3 >/dev/null 2>&1 || { echo >&2 "Require \"python3\" executable in MacOSX bash but it's not in the
#    PATH.  Aborting.";
# exit 1; } || python3 -V
#else
#    type python >/dev/null 2>&1 || { echo >&2 "Require \"python\" executable but it's not in the PATH.  Aborting.";
# exit 1; } || python -V
#    python_main_version_number=`python -V 2>&1 | sed 's/Python //g' | cut -d"." -f1 `
#    echo "python main version number captured: ${python_main_version_number}"
#    if [[ ! "${python_main_version_number}" == "3" ]] ; then echo -e "ERROR: Python 3 or up Expected in PATH; Aborting " ;
# exit 1 ; fi
#fi

type python >/dev/null 2>&1 || { echo >&2 "Require \"python\" executable but it's not in the PATH.  Aborting."; exit
1; } || python -V
python_main_version_number=`python -V 2>&1 | sed 's/Python //g' | cut -d"." -f1 `
if [[ ! "${python_main_version_number}" == "3" ]] ; then echo -e "ERROR: Python 3 or up Expected in PATH; Aborting " ; exit 1
fi
type vt >/dev/null 2>&1 || { echo >&2 "Require \"vt\" executable but it's not in the PATH.  Aborting."; exit 1; } ||
vt --version
type bcftools >/dev/null 2>&1 || { echo >&2 "Require \"bcftools\" executable but it's not in the PATH.  Aborting."; exit 1; }
if [[ $( echo "`bcftools --version-only  2>&1 | sed 's/+.*//g'` <  1.7 " | bc -l ) -eq 1  ]] ; then echo -e "ERROR: bcftools 1.7 or up Expected in PATH; Aborting " ; exit 1 ; fi

## init variables
init_some_vars
## get options
getOptions $@
## recap inputs
recap_input
## now we have all our inputs set up, let's process data ...
main
check_ev $? "main function "
exit