#!/usr/bin/env bash

function usage(){
	echo -e "\nUSAGE:"
	echo -e "`basename $0` \\
-d|--dir-work    	DIR_WORK directory where outputs will be written (needs to exist))  \\
-g|--ref-genome   	REFERENCE GENOME FASTA FILE  [Required] \\
-t|--toolname		Provide the toolname associated to the input vcf [REQUIRED]; see valid toolnames in prep_vcf_defaults.ini file, or use --list-valid-toolnames in in-line command \\
-o|--prepped-vcf-outfilename	Provide the name for the uptospecs vcf file that will be use as input for the vcfMerger2.0 tool \\
--make-bed-for-venn   enable making BED file for the Intervene python tool [default is disable] \\
--print-default-toolnames	Print default valid toolnames accepted so far (case insensitive) and exit \\
--vcf			vcf having all types of variants already (no need to concatenate) \\
--vcf-indels		tool's VCF with indels (.vcf or .vcf.gz) ; note: toDate, concerns strelka2 only \\
--vcf-snvs		tool's VCF with snvs (.vcf or .vcf.gz) ; note: toDate, concerns strelka2 only \\
--tumor-sname    	TUMOR SAMPLE NAME [Required]\\
--normal-sname   	NORMAL SAMPLE NAME [Required] \\
--contigs-file    	FILE_WITH_CONTIGS FORMATTED AS IT IS IN VCF HEADERS (Optional; depend on tool ; use the script 'convert_contig_list_from_bam_to_vcf_format.sh' located in utils directory to create the appropriate file ) [default is null ] \\
--do-not-normalize	disable normalization [default is enable] \\
--bam			BAM file to provide to generate intermediate contig file in case --contigs-file option is not provided but needed for current tool's vcf in process
--th-AR|--threshold-AR      AR value (float from 0.000001 to 1 ). Based on that value, the GT flag (genotype) will be assigned to 0/1 if below that threshold or 1/1 if equal or above that threshold [default value is 0.90 ] ;


WARNING: This script does NOT work on MacOS; Sorry. This is due to getOpt incompatibiity with macOS's unix-like system
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
	if [[ ! -e ${D} ]] ; then echo -e "DIR NOT FOUND << ${D} >> ; curdir = ${PWD}; Aborting!" ; fexit; fi ;
}

function checkFile(){
	local F=$1
	if [[ ! -e ${F} ]] ; then echo -e "FILE NOT FOUND << ${F} >> ; curdir = ${PWD}; Aborting!" ; fexit ; fi ;
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
	MAKE_BED_FOR_VENN="no"
	DELETE_TEMPS=0 ; ## 0 means keep-temps; 1 means delete temp files
}

function getOptions(){
# options may be followed by one colon to indicate they have a required argument
# NOTE: long option are not working on MacOS Sierra
if ! options=`getopt -o hd:b:g:o:t: -l help,dir-work:,ref-genome:,tumor-sname:,normal-sname:,vcf-indels:,vcf-snvs:,vcf:,toolname:,prepped-vcf-outfilename:,bam:,contigs-file:,print-valid-toolnames,do-not-normalize,make-bed-for-venn,delete-temps -- "$@" `
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
		-d|--dir-work) export DIR_OUTPUT=$( readlink -f $2) ; LI="${LI}\nDIR_WORK==\"${DIR_OUTPUT}\"" ;  shift ;;
		--tumor-sname) export TUMOR_SNAME=$2 ; LI="${LI}\nTUMOR_SNAME==\"${TUMOR_SNAME}\"" ; shift ;; ## TUMOR SAMPLE NAME as represented in SM tag in BAM
		--normal-sname) export NORMAL_SNAME=$2 ; LI="${LI}\nNORMAL_SNAME==\"${NORMAL_SNAME}\"";  shift ;; ## NORMAL SAMPLE NAME as represented in SM tag in BAM
		--vcf) export VCF_ALL_CALLS="$2" ; LI="${LI}\nVCF_ALL_CALLS==\"${VCF_ALL_CALLS}\"";  shift ;;
		--vcf-indels)	export VCF_INDELS_FILE="$2"  ; LI="${LI}\nVCF_INDELS_FILE==\"${VCF_INDELS_FILE}\"";  shift ;;
		--vcf-snvs) export VCF_SNVS_FILE="$2"  ; LI="${LI}\nVCF_SNVS_FILE==\"${VCF_SNVS_FILE}\"";  shift ;;
		--bam) export BAM_FILE=$2 ; LI="${LI}\nBAM_FILE==\"${BAM_FILE}\"";  shift ;;
		-g|--ref-genome) export REF_GENOME_FASTA="${2}" ;  LI="${LI}\nREF_GENOME_FASTA==\"${REF_GENOME_FASTA}\"";  shift ;; ## FASTA FILE ; .fai file is needed
		-t|--toolname) export TOOLNAME=$2 ; LI="${LI}\nTOOLNAME==\"${TOOLNAME}\"";  shift ;;
		--do-not-normalize) export NORMALIZE="no" ; LI="${LI}\nNORMALIZE==\"${NORMALIZE}\"" ;;
		--contigs-file) export CONTIGS_FILE="$2" ; LI="${LI}\nCONTIGS_FILE==\"${CONTIGS_FILE}\"";  shift ;; ## File containing the contigs in the same format as expected within a VCF file
		--print-valid-toolnames) echo ${VALID_TOOLNAMES} ; exit ;; ## print possible toolnames to be used with the vcfMerger prep module
		--delete-temps) export DELETE_TEMPS=1 ; LI="${LI}\nDELETE_TEMPS==\"${DELETE_TEMPS}\"" ;;
		-o|--prepped-vcf-outfilename) export VCF_FINAL_USER_GIVEN_NAME="$2" ; LI="${LI}\nVCF_FINAL_USER_GIVEN_NAME==\"${VCF_FINAL_USER_GIVEN_NAME}\"";  shift ;;
		--make-bed-for-venn) export MAKE_BED_FOR_VENN="yes" ; LI="${LI}\nMAKE_BED_FOR_VENN==\"${MAKE_BED_FOR_VENN}\"" ;;
		-h|--help) usage ; exit ;;
		(--) shift ;;
		(-*) echo -e "$0: error - unrecognized option $1\n\n" 1>&2   ; usage;  exit 1  ;;
		(*) break ; echo "$0: error --- unrecognized option $1" 1>&2 ; usage;  exit 1  ;;
		esac
		shift
	done
}

function delete_temporary_files(){
    local delete_temps=$1
#    local PATTERN="sname.vcf$|sname.prep.vcf$|sname.prep.norm.vcf$|sname.decomp.vcf$|sname.decomp.prep.vcf$|sname.decomp.prep.norm.vcf"
    local PATTERN="*sname.*.vcf"
    if [[ ${delete_temp} == 1 ]]
    then
        echo -e "file with pattern \"${PATTERN}\" will be deleted ... "
        rm $( find ${DIR_OUTPUT} -type f -name "${PATTERN}")
        if [[ $? -ne 0 ]]
        then
            echo "ERROR in deleting files with pattern ${PATTERN}"
        fi
    fi
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
	local VCF_OUT=$(basename ${VCF} ".vcf").sname.vcf
	echo " in ${funcname} :  ${VCF} ${TOOLNAME} ${NORMAL_SNAME} ${TUMOR_SNAME} " 1>&2

	echo -e "## Checking the Sample names columns and swapping them if necessary (we assume that the VCF is a SOMATIC calls vcf )" 1>&2
	COL10_VALUE=`grep -m1 -E "^#CHROM" ${VCF} | cut -f10`
	COL11_VALUE=`grep -m1 -E "^#CHROM" ${VCF} | cut -f11`

	if [[ ( "${COL10_VALUE}" == "NORMAL" && "${COL11_VALUE}" == "TUMOR" ) ]] ;
	then
		echo -e "\tsample name in column 10 is  NORMAL and column 11 is named TUMOR" 1>&2
		echo -e "\twe are updating the sample names appropriately here with the ones given by the user" 1>&2
		cat  ${VCF}| sed "/^#CHROM/ s/NORMAL/${NORMAL_SNAME}/ ; /^#CHROM/ s/TUMOR/${TUMOR_SNAME}/" > temp_sname_${VCF}
		mv temp_sname_${VCF} ${VCF_OUT}

	elif [[ "${COL11_VALUE}" == "NORMAL" && "${COL10_VALUE}" == "TUMOR"  ]] ;
	then
		echo -e "\tsample name in column 10 is  TUMOR and column 11 is named NORMAL" 1>&2
		echo -e "\tas we decided that column 10 should be NORMAL sample and column 11 the TUMOR one" 1>&2
		echo -e "\twe RENAME and SWAPPED the sample names." 1>&2
		cat ${VCF} | awk -v TUMORSNAME=${TUMOR_SNAME} -v NORMALSNAME=${NORMAL_SNAME} -F"\t" '{OFS="\t" ; if($1~/^##/){print ; next} else if ($1~/^#CHROM/){ sub("TUMOR",TUMORSNAME,$10) ; sub("NORMAL",NORMALSNAME,$11) ;tempCol=$10 ; $10=$11; $11=tempCol ; print } else { print } }' > temp_${TOOLNAME}.renamed_swapped_samples.vcf
		check_ev $? "swap column 10 and 11"  1>&2
		mv temp_${TOOLNAME}.renamed_swapped_samples.vcf ${VCF_OUT}

	elif [[ ( "${COL11_VALUE}" == "${TUMOR_SNAME}" && "${COL10_VALUE}" == "${NORMAL_SNAME}" ) || ( "${COL10_VALUE}" == "${TUMOR_SNAME}" && "${COL11_VALUE}" == "${NORMAL_SNAME}" )   ]] ;
	then
		echo -e "At least the sample name in column 10 is either Tumor or Normal; As we want normal sample in Column 10, we check if it is the case ..." 1>&2
		if [[ "${COL10_VALUE}" == "${TUMOR_SNAME}" ]]
		then
			echo "## we swap column 10 and 11 as we found that sample name in column 10 is the Tumor sample" 1>&2
			(grep -E "^##" ${VCF} ; grep -vE "^##" ${VCF} | awk -F"\t" '{OFS="\t" ; TEMP=$10 ; $10=$11 ; $11=TEMP ; print }' ) > temp_swap_samples_column.${TOOLNAME}.vcf
			if [[ $? -ne 0 ]] ; then echo -e "ERROR in swapping column; please check logs and your inputs to understand why it failed; Aborting LANCET post-processing; " ; fexit ; fi
			mv temp_swap_samples_column.${TOOLNAME}.vcf ${VCF_OUT}
		else
			echo -e "## We found that sample in column 10 is the NORMAL sample; we do not swap the columns" 1>&2
			cp ${VCF} ${VCF_OUT}
		fi
	else
		echo -e "ERROR: ${TOOLNAME}'s Sample name in VCF do NOT match 'NORMAL' or 'TUMOR' names OR any expected names already present in the VCF, sample name that was normally captured from SM tag in the BAM file;
		\nplease check your inputs; Aborting!" 1>&2
		fexit
	fi
	## return value which is the vcf filename
	echo "${VCF_OUT}"
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
			    echo -e "WARNING WARNING: neither --bam or --contigs-file got provided; We DO NOT garuantee that the contigs are present in LANCET vcfs in correct way; That might be the cause of Errors later on in the vcfMerge2 tool; Check that the First LANCET vcf has the correct expected contig names; WARNING WARNING" ;
				#echo -e "ERROR: at least one of the option --bam or --contigs-file has to be provided; Check your inputs; Aborting" ;
				#fexit
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

function modify_contig_info_in_lancet(){
	local VCF=$1
	## Compared with other tools, lancet does add the assembly name to the contig definition
	## here we remove the extra assembly info in order to get the same contig lines accross the tools
	fout_name=${VCF%.*}.contigs.vcf
	echo -e "## adding contigs to ${TOOLNAME}'s vcf header ..." 1>&2
	sed '/^##contig/s/,assembly.*>$/>/' ${VCF} > ${fout_name}
	check_ev $? "addContigs " 1>&2
	VCF=${fout_name}
	echo "${VCF}"
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
	mycmd="python -W ignore ${PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER} -i ${VCF} --normal_column 10 --tumor_column 11 --outfilename ${fout_name}"
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

function prepare_input_file_for_Venn(){
    local VCF="$1"
    local DIROUT=$2
    if [[ ! -e ${VCF} ]] ; then echo -e "ERROR: VCF NOT FOUND --> ${VCF}" ; fi
    local INPUT_FILE_FOR_VENN=$( basename ${VCF} ".vcf" ).intervene.bed
    #echo -e "grep -vE "^#" ${VCF} | awk -F"\t" '{OFS="_" ; print $1,$2,$4,$5 }' > ${INPUT_FILE_FOR_VENN} " 1>&2
    grep -vE "^#" ${VCF} | awk -F"\t" '{OFS="\t" ; print $1,$2,$2,$4,$5 }' | sort -k1,1V -k2,2n -k3,3n > ${DIROUT}/${INPUT_FILE_FOR_VENN}
    check_ev $? "prepare_input_file_for_Venn"
}


function prepare_input_file_for_Venn_SplitbyVariantType(){
    local VCF="$1"
    local DIROUT=$2
    if [[ ! -e ${VCF} ]] ; then echo -e "ERROR: VCF NOT FOUND --> ${VCF}" ; fi
    local INPUT_FILE_FOR_VENN_SNVS=$( basename ${VCF} ".vcf" ).intervene.snvs.bed
    local INPUT_FILE_FOR_VENN_INDELS=$( basename ${VCF} ".vcf" ).intervene.indels.bed
    #echo -e "grep -vE "^#" ${VCF} | awk -F"\t" '{OFS="_" ; print $1,$2,$4,$5 }' > ${INPUT_FILE_FOR_VENN} " 1>&2
    grep -vE "^#" ${VCF} | awk -F"\t" ' $4~/[ATCG]/ && $5~/[ATCG]/ && length($4)==1 && length($5)==1 || ( length($4)>1 && length($4)==length($5) ) {OFS="\t" ; print $1,$2,$2,$4,$5 }' | sort -k1,1V -k2,2n -k3,3n > ${DIROUT}/${INPUT_FILE_FOR_VENN_SNVS}
    check_ev $? "prepare_input_file_for_Venn_SNVS"
    grep -vE "^#" ${VCF} | awk -F"\t" ' length($4)>length($5) || length($5)>length($4) || $4=="." || $5=="."  {OFS="\t" ; print $1,$2,$2,$4,$5 }' | sort -k1,1V -k2,2n -k3,3n > ${DIROUT}/${INPUT_FILE_FOR_VENN_INDELS}
    check_ev $? "prepare_input_file_for_Venn_INDELS"
}

function final_msg(){
	local VCF=$1
	if [[ ${VCF_FINAL_USER_GIVEN_NAME} != "" ]] ;
	then
		VCF_FINAL=${DIR_OUPUT}/${VCF_FINAL_USER_GIVEN_NAME}
	else
		VCF_FINAL=${DIR_OUTPUT}/${TOOLNAME}.somatic.uts.vcf ; ## uts stands for up-to-specs for vcfMerger2
	fi
	cp ${VCF} ${VCF_FINAL}
	if [[ ${MAKE_BED_FOR_VENN} == "yes" ]]
	then
	    echo -e "preparing input file for intervene python module to make Venns" 1>&2
	    prepare_input_file_for_Venn "${VCF_FINAL}" "${DIR_OUTPUT}"
	fi

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
	VCF=$( modify_contig_info_in_lancet ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}

function process_octopus_vcf(){
	local VCF=$1
	VCF=$( check_and_update_sample_names ${VCF} )
	VCF=$( look_for_block_substitution_in_octopus ${VCF}) ## why do we put look for blocs before decompose? b/c we only use the first allele for collapsing block
	VCF=$( decompose ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}

function process_vardictjava_vcf(){
	local VCF=$1
	VCF=$( check_and_update_sample_names ${VCF} )
	#VCF=$( add_Contigs ${VCF} )
	VCF=$( make_vcf_upto_specs_for_VcfMerger ${VCF} )
	VCF=$( normalize_vcf ${VCF})
	final_msg ${VCF}
}

function run_tool(){
	local TOOLNAME=$( echo $1 | tr '[A-Z]' '[a-z]')
	local VCF=$2
	local DIR_PATH_TO_SCRIPTS="${DIR_PATH_TO_SCRIPTS}"

	case $TOOLNAME in
		strelka2|slk|strelka)
			PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/strelka2/strelka2.somatic.addFieldsForVcfMerger.py"
			process_strelka2_vcf ${VCF}
			;;
		mutect2|mtc|mutect)
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
			#get_contigs_file "${TOOLNAME}" "${BAM_FILE}" "${CONTIGS_FILE}"
			process_octopus_vcf ${VCF}
			;;
		vardict|vdt|vdj)
			PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/vardictjava/vardictjava.somatic.addFieldsForVcfMerger.py"
			#get_contigs_file "${TOOLNAME}" "${BAM_FILE}" "${CONTIGS_FILE}"
			process_vardictjava_vcf ${VCF}
			;;
		haplotypecaller|hc)
		    PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/haplotypecaller/haplotypecaller.germline.1s.addFieldsForVcfMerger.py"
		    process_haplotypecaller_vcf ${VCF}
		    ;;
        freebayes|fby|fbs|fb)
            PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/freebayes/freebayes.germline.1s.addFieldsForVcfMerger.py"
            process_freebayes_vcf ${VCF}
            ;;
        samtools|mpileup|st|mpl|mpp|mpup)
            PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/samtools/samtools.germline.1s.addFieldsForVcfMerger.py"
            process_samtools_mpileup_vcf ${VCF}
            ;;
        deepvariant|dpv)
            PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/deepvariant/deepvariant.germline.1s.addFieldsForVcfMerger.py"
            process_deepvariant_vcf ${VCF}
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
	checkDir ${DIR_OUTPUT}
	checkFile ${REF_GENOME_FASTA}


	## check if we deal with indels and snvs in separated vcf or in all-in-one vcf
	if [[ ${VCF_ALL_CALLS} != "" ]] ;
	then
		checkFile ${VCF_ALL_CALLS}
		VCF_ALL_CALLS=$(readlink -f ${VCF_ALL_CALLS})
		echo "FULL PATH to current VCF is: ${VCF_ALL_CALLS}"
		cd ${DIR_OUTPUT}
		if [[ ! -e $( basename ${VCF_ALL_CALLS}) ]]
		then
		    echo -e "CREATING SYMLINK in CURR DIR ${PWD}"
		    ln -sf ${VCF_ALL_CALLS} $(basename ${VCF_ALL_CALLS}) &>/dev/null ## we create a symlimk in current working directory (in case original vcf folder is not writable)
		fi
		VCF=$(basename ${VCF_ALL_CALLS}) ## make basename vcf the new VCF name
		echo "processing vcf:  ${VCF}" 1>&2
		run_tool ${TOOLNAME} ${VCF}
		delete_temporary_files ${DELETE_TEMPS}
	elif [[ ( ${VCF_SNVS_FILE} != "" && ${VCF_INDELS_FILE} != "" ) &&  ( -e ${VCF_SNVS_FILE} && -e ${VCF_INDELS_FILE} )  ]]
	then
		checkFile ${VCF_SNVS_FILE} ; checkFile ${VCF_INDELS_FILE}
		VCF=$( concatenate_snvs_indels ${TOOLNAME} ${VCF_SNVS_FILE} ${VCF_INDELS_FILE} )
		run_tool ${TOOLNAME} ${VCF} ${DIR_OUTPUT}
		delete_temporary_files ${DELETE_TEMPS}
	else
		echo -e "ERROR: Check your inputs ; VCF files information is missing or erroneous; Aborting!" ; 1>&2
		fexit
	fi

}

