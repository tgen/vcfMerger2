#!/usr/bin/env bash

## CONSTANT VARIABLE (modified accordingly)
DIR_PATH_TO_SCRIPTS=$(dirname $( dirname $(readlink -f ${BASH_SOURCE[0]}) )  )
echo -e "${DIR_PATH_TO_SCRIPTS}/prep_vcfs/prep_vcf_functions.sh" 1>&2
source ${DIR_PATH_TO_SCRIPTS}/prep_vcfs/prep_vcf_functions.sh  ## allows to load functions and reused them; below we write function we want to OVERWRITE from the sourced file
DIR_PATH_TO_SCRIPTS="$( dirname $(readlink -f ${BASH_SOURCE[0]}) )"


function init_some_vars(){
	LI="RECAP_INPUTS_USED:"
	TOOLNAME=""
	VCF_ALL_CALLS=""
	NORMAL_SNAME=""
	TUMOR_SNAME=""
	GERMLINE_SNAMES=""
	VCF_INDELS_FILE=""
	VCF_SNVS_FILE=""
	BAM_FILE=""
	NORMALIZE="yes"
	TARGETS_BED_FILE_FTT=""
	CONTIGS_FILE=""
	VCF_FINAL_USER_GIVEN_NAME=""
	TH_AR=""
	MAKE_BED_FOR_VENN="no"
}

function getOptions(){
# options may be followed by one colon to indicate they have a required argument
# NOTE: long option are not working on MacOS Sierra
if ! options=`getopt -o hd:b:g:o:t: -l help,dir-work:,ref-genome:,germline-snames:,vcf:,toolname:,prepped-vcf-outfilename:,bam:,contigs-file:,print-valid-toolnames,do-not-normalize,make-bed-for-venn -- "$@" `
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
		--vcf) export VCF_ALL_CALLS="$2" ; LI="${LI}\nVCF_ALL_CALLS==\"${VCF_ALL_CALLS}\"";  shift ;;
		--germline-snames) export GERMLINE_SNAMES="$2" ; LI="${LI}\nGERMLINE_SNAMES==\"${VCF_ALL_CALLS}\"";  shift ;;
		--bam) export BAM_FILE=$2 ; LI="${LI}\nBAM_FILE==\"${BAM_FILE}\"";  shift ;;
		-g|--ref-genome) export REF_GENOME_FASTA="${2}" ;  LI="${LI}\nREF_GENOME_FASTA==\"${REF_GENOME_FASTA}\"";  shift ;; ## FASTA FILE ; .fai file is needed
		-t|--toolname) export TOOLNAME=$2 ; LI="${LI}\nTOOLNAME==\"${TOOLNAME}\"";  shift ;;
		--do-not-normalize) export NORMALIZE="no" ; LI="${LI}\nNORMALIZE==\"${NORMALIZE}\"" ;;
		--contigs-file) export CONTIGS_FILE="$2" ; LI="${LI}\nCONTIGS_FILE==\"${CONTIGS_FILE}\"";  shift ;; ## File containing the contigs in the same format as expected within a VCF file
		--print-valid-toolnames) echo ${VALID_TOOLNAMES} ; exit ;; ## print possible toolnames to be used with the
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


function check_and_update_sample_names(){
	##@@@@@@@@@@@@@@@@@@@@@@@
	## CHECK SAMPLES NAMES
	##@@@@@@@@@@@@@@@@@@@@@@@
	## so far the sample name should be NORMAL, TUMOR or the SM tag from BAM file ;
	## if not ##TODO revise the way we handle the sample names to make it much more generic;
	## will need discussion to do so, and will need to capture all the use cases possible ;

	local VCF=$1
	local LSNAMES=( $( echo -e "$2" | sed 's/\|/ /g')  ) ### list of samples in order the user want them to be
	local VCF_OUT=$(basename ${VCF} ".vcf").sname.vcf
	SNAMES_IN_VCF=( $( zcat -f ${VCF} | grep -m 1 "#CHROM"  | cut -f10- | sed 's/\t/ /' )  )
	HEADERLINE_VCF="$( zcat -f ${VCF} | grep -m 1 "#CHROM"  | cut -f1-9 | sed 's/\t/ /'  )"

    echo " in ${funcname} :  ${VCF} ${TOOLNAME} ${SNAMES_IN_VCF[@]} " 1>&2
    if [[ ${#SNAMES_IN_VCF[@]} -ne ${#LSNAMES[@]} ]] ;
    then
        echo -e "ERROR: number of sample name(s) given by user and number of sample(s) found in VCF are DIFFERENT ; Aborting ;"
        exit -1
    fi
    if [[ ${#SNAMES_IN_VCF[@]} -eq 1  ]] ;
    then
        echo -e "Found only 1 Sample in List, therefore checking if name matches user-given name ; if not we update name" ;
        for I in `seq 0 $((${#SNAMES_IN_VCF[@]}-1))`
        do
            if [[ ${SNAMES_IN_VCF[I]} -ne ${LSNAMES[I]} ]] ;
            then
                ${SNAMES_IN_VCF[I]}=${LSNAMES[I]}
            fi
        done
        zcat -f ${VCF} | sed "/^##CHROM/s/\t$( echo -e "${SNAMES_IN_VCF[@]}" | sed 's/ \+/\t/g ; s/\t$// ') /\t$( echo -e "${LSNAMES[@]}" | sed 's/ \+/\t/g ; s/\t$// ')/" > ${VCF_OUT}
        echo ${VCF_OUT}
    elif [[ ${#SNAMES_IN_VCF[@]} -gt 1  ]] ;
    then
        ##TODO: check and manage more than one sample in Germline VCF; Swapping column around might be necessary if sample do not output the sample in the same order (mostly becaue they do not use the same name either)
        echo -e "Found more than 1 Sample in List, therefore checking if name matches user-given name ;"
        echo -e "if name is not present in vcf, we raise an exception or we assume the user-given list is for renaming the sample names in vcf and that list is in correct order; Let's keep that in min for a future #TODO" ;
        echo -e "Normally the calls should have been done on the same BAM file, so if variant caller follow specs, the same name should have been used for the Sample Names."
        echo -e "but as usual, there is probably no harmony ..."

        ## HERE we ASSUME that the user-given list of germline sample names is for renaming only and the order of the sample is ALREADY the same for ALL vcf
#        for I in `seq 0 $((${#SNAMES_IN_VCF[@]}-1))`
#        do
#            CURNAME=${SNAMES_IN_VCF[I]}
#            SN_IN_VCF=false
#            POSITION_J=0
#            for J in `seq 0 $((${#LSNAMES[@]}-1))`
#            do
#                if [[ ${SNAMES_IN_VCF[I]} -eq ${LSNAMES[J]} ]] ;
#                then
#                    SN_IN_VCF=true
#                    POSITION_J=${J}
#                    break
#                fi
#            done
#            if [[ ${SN_IN_VCF} == "true" ]] ;
#            then
#                if [[ ${I} == ${POSITION_J} ]] ; then echo -e "sample names are in same column ; no modification of headerline line in vcf necessary" ;
#                else
#                    echo -e "we need to update the column by moving it to position $I"
#                    ##TODO
#                fi
#            fi
#        done
        zcat -f ${VCF} > ${VCF_OUT}
        echo ${VCF_OUT}
    fi

}


function make_vcf_upto_specs_for_VcfMerger_Germline(){
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	## updating ${TOOLNAME} VCF to specs for vcfMerger2_Germline
	##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	local VCF="$1"
	#local LSNAMES="${2}"  ## must be in final wanted order
	fout_name=${VCF%.*}.prep.vcf
	echo -e "## prep vcf for vcfMerger2 Germline ..." 1>&2
	echo -e "## VCF == ${VCF}" 1>&2
	echo -e "## outVCF == ${fout_name}" 1>&2

	mycmd="python ${PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER} -i ${VCF} --outfilename ${fout_name}"
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



function process_haplotypecaller_vcf(){
    local VCF=${1}
    VCF=$( check_and_update_sample_names ${VCF} )
    VCF=$( make_vcf_upto_specs_for_VcfMerger_Germline ${VCF}  )
}

function process_deepvariant_vcf(){
    local VCF=${1}
    VCF=$( check_and_update_sample_names ${VCF} )
    VCF=$( make_vcf_upto_specs_for_VcfMerger_Germline ${VCF}  )
}

function process_freebayes_vcf(){
    local VCF=${1}
    VCF=$( check_and_update_sample_names ${VCF} )
    VCF=$( make_vcf_upto_specs_for_VcfMerger_Germline ${VCF}  )
}

function process_samtools_mpileup_vcf(){
    local VCF=${1}
    VCF=$( check_and_update_sample_names ${VCF} )
    VCF=$( make_vcf_upto_specs_for_VcfMerger_Germline ${VCF}  )
}

function process_octopus_vcf(){
    local VCF=${1}
    VCF=$( check_and_update_sample_names ${VCF} )
    VCF=$( make_vcf_upto_specs_for_VcfMerger_Germline ${VCF}  )
}

function process_strelka2_vcf(){
    local VCF=${1}
    VCF=$( check_and_update_sample_names ${VCF} )
    VCF=$( make_vcf_upto_specs_for_VcfMerger_Germline ${VCF}  )
}


function run_tool(){
    local TOOLNAME=$( echo $1 | tr '[A-Z]' '[a-z]' | tr ' ' '_' )
    local VCF="$2"
	local DIR_PATH_TO_SCRIPTS="${DIR_PATH_TO_SCRIPTS}"

    case $TOOLNAME in
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
        deepvariant|dv|dvt)
            PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/deepvariant/deepvariant.germline.1s.addFieldsForVcfMerger.py"
            process_deepvariant_vcf ${VCF}
        ;;
        octopus|oct)
            PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/octopus/octopus.germline.1s.addFieldsForVcfMerger.py"
            process_octopus_vcf ${VCF}
        ;;
        strelka2|slk)
            PYTHON_SCRIPT_PREP_VCF_FOR_VCFMERGER="${DIR_PATH_TO_SCRIPTS}/octopus/strelka2.germline.1s.addFieldsForVcfMerger.py"
            process_strelka2_vcf ${VCF}
        ;;

		(*) echo -e "\nERROR: unrecognized toolname:  $1\nERROR: valid toolnames are << ${VALID_TOOLNAMES} >>" 1>&2 ; fexit  ;;
	esac
	echo "${SCRIPT_FP}" ## we return the full path of the script
}

function main(){

	##@@@@@@@@@@@@@@@@@##
	##  check inputs   ##
	##@@@@@@@@@@@@@@@@@##
#	echo -e "## Checking inputs ..."  1>&2
#	if [[ ${TOOLNAME} == ""  ]] ; then echo -e "ERROR: --toolname has to be provided ; Aborting." 1>&2 ; fexit ; fi
#	if [[ ${REF_GENOME_FASTA} == ""  ]] ; then echo -e "ERROR: reference genome option is required ; here we have a missing value; provide --ref-genome ; Aborting." 1>&2 ; fexit ; fi
#	if [[ ${GERMLINE_SNAMES} == ""  ]] ; then echo -e "ERROR: GERMLINE SAMPLE NAMES MUST be provided in the same order as present in all VCF files; REQUIREMENT: \
#	if germline VCFs have more than sample after column 9 in vcf, the order of these sample (aka column order) MUST be IDENTICAL even though their name might not be \
#	 ; Missing Values Found ; Check your inputs;  Aborting."  1>&2 ; fexit ; fi


	## check files and folders if exist
	checkDir ${DIR_WORK}
	checkFile ${VCF}
#	checkFile ${REF_GENOME_FASTA}

	cd ${DIR_WORK}

	## check if we deal with indels and snvs in separated vcf or in all-in-one vcf
	if [[ ${VCF_ALL_CALLS} != "" ]] ;
	then
		checkFile ${VCF_ALL_CALLS}
		if [[ ! -e $( basename ${VCF_ALL_CALLS}) ]]
		then
		    echo -e "CREATING SYMLINK in CURR DIR ${PWD}"
		    ln -sf ${VCF_ALL_CALLS} $(basename ${VCF_ALL_CALLS}) &>/dev/null ## we create a symlimk in current working directory (in case original vcf folder is not writable)
		fi
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


