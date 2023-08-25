#!/usr/bin/env bash

##@@@@@@@@@@@@@@@
##FUNCTIONS
##@@@@@@@@@@@@@@@
function check_ev(){
    EV=$1 ; MSG=$2
    RED="\e[31m"
    GREEN="\e[32m"
    ENDCOLOR="\e[0m"
    if [ ${EV} -ne 0 ] ;
    then
        echo -e "${RED}ERROR with: ${MSG}${ENDCOLOR}" ; echo "ERROR__ABORTING ! "  1>&2 ; exit 1 ;
    else
        echo -e "${MSG}: ${GREEN}OK! ${ENDCOLOR}" 1>&2 ;
    fi
}

function check_file(){
    RED="\e[31m"
    GREEN="\e[32m"
    ENDCOLOR="\e[0m"
    ALL_OK=1
    for F in $*
    do
        if [[ ! -e $F ]] ; then echo -e "${RED}ERROR: FNF ${ENDCOLOR} for << ${F} >> "  1>&2 ; ALL_OK=0; fi
        echo -e "File ${GREEN}FOUND${ENDCOLOR}:  ${F}"  1>&2
    done
    if [[ ${ALL_OK} -eq 0 ]] ; then echo -e "${RED}File(s) NOT FOUND; Aborting${ENDCOLOR}; Check your inputs" ; exit 1 ; fi
}

function usage(){
    echo -e "\nUSAGE: $0  1__VCF_SNV_FILE[FILE] 2__SNAME_in_10th_column_VCF[STRING] 3__outFileName.vcf.gz[FILE_PATHNAME/FULL/REL/BASENAME] "
    echo ""
}


VCF_SNV=$1

check_file "${VCF_SNV}"  1>&2

ENOA=1
if [[ $# -ne ${ENOA} ]] ; then echo -e "Expected ${ENOA} args found $# ; Aborting;" ; usage ; exit 1 ; fi

VCF_PREPPED_OUTFILENAME=${VCF_SNV%.*}.prep.vcf

## Command for preprocessing SNV calls VCF:
bcftools view --threads 2 -O v "${VCF_SNV}"  | \
bcftools +fill-tags - -- -t 'FORMAT/AR=AF' | \
bcftools view -O v -o "${VCF_PREPPED_OUTFILENAME}"
check_ev $? "bcftools snv prep"  1>&2

echo ${VCF_PREPPED_OUTFILENAME}