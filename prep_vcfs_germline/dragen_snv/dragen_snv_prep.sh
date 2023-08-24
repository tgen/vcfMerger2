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
        echo -e "${RED}ERROR with: ${MSG}${ENDCOLOR}" ; echo "ERROR__ABORTING ! " ; exit 1 ;
    else
        echo -e "${MSG}: ${GREEN}OK! ${ENDCOLOR}" ;
    fi
}

function check_file(){
    RED="\e[31m"
    GREEN="\e[32m"
    ENDCOLOR="\e[0m"
    ALL_OK=1
    for F in $*
    do
        if [[ ! -e $F ]] ; then echo -e "${RED}ERROR: FNF ${ENDCOLOR} for << ${F} >> " ; ALL_OK=0; fi
        echo -e "File ${GREEN}FOUND${ENDCOLOR}:  ${F}"
    done
    if [[ ${ALL_OK} -eq 0 ]] ; then echo -e "${RED}File(s) NOT FOUND; Aborting${ENDCOLOR}; Check your inputs" ; exit 1 ; fi
}

function usage(){
    echo -e "\nUSAGE: $0  1__REGIONS_BED_FILE[FILE]  2__VCF_SNV_FILE[FILE] 3__SNAME_in_10th_column_VCF[STRING] 4__outFileName.vcf.gz[FILE_PATHNAME/FULL/REL/BASENAME]"
    echo ""
}


BED_REGIONS_FILE=$1
VCF_SNV=$2
SNAME=$3
OUTFILENAME_VCF=$4

ENOA=4
if [[ $# -ne ${ENOA} ]] ; then echo -e "Expected ${ENOA} args found $# ; Aborting;" ; usage ; exit 1 ; fi

VCF_PREPPED_OUTFILENAME=${OUTFILENAME_VCF}

## Command for preprocessing SNV calls VCF:
bcftools filter --regions-file "${BED_REGIONS_FILE}"  -i 'TYPE!="snp"' --threads 2 "${VCF_SNV}"  | \
bcftools +fill-tags - -- -t 'FORMAT/AR=AF' | \
bcftools view -O z -o "${VCF_PREPPED_OUTFILENAME}"
check_ev $? "bcftools snv prep"

bcftools index --force "${VCF_PREPPED_OUTFILENAME}"
check_ev $? "bcftools snv prep index"