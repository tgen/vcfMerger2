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

## python script that MUST be in the PATH or we should the path location
PYTHON_SCRIPT_DEDUP_DRAGEN_SV_CALLS="/home/clegendre/qscripts/gits/tcl_check_duppos_in_vcf/main_dedup_sv.py"

BED_REGIONS_FILE=$1
VCF_SV=$2
SNAME=$3
OUTFILENAME_VCF=$4

ENOA=4
if [[ $# -ne ${ENOA} ]] ; then echo -e "Expected ${ENOA} args found $# ; Aborting;" ; usage ; exit 1 ; fi

VCF_PREPPED_OUTFILENAME=${OUTFILENAME_VCF}

## Command for preprocessing SV calls VCF:
## Preprocessing SV VCF
## Adding GT to all the filtered calls
## Determining the GT value from the threshold defined at the beginning by the user or based on the default value (see --threshold_GT option in vcfMerger2.py)
## Adding with DP, AR, AD to FORMAT using the given PR and SR values


( bcftools view --threads 2 "${VCF_SV}" -h | \
grep -vE "#CHROM" ; echo -e '##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">' ; \
bcftools view --threads 2 "${VCF_SV}" -h | \
grep -P "#CHROM\tPOS"  ; \
bcftools filter --threads 2 -i 'TYPE!="snp" || TYPE=="indel" || TYPE=="other" ' --regions-file "${BED_REGIONS_FILE}" "${VCF_SV}" | \
bcftools view -H | \
awk '{FS=OFS="\t" ; $9="GT:"$9 ; $10="./.:"$10 ; print }' ) | \
bcftools +fill-tags - -- -t 'FORMAT/DP=int(smpl_sum(PR+SR)),FORMAT/AR=float((PR + SR)/smpl_sum(PR+SR)),FORMAT/AD=int(PR+SR)' | \
bcftools +setGT -O z -o "${VCF_SV_TEMPFILE}"  - -- -t a -n m
check_ev $? "bcftools SV prep"

## indexing newly created vcf.gz file
bcftools index --force "${VCF_SV_TEMPFILE}"
check_ev $? "bcftools SV prep index"


##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## DEDUP STEP for SV calls
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## check now if the sv vcf file has lines of variants with the same position
## for instance, we can have the same position 28034157 position called twice by the vairant caller, but one is defined has a DRAGEN::INS and the other has DRAGEN::DUP; both of them refer has an insertion; The goal is to make sure that these teo lines do not represent the exact same call;
COUNT_DUP_POS=$( bcftools query -f '%CHROM:%POS\n' ${VCF_SV_TEMPFILE} | sort --parallel=2 | uniq -d | wc -l )
VCF_SV_TEMPFILE_DEDUP=${VCF_SV_TEMPFILE/.vcf.gz/.dedup.vcf}
if [[ ${COUNT_DUP_POS} -ne 0 ]]
then
    echo -e "more than one variant at same position have been found ... let's process the VCF to determine if the same variant has been called twice or if the variants found at the same positionare in fact two different variants ..."
    echo -e "running python script tcl_check_duppos_in_VCF"

    python3 ${PYTHON_SCRIPT_DEDUP_DRAGEN_SV_CALLS} -i "${VCF_SV_TEMPFILE}" -o "${VCF_SV_TEMPFILE_DEDUP}"
    check_ev $? "python dedup SV VCF"

    echo -e "converting vcf to vcf.gz ......"
    mycmd="bcftools view --threads 2 -O z -o ${VCF_SV_TEMPFILE_DEDUP}.gz ${VCF_SV_TEMPFILE_DEDUP}"
    echo "${mycmd}"
    eval ${mycmd}
    check_ev $? "conversion vcf to vcf.gz"
    bcftools index --force --threads 2 --csi "${VCF_SV_TEMPFILE_DEDUP}.gz"
    check_ev $? "bcftools SV dedup index"
else
    echo -e "No duplicated positions found in VCF. We did not parse or update the VCF"
    echo -e "for consistency, we are creating the expected out file with dedup information here by copying as is the file << ${VCF_SV_TEMPFILE} >> to << ${VCF_SV_TEMPFILE_DEDUP} >> "

    cp "${VCF_SV_TEMPFILE}" "${VCF_SV_TEMPFILE_DEDUP}.gz"
    check_ev $? "copy VCF"

    bcftools index --force --threads 2 --csi "${VCF_SV_TEMPFILE_DEDUP}.gz"
    check_ev $? "bcftools SV dedup index"
fi