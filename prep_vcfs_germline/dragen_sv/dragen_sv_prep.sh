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
        echo -e "${RED}ERROR with: ${MSG}${ENDCOLOR}" ; echo "ERROR__ABORTING ! " 1>&2 ; exit 1 ;
    else
        echo -e "${MSG}: ${GREEN}OK! ${ENDCOLOR}"  1>&2;
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
    if [[ ${ALL_OK} -eq 0 ]] ; then echo -e "${RED}File(s) NOT FOUND; Aborting${ENDCOLOR}; Check your inputs"  1>&2 ; exit 1 ; fi
}

function usage(){
    echo -e "\nUSAGE: $0  1__REGIONS_BED_FILE[FILE]  2__VCF_SNV_FILE[FILE] 3__SNAME_in_10th_column_VCF[STRING] 4__outFileName.vcf.gz[FILE_PATHNAME/FULL/REL/BASENAME]"
    echo ""
}

## python script that MUST be in the PATH or we should the path location
PYTHON_SCRIPT_DEDUP_DRAGEN_SV_CALLS="main_dedup_sv.py"

VCF_SV=$1

ENOA=1
if [[ $# -ne ${ENOA} ]] ; then echo -e "Expected ${ENOA} args found $# ; Aborting;" ; usage ; exit 1 ; fi

VCF_PREPPED_OUTFILENAME=${VCF_SV%.*}.prep.vcf

## Command for preprocessing SV calls VCF:
## Preprocessing SV VCF
## Adding GT to all the filtered calls
## Determining the GT value from the threshold defined at the beginning by the user or based on the default value (see --threshold_GT option in vcfMerger2.py)
## Adding with DP, AR, AD to FORMAT using the given PR and SR values

## NOTE: HARDCODED VALUE for THRESHOLD GT: 0.9
## Why all these lines? Because the smpl_sum() function did not have the expected behavior when trying to calculate the DP and later on the AR tags particularly if the SR tag is
## absent from the record
## so we tricked the system and set manually the DP1 and DP2 from PR first, then SR respectively.
## then we set the DP value DP1+DP2
## and finally we are able to calculate the AR;
## Then to set the GT with using the setGT command, we need VAF and collect it from AR we calculted earlier
## This can probably be shorten to less lines with future bcftools version (got tested with bcftools v1.16, v1.17 & v1.18)
( bcftools view --threads 2 "${VCF_SV}" -h | \
grep -vE "#CHROM" ; echo -e '##FORMAT=<ID=GT,Number=.,Type=String,Description="Genotype">' ; \
bcftools view --threads 2 "${VCF_SV}" -h | \
grep -P "#CHROM\tPOS"  ; \
bcftools view -H "${VCF_SV}" | \
awk '{FS=OFS="\t" ; $9="GT:"$9 ; $10="./.:"$10 ; print }' | \
awk '{FS=OFS="\t" ; if($9 !~ /SR/ ){ $9=$9":SR" ; $10=$10":0,0"} ; print }' ) | \
bcftools +fill-tags - -- -t  'FORMAT/DP1=int(smpl_sum(FORMAT/PR) )' | \
bcftools +fill-tags - -- -t  'FORMAT/DP2=int(smpl_sum(FORMAT/SR) )' | \
bcftools +fill-tags - -- -t  'FORMAT/DP=int(DP1+DP2)' | \
sed 's/Description="Added by +fill-tags expression FORMAT\/DP=int(DP1+DP2)"/Description="Variant depth calculated using PR, SR or PR+SR when tags available"/  ; s/##FORMAT=<ID=DP,Number=.,
/##FORMAT=<ID=DP,Number=1,/' | \
bcftools +fill-tags - -- -t  'FORMAT/AR=float((PR[0:1]+SR[0:1])/DP)' | \
sed 's/Description="Added by +fill-tags expression FORMAT\/AR=float((PR[0:1]/Description="Variant Allelic Ratio FORMAT\/AR=float((PR[0:1]/ ; s/##FORMAT=<ID=AR,Number=.,
/##FORMAT=<ID=AR,Number=1,/' | \
bcftools +fill-tags - -- -t  'FORMAT/AD=int(PR+SR)' | \
sed 's/Description="Added by +fill-tags expression FORMAT\/AD=int(PR/Description="Variant Allelic Depth FORMAT\/AD=int(PR/' | \
bcftools +fill-tags - -- -t  'FORMAT/VAF=AR' | \
sed 's/Description="Added by +fill-tags expression FORMAT\/VAF=AR/Description="Variant Allele Frequency FORMAT\/VAF=AR/ ; s/##FORMAT=<ID=VAF,Number=.,
/##FORMAT=<ID=VAF,Number=1,/' | \
bcftools +setGT - -- -t ./. -n c:'1/1' | \
bcftools +setGT - -- -t q  -n c:'0/1' -i "FMT/VAF<0.90"  | \
bcftools annotate -x FORMAT/DP1,FORMAT/DP2 -O v -o "${VCF_PREPPED_OUTFILENAME}"  ;

check_ev $? "bcftools SV prep"  1>&2
## TODO: if sed lines slow down the process, we can later have only one sed line using `-e` and replace all in one command ;
## But this means we will have to be more specific in the search and replace of the strings

##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## DEDUP STEP for SV calls
##@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## check now if the sv vcf file has lines of variants with the same position
## for instance, we can have the same position 28034157 position called twice by the vairant caller, but one is defined has a DRAGEN::INS and the other has DRAGEN::DUP; both of them refer has an insertion; The goal is to make sure that these teo lines do not represent the exact same call;
COUNT_DUP_POS=$( bcftools query -f '%CHROM:%POS\n' ${VCF_PREPPED_OUTFILENAME} | sort --parallel=2 | uniq -d | wc -l )

VCF_PREPPED_OUTFILENAME_DEDUP=${VCF_PREPPED_OUTFILENAME/.vcf/.dedup.vcf}

if [[ ${COUNT_DUP_POS} -ne 0 ]]
then
    echo -e "more than one variant at same position have been found ... let's process the VCF to determine if the same variant has been called twice or if the variants found at
    the same position are in fact two different variants ..."  1>&2
    echo -e "running python script tcl_check_duppos_in_VCF" 1>&2

    ${PYTHON_SCRIPT_DEDUP_DRAGEN_SV_CALLS} -i "${VCF_PREPPED_OUTFILENAME}" -o "${VCF_PREPPED_OUTFILENAME_DEDUP}"
    check_ev $? "python dedup SV VCF"
else
    echo -e "No duplicated positions has been found in VCF. We did not parse or update the VCF"  1>&2
    echo -e "for consistency, we are creating the expected out file here by copying as is the file << ${VCF_PREPPED_OUTFILENAME} >> to << ${VCF_PREPPED_OUTFILENAME_DEDUP} >> "  1>&2

    cp "${VCF_PREPPED_OUTFILENAME}" "${VCF_PREPPED_OUTFILENAME_DEDUP}"
    check_ev $? "copy VCF"  1>&2

fi

echo ${VCF_PREPPED_OUTFILENAME_DEDUP}