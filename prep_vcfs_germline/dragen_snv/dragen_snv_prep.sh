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

function check_given_value(){
    VARNAME=$1
    VARVALUE=$2

    if [[ "${VARVALUE}" != "yes" && ${VARVALUE} != "no" ]]
    then
        echo -e "ERROR: Expected values are 'yes' or 'no' (case insensitive) for the variable ${VARNAME}; Aborting;"
        exit 1
    fi

}

function usage(){
    echo -e "USAGE:   $0  1__ADD_SVLEN?[STRING:valuespossible:yes/no] 2__KEEP_TEMP?[String:valuespossible:yes/no]  3__VCF_SNV_FILE[FILE]  "
    echo -e "Example: $0  yes no  \${MY_DRAGEN_SNV_VCF} "
    echo -e "Note: the output file, i.e. prepped vcf is an uncompressed vcf."
}

ADD_SVLEN_TO_VCF=$1
KEEP_TEMP_FILE=$2
VCF_SNV=$3
#VCF_PREPPED_OUTFILENAME=$4

check_file "${VCF_SNV}"  1>&2

ENOA=3
if [[ $# -ne ${ENOA} ]] ; then echo -e "ERROR:  Expected ${ENOA} args found $# ; Aborting;" ; usage ; exit 1 ; fi

VCF_PREPPED_OUTFILENAME=${VCF_SNV%.*}.prep.vcf

set -euo pipefail

ADD_SVLEN_TO_VCF=$( echo "$1" | tr '[A-Z]' '[a-z]')
KEEP_TEMP_FILE=$( echo "$2" | tr '[A-Z]' '[a-z]')
check_given_value "\$1_ADD_SVLEN_TO_VCF" ${ADD_SVLEN_TO_VCF}
check_given_value "\$2_KEEP_TEMP_FILE" ${KEEP_TEMP_FILE}


if [[ ${ADD_SVLEN_TO_VCF} == "yes" ]]
then
	# Commands for preprocessing DRAGEN SNV calls VCF:
	# Adding SVLEN implies that the ALT sequence is given in its entirety ; The script is expecting this format; if ALT are represented by something else such as ...
	# let's say <INS>, that case has not been taken into account here since in ALL the dragen_snv vcf created this case has never been seen or encountered.
	echo -e "running bcftools query ..."   1>&2
	bcftools query -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\n' "${VCF_SNV}" | awk '{FS=OFS="\t" ; L=length($5)-length($4) ; if(length($5)==length($4)){L=1} ; print $0,L}' | bgzip -c > "${VCF_SNV}.temp_for_anno_svlen.tsv.gz"
	check_ev $? "bcftools query to capture svlen"  1>&2

	tabix -s 1 -b 2 -e 3 ${VCF_SNV}.temp_for_anno_svlen.tsv.gz
	check_ev $? "tabix on temp file for anno"  1>&2
	# making the annotation file for the header in vcf
	echo -e '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' > ${VCF_SNV}.annots.hdr

	echo -e "running bcftools annotate ..."   1>&2
	bcftools annotate -a ${VCF_SNV}.temp_for_anno_svlen.tsv.gz  -h ${VCF_SNV}.annots.hdr -c CHROM,POS,POS,REF,ALT,SVLEN -O v -o "${VCF_SNV}.anno.temp.vcf" "${VCF_SNV}"
	check_ev $? "bcftools snv prep"  1>&2

	echo -e "running bcftools fill-tags ..."   1>&2
	bcftools +fill-tags "${VCF_SNV}.anno.temp.vcf" -O v -o "${VCF_PREPPED_OUTFILENAME}" -- -t 'FORMAT/AR=AF'
	check_ev $? "bcftools snv prep"  1>&2

	if [[ ${KEEP_TEMP_FILE} == "no" ]]
	then
			echo -e "removing temp files ..." 1>&2
			rm "${VCF_SNV}.anno.temp.vcf" "${VCF_SNV}.annots.hdr" "${VCF_SNV}.temp_for_anno_svlen.tsv.gz" "${VCF_SNV}.temp_for_anno_svlen.tsv.gz.tbi" 1>&2 2> /dev/null
	fi
else
    echo -e "skipped adding svlen"  1>&2
	echo -e "running bcftools fill-tags ..."  1>&2
	bcftools +fill-tags "${VCF_SNV}" -O v -o "${VCF_PREPPED_OUTFILENAME}" -- -t 'FORMAT/AR=AF'
	check_ev $? "bcftools snv prep"  1>&2
fi

## return value of this script is an uncompressed VCF file:
echo ${VCF_PREPPED_OUTFILENAME}
