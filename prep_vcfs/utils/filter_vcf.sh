#!/usr/bin/env bash

SNPSIFT_JAR_PATH=$1
CODE_EXT_FILT=$2  # can be either "filt" or "pass"
STRING_FOR_SNPSIFT_FILTERING="$3"  ## double quote mandatory here
VCF=$4

if [[ ${VCF##*.} == "vcf" ]] ; then
    VCF_OUT=$(basename ${VCF} ".vcf").${CODE_EXT_FILT}.vcf
elif [[ ${VCF##*.} == "gz" ]] ; then
    VCF_OUT=$(basename ${VCF} ".vcf.gz").${CODE_EXT_FILT}.vcf
else
    echo -e "ERROR: Unexpected file extension. File Extension is neither .vcf or .vcf.gz; \
    Please Comply to the expected file extension of the VCF file; Aborting."
   exit 1
fi

set -eu

zcat -f ${VCF} | java -jar ${SNPSIFT_JAR_PATH} filter " ${STRING_FOR_SNPSIFT_FILTERING} " > ${VCF_OUT}
if [[ $? -ne 0 ]] ;
then
    echo -e "ERROR: Filtering VCF step FAILED;\n Command  << zcat -f ${VCF} | java -jar ${SNPSIFT_JAR_PATH} filter ' ${STRING_FOR_SNPSIFT_FILTERING} ' > ${VCF_OUT} >> FAILED"
    exit 1
fi